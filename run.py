import os
from os.path import basename
import sys
import glob
import shutil
import pandas as pd
import json

from pandas.core import base

from snakemake import snakemake
from snakemake.io import load_configfile

from scripts.utilize import bark_notification, feishu_notification
from scripts.utilize import build_metadata_table, build_sample_table
from scripts.utilize import table_to_sql, table_from_sql, update_status

# check config
def check_config(config):
        
    if config['mail']:
        if 'sender' not in config.keys() or config['sender'] is None:
            print("sender should not be empty", file = sys.stderr)
            return False

        if  'sender_password' not in config.keys() or config['sender_password'] is None:
            print("sender_password should not be empty", file = sys.stderr)
            return False

        if  'mail_to' not in config.keys() or config['mail_to'] is None:
            print("mail_to should not be empty", file = sys.stderr)
            return False

    if config['bark']:
        if 'bark_api' not in config.keys() or config['bark_api'] is None:
            print("bark_api should not be empty", file = sys.stderr)
            return False

    return True

# remove the depulication samples aside loop
def remove_duplication(sample_files, dup_file = ".file_duplication.json", 
                       duplication_dir="duplication"):
    dup_dict = {}
    if os.path.exists(dup_file):
        with open(dup_file, "r") as f:
            dup_dict = json.load(f)

    for file in sample_files:
        df = pd.read_csv(file, sep = "\t")
        dict_key = "_".join(sorted(df['GSM'].to_list()))

        if dict_key in dup_dict.keys() and \
            os.path.basename(file) != os.path.basename(dup_dict[dict_key]):
            if os.path.isfile(os.path.join(duplication_dir, os.path.basename(file))):
                shutil.copy2(file, duplication_dir)
                os.unlink(file)
            sample_files.remove(file)
        else:
            dup_dict[dict_key] = file
    
    if len(dup_dict) > 0 :
        with open(dup_file, "w") as f:
            json.dump(dup_dict, f)

    return sample_files

def get_snakefile(root_dir = ".", file = "Snakefile"):
    sf = os.path.join(root_dir, file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file;  tried %s" %sf)
    return sf

# hard decode total download speed
def run_snakemake(snakefile, configfiles, cores, unlock=False):
    status = snakemake(snakefile= snakefile, 
                      configfiles = configfiles, 
                      cores = cores,
                      resources= {"rx": 80, 'limit_dump': 2, 'limit_merge': 2},
                      force_incomplete = True,
                      latency_wait = 5,
                      attempt = 3,
					  unlock=unlock)
    return status

# 
def main(root_dir, args):
    if len(args) < 3:
        print("python run.py unfinished config.yaml cores")
        sys.exit(1)

    unfinished_dir = args[1]
    config_file    = [args[2]]
    cores = 79 if len(args) <= 3 else int(args[3])
    parallel = 10 if len(args) <= 4 else int(args[4])
    # load config
    config = load_configfile(config_file[0])
    # set defualt parameter
    if 'mail' not in config.keys():
        config['mail'] = False
    if 'bark' not in config.keys():
        config['bark'] = False
    # check_config
    if not check_config(config):
        sys.exit(1)

    finished_dir = "finished"
    duplication_dir = "duplication"
    metdata_dir = "metadata"
    

    if not os.path.exists(metdata_dir):
        os.makedirs(metdata_dir)
    # finished dir
    if not os.path.exists(finished_dir):
        os.makedirs(finished_dir)
    # duplication metadir
    if not os.path.exists(duplication_dir):
        os.makedirs(duplication_dir)
    
    
    dup_file = ".file_duplication.json"
    # preprocess
    sample_files = glob.glob( os.path.join(unfinished_dir,  "*.txt") )
    # remove duplication before running
    sample_files = remove_duplication(sample_files, dup_file, duplication_dir)
    
    # select unfinished files
    db = "meta_info.sqlite3"
    if not os.path.exists(db):
        # create database if not exists 
        df = build_metadata_table(sample_files)
        df2 = build_sample_table(df, unfinished_dir)
        table_to_sql(df, table_name="meta", db=db)
        table_to_sql(df2, table_name="sample", db=db)
    else:
        # otherwise, upgrade the db
        df = build_metadata_table(sample_files, table_name="meta", db=db)
        if df.shape[0] > 0 and df.shape[1] > 0:
            table_to_sql(df, table_name="meta", db=db)
            df2 = build_sample_table(df, unfinished_dir)
            table_to_sql(df2, table_name="sample", db=db)

    # get the unfinished meta files
    df = table_from_sql( table_name = "meta", db = db )
    sample_files = df.loc[df['status'] == 0, 'meta_file'].to_list()
    sample_files = [os.path.join(unfinished_dir, f) for f in sample_files ]
    
    sf = get_snakefile(root_dir, args[5] if len(args) > 5 else "Snakefile")
    
    #todo_files = []

    while len(sample_files) > 0 :
        for i in range(min(parallel, len(sample_files))):
            file = sample_files.pop()
            shutil.move(file, metdata_dir)
            #todo_files.append(os.path.join( metdata_dir ,os.path.basename(file))  )

        run_snakemake(sf, config_file, cores, unlock=True )
        status = run_snakemake(sf, config_file, cores )
        # move the finished file to finished

        if status:
            contents="snakemake run successfully"
            if config['bark']:
                bark_notification(config['bark_api'], contents)
            if config['feishu']:
                feishu_notification(config['feishu_api'], contents)
            finished_file = glob.glob( os.path.join(metdata_dir,  "*.txt") )
            for _ in range(len(finished_file)):
                file = finished_file.pop()
                update_status(file, table_name="meta", db = db)
                
                if os.path.isfile(os.path.join(finished_dir, os.path.basename(file))):
                   shutil.copy2(file, finished_dir)
                   os.unlink(file)
                else:
                    shutil.move(file, finished_dir)
        else:
            contents = "snakemake run failed"
            if config['bark']:
                bark_notification(config['bark_api'], contents)
            if config['feishu']:
                feishu_notification(config['feishu_api'], contents)  
            broken_file = glob.glob( os.path.join(metdata_dir,  "*.txt") )         
            for _ in range(len(broken_file)):
                file = broken_file.pop()
                shutil.move(file, "unfinished")
            #sys.exit(1)

if __name__ == '__main__':
    root_dir = os.path.dirname(os.path.abspath(__file__))
    args = sys.argv
    main(root_dir, args) 
