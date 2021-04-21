import os
from os.path import basename
import sys
import glob
import shutil
import pandas as pd
import json

from urllib import request, parse
from urllib.request import urlopen
from urllib.parse import quote
from urllib.parse import urljoin

from pandas.core import base

from snakemake import snakemake
from snakemake.io import load_configfile

# arguments order
# 


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

# remove the finished samples
def remove_finished(sample_files, finished_file, finished_dir):
    if not os.path.exists(finished_file):
        return sample_files
    
    with open(finished_file, "r") as f:
        file_dict = json.load(f)
    
    files_set = set([ os.path.basename(f) for f in file_dict.values() ])
    
    for file in sample_files:
        if os.path.basename(file) in files_set:
            shutil.move(file, finished_dir)
            sample_files.remove(file)
    return sample_files

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
            shutil.move(file, duplication_dir)
            sample_files.remove(file)
        else:
            dup_dict[dict_key] = file
    
    if len(dup_dict) > 0 :
        with open(dup_file, "w") as f:
            json.dump(dup_dict, f)

    return sample_files


# check whether new job  was finished or not ?
def check_finished(file, finished_dict):
    if not os.path.exists(finished_dict):
        return False
    with open(finished_dict, "r") as f:
        file_dict = json.load(f)
    
    if file in file_dict.keys():
        if os.path.basename(file) == os.path.basename(file_dict[file]):
            return True
    return False

def check_duplication(file, dup_file):
    dup_dict = {}
    if os.path.exists(dup_file):
        with open(dup_file, "r") as f:
            dup_dict = json.load(f)
    
    df = pd.read_csv(file, sep = "\t")
    dict_key = "_".join(sorted(df['GSM'].to_list()))

    if dict_key not in dup_dict.keys():
        dup_dict[dict_key] = file
        with open(dup_file, "w") as f:
            json.dump(dup_dict, f)
        return False
    if os.path.basename(file) == os.path.basename(dup_dict[dict_key]):
        return False
    return True


def bark_notification(api, contents):
    base_url = api
    content = quote(contents)
    full_url = urljoin(base_url,  content)
    urlopen(full_url)

def feishu_notification(api,contents):
    req =  request.Request(api, method="POST") # this will make the method "POST"
    req.add_header('Content-Type', 'application/json')
    data_dict = {
        "msg_type": "text",
        "content": {"text": "进展报告: " + contents}
    }
    data = json.dumps(data_dict).encode()
    resp = urlopen(req, data = data)
    return  resp


def get_snakefile(root_dir = ".", file = "Snakefile"):
    sf = os.path.join(root_dir, file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file;  tried %s" %sf)
    return sf

# hard decode total download speed
def run_snakemake(Snakefile, configfiles, cores):
    status = snakemake(snakefile= Snakefile, 
                      configfiles = configfiles, 
                      cores = cores,
                      resources= {"rx": 80},
                      force_incomplete = True,
                      latency_wait = 5,
                      attempt = 3)
    return status

# 
def main(root_dir, args):
    if len(args) < 3:
        print("python run.py unfinished config.yaml cores")
        sys.exit(1)

    unfinished_dir = args[1]
    config_file    = [args[2]]
    cores = 79 if len(args) <= 3 else int(args[3])
    parallel = 5
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
    
    # file_dcit record the finished file which generate by snakemake workflow
    file_dict = "file_dict.json"
    # file record the duplication
    dup_file = ".file_duplication.json"
    # preprocess
    sample_files = glob.glob( os.path.join(unfinished_dir,  "*.txt") )

    # remove finished samples
    sample_files = remove_finished(sample_files, file_dict, finished_dir)
    # remove duplication before running
    sample_files = remove_duplication(sample_files, dup_file, duplication_dir)

    sf = get_snakefile(root_dir, "Snakefile")
    
    #todo_files = []

    while len(sample_files) > 0 :
        for i in range(min(parallel, len(sample_files))):
            file = sample_files.pop()
            if check_finished(file, finished_dict=file_dict):
                shutil.move(file, finished_dir)
                #todo_files.append(os.path.join( finished_dir ,os.path.basename(file)) )
            elif check_duplication(file, dup_file=dup_file):
                shutil.move(file, duplication_dir)
                #todo_files.append(os.path.join( duplication_dir ,os.path.basename(file)) )
            else:
                shutil.move(file, metdata_dir)
                #todo_files.append(os.path.join( metdata_dir ,os.path.basename(file))  )

        status = run_snakemake(sf, config_file, cores )
        # move the finished file to finished

        if status:
            contents="snakemake run successfully"
            if config['bark']:
                bark_notification(config['bark_api'], contents)
            if config['feishu']:
                feishu_notification(config['feishu_api'], contents)
            finished_file = glob.glob( os.path.join(metdata_dir,  "*.txt") )
            for _ in range(min(parallel, len(finished_file))):
                file = finished_file.pop()
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
