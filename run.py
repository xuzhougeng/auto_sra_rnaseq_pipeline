import os
import sys
import glob
import shutil
import pandas as pd
import json

from snakemake import snakemake

# arguments order
# 

def get_snakefile(root_dir = ".", file = "Snakefile"):
    sf = os.path.join(root_dir, file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file;  tried %s" %sf)
    return sf

# check the duplication metadata file
def check_duplication(file, file_dict, finished_dir):
    sample_files = glob.glob( os.path.join(finished_dir,  "*.txt") )
    if file in sample_files:
        sample_files.remove(file)
    existed_files = set()
    if os.path.exists(file_dict):
        with open("file_dict.json", "r") as f:
            file_dict = json.load(f)
        existed_files = set(file_dict.keys())
    
    if len(sample_files) > 0 :
        for f in sample_files:
            df = pd.read_csv(f, sep = "\t")
            dict_key = "_".join(sorted(df['GSM'].to_list()))
            existed_files.add(dict_key)
    
    df = pd.read_csv(file, sep = "\t")
    dict_key = "_".join(sorted(df['GSM'].to_list()))
    if dict_key in existed_files:
        return True
    else:
        False

def run_snakemake(Snakefile, configfiles, cores):
    status = snakemake(snakefile= Snakefile, 
                      configfiles = configfiles, 
                      cores = cores,
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

    finished_dir = "finished"

    if not os.path.exists("metadata"):
        os.makedirs("metadata")
    # finished dir
    if not os.path.exists(finished_dir):
        os.makedirs(finished_dir)
    # duplication metadir
    if not os.path.exists("duplication"):
        os.makedirs("duplication")
    
    file_dict = "file_dict.json"

    sf = get_snakefile(root_dir, "Snakefile")
    sample_files = glob.glob( os.path.join(unfinished_dir,  "*.txt") )
    todo_files = []

    while len(sample_files) > 0 :
        for i in range(min(10, len(sample_files))):
            file = sample_files.pop()
            if check_duplication(file, file_dict, finished_dir):
                shutil.move(file, "duplication")
            elif check_duplication(file, file_dict, unfinished_dir):
                shutil.move(file, "duplication")
            else:
                todo_files.append(file)
                shutil.move(file, "metadata")
        status = run_snakemake(sf, config_file, cores )
        # move the finished file to finished
        if status:
            finished_file = glob.glob( os.path.join("metadata",  "*.txt") )
            for i in range(min(parallel, len(finished_file))):
                file = finished_file.pop()
                shutil.move(file, finished_dir)
        else:
            for i in range(len(todo_files)):
                file = todo_files.pop()
                shutil.move(file, "unfinished")
            sys.exit(1)

if __name__ == '__main__':
    root_dir = os.path.dirname(os.path.abspath(__file__))
    args = sys.argv
    main(root_dir, args) 


