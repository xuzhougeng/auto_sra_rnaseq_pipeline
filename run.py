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

from concurrent.futures import ThreadPoolExecutor, as_completed

def process_sample_file(sample_file, metdata_dir, sf, config_file, cores, config):
    shutil.move(sample_file, metdata_dir)
    run_snakemake(sf, config_file, cores, unlock=True)
    status = run_snakemake(sf, config_file, cores)
    if status:
        contents = "snakemake run successfully"
        if config['bark']:
            bark_notification(config['bark_api'], contents)
        if config['feishu']:
            feishu_notification(config['feishu_api'], contents)

        if os.path.isfile(os.path.join(finished_dir, os.path.basename(sample_file))):
            shutil.copy2(sample_file, finished_dir)
            os.unlink(sample_file)
        else:
            shutil.move(sample_file, finished_dir)
    else:
        contents = "snakemake run failed"
        if config['bark']:
            bark_notification(config['bark_api'], contents)
        if config['feishu']:
            feishu_notification(config['feishu_api'], contents)
        shutil.move(sample_file, "unfinished")

def main(root_dir, args):
    import yaml
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
    metadata_files = glob.glob(os.path.join(unfinished_dir, "*.txt"))
    sf = get_snakefile(root_dir, args[5] if len(args) > 5 else "Snakefile")

    # Read the base config from config.yaml
    with open(config_file[0], 'r') as stream:
        base_config = yaml.safe_load(stream)

    task_dict = {}
    for metadata_file in metadata_files:
        task_name = os.path.basename(metadata_file)
        task_config = base_config.copy()
        task_config['metadata'] = metadata_file
        task_dict[task_name] = task_config

    def submit_task(task_name, task_config):
        config_path = os.path.join(metdata_dir, f"{task_name}_config.yaml")
        with open(config_path, 'w') as config_out:
            yaml.dump(task_config, config_out)
        run_snakemake(sf, [config_path], cores, unlock=True)
        status = run_snakemake(sf, [config_path], cores)
        if status:
            contents = "snakemake run successfully"
            if config['bark']:
                bark_notification(config['bark_api'], contents)
            if config['feishu']:
                feishu_notification(config['feishu_api'], contents)
            if os.path.isfile(os.path.join(finished_dir, task_name)):
                shutil.copy2(metadata_file, finished_dir)
                os.unlink(metadata_file)
            else:
                shutil.move(metadata_file, finished_dir)
        else:
            contents = "snakemake run failed"
            if config['bark']:
                bark_notification(config['bark_api'], contents)
            if config['feishu']:
                feishu_notification(config['feishu_api'], contents)
            shutil.move(metadata_file, "unfinished")

    with ThreadPoolExecutor(max_workers=parallel) as executor:
        future_to_task = {executor.submit(submit_task, task_name, task_config): task_name for task_name, task_config in task_dict.items()}
        for future in as_completed(future_to_task):
            task_name = future_to_task[future]
            try:
                future.result()
            except Exception as exc:
                print(f'{task_name} generated an exception: {exc}')

if __name__ == '__main__':
    root_dir = os.path.dirname(os.path.abspath(__file__))
    args = sys.argv
    main(root_dir, args) 
