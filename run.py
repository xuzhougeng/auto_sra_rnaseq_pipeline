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

def process_sample_file(sample_file, metadata_dir, sf, config_file, cores, config):
    # Ensure the sample file is moved to the metadata directory
    shutil.move(sample_file, os.path.join(metadata_dir, os.path.basename(sample_file)))
    
    # Update the config with the new metadata file path
    config['metadata'] = os.path.join(metadata_dir, os.path.basename(sample_file))
    
    # Run Snakemake to unlock any potential issues
    run_snakemake(sf, [config_file], cores, unlock=True)
    
    # Execute Snakemake with the provided configuration
    status = run_snakemake(sf, [config_file], cores)
    
    # Check the status and handle notifications accordingly
    if status:
        contents = "snakemake run successfully"
        if config.get('bark'):
            bark_notification(config['bark_api'], contents)
        if config.get('feishu'):
            feishu_notification(config['feishu_api'], contents)

        finished_sample_file = os.path.join("finished", os.path.basename(sample_file))
        shutil.move(config['metadata'], finished_sample_file)
    else:
        contents = "snakemake run failed"
        if config.get('bark'):
            bark_notification(config['bark_api'], contents)
        if config.get('feishu'):
            feishu_notification(config['feishu_api'], contents)
        shutil.move(config['metadata'], os.path.join("unfinished", os.path.basename(sample_file)))

def main(root_dir, args):
    if len(args) < 3:
        print("python run.py unfinished config.yaml cores")
        sys.exit(1)

    unfinished_dir = args[1]
    config_file_path = args[2]
    cores = 79 if len(args) <= 3 else int(args[3])
    parallel = 10 if len(args) <= 4 else int(args[4])

    # Load base config
    base_config = load_configfile(config_file_path)

    # Prepare directories
    finished_dir = "finished"
    duplication_dir = "duplication"
    metadata_dir = "metadata"
    for dir_path in [finished_dir, duplication_dir, metadata_dir]:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

    # Get metadata files
    metadata_files = glob.glob(os.path.join(unfinished_dir, "*.txt"))

    # Define task dictionary
    task_dict = {}
    for metadata_file in metadata_files:
        config = base_config.copy()
        config['metadata'] = metadata_file
        task_dict[os.path.basename(metadata_file)] = config

    # Get Snakefile path
    sf = get_snakefile(root_dir, args[5] if len(args) > 5 else "Snakefile")

    # Process tasks
    with ThreadPoolExecutor(max_workers=parallel) as executor:
        future_to_task = {executor.submit(process_sample_file, metadata_file, metadata_dir, sf, [config_file_path], cores, config): metadata_file for metadata_file, config in task_dict.items()}
        for future in as_completed(future_to_task):
            metadata_file = future_to_task[future]
            try:
                future.result()
            except Exception as exc:
                print(f'{metadata_file} generated an exception: {exc}')

if __name__ == '__main__':
    root_dir = os.path.dirname(os.path.abspath(__file__))
    args = sys.argv
    main(root_dir, args) 
