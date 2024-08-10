import os
import sys
import glob
import shutil
import subprocess
import time

import yaml

from scripts.utilize import bark_notification, feishu_notification
from concurrent.futures import ThreadPoolExecutor, TimeoutError


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
    cmd = [
        "snakemake",
        "-s", snakefile,
        "--configfile", configfiles[0],
        "--cores", str(cores),
        "--resources", "rx=80", "limit_dump=2", "limit_merge=2",
        "--rerun-incomplete",
        "--latency-wait", "5",
        "--restart-times", "3"
    ]
    
    if unlock:
        cmd.append("--unlock")
    
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    def run_process():
        stdout, stderr = process.communicate()
        return process.returncode, stdout, stderr
    
    with ThreadPoolExecutor() as executor:
        future = executor.submit(run_process)
        try:
            returncode, stdout, stderr = future.result(timeout=3600)  # 1 hour timeout
            if returncode == 0:
                return True
            else:
                print(f"Snakemake command failed with return code {returncode}", file=sys.stderr)
                print(f"Stdout: {stdout}", file=sys.stderr)
                print(f"Stderr: {stderr}", file=sys.stderr)
                return False
        except TimeoutError:
            process.kill()
            print("Snakemake execution timed out after 1 hour", file=sys.stderr)
            return False

from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

def process_sample_file(args):
    metadata_file, metadata_dir, sf, config_file, cores = args
    try:
        # Run Snakemake to unlock any potential issues
        run_snakemake(sf, [config_file], cores, unlock=True)
        
        # Execute Snakemake with the provided configuration
        start_time = time.time()
        status = run_snakemake(sf, [config_file], cores)
        end_time = time.time()
        
        # Load the config to access notification settings
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # Check the status and handle notifications accordingly
        if status:
            contents = f"Snakemake run successfully for {metadata_file}"
            # Move the metadata file to finished directory
            shutil.move(os.path.join(metadata_dir, metadata_file), os.path.join("finished", metadata_file))
        else:
            if end_time - start_time >= 3600:
                contents = f"Snakemake run timed out after 1 hour for {metadata_file}"
            else:
                contents = f"Snakemake run failed for {metadata_file}"
        
        # Send notifications
        if config.get('bark'):
            bark_notification(config['bark_api'], contents)
        if config.get('feishu'):
            feishu_notification(config['feishu_api'], contents)
        
        return contents
    except Exception as e:
        error_message = f"Error processing {metadata_file}: {str(e)}"
        print(error_message, file=sys.stderr)
        return error_message

def main(root_dir, args):
    usage = """
    Usage: python run.py <unfinished_dir> <config.yaml> [cores] [parallel] [Snakefile]
    
    Arguments:
    unfinished_dir   : Directory containing unfinished metadata files
    config.yaml      : Path to the configuration file
    cores            : Number of cores to use (default: 79)
    parallel         : Number of parallel tasks (default: 10)
    Snakefile        : Path to custom Snakefile (default: 'Snakefile' in root directory)
    
    Example: python run.py unfinished config.yaml 40 5
    """

    if len(args) < 2 or args[1] in ['-h', '--help']:
        print(usage)
        sys.exit(0)

    if len(args) < 3:
        print("Error: Not enough arguments.")
        print(usage)
        sys.exit(1)

    unfinished_dir = args[1]
    config_file_path = args[2]
    cores = 79 if len(args) <= 3 else int(args[3])
    parallel = 10 if len(args) <= 4 else int(args[4])

    print(f"Running with the following parameters:")
    print(f"Unfinished directory: {unfinished_dir}")
    print(f"Config file: {config_file_path}")
    print(f"Cores: {cores}")
    print(f"Parallel tasks: {parallel}")
    print(f"Snakefile: {args[5] if len(args) > 5 else 'Default'}")

    # Load base config
    with open(config_file_path, 'r') as config_file:
        base_config = yaml.safe_load(config_file)

    # Prepare directories
    finished_dir = "finished"
    duplication_dir = "duplication"
    metadata_dir = "metadata"
    temp_config_dir = "temp_configs"
    for dir_path in [finished_dir, duplication_dir, metadata_dir, temp_config_dir]:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

    # Get metadata files
    metadata_files = glob.glob(os.path.join(unfinished_dir, "*.txt"))

    # Define task dictionary
    task_dict = {}
    for metadata_file in metadata_files:
        config = base_config.copy()
        config['metadata'] = metadata_file
        
        # Create temporary YAML file
        temp_config_file = os.path.join(temp_config_dir, f"{os.path.basename(metadata_file)}.yaml")
        with open(temp_config_file, 'w') as temp_file:
            yaml.dump(config, temp_file)
        
        task_dict[os.path.basename(metadata_file)] = temp_config_file

    # Get Snakefile path
    sf = get_snakefile(root_dir, args[5] if len(args) > 5 else "Snakefile")

    # Process tasks
    with ProcessPoolExecutor(max_workers=parallel) as executor:
        tasks = [(metadata_file, unfinished_dir, sf, temp_config_file, cores) for metadata_file, temp_config_file in task_dict.items()]
        future_to_task = {executor.submit(process_sample_file, task): task[0] for task in tasks}
        for future in as_completed(future_to_task):
            metadata_file = future_to_task[future]
            try:
                result = future.result()
                print(f'Task for {metadata_file} completed with result: {result}')
            except Exception as exc:
                print(f'Task for {metadata_file} generated an exception: {exc}')

    # Clean up temporary config files
    shutil.rmtree(temp_config_dir)

if __name__ == '__main__':
    root_dir = os.path.dirname(os.path.abspath(__file__))
    args = sys.argv
    main(root_dir, args) 
