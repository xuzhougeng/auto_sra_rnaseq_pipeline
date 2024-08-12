import os
import sys
import glob
import shutil
import subprocess
import argparse

import yaml

from scripts.utilize import bark_notification, feishu_notification


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
def run_snakemake(snakefile, configfiles, cores, unlock=False, timeout=None, executor=None, executor_profile_path=None):
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
    
    if executor:
        cmd.extend(["--executor", executor])
        if executor_profile_path:
            cmd.extend(["--profile", executor_profile_path])

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=timeout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Snakemake command failed: {e}", file=sys.stderr)
        print(f"Stdout: {e.stdout}", file=sys.stderr)
        print(f"Stderr: {e.stderr}", file=sys.stderr)
        return False
    except subprocess.TimeoutExpired:
        print(f"Snakemake command timed out after {timeout} seconds", file=sys.stderr)
        return False

from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

def process_sample_file(metadata_file, metadata_dir, sf, config_file, cores, executor, executor_profile_path, timeout):
    try:
        run_snakemake(sf, [config_file], cores, unlock=True, timeout=timeout, executor=executor, executor_profile_path=executor_profile_path)
        status = run_snakemake(sf, [config_file], cores, timeout=timeout, executor=executor, executor_profile_path=executor_profile_path)

        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        if status:
            contents = f"Snakemake run successfully for {metadata_file}"
            shutil.move(os.path.join(metadata_dir, metadata_file), os.path.join("finished", metadata_file))
        else:
            contents = f"Snakemake run failed or timed out for {metadata_file}"

        if config.get('bark'):
            bark_notification(config['bark_api'], contents)
        if config.get('feishu'):
            feishu_notification(config['feishu_api'], contents)
        
        return contents
    except Exception as e:
        error_message = f"Error processing {metadata_file}: {str(e)}"
        print(error_message, file=sys.stderr)
        return error_message

def main():
    parser = argparse.ArgumentParser(description="Process Snakemake metadata files.")
    parser.add_argument('unfinished_dir', help="Directory containing unfinished metadata files")
    parser.add_argument('config_file', help="Path to the configuration file")
    parser.add_argument('--cores', type=int, default=79, help="Number of cores to use (default: 79)")
    parser.add_argument('--parallel', type=int, default=10, help="Number of parallel tasks (default: 10)")
    parser.add_argument('--snakefile', default='Snakefile', help="Path to custom Snakefile (default: 'Snakefile' in root directory)")
    parser.add_argument('--root_dir', default='.', help="Root directory for the Snakefile (default: current directory)")
    parser.add_argument('--executor', default=None, help="Specify a Snakemake executor (default: None)")
    parser.add_argument('--executor_profile_path', default=None, help="Path to the executor profile file (default: None)")

    args = parser.parse_args()

    print(f"Running with the following parameters:")
    print(f"Unfinished directory: {args.unfinished_dir}")
    print(f"Config file: {args.config_file}")
    print(f"Cores: {args.cores}")
    print(f"Parallel tasks: {args.parallel}")
    print(f"Snakefile: {args.snakefile}")
    print(f"Executor: {args.executor}")
    print(f"Executor profile path: {args.executor_profile_path}")

    with open(args.config_file, 'r') as config_file:
        base_config = yaml.safe_load(config_file)

    finished_dir = "finished"
    duplication_dir = "duplication"
    metadata_dir = "metadata"
    temp_config_dir = "temp_configs"
    for dir_path in [finished_dir, duplication_dir, metadata_dir, temp_config_dir]:
        os.makedirs(dir_path, exist_ok=True)

    metadata_files = glob.glob(os.path.join(args.unfinished_dir, "*.txt"))

    task_dict = {}
    for metadata_file in metadata_files:
        config = base_config.copy()
        config['metadata'] = metadata_file
        
        temp_config_file = os.path.join(temp_config_dir, f"{os.path.basename(metadata_file)}.yaml")
        with open(temp_config_file, 'w') as temp_file:
            yaml.dump(config, temp_file)
        
        task_dict[os.path.basename(metadata_file)] = temp_config_file

    sf = get_snakefile(args.root_dir, args.snakefile)

    with ProcessPoolExecutor(max_workers=args.parallel) as executor:
        tasks = [
            (metadata_file, args.unfinished_dir, sf, temp_config_file, args.cores,  args.executor, args.executor_profile_path) 
            for metadata_file, temp_config_file in task_dict.items()
        ]
        future_to_task = {executor.submit(process_sample_file, *task): task[0] for task in tasks}
        for future in as_completed(future_to_task):
            metadata_file = future_to_task[future]
            try:
                result = future.result()
                print(f'Task for {metadata_file} completed with result: {result}')
            except Exception as exc:
                print(f'Task for {metadata_file} generated an exception: {exc}')

    shutil.rmtree(temp_config_dir)

if __name__ == '__main__':
    main()
