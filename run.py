import os
import sys
import glob
import shutil
import subprocess
import argparse
import pandas as pd

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

    print(f"Running command: {' '.join(cmd)}")
    
    try:
        # 不捕获输出，直接显示到终端
        result = subprocess.run(cmd, check=True, timeout=timeout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Snakemake command failed with return code: {e.returncode}", file=sys.stderr)
        return False
    except subprocess.TimeoutExpired:
        print(f"Snakemake command timed out after {timeout} seconds", file=sys.stderr)
        return False



def process_sample_file(metadata_file, metadata_dir, sf, config_file, cores, executor, executor_profile_path, timeout):
    try:
        print(f"Unlocking workflow for {metadata_file}...")
        run_snakemake(sf, [config_file], cores, unlock=True, timeout=timeout, executor=executor, executor_profile_path=executor_profile_path)
        
        print(f"Starting Snakemake execution for {metadata_file}...")
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

def check_sra_files(metadata_file, sra_dir):
    """检查metadata文件中所有SRR对应的SRA文件是否存在"""
    try:
        df = pd.read_csv(metadata_file, sep='\t')
        missing_files = []
        
        for _, row in df.iterrows():
            srr_list = row['SRR'].strip().split(',')
            
            for srr in srr_list:
                srr = srr.strip()
                sra_file = os.path.join(sra_dir, srr, f"{srr}.sra")
                if not os.path.exists(sra_file):
                    missing_files.append(sra_file)
        
        return missing_files
    except Exception as e:
        print(f"Error checking SRA files for {metadata_file}: {e}", file=sys.stderr)
        return ["Error reading metadata file"]

def main():
    parser = argparse.ArgumentParser(description="Process Snakemake metadata files.")
    parser.add_argument('unfinished_dir', help="Directory containing unfinished metadata files")
    parser.add_argument('config_file', help="Path to the configuration file")
    parser.add_argument('--cores', type=int, default=79, help="Number of cores to use (default: 79)")
    parser.add_argument('--snakefile', default='Snakefile', help="Path to custom Snakefile (default: 'Snakefile' in root directory)")
    parser.add_argument('--root_dir', default='.', help="Root directory for the Snakefile (default: current directory)")
    parser.add_argument('--executor', default=None, help="Specify a Snakemake executor (default: None)")
    parser.add_argument('--executor_profile_path', default=None, help="Path to the executor profile file (default: None)")
    parser.add_argument('--timeout', type=int, default=None, help="Timeout for each Snakemake run in seconds (default: None)")
    parser.add_argument('--sra_dir', default='sra', help="Directory containing SRA files (default: 'sra')")

    args = parser.parse_args()

    print(f"Running with the following parameters:")
    print(f"Unfinished directory: {args.unfinished_dir}")
    print(f"Config file: {args.config_file}")
    print(f"Cores: {args.cores}")
    print(f"Snakefile: {args.snakefile}")
    print(f"Executor: {args.executor}")
    print(f"Executor profile path: {args.executor_profile_path}")
    print(f"Timeout: {args.timeout}")
    print(f"SRA directory: {args.sra_dir}")

    with open(args.config_file, 'r') as config_file:
        base_config = yaml.safe_load(config_file)

    finished_dir = "finished"
    duplication_dir = "duplication"
    metadata_dir = "metadata"
    temp_config_dir = "temp_configs"
    for dir_path in [finished_dir, duplication_dir, metadata_dir, temp_config_dir]:
        os.makedirs(dir_path, exist_ok=True)

    metadata_files = glob.glob(os.path.join(args.unfinished_dir, "*.txt"))

    sf = get_snakefile(args.root_dir, args.snakefile)

    # 统计变量
    total_tasks = len(metadata_files)
    processed_tasks = 0
    skipped_tasks = 0
    skipped_files = []

    # 串行处理每个metadata文件
    for metadata_file in metadata_files:
        print(f"Checking SRA files for {metadata_file}...")
        
        # 检查SRA文件是否存在
        missing_files = check_sra_files(metadata_file, args.sra_dir)
        
        if missing_files:
            print(f"Skipping {metadata_file} - Missing SRA files:")
            for missing_file in missing_files:
                print(f"  - {missing_file}")
            
            skipped_tasks += 1
            skipped_files.append(os.path.basename(metadata_file))
            continue
        
        print(f"All SRA files found for {metadata_file}. Processing...")
        
        config = base_config.copy()
        config['metadata'] = metadata_file
        
        temp_config_file = os.path.join(temp_config_dir, f"{os.path.basename(metadata_file)}.yaml")
        with open(temp_config_file, 'w') as temp_file:
            yaml.dump(config, temp_file)
        
        result = process_sample_file(
            os.path.basename(metadata_file), 
            args.unfinished_dir, 
            sf, 
            temp_config_file, 
            args.cores, 
            args.executor, 
            args.executor_profile_path, 
            args.timeout
        )
        print(f'Task for {metadata_file} completed with result: {result}')
        processed_tasks += 1

    shutil.rmtree(temp_config_dir)

    # 输出统计信息
    print(f"\n=== Processing Summary ===")
    print(f"Total tasks: {total_tasks}")
    print(f"Processed tasks: {processed_tasks}")
    print(f"Skipped tasks: {skipped_tasks}")
    
    if skipped_tasks > 0:
        print(f"\nSkipped files due to missing SRA files:")
        for skipped_file in skipped_files:
            print(f"  - {skipped_file}")

if __name__ == '__main__':
    main()
