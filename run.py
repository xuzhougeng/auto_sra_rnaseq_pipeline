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
        # 捕获stderr来检查错误类型，但stdout仍然实时显示
        result = subprocess.run(cmd, check=True, timeout=timeout, stderr=subprocess.PIPE, text=True)
        return True, ""
    except subprocess.CalledProcessError as e:
        print(f"Snakemake command failed with return code: {e.returncode}", file=sys.stderr)
        stderr_output = e.stderr if e.stderr else ""
        return False, stderr_output
    except subprocess.TimeoutExpired:
        print(f"Snakemake command timed out after {timeout} seconds", file=sys.stderr)
        return False, "Timeout"

def process_sample_file(metadata_file, metadata_dir, 
                        finished_dir, failed_dir, 
                        sf, config_file, cores, executor, executor_profile_path, timeout):
    try:
        print(f"Unlocking workflow for {metadata_file}...")
        unlock_status, unlock_error = run_snakemake(sf, [config_file], cores, unlock=True, timeout=timeout, executor=executor, executor_profile_path=executor_profile_path)
        
        print(f"Starting Snakemake execution for {metadata_file}...")
        status, error_output = run_snakemake(sf, [config_file], cores, timeout=timeout, executor=executor, executor_profile_path=executor_profile_path)

        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # 检查是否有MissingOutputException错误
        has_missing_output_error = "MissingOutputException" in error_output
        
        if status:
            contents = f"Snakemake run successfully for {metadata_file}"
            shutil.move(os.path.join(metadata_dir, metadata_file), os.path.join(finished_dir, metadata_file))
        else:
            if has_missing_output_error:
                contents = f"Snakemake run failed for {metadata_file} - MissingOutputException (filesystem latency issue)"
            else:
                contents = f"Snakemake run failed or timed out for {metadata_file}"
            shutil.move(os.path.join(metadata_dir, metadata_file), os.path.join(failed_dir, metadata_file))

        if config.get('bark'):
            bark_notification(config['bark_api'], contents)
        if config.get('feishu'):
            feishu_notification(config['feishu_api'], contents)
        
        return contents, has_missing_output_error
    except Exception as e:
        error_message = f"Error processing {metadata_file}: {str(e)}"
        print(error_message, file=sys.stderr)
        return error_message, False

def validate_metadata_file(metadata_file):
    """验证metadata文件格式是否符合要求"""
    try:
        df = pd.read_csv(metadata_file, sep='\t')
        errors = []
        
        # 检查必需的列是否存在
        required_columns = ['SRR', 'paired', 'GSM', 'GSE']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            errors.append(f"Missing required columns: {', '.join(missing_columns)}")
            return errors
        
        # 检查每一行的SRR和paired字段
        for idx, row in df.iterrows():
            row_num = idx + 2  # +2 because pandas is 0-indexed and we have header
            
            # 检查SRR字段
            srr_value = row['SRR']
            if pd.isna(srr_value) or str(srr_value).strip().upper() in ['NA', 'NAN', '']:
                errors.append(f"Row {row_num}: SRR field is empty or NA")
            
            # 检查paired字段
            paired_value = row['paired']
            if pd.isna(paired_value) or str(paired_value).strip().upper() in ['NA', 'NAN', '']:
                errors.append(f"Row {row_num}: paired field is empty or NA")
            elif str(paired_value).strip().upper() not in ['PAIRED', 'SINGLE']:
                errors.append(f"Row {row_num}: paired field must be 'PAIRED' or 'SINGLE', got '{paired_value}'")
        
        return errors
    except Exception as e:
        return [f"Error reading metadata file: {str(e)}"]

def check_sra_files(metadata_file, sra_dir):
    """检查metadata文件中所有SRR对应的SRA文件是否存在"""
    try:
        df = pd.read_csv(metadata_file, sep='\t')
        missing_files = []
        
        for _, row in df.iterrows():
            # 检查SRR值是否为NaN或NA
            srr_value = row['SRR']
            if pd.isna(srr_value) or str(srr_value).strip().upper() in ['NA', 'NAN', '']:
                continue  # 跳过无效的SRR值（这些应该已经在validate_metadata_file中被捕获）
                
            srr_list = str(srr_value).strip().split(',')
            
            for srr in srr_list:
                srr = srr.strip()
                if srr and srr.upper() not in ['NA', 'NAN']:  # 确保SRR值有效
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
    failed_dir = "failed"
    metadata_dir = "metadata"
    temp_config_dir = "temp_configs"
    for dir_path in [finished_dir, failed_dir, metadata_dir, temp_config_dir]:
        os.makedirs(dir_path, exist_ok=True)

    metadata_files = glob.glob(os.path.join(args.unfinished_dir, "*.txt"))

    sf = get_snakefile(args.root_dir, args.snakefile)

    # 统计变量
    total_tasks = len(metadata_files)
    processed_tasks = 0
    skipped_tasks = 0
    failed_tasks = 0
    missing_output_tasks = 0
    skipped_files = []
    failed_files = []
    missing_output_files = []

    # 串行处理每个metadata文件
    for metadata_file in metadata_files:
        print(f"Validating metadata file {metadata_file}...")
        
        # 首先验证metadata文件格式
        validation_errors = validate_metadata_file(metadata_file)
        
        if validation_errors:
            print(f"Skipping {metadata_file} - Metadata validation failed:")
            for error in validation_errors:
                print(f"  - {error}")
            
            skipped_tasks += 1
            skipped_files.append(os.path.basename(metadata_file))
            continue
        
        print(f"Metadata file {metadata_file} is valid. Checking SRA files...")
        
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
        
        result, has_missing_output = process_sample_file(
            os.path.basename(metadata_file), 
            args.unfinished_dir, 
            finished_dir,
            failed_dir,
            sf, 
            temp_config_file, 
            args.cores, 
            args.executor, 
            args.executor_profile_path, 
            args.timeout
        )
        
        print(f'Task for {metadata_file} completed with result: {result}')
        
        if "successfully" in result:
            processed_tasks += 1
        else:
            failed_tasks += 1
            failed_files.append(os.path.basename(metadata_file))
            
            if has_missing_output:
                missing_output_tasks += 1
                missing_output_files.append(os.path.basename(metadata_file))

    shutil.rmtree(temp_config_dir)

    # 输出统计信息
    print(f"\n=== Processing Summary ===")
    print(f"Total tasks: {total_tasks}")
    print(f"Successfully processed: {processed_tasks}")
    print(f"Failed tasks: {failed_tasks}")
    print(f"Skipped tasks: {skipped_tasks}")
    
    if missing_output_tasks > 0:
        print(f"\nTasks with MissingOutputException errors: {missing_output_tasks}")
        for missing_output_file in missing_output_files:
            print(f"  - {missing_output_file}")
    
    if failed_tasks > 0:
        print(f"\nAll failed files:")
        for failed_file in failed_files:
            print(f"  - {failed_file}")
    
    if skipped_tasks > 0:
        print(f"\nSkipped files due to missing SRA files:")
        for skipped_file in skipped_files:
            print(f"  - {skipped_file}")

if __name__ == '__main__':
    main()
