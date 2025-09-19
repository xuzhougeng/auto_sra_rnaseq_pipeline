import os
import sys
import glob
import json
import pandas as pd
import numpy as np
from os.path import basename, join
from collections import deque
#print(workflow.scheduler_type, file = sys.stderr)

# Auxiliary functions

def send_mail(subject, content, 
              sender, sender_passwd, 
              smtp_server = 'smtp.qq.com',
              msg_to="xuzhougeng@163.com"):
    import smtplib
    from email.mime.text import MIMEText

    msg = MIMEText(content, 'plain', 'utf-8')
    
    msg['Subject'] = subject
    msg['From'] = sender
    msg['To'] = msg_to
    # Send the email via our own SMTP server.
    try:
        client = smtplib.SMTP_SSL(smtp_server, smtplib.SMTP_SSL_PORT)
        print("connecting mail server successfully")
        client.login(sender, sender_passwd)
        print("loging mail server successfully")
        client.sendmail(sender, msg_to, msg.as_string())
        print("sending mail successfully")
    except smtplib.SMTPException as e:
        print("unable to send mail, check your SMTP server.", file = sys.stderr)
    finally:
        client.quit()


# get the scripts dir

configfile: "config.yaml"

root_dir = os.path.dirname(os.path.abspath(workflow.snakefile))
script_dir = os.path.join(root_dir, "scripts")
if not os.path.exists(script_dir):
    sys.exit("Unable to locate the Snakemake workflow file;  tried %s" %script_dir)


# filename of metadata is DB_ID
# column of metadata_file is 'GSE     GSM     gene    method  celline group   group_name      type    platform        SRR     paired'
# separated by '\t'
metadata_file = config['metadata']
metadata_df = pd.read_csv(metadata_file, sep = "\t")
DB_ID = os.path.basename(metadata_file).replace(".txt", "")


# 中间文件: SRA的文件名
sra_files = []      
for sra_id in metadata_df['SRR'].to_list():
    sra_files.append("sra/{sra}/{sra}.sra".format(sra=sra_id))

# 中间文件: 记录COUNT_FILE
count_files = []
for sample in metadata_df['GSM'].to_list():
    count_files.append("02_read_align/{sample}_ReadsPerGene.out.tab".format(sample=sample))

# 中间文件: 记录bigwig的文件名
bigwig_files = []
for sample in metadata_df['GSM'].to_list():
    bigwig_files.append("04_bigwig/{sample}.bw".format(sample=sample))


# 输出文件: counts_file和deseq_file
expr_matrix_file = f"03_merged_counts/{DB_ID}.tsv"
deseq_file  = f"05_DGE_analysis/{DB_ID}.Rds"


# functions for generate input of Snakemake rules
def get_input_data(wildcards):
    df = metadata_df
    sample = wildcards.sample
    paired = df.loc[df['GSM'] == sample, 'paired'].tolist()[0]  == "PAIRED"
    if paired:
        return [ os.path.join("01_clean_data", sample + '_R1.fq'), 
        os.path.join("01_clean_data",sample + '_R2.fq')]
    else:
        return [  os.path.join("01_clean_data", sample + '.fq') ]


# SRA数据路径配置
sra_data_path = config.get('sra_data_path', 'sra')

localrules: all, get_sra, merge_data, merge_R1_data, merge_R2_data, data_conversion_single, data_conversion_pair

rule all:
    input:
        deseq_file,
        expr_matrix_file,
        bigwig_files,
        count_files

# 直接使用已有的SRA数据，不进行下载
rule get_sra:
    output: 
        "sra/{sra}/{sra}.sra"
    run:
        # 构建SRA文件路径
        sra_file = os.path.join(sra_data_path, wildcards.sra, f"{wildcards.sra}.sra")
        
        # 检查文件是否存在
        if not os.path.exists(sra_file):
            raise FileNotFoundError(f"SRA file not found: {sra_file}")
        
        # 创建输出目录
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        
        # 如果源文件和目标文件不是同一个路径，创建符号链接
        if os.path.abspath(sra_file) != os.path.abspath(output[0]):
            if os.path.exists(output[0]):
                os.remove(output[0])
            os.symlink(os.path.abspath(sra_file), output[0])
            print(f"Created symlink: {output[0]} -> {sra_file}")
        else:
            print(f"Using existing SRA file: {output[0]}")



include: "rules/single_end_process.smk" # single end
include: "rules/paired_end_process.smk" # pair end

# alignment
rule align_and_count:
    input:
        fastq = get_input_data
    wildcard_constraints:
        sample="[A-Za-z0-9]+"
    params:
        index = config['index'],
        prefix = lambda wildcards : wildcards.sample,
        GTF = config['GTF']
    output: 
        bam = temp("02_read_align/{sample}_Aligned.sortedByCoord.out.bam"),
        counts = "02_read_align/{sample}_ReadsPerGene.out.tab"
    priority: 40
    resources:
        mem_mb = lambda wildcards, attempt: 10000 if attempt == 1 else 60000 * (attempt-1)
    threads: config['star_threads']
    benchmark:
        "benchmark/STAR/{sample}.tsv"
    conda:
        "envs/align.yaml"
    log: os.path.join( "log", '{sample}_Log.final.out'),
    shell:"""
        STAR \
    	--genomeDir {params.index} \
    	--runThreadN {threads} \
    	--readFilesIn {input.fastq} \
    	--outFileNamePrefix 02_read_align/{params.prefix}_ \
    	--outSAMtype BAM SortedByCoordinate \
    	--outBAMsortingThreadN 10 \
        --limitBAMsortRAM $(({resources.mem_mb} * 1000000)) \
        --quantMode GeneCounts --sjdbGTFfile {params.GTF} \
        --outTmpKeep None && \
        mv 02_read_align/{params.prefix}_Log.final.out {log}
    """
    # delete --readFilesCommand zcat \

rule build_bam_index:
    input: "02_read_align/{sample}_Aligned.sortedByCoord.out.bam"
    priority: 40
    output: temp("02_read_align/{sample}_Aligned.sortedByCoord.out.bam.bai")
    conda:
        "envs/align.yaml"
    threads: 4
    shell:"samtools index -@ {threads} {input}"

rule bamtobw:
    input: 
        bam = "02_read_align/{sample}_Aligned.sortedByCoord.out.bam",
        bai = "02_read_align/{sample}_Aligned.sortedByCoord.out.bam.bai"
    output: "04_bigwig/{sample}.bw"
    priority: 45
    params:
        bs = "50",
        gs = "2913022398",
        norm = "BPM"
    conda:
        "envs/align.yaml"
    threads: 10
    shell:"""
    bamCoverage -p {threads} \
        --binSize {params.bs} \
        --effectiveGenomeSize {params.gs} \
        --normalizeUsing {params.norm} \
        -b {input.bam} -o {output}
    """

rule combine_count:
    input: count_files
    output: "03_merged_counts/{DB_ID}.tsv"
    priority: 35
    run:
        #if not os.path.exists("03_merged_counts"):
            #os.mkir("03_merged_counts")
        df = pd.read_csv(input[0],header=None, sep="\t",index_col=0 )
        df = df.iloc[:,[0]]
        rename_dict = {1: os.path.basename(input[0]).replace('_ReadsPerGene.out.tab','')}
        df  =  df.rename(columns=rename_dict) 
        for file in input[1:]:
            df2 = pd.read_csv(file,header=None, sep="\t",index_col=0 )
            df2 = df2.iloc[:,[0]]
            rename_dict = {1: os.path.basename(file).replace('_ReadsPerGene.out.tab','')}
            df2  =  df2.rename(columns=rename_dict) 
            df = df2.merge(df,left_index=True, right_index=True)

        df.to_csv(output[0], sep='\t', encoding='utf-8')


rule DGE_analysis:
    input: 
        expr_matrix_file,
        metadata_file
    output: "05_DGE_analysis/{DB_ID}.Rds"
    priority: 35
    params:
        script_dir = script_dir,
    shell:"""
    Rscript {params.script_dir}/DESeq2_diff.R "{input[0]}" "{input[1]}" "{output}"
    """

# 清理STAR共享内存索引

onsuccess:
    print("Deleting the unnessary file")
    from shutil import rmtree
    if os.path.exists("02_read_align"):
        rmtree("02_read_align")
    
    contents = "snakemake run successful\nFollowing jobs fininished:\n "+ "\n".join(deseq_file)
    if config['mail']:
    #if len(config["sender"]) > 0 and len(config["sender_password"]) > 0:
        send_mail(subject = "snakemake run successful", content = contents, 
            sender = config["sender"],
            sender_passwd = config["sender_password"],
            msg_to=config['mail_to'])

onerror:
    contents = open(log, "r").read()
    if config['mail']:
    #if len(config["sender"]) > 0 and len(config["sender_password"]) > 0:
        send_mail(subject = "snakemake run failed", content = "snakemake run failed" + contents, 
            sender = config["sender"],
            sender_passwd = config["sender_password"],
            msg_to=config['mail_to'])
