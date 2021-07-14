import os
import sys
import glob
import json
import pandas as pd
import numpy as np
from collections import deque
from .scripts.utilize import table_from_sql

#print(workflow.scheduler_type, file = sys.stderr)

# Auxiliary function

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

# define hash function
def myhash(string, size=8):
    import math
    string = string.replace("GSM", "")
    string = string.replace("_","0")[0:100]
    hash_value = int(string)<< 4
    hash_value = int( abs(hash_value) / math.pow(10, size) )
    hash_value = str(abs(hash_value))[0:8]

    return int(hash_value)

# get the input data of R1 and R2 or single
# get input GSM data
def get_input_data(wildcards):
    
    df = sample_df
    sample = wildcards.sample
    paired = df.loc[df['GSM'] == sample, 'paired'].tolist()[0]  == "PAIRED"
    if paired:
        return [ os.path.join("01_clean_data", sample + '_R1.fq.gz'), 
        os.path.join("01_clean_data",sample + '_R2.fq.gz')]
    else:
        return [  os.path.join("01_clean_data", sample + '.fq.gz') ]

# get GSM ID
def get_counts_file(wildcards):
    number  = wildcards.number
    GSE_ID = wildcards.GSE_ID
    gene   = wildcards.gene
    df = metadata_dict["{}_{}_{}".format(GSE_ID, gene, number)]
    samples =  df['GSM'].to_list()
    count_files = ["02_read_align/{sample}_ReadsPerGene.out.tab".format(sample=sample) for sample in samples]
    return  count_files

# get the count and meta
def get_counts_and_meta(wildcards):
    number  = wildcards.number
    GSE_ID = wildcards.GSE_ID
    gene   = wildcards.gene
    df = metadata_dict["{}_{}_{}".format(GSE_ID, gene, number)]
    dict_key = "_".join(sorted(df['GSM'].to_list() ) )
    meta_file =  file_dict[dict_key]
    counts_file  =  "03_merged_counts/{}_{}_{}.tsv".format(GSE_ID, gene, number)
    return [counts_file, meta_file]

# get the dictionary
def get_dict(wildcards):
    number  = wildcards.number
    GSE_ID = wildcards.GSE_ID
    gene   = wildcards.gene
    df = metadata_dict["{}_{}_{}".format(GSE_ID, gene, number)]
    dict_key = "_".join(sorted(df['GSM'].to_list() ) )
    return dict_key



root_dir = os.path.dirname(os.path.abspath(workflow.snakefile))
script_dir = os.path.join(root_dir, "scripts")
if not os.path.exists(script_dir):
    sys.exit("Unable to locate the Snakemake workflow file;  tried %s" %script_dir)

# get the scripts dir

configfile: "config.yaml"

# get metadata file directory
metadata = config['metadata']

# retrieve the data from database
db = "meta_info.sqlite3"

sample_files = glob.glob( os.path.join(metadata,  "*.txt") )
sample_file_names = [ basename(f) for f in sample_files ]

metadata_df = table_from_sql(table_name = "meta", db=db)
metadata_df = metadata_df.loc[metadata_df['meta_file'].isin(sample_file_names), ]

accessions = metadata_df['accesion']
counts_file = metadata_df['count_file']
deseq_file  = metadata_df['deseq_file']

# 记录所有的元信息, 用于后续查询
sample_df = table_from_sql(table_name = "sample", db=db)
sample_df = sample_df.loc[sample_df['accesion'].isin(accessions), ]


# 如果counts_file 没有内容，则直接退出
if len(counts_file) < 0:
    print("no job to do")
    os._exit(0)


samples =  sample_df['GSM'].to_list()

bigwig_files = expand('04_bigwig/{sample}.bw', sample = samples)


localrules: all, data_downloader
rule all:
    input:
        deseq_file,
        counts_file,
        bigwig_files
        
# download data from NCBI
rule data_downloader:
    priority: 5
    params: 
        sra_id = lambda wildcards: wildcards.sra,
        maxsize = "100G"
    output: temp("sra/{sra}/{sra}.sra")
    threads: config['download_threads']
    resources:
        rx = 40
    conda:
        "envs/download.yaml"
    benchmark:
        "benchmark/download/{sra}.tsv"
    shell:"""
    if [ -f {output} ] ;then \
        echo "{params.sra_id} has beed downloaded" ;\   
    elif [ -f sra/{params.sra_id}.sra.lock ] ; then \
        rm -f sra/{params.sra_id}.sra.lock sra/{params.sra_id}.sra.prf sra/{params.sra_id}.sra.tmp && \
        prefetch --max-size {params.maxsize} -O sra {params.sra_id} && mv sra/{params.sra_id}.sra {output} ;\
    elif [ ! -f {output} ] ;then \
        prefetch --max-size {params.maxsize} -O sra {params.sra_id} && mv sra/{params.sra_id}.sra {output} ;\
    else  \
        exit 1 ;\
    fi
    """

include: "rules/single_end_process.smk" # single end
include: "rules/paired_end_process.smk" # pair end

# alignment
rule align_and_count:
    input:
        get_input_data
    wildcard_constraints:
        sample="[A-Za-z0-9]+"
    params:
        index = config['index'],
        prefix = lambda wildcards : wildcards.sample,
        GTF = config['GTF']
    output: 
        bam = temp("02_read_align/{sample}_Aligned.sortedByCoord.out.bam"),
        counts = temp("02_read_align/{sample}_ReadsPerGene.out.tab")
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
    	--readFilesIn {input} \
    	--readFilesCommand zcat \
    	--outFileNamePrefix 02_read_align/{params.prefix}_ \
    	--outSAMtype BAM SortedByCoordinate \
    	--outBAMsortingThreadN 10 \
        --limitBAMsortRAM $(({resources.mem_mb} * 1000000)) \
        --quantMode GeneCounts --sjdbGTFfile {params.GTF} \
        --outTmpKeep None && \
        mv 02_read_align/{params.prefix}_Log.final.out {log}
    """
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
    input: get_counts_file
    output: "03_merged_counts/{GSE_ID}_{gene}_{number}.tsv"
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
        get_counts_and_meta
    output: "05_DGE_analysis/{GSE_ID}_{gene}_{number}.Rds"
    priority: 35
    params:
        script_dir = script_dir,
        #dict_key = lambda wildcards : get_dict(wildcards)
    shell:"""
    Rscript {params.script_dir}/DESeq2_diff.R "{input[0]}" "{input[1]}" "{output}"
    """
    # shell:"""
    # Rscript {params.script_dir}/DESeq2_diff.R "{input[0]}" "{input[1]}" "{output}" &&
    # python {params.script_dir}/update_json.py {params.dict_key} "{input[1]}"
    # """

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