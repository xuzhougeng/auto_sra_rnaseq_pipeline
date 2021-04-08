import os
import sys
import glob
import json
import pandas as pd
import numpy as np
from collections import deque


configfile: "config.yaml"

# get metadata file directory
metadata = config['metadata']

file_dict = {}
# file_dict.json:  key为GSE号的排序后组合, value为count+原始文件名
# 每次处理新的文件的时候, 通过对key的查找，可以用来确定该文件是否已经处理过

if os.path.exists("file_dict.json"):
    with open("file_dict.json", "r") as f:
        file_dict = json.load(f)


# 新的文件

counter = len(file_dict) # 已有的文件数

counts_file = []    # 记录输出的tsv文件名
metadata_dict = {}  # 字典, key为tsv文件名, value为对应的DataFrame

sample_files = glob.glob( os.path.join(metadata,  "*.txt") )
for file in sample_files:
    df = pd.read_csv(file, sep = "\t")
    dict_key = "_".join(sorted(df['GSM'].to_list()))
    # 判断当前文件是否存在重复或是否已经处理
    if dict_key not in file_dict:
        file_dict[dict_key] = [ counter, file ]
        
        GSE_ID = np.unique(df['GSE'])[0]
        gene   = np.unique(df['gene'])[0]
        file_name = "03_merged_counts/dataset{}_{}_{}.tsv".format(counter, GSE_ID, gene)
        
        counts_file.append(file_name)
        metadata_dict["{}_{}_{}".format(counter, GSE_ID, gene)] = df
        counter += 1


# 记录所有元信息, 用于后续查询
# 如果counts_file 没有内容，则直接退出
if len(counts_file) > 0:
    metadata_df = pd.concat(metadata_dict.values(), ignore_index=True)
    rep_len = list(map(len, metadata_dict.values()))
    metadata_df['key'] = np.repeat(list(metadata_dict.keys()), rep_len )
else:
    print("no job to do")
    sys.exit(0)


samples =  metadata_df['GSM'].to_list()

bigwig_files = expand('04_bigwig/{sample}.bw', sample = samples)

# get the input data of R1 and R2 or single

# get input GSM data
def get_input_data(wildcards):
    
    df = metadata_df
    sample = wildcards.sample
    paired = df.loc[df['GSM'] == sample, 'paired'].tolist()[0]  == "PAIRED"
    if paired:
        return [ os.path.join("01_clean_data", sample + '_R1.fq.gz'), 
        os.path.join("01_clean_data",sample + '_R2.fq.gz')]
    else:
        return [  os.path.join("01_clean_data", sample + '.fq.gz')]

# get GSM ID
def get_counts_file(wildcards):
    number  = wildcards.number
    GSE_ID = wildcards.GSE_ID
    gene   = wildcards.gene
    df = metadata_dict["{}_{}_{}".format(number, GSE_ID, gene)]
    samples =  df['GSM'].to_list()
    count_files = ["02_read_align/{sample}_ReadsPerGene.out.tab".format(sample=sample) for sample in samples]
    return  count_files

    
localrules: all, data_downloader
rule all:
    input:
        counts_file,
        bigwig_files


# download data from NCBI
rule data_downloader:
    params: 
        sra_id = lambda wildcards: wildcards.sra
    output: temp("sra/{sra}/{sra}.sra")
    threads: config['download_threads']
    conda:
        "envs/download.yaml"
    shell:"""
    prefetch -O sra {params.sra_id} && \
        [[ ! -f {output} ]] && \
        mv sra/{params.sra_id}.sra {output} || \
        echo "{params.sra_id} finished download"
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
        counts = "02_read_align/{sample}_ReadsPerGene.out.tab"
    priority: 10
    threads: config['star_threads']
    conda:
        "envs/align.yaml"
    shell:"""
        STAR \
    	--genomeDir {params.index} \
    	--runThreadN {threads} \
    	--readFilesIn {input} \
    	--readFilesCommand zcat \
    	--outFileNamePrefix 02_read_align/{params.prefix}_ \
    	--outSAMtype BAM SortedByCoordinate \
    	--outBAMsortingThreadN 10 \
        --quantMode GeneCounts --sjdbGTFfile {params.GTF} &&
        rm -rf 02_read_align/{params.prefix}__STARgenome
    """
rule build_bam_index:
    input: "02_read_align/{sample}_Aligned.sortedByCoord.out.bam"
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
    output: "03_merged_counts/dataset{number}_{GSE_ID}_{gene}.tsv"
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


onsuccess:
    # when 
    with open("file_dict.json", "w") as f:
        json.dump(file_dict, f)

onerror:
    print("An error occurred")


