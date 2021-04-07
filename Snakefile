import os
import pandas as pd
import numpy as np


"""
SRA数据自动处理流程, 要求，必须有SRR和paired两列
如果SRR存在多个, 在解压缩完成之后会进行合并

requirement:
- sratoolkit
- fastp
- star
"""
configfile: "config.yaml"

sample_file = config['sample']

df = pd.read_csv(sample_file, sep = "\t")

samples = df['GSM'].to_list()

result_files = expand('02_read_align/{sample}_Aligned.sortedByCoord.out.bam', sample = samples)
all_counts   = expand('02_read_align/{sample}_ReadsPerGene.out.tab', sample = samples)

# get the input data of R1 and R2 or single
def get_merged_input_data(wildcards):
    sample = wildcards.sample

    SRR_ID = df.loc[df['GSM'] == sample, "SRR"].tolist()[0]
    SRR_ID = SRR_ID.split(",")
    print(SRR_ID)

    return [ "sra/{x}.fastq".format(x=x) for x in SRR_ID ]

def get_merged_input_data_R1(wildcards):
    sample = wildcards.sample

    SRR_ID = df.loc[df['GSM'] == sample, "SRR"].tolist()[0]
    SRR_ID = SRR_ID.split(",")

    return [ "sra/{x}_1.fastq".format(x=x) for x in SRR_ID ]

def get_merged_input_data_R2(wildcards):
    sample = wildcards.sample

    SRR_ID = df.loc[df['GSM'] == sample, "SRR"].tolist()[0]
    SRR_ID = SRR_ID.split(",")

    return [ "sra/{x}_2.fastq".format(x=x) for x in SRR_ID ]

# get input GSM data
def get_input_data(wildcards):
    sample = wildcards.sample
    paired = df.loc[df['GSM'] == sample, 'paired'].tolist()[0]  == "PAIRED"
    if paired:
        return [ os.path.join("01_clean_data", sample + '_R1.fq.gz'), 
        os.path.join("01_clean_data",sample + '_R2.fq.gz')]
    else:
        return [  os.path.join("01_clean_data", sample + '.fq.gz')]

rule all:
    input:
        result_files


# download data from NCBI
rule data_downloader:
    params: 
        sra_id = lambda wildcards: wildcards.sra
    output: temp("sra/{sra}.sra")
    conda:
        "envs/download.yaml"
    shell:"""
    prefetch -O {input} -O sra 
    """

# covert the fastq
# for single
rule data_conversion_single:
    input: "sra/{sra}.sra"
    wildcard_constraints:
        sra="[A-Za-z0-9]+"
    output: temp("sra/{sra}.fastq")
    conda:
        "envs/download.yaml"
    shell:"fastq-dump {input} -O sra" 

rule merge_data:
    input: get_merged_input_data
    output: "00_raw_data/{sample}.fq.gz"
    shell: "cat {input} | pigz > {output}"

rule data_clean_single:
    input: "00_raw_data/{sample}.fq.gz"
    wildcard_constraints:
        sample="[A-Za-z0-9]+"
    params:
        json = lambda wildcards : os.path.join( "qc",  wildcards.sample + '.json'),
        html = lambda wildcards : os.path.join( "qc",  wildcards.sample + '.html')
    output: "01_clean_data/{sample}.fq.gz"
    conda:
        "envs/preprocess.yaml"
    shell:"""
    fastp -w {threads} -i {input}  -o {output}  \
		-j {params.json} -h {params.html}
    """

# for pair-end
rule data_conversion_pair:
    input: "sra/{sra}.sra"
    wildcard_constraints:
        sra="[A-Za-z0-9]+"
    output:
        temp("sra/{sra}_1.fastq"),
        temp("sra/{sra}_2.fastq") 
    conda:
        "envs/download.yaml"
    shell:"fastq-dump --split-files {input} -O sra" 

rule merge_R1_data:
    input: get_merged_input_data_R1
    output: "00_raw_data/{sample}_R1.fq.gz"
    shell: "cat {input} | pigz > {output}"

rule merge_R2_data:
    input: get_merged_input_data_R2
    output: "00_raw_data/{sample}_R2.fq.gz"
    shell: "cat {input} | pigz > {output}"

rule data_clean_pair:
    input:
        r1 = "00_raw_data/{sample}_R1.fq.gz",
        r2 = "00_raw_data/{sample}_R2.fq.gz"
    wildcard_constraints:
        sample="[A-Za-z0-9]+"
    params:
        json = lambda wildcards : os.path.join( "qc",  wildcards.sample + '.json'),
        html = lambda wildcards : os.path.join( "qc",  wildcards.sample + '.html'),
    output:
        r1 = "01_clean_data/{sample}_R1.fq.gz",
        r2 = "01_clean_data/{sample}_R2.fq.gz" 
    conda:
        "envs/preprocess.yaml"
    threads: 8
    shell:"""
	fastp -w {threads} -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
		-j {params.json} -h {params.html}
    """

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
        bam = "02_read_align/{sample}_Aligned.sortedByCoord.out.bam"
        count = "02_read_align/{sample}_ReadsPerGene.out.tab"
    threads: 40
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
        --quantMode GeneCounts --sjdbGTFfile {params.GTF}
    """

# rule merge count
rule combine_count:
    input: all_counts
    output: "expreSet.tsv"
