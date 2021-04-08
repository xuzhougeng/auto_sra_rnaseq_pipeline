def get_merged_input_data_R1(wildcards):
    df = metadata_df
    sample = wildcards.sample

    SRR_ID = df.loc[df['GSM'] == sample, "SRR"].tolist()[0]
    SRR_ID = SRR_ID.split(",")

    return [ "sra/{x}_1.fastq".format(x=x) for x in SRR_ID ]

def get_merged_input_data_R2(wildcards):
    df = metadata_df
    sample = wildcards.sample

    SRR_ID = df.loc[df['GSM'] == sample, "SRR"].tolist()[0]
    SRR_ID = SRR_ID.split(",")

    return [ "sra/{x}_2.fastq".format(x=x) for x in SRR_ID ]



rule data_conversion_pair:
    input: "sra/{sra}/{sra}.sra"
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
    threads: config['pigz_threads']
    output: "00_raw_data/{sample}_R1.fq.gz"
    shell: "cat {input} | pigz -p {threads} > {output}"

rule merge_R2_data:
    input: get_merged_input_data_R2
    threads: config['pigz_threads']
    output: "00_raw_data/{sample}_R2.fq.gz"
    shell: "cat {input} | pigz -p {threads} > {output}"

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
    threads: config['fastp_threads']
    shell:"""
	fastp -w {threads} -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
		-j {params.json} -h {params.html}
    """

