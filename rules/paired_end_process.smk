def get_merged_input_data_R1(wildcards):
    df = metadata_df
    sample = wildcards.sample
    separator = config.get('srr_separator', ',')

    SRR_ID = df.loc[df['GSM'] == sample, "SRR"].tolist()[0]
    SRR_ID = SRR_ID.split(separator)

    return [ "sra/{x}_1.fastq".format(x=x) for x in SRR_ID ]

def get_merged_input_data_R2(wildcards):
    df = metadata_df
    sample = wildcards.sample
    separator = config.get('srr_separator', ',')

    SRR_ID = df.loc[df['GSM'] == sample, "SRR"].tolist()[0]
    SRR_ID = SRR_ID.split(separator)

    return [ "sra/{x}_2.fastq".format(x=x) for x in SRR_ID ]



rule data_conversion_pair:
    input: "sra/{sra}/{sra}.sra"
    priority: 10
    wildcard_constraints:
        sra="[A-Za-z0-9]+"
    output:
        temp("sra/{sra}_1.fastq"),
        temp("sra/{sra}_2.fastq") 
    conda:
        "envs/download.yaml"
    resources:
        limit_dump = 1
    shell:"fasterq-dump sra/{wildcards.sra} -O sra" 
    #shell:"fastq-dump --split-files sra/{wildcards.sra} -O sra" 

rule merge_R1_data:
    input: get_merged_input_data_R1
    priority: 20
    #threads: config['pigz_threads']
    resources:
        limit_merge = 1
    output: temp("00_raw_data/{sample}_R1.fq")
    #shell: "cat {input} | pigz -p {threads} > {output}"
    shell: "cat {input}  > {output}"

rule merge_R2_data:
    input: get_merged_input_data_R2
    priority: 20
    #threads: config['pigz_threads']
    resources:
        limit_merge = 1
    output: temp("00_raw_data/{sample}_R2.fq")
    #shell: "cat {input} | pigz -p {threads} > {output}"
    shell: "cat {input}  > {output}"

rule data_clean_pair:
    input:
        r1 = "00_raw_data/{sample}_R1.fq",
        r2 = "00_raw_data/{sample}_R2.fq"
    priority: 30
    wildcard_constraints:
        sample="[A-Za-z0-9]+"
    output:
        r1 = temp("01_clean_data/{sample}_R1.fq"),
        r2 = temp("01_clean_data/{sample}_R2.fq") 
    conda:
        "envs/preprocess.yaml"
    threads: config['fastp_threads']
    log:
        json = os.path.join( "log", '{sample}.json'),
        html = os.path.join( "log", '{sample}.html'),
        logs = os.path.join( "log", "{sample}_fastp.log")
    shell:"""
	fastp -w {threads} -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
		-j {log.json} -h {log.html} &> {log.logs}
    """

