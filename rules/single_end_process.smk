def get_merged_input_data(wildcards):
    df = metadata_df
    sample = wildcards.sample

    SRR_ID = df.loc[df['GSM'] == sample, "SRR"].tolist()[0]
    SRR_ID = SRR_ID.split(",")
    #print(SRR_ID)

    return [ "sra/{x}.fastq".format(x=x) for x in SRR_ID ]


rule data_conversion_single:
    input: "sra/{sra}/{sra}.sra"
    priority: 10
    wildcard_constraints:
        sra="[A-Za-z0-9]+"
    output: temp("sra/{sra}.fastq")
    conda:
        "envs/download.yaml"
    resources:
        limit_dump = 1
    shell:"fasterq-dump sra/{wildcards.sra} -O sra" 
    #shell:"fastq-dump --split-files sra/{wildcards.sra} -O sra" 

rule merge_data:
    input: get_merged_input_data
    priority: 20
    resources:
        limit_merge = 1
    output: temp("00_raw_data/{sample}.fq")
    #shell: "cat {input} | pigz > {output}"
    shell: "cat {input}  > {output}"

rule data_clean_single:
    input: "00_raw_data/{sample}.fq"
    priority: 30
    wildcard_constraints:
        sample="[A-Za-z0-9]+"
    output: temp("01_clean_data/{sample}.fq")
    threads: config['fastp_threads']
    log:
        json = os.path.join( "log", '{sample}.json'),
        html = os.path.join( "log", '{sample}.html'),
        logs = os.path.join( "log", "{sample}_fastp.log")
    conda:
        "envs/preprocess.yaml"
    shell:"""
    fastp -w {threads} -i {input}  -o {output}  \
		-j {log.json} -h {log.html} &> {log.logs}
    """

