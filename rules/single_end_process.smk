def get_merged_input_data(wildcards):
    df = metadata_df
    sample = wildcards.sample

    SRR_ID = df.loc[df['GSM'] == sample, "SRR"].tolist()[0]
    SRR_ID = SRR_ID.split(",")
    print(SRR_ID)

    return [ "sra/{x}.fastq".format(x=x) for x in SRR_ID ]


rule data_conversion_single:
    input: "sra/{sra}/{sra}.sra"
    wildcard_constraints:
        sra="[A-Za-z0-9]+"
    output: temp("sra/{sra}.fastq")
    conda:
        "envs/download.yaml"
    shell:"fastq-dump {input} -O sra" 

rule merge_data:
    input: get_merged_input_data
    output: temp("00_raw_data/{sample}.fq.gz")
    shell: "cat {input} | pigz > {output}"

rule data_clean_single:
    input: temp("00_raw_data/{sample}.fq.gz")
    wildcard_constraints:
        sample="[A-Za-z0-9]+"
    params:
        json = lambda wildcards : os.path.join( "qc",  wildcards.sample + '.json'),
        html = lambda wildcards : os.path.join( "qc",  wildcards.sample + '.html')
    output: temp("01_clean_data/{sample}.fq.gz")
    conda:
        "envs/preprocess.yaml"
    shell:"""
    fastp -w {threads} -i {input}  -o {output}  \
		-j {params.json} -h {params.html}
    """

