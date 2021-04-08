
# RNA-seq automation process pipeline 


## Usage

First, clone this repos to local with gti

```bash
git clone https://github.com/xuzhougeng/auto_sra_rnaseq_pipeline.git
# or mirror
git clone https://gitclone.com/github.com/xuzhougeng/auto_sra_rnaseq_pipeline.git
```

Then, install the requirement with mannual or with bioconda

- pigz: parallel gzip
- snakemake
- sra-tool
- fastp
- star

For example, we use bioconda to create a new environment for this pipeline

```bash
conda create -n rna_seq snakemake sra-tools fastp star
# activate the environment
conda activate rna_seq
```

Next, preprare the STAR index and annotation file in gtf format for your genome.

```bash
reference=/path/to/your/genome
index=/path/of/index/directory
STAR \
    --runThreadN 50 \
    --runMode genomeGenerate \
    --genomeDir $index \
    --genomeFastaFiles ${reference}
```

Create a project directory and copy the config.yaml in this repo

```bash
mkdir results
cp /path/to/config.yaml results/
cd results
```

Run this pipepline

```bash
snakemake --restart-times 10 -j 120  --configfile config.yaml -s /path/to/auto_sra_rnaseq_pipeline/Snakefile  -n
```

`-restart-times` will retry 10 times for failer jobs, if your know some rule will stop because of network


