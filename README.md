# RNA-seq Automation Pipeline

[中文版本](./README_zh.md)

## Usage

Clone this repository locally:

```bash
git clone https://github.com/xuzhougeng/auto_sra_rnaseq_pipeline.git
```

### Environment Setup

#### Snakemake

```bash
conda create -n snakemake -c conda-forge python=3.11 -y
conda run -n snakemake python -m \
pip install snakemake==8.16 pandas -i https://pypi.tuna.tsinghua.edu.cn/simple
```

#### R

```bash
conda create -c conda-forge -n r432 r-base==4.3.2 -y && \
conda install -n r432 -y -c conda-forge -c bioconda bioconductor-deseq2=1.42.0 r-ashr=2.2.63 r-data.table
conda install -n r432 -y -c conda-forge -c bioconda bioconductor-genomeinfodbdata=1.2.11 --force-reinstall
```

#### gpsa Environment

```bash
# micromamba
conda create -n gpsa -c conda-forge -c bioconda \
    star=2.7.1a samtools fastp sra-tools deeptools pigz -y
```

### Build the STAR Index

Prepare a STAR index and the corresponding genome annotation (GTF):

```bash
reference=/path/to/your/genome
index=/path/of/index/directory
STAR \
    --runThreadN 50 \
    --runMode genomeGenerate \
    --genomeDir $index \
    --genomeFastaFiles ${reference}
```

### Run the Pipeline

Create a project directory and copy `config.yaml` from this repository into it:

```bash
mkdir results
cp /path/to/config.yaml results/
cd results
```

Update `config.yaml` with your settings. The `metadata` entry must point to the directory that stores metadata files, and every metadata file must end with `.txt`.

### Metadata File Requirements

Metadata files must satisfy the following rules or they will be skipped:

**Required columns**
- GSM, GSE, gene, SRR, paired

**Data integrity**
- `SRR` cannot be empty or `NA`
- `paired` cannot be empty or `NA`, and must be either `PAIRED` or `SINGLE`

**Example**
```
GSE     GSM     gene    method  celline group   group_name      type    platform        SRR     paired
GSE251750       GSM7987315      SMCHD1  ko      LHCN-M2 control LHCN-M2 cells, wildtype RNA     SRA     GPL18573        SRR12345        PAIRED
```

The pipeline validates metadata files before they are processed. Files with problems are skipped and the specific error will be reported.

Run the pipeline:

```bash
export PATH=~/micromamba/envs/r432/bin/:~/micromamba/envs/gpsa/bin/:$PATH;
ulimit -n 10240

python3 /path/to/auto_sra_rnaseq_pipeline/run.py --cores 120 unfinished config.yaml --SNake

# unfinished is the directory that stores pending tasks
# config.yaml is your configuration file
# 79 represents the number of tasks
```

If a run fails and suggests using `--unlock`, execute the following command:

```bash
snakemake --configfile config.yaml -s auto_sra_rnaseq_pipeline/Snakefile --unlock
```

If you stop the process with `kill` or `Ctrl+C`, files already moved into the metadata directory will not be moved back to `unfinished`, so you'll need to move them manually.

## Configuration File Notes

Thread counts per rule:
- `pigz_threads`: 10
- `fastp_threads`: 8
- `star_threads`: 20

SRA data path configuration:
- `sra_data_path`: `"sra"`  # directory containing downloaded SRA files

Notification settings:

Email (currently only QQ Mail POP3 is supported). Set to `False` to disable.
- `mail`: False
- `sender`:
- `sender_password`:
- `mail_to`:

iOS Bark notifications
- `bark`: False
- `bark_api`: `#"https://api.day.app/xxxxxxxxxxxxxxxxx"`

Feishu notifications
- `feishu`: False
- `feishu_api`:
