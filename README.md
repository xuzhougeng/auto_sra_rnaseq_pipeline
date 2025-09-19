# RNA-seq Automation Pipeline

[中文版本](./README_zh.md)

## Pipeline Overview

This Snakemake workflow automates SRA-powered RNA-seq processing from raw archives to differential expression results. The table below highlights the main steps, tools, and notable parameters.

| Stage | Tool (version) | Purpose | Key parameters |
| --- | --- | --- | --- |
| Index preload | STAR 2.7.1a | Load genome index into shared memory before alignment | `--genomeLoad LoadAndExit` |
| SRA access | symlink + fasterq-dump (sra-tools 3.1.1) | Reuse locally downloaded `.sra` files and convert to FASTQ | `fasterq-dump sra/<SRR> -O sra` |
| Read merging | coreutils `cat` | Combine multiple SRR runs per sample | Controlled by `srr_separator` (default `,`) |
| Read QC | fastp | Trim/filter reads and generate QC reports | `-w {fastp_threads}` with JSON/HTML logs in `log/` |
| Alignment & counting | STAR 2.7.1a | Produce coordinate-sorted BAM and gene counts | `--quantMode GeneCounts`, `--outSAMtype BAM SortedByCoordinate`, threads from `star_threads` |
| BAM indexing | samtools | Index BAM files | `samtools index -@ 4` |
| Signal tracks | bamCoverage (deepTools) | Generate normalized bigWig coverage files | `--binSize 50`, `--normalizeUsing BPM`, `--effectiveGenomeSize 2913022398` |
| Count matrix | pandas | Merge `ReadsPerGene` tables into a single matrix | Outputs `03_merged_counts/<DB_ID>.tsv` |
| Differential analysis | DESeq2 1.42.0 + ashr 2.2.63 | Shrinkage-based differential expression | `lfcShrink(..., type="ashr")` with groups from metadata |
| Cleanup & notify | STAR, optional email/Bark/Feishu | Release shared-memory index and send run status | Toggle via `mail`, `bark`, `feishu` flags |

## Configuration
See detailed settings and examples in `doc/config.md`.


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

## Citation
If this project helps your research, please cite:

Shipeng Guo, Zhougeng Xu, Xiangjun Dong, Dongjie Hu, Yanshuang Jiang, Qunxian Wang, Jie Zhang, Qian Zhou, Shengchun Liu, Weihong Song, GPSAdb: a comprehensive web resource for interactive exploration of genetic perturbation RNA-seq datasets, Nucleic Acids Research, Volume 51, Issue D1, 6 January 2023, Pages D964–D968, https://doi.org/10.1093/nar/gkac1066
