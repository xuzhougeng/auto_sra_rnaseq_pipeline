\n+# Run the Pipeline Using the D21122 Example

This tutorial uses the sample metadata file `doc/D21122.txt` to demonstrate the complete process: preparing environments, configuring parameters, preparing SRA data, running the workflow, and obtaining results.

[中文版本](./tutorial.md)

## Environment Setup (Software Installation)

The examples below use conda/micromamba to create three environments and expose them via `PATH` for the workflow:

1) Snakemake environment (for workflow orchestration)

```bash
conda create -n snakemake -c conda-forge python=3.11 -y
conda run -n snakemake python -m \
  pip install snakemake==8.16 pandas -i https://pypi.tuna.tsinghua.edu.cn/simple
```

2) R environment (for differential expression analysis)

```bash
conda create -c conda-forge -n r432 r-base==4.3.2 -y && \
conda install -n r432 -y -c conda-forge -c bioconda \
  bioconductor-deseq2=1.42.0 r-ashr=2.2.63 r-data.table
conda install -n r432 -y -c conda-forge -c bioconda \
  bioconductor-genomeinfodbdata=1.2.11 --force-reinstall
```

3) Alignment and data processing environment (STAR, samtools, sra-tools, fastp, deepTools, etc.)

```bash
conda create -n gpsa -c conda-forge -c bioconda \
  star=2.7.1a samtools fastp sra-tools deeptools pigz -y
```

4) Add required executables to `PATH` before running (replace with your actual paths):

```bash
export PATH=~/micromamba/envs/r432/bin/:~/micromamba/envs/gpsa/bin/:$PATH
```

Recommended: check key tool versions

```bash
snakemake --version
STAR --version
Rscript --version
fasterq-dump --version
samtools --version
```

5) Pre-build the STAR genome index and prepare a GTF annotation (example):

```bash
reference=/path/to/your/genome.fa
index=/path/to/star/index
STAR \
  --runThreadN 50 \
  --runMode genomeGenerate \
  --genomeDir "$index" \
  --genomeFastaFiles "$reference"
```

See more configuration details in `doc/config.md`.

## Step 1: Prepare Working Directory and Config

1) Create a project directory anywhere, for example:

```bash
mkdir -p ~/rna_seq_results/{unfinished,finished,failed}
cd ~/rna_seq_results
```

2) Copy the example config and edit it:

```bash
cp /path/to/auto_sra_rnaseq_pipeline/config.yaml.example ./config.yaml
```

Open `config.yaml` and ensure/modify at least these keys:

- `index`: your STAR index directory, e.g. `/data/reference/genome/GRCh38/STAR`
- `GTF`: path to the reference genome GTF
- `sra_data_path`: directory containing pre-downloaded SRA files (default `sra`)
- `star_threads` / `fastp_threads`: adjust based on your machine

You do not need to set `metadata` in `config.yaml`; `run.py` will write it dynamically for each task.

## Step 2: Place the Sample Metadata File

Copy the repository’s sample metadata into the `unfinished` directory:

```bash
cp /path/to/auto_sra_rnaseq_pipeline/doc/D21122.txt unfinished/
```

Notes:
- `D21122.txt` already contains the required fields `GSE/GSM/SRR/paired`.
- `<DB_ID>` in output filenames is derived from the metadata filename (without extension). In this example, the final merged counts and DGE results will be:
  - `03_merged_counts/D21122.tsv`
  - `05_DGE_analysis/D21122.Rds`

## Step 3: Prepare Raw SRA Data

This workflow does not download `.sra` files online; instead, it reuses locally downloaded SRA archives. Please use the NCBI sra-tools `prefetch` in advance and place files under the directory specified by `sra_data_path` in `config.yaml` (recommended: `sra`).

1) Extract SRR accessions from the metadata:

```bash
cut -f 11 doc/D21122.txt | tail -n +2 | tr ',' '\n' | sort -u > srr.list
```

2) Download each SRR into `sra/` (this produces `sra/SRRXXXX/SRRXXXX.sra`, which matches workflow expectations):

```bash
mkdir -p sra
while read srr; do
  prefetch -O sra "$srr"
done < srr.list
```

If you already have `.sra` files elsewhere, create symlinks as `sra/<SRR>/<SRR>.sra`.

## Step 4: Run the Workflow (Recommended: `run.py`)

`run.py` will automatically:
- Check the Snakemake environment
- Validate metadata format and existence of SRA files
- Perform an initial `--unlock`
- Write a per-metadata temporary config and run Snakemake
- Move metadata files to `finished` or `failed` based on the run result
- Optionally send notifications

Run in the project directory:

```bash
python3 /path/to/auto_sra_rnaseq_pipeline/run.py \
  --cores 20 \
  unfinished \
  config.yaml \
  --snakefile Snakefile \
  --sra_dir sra
```

Optional: using the Slurm executor on a cluster (the repo includes an example at `slurm/config.yaml`):

```bash
python3 /path/to/auto_sra_rnaseq_pipeline/run.py \
  --cores 100 \
  unfinished \
  config.yaml \
  --snakefile Snakefile \
  --executor slurm \
  --executor_profile_path /path/to/auto_sra_rnaseq_pipeline/slurm
```

## Step 5: Inspect Results and Outputs

If the workflow completes successfully, key outputs include:

- Expression matrix: `03_merged_counts/D21122.tsv`
- Per-sample bigWig files: `04_bigwig/<GSM>.bw`
- Differential analysis result (R object): `05_DGE_analysis/D21122.Rds`

Note: `05_DGE_analysis/D21122.Rds` is written by `scripts/DESeq2_diff.R` using `save()` (not `saveRDS()`), so in R use `load()`:

```r
load('05_DGE_analysis/D21122.Rds')
names(kosadata)
head(kosadata$diffResults)
```

QC and run logs are also written under `log/` and `logs/`.

## FAQ

- MissingOutputException or errors due to filesystem latency: run `--unlock` as suggested, or remove incomplete intermediates before retrying; `run.py` attempts `--unlock` automatically before each run.
- SRA files not found: ensure `config.yaml:sra_data_path` and `--sra_dir` contain `sra/<SRR>/<SRR>.sra`; using `prefetch -O sra` produces exactly this structure.
- Metadata validation failed: you must provide `GSM/GSE/SRR/paired`; `paired` must be either `PAIRED` or `SINGLE`. For multiple SRRs, separate with commas.

This completes an end-to-end demonstration based on `doc/D21122.txt`.

