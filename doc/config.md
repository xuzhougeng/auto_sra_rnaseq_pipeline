# Configuration Guide

This document centralizes all configuration options used by the RNA-seq automation pipeline. Start by copying `config.yaml.example` to your working directory as `config.yaml`, then update fields as needed.

## Example

```yaml
index: "/data/reference/genome/GRCh38/STAR"
GTF: "/data/reference/genome/GRCh38/Homo_sapiens.GRCh38.95.sort.gtf"
metadata: "metadata"               # Overridden by run.py per-task

# Path to pre-downloaded SRA archives
sra_data_path: "sra"

# Separator for multiple SRR accessions in one row
srr_separator: ","

# Notification settings (optional)
mail: False                         # Email via QQ POP3 only
sender:                             # e.g., "example@qq.com"
sender_password:                    # QQ mail POP3/SMTP auth code
mail_to:                            # Recipient email

bark: False
bark_api:                           # e.g., "https://api.day.app/xxxxxxxxxxxxxxxxx"

feishu: False
feishu_api:

# Per-rule threads
pigz_threads: 10
fastp_threads: 8
star_threads: 20
```

## Field Notes

- `index`, `GTF`
  - Paths to the STAR genome index directory and the annotation GTF file.
  - Build the index ahead of time with STAR `--runMode genomeGenerate`.

- `metadata`
  - When using `run.py`, this value is set automatically per dataset by writing a temporary config; you usually do not need to change it manually.
  - If running Snakemake directly, set `metadata` to the full path of your tab-delimited metadata file ending with `.txt`.

- `sra_data_path`
  - Directory containing pre-downloaded `.sra` archives. The workflow expects files at:
    - `<sra_data_path>/<SRR>/<SRR>.sra`
  - Using NCBI `prefetch -O sra SRRXXXX` will create exactly this structure inside `sra/`.
  - If your files live elsewhere, you can keep them there; the workflow will create symlinks into the working directory.

- `srr_separator`
  - Character used to split multiple SRR accessions in one metadata row (defaults to comma).

- Notifications: `mail`, `bark`, `feishu`
  - Optional run-completion notifications. Set the toggle to `True` and provide the corresponding credentials/API URL.
  - Email currently supports QQ Mail POP3/SMTP only.

- Threads: `pigz_threads`, `fastp_threads`, `star_threads`
  - Per-rule thread settings used by the Snakemake rules for compression, read QC, and alignment.
  - Adjust based on your hardware and job scheduler constraints.

## Advanced: Executors and Profiles

When running on a cluster, you can pass an executor and a profile via `run.py`, for example with Slurm:

```bash
python3 run.py \
  --cores 100 \
  unfinished \
  config.yaml \
  --snakefile Snakefile \
  --executor slurm \
  --executor_profile_path slurm
```

See `slurm/config.yaml` for a sample profile.

## Environments and Tool Versions

Conda environments for rules are defined in:
- `envs/align.yaml` (STAR 2.7.1a, samtools, deepTools)
- `envs/download.yaml` (sra-tools 3.1.1)
- `rules/envs/preprocess.yaml` (fastp)

Refer to these files for exact versions used at runtime.

