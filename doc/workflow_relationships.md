# Snakemake Workflow Structure

This document summarizes how the two primary `Snakefile`s in the repository coordinate the available rules, clarifying data flow and module responsibilities.

## 1. Top-Level Snakefile (`Snakefile`)

The default workflow processes generic SRA-based RNA-seq datasets. Key characteristics:

- **Configuration loading**: reads metadata paths, STAR index locations, thread counts, and notification toggles from `config.yaml`.
- **Modular rules**: uses `include` statements for `rules/single_end_process.smk` and `rules/paired_end_process.smk`, which encapsulate preprocessing logic for single-end and paired-end libraries respectively.
- **Local execution**: declares `localrules` (e.g., `all`, `get_sra`, `preload_star_index`) to force selected steps to run on the submission host.

### 1.1 Data Flow Overview

```
get_sra
  └─> data_conversion_(single|pair)
        └─> merge_(data|R1_data|R2_data)
              └─> data_clean_(single|pair)
                    └─> align_and_count
                          ├─> build_bam_index
                          ├─> bamtobw
                          └─> combine_count ──> DGE_analysis ──> cleanup_star_index
```

- `preload_star_index` and `cleanup_star_index` load and release the STAR genome index in shared memory before and after alignment.
- `onsuccess` and `onerror` hooks remove temporary artifacts and optionally send email/Bark/Feishu notifications depending on configuration.

### 1.2 Key Rules and Dependencies

| Rule | Inputs | Outputs | Notes |
| --- | --- | --- | --- |
| `get_sra` | Pre-downloaded `.sra` | `sra/<SRR>/<SRR>.sra` | Verifies source files and creates symlinks in the working directory. |
| `data_conversion_single` / `data_conversion_pair` | `.sra` | `.fastq` / `.fastq` (R1/R2) | Calls `fasterq-dump` to unpack SRA archives. |
| `merge_data` / `merge_R1_data` / `merge_R2_data` | Multiple `.fastq` | `00_raw_data/*.fq` | Concatenates all SRR runs belonging to the same sample. |
| `data_clean_single` / `data_clean_pair` | Raw `.fq` | `01_clean_data/*.fq` | Uses `fastp` for trimming, filtering, and QC reporting. |
| `align_and_count` | Clean `.fq`, STAR index, GTF | `Aligned.sortedByCoord.out.bam`, `ReadsPerGene.out.tab` | Executes STAR to produce sorted BAM files and gene-level counts. |
| `build_bam_index` | BAM | BAM index | Runs `samtools index`. |
| `bamtobw` | BAM + BAI | `04_bigwig/*.bw` | Uses `bamCoverage` to create normalized coverage tracks. |
| `combine_count` | `ReadsPerGene` tables | `03_merged_counts/<DB_ID>.tsv` | Aggregates counts with `pandas` into a single expression matrix. |
| `DGE_analysis` | Expression matrix + metadata | `05_DGE_analysis/<DB_ID>.Rds` | Calls `scripts/DESeq2_diff.R` for differential expression. |
| `cleanup_star_index` | DGE result | `logs/star_index_cleaned.flag` | Frees the STAR index from shared memory. |

## 2. Single-End and Paired-End Rule Modules

- `rules/single_end_process.smk`
  - Defines rules such as `data_conversion_single`, `merge_data`, and `data_clean_single` for single-end datasets.
  - `get_merged_input_data` parses the `SRR` column in metadata, returning the list of FASTQ files per sample.

- `rules/paired_end_process.smk`
  - Provides `data_conversion_pair`, `merge_R1_data`, `merge_R2_data`, and `data_clean_pair` for paired-end libraries.
  - Uses `srr_separator` (comma by default) to split multiple SRR accessions per sample, constructing R1/R2 input lists separately.

Both modules hand their cleaned FASTQ outputs to the `get_input_data` helper defined in `Snakefile`, which feeds `align_and_count` so that alignment is shared between single-end and paired-end branches.

## 3. `Snakefile_ENCODE`

`Snakefile_ENCODE` adapts the workflow to ENCODE-style metadata:

- Expects columns such as `sample`, `R1_file_accession`, and `R2_file_accession`.
- Assumes paired-end reads (values `paired-ended`) and consumes pre-downloaded `00_raw_data/*.fastq.gz` files instead of running SRA extraction.
- Reuses the `fastp`, STAR, and `bamCoverage` rules found in the main workflow but omits STAR index preload/cleanup hooks.

## 4. Extension and Maintenance Tips

- To add new processing stages, create dedicated `.smk` modules under `rules/` and import them from the main `Snakefile` via `include`.
- If new input types (e.g., CRAM) are introduced, mirror the single/paired pattern with appropriate conversion and cleaning rules.
- Conda environments used by the rules live in `envs/` and `rules/envs/`; update those definitions—and the README tool versions—whenever dependencies change.

