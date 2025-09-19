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

## Documentation
- Quickstart tutorial: `doc/tutorial.md`
- Configuration guide: `doc/config.md`
- Workflow structure: `doc/workflow_relationships.md`




## Citation
If this project helps your research, please cite:

Shipeng Guo, Zhougeng Xu, Xiangjun Dong, Dongjie Hu, Yanshuang Jiang, Qunxian Wang, Jie Zhang, Qian Zhou, Shengchun Liu, Weihong Song, GPSAdb: a comprehensive web resource for interactive exploration of genetic perturbation RNA-seq datasets, Nucleic Acids Research, Volume 51, Issue D1, 6 January 2023, Pages D964–D968, https://doi.org/10.1093/nar/gkac1066
