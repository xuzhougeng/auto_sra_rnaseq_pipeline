

# 使用 D21122 示例运行本流程

[English Version](./tutorial_en.md)

本文以仓库内的示例元数据文件 `doc/D21122.txt` 为例，演示如何从准备环境、配置参数、准备 SRA 数据到运行流程并获取结果的完整步骤。

## 环境准备（软件安装）

以下示例以 conda/micromamba 为例，创建 3 个环境，并通过 PATH 暴露给流程使用：

1) Snakemake 环境（用于调度流程）

```bash
conda create -n snakemake -c conda-forge python=3.11 -y
conda run -n snakemake python -m \\
  pip install snakemake==8.16 pandas -i https://pypi.tuna.tsinghua.edu.cn/simple
```

2) R 环境（用于差异分析）

```bash
conda create -c conda-forge -n r432 r-base==4.3.2 -y && \\
conda install -n r432 -y -c conda-forge -c bioconda \\
  bioconductor-deseq2=1.42.0 r-ashr=2.2.63 r-data.table
conda install -n r432 -y -c conda-forge -c bioconda \\
  bioconductor-genomeinfodbdata=1.2.11 --force-reinstall
```

3) 对齐与数据处理环境（STAR、samtools、sra-tools、fastp、deepTools 等）

```bash
conda create -n gpsa -c conda-forge -c bioconda \\
  star=2.7.1a samtools fastp sra-tools deeptools pigz -y
```

4) 运行前将所需可执行文件加入 PATH（按实际安装路径替换）：

```bash
export PATH=~/micromamba/envs/r432/bin/:~/micromamba/envs/gpsa/bin/:$PATH
```

建议检查关键工具版本：

```bash
snakemake --version
STAR --version
Rscript --version
fasterq-dump --version
samtools --version
```

5) 预先构建 STAR 索引与准备 GTF 注释（示例）：

```bash
reference=/path/to/your/genome.fa
index=/path/to/star/index
STAR \\
  --runThreadN 50 \\
  --runMode genomeGenerate \\
  --genomeDir "$index" \\
  --genomeFastaFiles "$reference"
```

更多配置项说明见 `doc/config.md`。

## 第一步：准备工作目录与配置

1) 在任意位置创建一个项目目录，例如：

```bash
mkdir -p ~/rna_seq_results/{unfinished,finished,failed}
cd ~/rna_seq_results
```

2) 复制示例配置并编辑：

```bash
cp /path/to/auto_sra_rnaseq_pipeline/config.yaml.example ./config.yaml
```

用编辑器打开 `config.yaml`，至少确认与修改以下键：

- `index`: 你的 STAR 索引目录，例如 `/data/reference/genome/GRCh38/STAR`
- `GTF`: 参考基因组的 GTF 文件路径
- `sra_data_path`: 预下载的 SRA 文件所在目录（默认 `sra`）
- `star_threads` / `fastp_threads`: 根据机器资源调整

无需在 `config.yaml` 中手动设置 `metadata`，后续用 `run.py` 会为每个任务动态写入。

## 第二步：放置示例元数据文件

将仓库中的示例元数据复制到 `unfinished` 目录：

```bash
cp /path/to/auto_sra_rnaseq_pipeline/doc/D21122.txt unfinished/
```

说明：
- `D21122.txt` 已包含 `GSE/GSM/SRR/paired` 等必需字段；
- 输出文件名中的 `<DB_ID>` 会取自元数据文件名（去掉扩展名）。因此本例最终的合并计数与差异分析将分别为：
  - `03_merged_counts/D21122.tsv`
  - `05_DGE_analysis/D21122.Rds`

## 第三步：准备 SRA 原始数据

本流程不负责在线下载 `.sra`，而是复用本地已下载的 SRA。请使用 NCBI sra-tools 的 `prefetch` 先行下载，并存放到 `config.yaml` 指定的 `sra_data_path` 目录（推荐直接使用 `sra`）。

1) 从元数据中提取 SRR 列表：

```bash
cut -f 11 doc/D21122.txt | tail -n +2 | tr ',' '\n' | sort -u > srr.list
```

2) 逐个下载到 `sra/` 目录（会自动形成 `sra/SRRXXXX/SRRXXXX.sra` 结构，正好符合流程期望）：

```bash
mkdir -p sra
while read srr; do
  prefetch -O sra "$srr"
done < srr.list
```

若已在其他位置下载好 `.sra` 文件，可创建同名符号链接到 `sra/<SRR>/<SRR>.sra`。

## 第四步：运行流程（推荐用 run.py）

`run.py` 会自动完成以下工作：
- 检查 Snakemake 环境；
- 校验元数据格式与 SRA 文件是否存在；
- 先执行一次 `--unlock`；
- 为每个元数据写入临时配置并运行 Snakemake；
- 按运行结果将元数据移动到 `finished` 或 `failed` 目录；
- 可选发送通知。

在项目目录中执行：

```bash
python3 /path/to/auto_sra_rnaseq_pipeline/run.py \
  --cores 20 \
  unfinished \
  config.yaml \
  --snakefile Snakefile \
  --sra_dir sra
```

可选：如在集群上使用 Slurm 执行器（仓库已自带 `slurm/config.yaml` 示例）：

```bash
python3 /path/to/auto_sra_rnaseq_pipeline/run.py \
  --cores 100 \
  unfinished \
  config.yaml \
  --snakefile Snakefile \
  --executor slurm \
  --executor_profile_path /path/to/auto_sra_rnaseq_pipeline/slurm
```

## 第五步：查看结果与产物

流程成功后，关键产物包括：

- 表达矩阵：`03_merged_counts/D21122.tsv`
- 每个样本的 bigWig：`04_bigwig/<GSM>.bw`
- 差异分析结果（R 格式对象）：`05_DGE_analysis/D21122.Rds`

注意：`05_DGE_analysis/D21122.Rds` 由脚本 `scripts/DESeq2_diff.R` 生成，使用的是 `save()` 而非 `saveRDS()`，因此在 R 中需用 `load()` 读取：

```r
load('05_DGE_analysis/D21122.Rds')
names(kosadata)
head(kosadata$diffResults)
```

同时会在 `log/` 与 `logs/` 目录下生成 QC 与运行日志。

## 常见问题（FAQ）

- MissingOutputException 或文件系统延迟导致报错：可按提示执行 `--unlock`，或在下一次运行前删除未完整的中间文件；`run.py` 已默认在每次开跑前尝试 `--unlock`。
- 找不到 SRA 文件：确保 `config.yaml:sra_data_path` 与 `--sra_dir` 指向的目录下存在 `sra/<SRR>/<SRR>.sra`；使用 `prefetch -O sra` 可直接得到该结构。
- 元数据校验失败：必须提供 `GSM/GSE/SRR/paired`，其中 `paired` 仅允许 `PAIRED` 或 `SINGLE`；若有多条 SRR 用英文逗号分隔。

至此，即可基于 `doc/D21122.txt` 完成一次端到端的流程演示。
