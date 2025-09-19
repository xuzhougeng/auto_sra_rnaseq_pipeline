
# RNA-seq automation process pipeline 

[English Version](./README.md)

## 使用方法

将这个仓库克隆到本地

```bash
git clone https://github.com/xuzhougeng/auto_sra_rnaseq_pipeline.git
```

### 配置环境

snkeamek

```
conda create -n snakemake -c conda-forge python=3.11 -y
conda run -n snakemake python -m \
pip install snakemake==8.16 pandas  -i https://pypi.tuna.tsinghua.edu.cn/simple
```

R

```
conda create -c conda-forge -n r432 r-base==4.3.2 -y && \
conda install -n r432 -y -c conda-forge -c bioconda bioconductor-deseq2=1.42.0 r-ashr=2.2.63 r-data.table
conda install -n r432 -y -c conda-forge -c bioconda bioconductor-genomeinfodbdata=1.2.11 --force-reinstall
```

gpsa

```
# micromamba
conda create -n gpsa -c conda-forge -c bioconda \
 star=2.7.1a samtools fastp sra-tools deeptools pigz -y
``` 


### 构建索引

我们需要构建STAR索引，以及准备参考基因组对应的GTF

```bash
reference=/path/to/your/genome
index=/path/of/index/directory
STAR \
    --runThreadN 50 \
    --runMode genomeGenerate \
    --genomeDir $index \
    --genomeFastaFiles ${reference}
```

### 执行

创建一个项目文件夹，然后将我们仓库中的config.yaml 复制到该目录下，

```bash
mkdir results
cp /path/to/config.yaml results/
cd results
```

注意修改config.yaml的配置信息，其中metadata指的是存放metadata文件的目录，metadata文件必须以.txt结尾，否则不识别。

### metadata文件格式要求

metadata文件必须满足以下要求，否则程序会跳过处理：

**必需列：**
- GSM, GSE, gene, SRR, paired

**数据完整性要求：**
- SRR字段不能为空或NA值
- paired字段不能为空或NA值，且必须是'PAIRED'或'SINGLE'

**示例：**
```
GSE     GSM     gene    method  celline group   group_name      type    platform        SRR     paired
GSE251750       GSM7987315      SMCHD1  ko      LHCN-M2 control LHCN-M2 cells, wildtype RNA     SRA     GPL18573        SRR12345        PAIRED
```

程序会在处理前自动验证metadata文件格式，发现问题的文件会被跳过并显示具体错误信息。

运行方法

```bash
export PATH=~/micromamba/envs/r432/bin/:~/micromamba/envs/gpsa/bin/:$PATH;
ulimit -n 10240

python3 /path/to/auto_sra_rnaseq_pipeline/run.py --cores 120  unfinished config.yaml --SNake

# unfinished指的是未完成任务的位置
# config.yaml是你的配置文件
# 79 表示任务数
```

如果任务失败，再次运行提醒中有 --unlock, 需要运行如下的代码

```bash
snakemake --configfile config.yaml -s auto_sra_rnaseq_pipeline/Snakefile  --unlock
```

如果通过Kill或者ctrl+C的方法停止已经运行的进程，已经移动到metadata中的文件不会移动回unfinished中，需要我们自己动手移动。

## 配置
完整配置说明与示例请参考 `doc/config.md`。

## 引用
如果本项目对你有所帮助，请引用：

Shipeng Guo, Zhougeng Xu, Xiangjun Dong, Dongjie Hu, Yanshuang Jiang, Qunxian Wang, Jie Zhang, Qian Zhou, Shengchun Liu, Weihong Song, GPSAdb: a comprehensive web resource for interactive exploration of genetic perturbation RNA-seq datasets, Nucleic Acids Research, Volume 51, Issue D1, 6 January 2023, Pages D964–D968, https://doi.org/10.1093/nar/gkac1066
