
# RNA-seq automation process pipeline 

[English Vesion](./README_en.md)

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

注意修改config.yaml的配置信息，其中metadata指的是存放metadata文件的目录，metadata文件必须以.txt结尾，否则不识别

另外metadata必须包括如下列， GSM, GSE, gene, SRR, 否则程序运行绝对失败

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


## 配置文件说明


如下参数控制不同规则的运行所需要的线程数

- pigz_threads: 10
- fastp_threads: 8
- star_threads: 20

SRA数据路径配置：
- sra_data_path: "sra"  # 指定已下载SRA文件的存储路径

如下参数和任务完成后的提醒有关

邮件提醒（目前只支持qq邮箱的pop3协议）,设置为False表示不启用，

- mail: False
- sender: 
- sender_password:  
- mail_to: 

IOS bark提醒

- bark: False
- bark_api: #"https://api.day.app/xxxxxxxxxxxxxxxxx"

Feishu提醒

- feishu: False
- feishu_api: 
