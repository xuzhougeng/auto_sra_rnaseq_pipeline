
# RNA-seq automation process pipeline 

[English Vesion](./README_en.md)

## 使用方法

将这个仓库克隆到本地

```bash
git clone https://github.com/xuzhougeng/auto_sra_rnaseq_pipeline.git
# or mirror
git clone https://gitclone.com/github.com/xuzhougeng/auto_sra_rnaseq_pipeline.git
```

在服务器上安装依赖如下的依赖环境

- [pigz](https://zlib.net/pigz/): 并行化文件压缩
- snakemake
- sra-tools=2.10.8
- fastp
- samtools
- star: 如果需要在多台服务器运行该流程，需要确保star的版本一致
- deeptools
- R/DESeq2, data.table, ashr

我们可以使用bioconda来安装相关环境

```bash
conda create -n rna_seq snakemake sra-tools>2.10.0 fastp star deeptools
# activate the environment
conda activate rna_seq
```

sra-tools的版本必须大于2.10.0, 因为之前版本会直接在目录下输出sra,但是新版本会新建一个以SRA编号
命名的文件夹，然后在该目录输出sra文件

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
python3 /path/to/auto_sra_rnaseq_pipeline/run.py unfinished config.yaml 79
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


如下参数控制不同规则的运行所需要的规则数

- download_threads: 1
- pigz_threads: 10
- fastp_threads: 8
- star_threads: 20

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
