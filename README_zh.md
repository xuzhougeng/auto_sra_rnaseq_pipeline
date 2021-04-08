
# RNA-seq automation process pipeline 


## 使用方法

将这个仓库克隆到本地

```bash
git clone https://github.com/xuzhougeng/auto_sra_rnaseq_pipeline.git
# or mirror
git clone https://gitclone.com/github.com/xuzhougeng/auto_sra_rnaseq_pipeline.git
```

在服务器上安装依赖如下的依赖环境

- pigz: parallel gzip
- snakemake
- sra-tool
- fastp
- star

我们可以使用bioconda来安装相关环境

```bash
conda create -n rna_seq snakemake sra-tools fastp star
# activate the environment
conda activate rna_seq
```

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

```bash
snakemake --restart-times 10 -j 120  --configfile config.yaml -s /path/to/auto_sra_rnaseq_pipeline/Snakefile 
```

`-restart-times 10`表示在失败的时候会重新尝试10次 


