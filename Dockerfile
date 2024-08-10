# 基础镜像
FROM ubuntu:22.04

# 设置环境变量
ENV DEBIAN_FRONTEND=noninteractive

# 更新APT并安装必要的软件包
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    bzip2 \
    ca-certificates \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# 安装Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

# 添加conda路径到环境变量
ENV PATH="/opt/conda/bin:$PATH"

# 添加.condarc配置文件
RUN echo 'channels:' > ~/.condarc && \
    echo '  - bioconda' >> ~/.condarc && \
    echo '  - conda-forge' >> ~/.condarc && \
    echo '  - nodefaults' >> ~/.condarc && \
    echo 'custom_channels:' >> ~/.condarc && \
    echo '  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud' >> ~/.condarc && \
    echo '  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud' >> ~/.condarc && \
    echo 'show_channel_urls: true' >> ~/.condarc && \
    echo 'channel_priority: strict' >> ~/.condarc

# 创建并激活Snakemake环境
RUN conda create -n snakemake -c conda-forge python=3.11 && \
    /opt/conda/bin/conda activate snakemake && \
    pip install snakemake==8.16 pandas -i https://pypi.tuna.tsinghua.edu.cn/simple

# 创建并配置GTBA环境，安装指定的软件包
RUN conda create -n gtba -c conda-forge -c bioconda \
    star=2.7.1a samtools fastp sra-tools deeptools pigz

# 创建并配置R环境，安装指定的软件包
RUN conda create -n r432 r-base=4.3.2 && \
    conda install -n r432 bioconductor-deseq2=1.42.0 r-ashr=2.2.63 r-data.table

# 设置R的镜像加速器
RUN echo 'options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))' >> ~/.Rprofile && \
    echo 'options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")' >> ~/.Rprofile

# 切换到工作目录
WORKDIR /workspace

# 默认启动命令
CMD ["/bin/bash"]