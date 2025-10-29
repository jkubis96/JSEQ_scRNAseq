# Dockerfile dla JSEQ_scRNAseq pipeline
FROM ubuntu:20.04

WORKDIR /app

RUN apt-get update && apt-get install -y locales \
    && localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8 \
    && rm -rf /var/lib/apt/lists/*
ENV LANG en_US.utf8

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    git \
    python3.8 \
    python3-pip \
    software-properties-common \
    libtbb-dev \
    r-base=3.6.3-2 \
    curl \
    openssl \
    libcurl4-openssl-dev \
    libhdf5-dev \
    libhdf5-serial-dev \
    h5utils \
    hdf5-tools \
    hdf5-helpers \
    unzip \
    default-jdk \
    wget \
    samtools \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir \
    pysam==0.16.0.1 \
    biopython==1.78 \
    umi_tools==1.0.1 \
    numba \
    umap-learn==0.5.1 \
    gdown



RUN Rscript -e "install.packages(c('httr','leiden','igraph','readxl','pheatmap','matrix','tidyverse','doparallel','dosnow','stringr','BiocManager','plotly','gridExtra','Seurat','metap','viridis','ape'), repos='https://cran.rstudio.com/')" \
    && Rscript -e "if (!requireNamespace('devtools', quietly = TRUE)) install.packages('devtools', repos='https://cran.rstudio.com/')" \
    && Rscript -e "devtools::install_url('https://github.com/jkubis96/GTF-tool/raw/refs/heads/main/packages/GTF.tool_0.1.2.tar.gz', dependencies = TRUE)" \
	&& Rscript -e "devtools::install_url('https://github.com/jkubis96/CSSG/raw/refs/heads/main/packages/CSSG.toolkit_0.1.0.tar.gz', dependencies = TRUE)"

RUN cd /app/JSEQ_scRNAseq/setup \
    && gdown 1ndAFxTqHUFjhfBEiFuVs-D1SMKBmhfyI \
    && dpkg -i rna-star_2.7.3a+dfsg-1build2_amd64.deb \
    && rm rna-star_2.7.3a+dfsg-1build2_amd64.deb \
    && gdown 1nQzT2deG9l0Ho_Nj0splNv9kZIIh-gYv \
    && dpkg -i fastp_0.20.0+dfsg-1build1_amd64.deb \
    && rm fastp_0.20.0+dfsg-1build1_amd64.deb \
    && gdown 1deqNjK2Ix_O0yPQTnXqD6ShX2WYH5PAz \
    && unzip Drop-seq_tools-2.4.0.zip \
    && rm Drop-seq_tools-2.4.0.zip \
    && mv Drop-seq_tools-2.4.0 DropSeq \
    && chmod -R +x DropSeq

