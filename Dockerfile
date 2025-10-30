# Dockerfile dla JSEQ_scRNAseq pipeline
FROM ubuntu:20.04

WORKDIR /app

RUN apt-get update && apt-get install -y locales \
    && localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8 \
    && rm -rf /var/lib/apt/lists/*
ENV LANG en_US.utf8
ENV DEBIAN_FRONTEND=noninteractive


RUN apt-get update &&\
	apt-get install -y \
    samtools \
	default-jdk \
	wget \
    python3.8 \
    python3-pip 

RUN pip3 install --no-cache-dir \
    numpy==1.21.6 \
    numba==0.53.1 \
    pysam==0.16.0.1 \
    biopython==1.78 \
    umi_tools==1.0.1 \
    umap-learn==0.5.1 \
    gdown

RUN apt-get update && \
    apt-get install -y software-properties-common curl openssl \
                       libcurl4-openssl-dev libhdf5-dev libhdf5-serial-dev \
                       h5utils hdf5-tools hdf5-helpers unzip libtbb-dev && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository "deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/" && \
    apt-get update && \
    apt-get install -y r-base=3.6.3-2

RUN apt-get install -y r-cran-httr=1.4.1-1ubuntu1 \
                       r-cran-leiden=0.3.3+dfsg-1 \
                       r-cran-igraph=1.2.4.2-2build1 \
                       r-cran-readxl=1.3.1-2build1 \
                       r-cran-pheatmap=1.0.12-1 \
                       r-cran-matrix=1.2-18-1 \
                       r-cran-tidyverse=1.3.0-1 \
                       r-cran-doparallel=1.0.15-1 \
                       r-cran-dosnow=1.0.18-1 \
                       r-cran-stringr=1.4.0-1 \
                       r-cran-biocmanager \
                       r-cran-plotly \
                       r-cran-gridextra \
                       r-cran-seurat=3.1.3-1 \
                       r-cran-metap \
                       r-cran-viridis \
                       r-cran-ape && \
					   apt-get update && \
					   apt-get clean && \
					   rm -rf /var/lib/apt/lists/*


					 

RUN R -e "Sys.setenv(R_INSTALL_STAGED = FALSE); \
          BiocManager::install('MAST'); \
          install.packages(c('umap', 'remotes')); \
          remotes::install_version('textclean', version='0.9.3'); \
          remotes::install_url('https://github.com/jkubis96/GTF-tool/raw/refs/heads/main/packages/GTF.tool_0.1.2.tar.gz', dependencies=TRUE); \
          remotes::install_url('https://github.com/jkubis96/CSSG/raw/refs/heads/main/packages/CSSG.toolkit_0.1.0.tar.gz', dependencies=TRUE)"


RUN mkdir -p /tools \
    && cd /tools \
    && gdown 1ndAFxTqHUFjhfBEiFuVs-D1SMKBmhfyI \
    && dpkg -i rna-star_2.7.3a+dfsg-1build2_amd64.deb || apt-get install -f -y \
    && rm rna-star_2.7.3a+dfsg-1build2_amd64.deb \
    && gdown 1nQzT2deG9l0Ho_Nj0splNv9kZIIh-gYv \
    && dpkg -i fastp_0.20.0+dfsg-1build1_amd64.deb \
    && rm fastp_0.20.0+dfsg-1build1_amd64.deb \
    && gdown 1deqNjK2Ix_O0yPQTnXqD6ShX2WYH5PAz \
    && unzip Drop-seq_tools-2.4.0.zip \
    && rm Drop-seq_tools-2.4.0.zip \
    && mv Drop-seq_tools-2.4.0 DropSeq \
    && chmod -R +x DropSeq


