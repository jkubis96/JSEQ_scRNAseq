#Dockerfile for JSEQ_scRNAseq pipeline
FROM ubuntu:20.04

WORKDIR /app
RUN apt-get update
RUN apt-get update && apt-get install -y locales && rm -rf /var/lib/apt/lists/* \
    && localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8
ENV LANG en_US.utf8
RUN apt-get update
RUN apt-get install -y sudo
RUN sudo apt-get install -y git
RUN git clone https://github.com/jkubis96/JSEQ_scRNAseq.git

RUN sudo apt-get update

RUN sudo apt -y install python3.8
RUN sudo apt -y install python3-pip
RUN sudo apt-get update
RUN pip3 install pysam==0.16.0.1
RUN pip3 install biopython==1.78
RUN pip3 install umi_tools==1.0.1
RUN pip3 install numba
RUN pip3 install umap-learn==0.5.1



RUN sudo DEBIAN_FRONTEND=noninteractive apt-get install -y software-properties-common
RUN sudo apt-get update
RUN sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
RUN sudo apt-get update
RUN sudo apt-get -y install r-base=3.6.3-2
RUN sudo apt-get install -y curl
RUN sudo apt-get -y install openssl
RUN sudo apt-get -y install libcurl4-openssl-dev
RUN sudo apt-get -y install libhdf5-dev
RUN sudo apt-get -y install libhdf5-serial-dev
RUN sudo apt-get -y install h5utils
RUN sudo apt-get -y install hdf5-tools
RUN sudo apt-get -y install hdf5-helpers
RUN sudo apt-get -y install r-cran-httr=1.4.1-1ubuntu1
RUN sudo apt-get -y install r-cran-leiden=0.3.3+dfsg-1
RUN sudo apt-get -y install r-cran-igraph=1.2.4.2-2build1
RUN sudo apt-get -y install r-cran-readxl=1.3.1-2build1
RUN sudo apt-get -y install r-cran-pheatmap=1.0.12-1
RUN sudo apt-get -y install r-cran-matrix=1.2-18-1
RUN sudo apt-get -y install r-cran-tidyverse=1.3.0-1
RUN sudo apt-get -y install r-cran-doparallel=1.0.15-1
RUN sudo apt-get -y install r-cran-dosnow=1.0.18-1
RUN sudo apt-get -y install r-cran-stringr=1.4.0-1
RUN sudo apt-get -y install r-cran-seurat
RUN sudo apt-get -y install r-cran-biocmanager
RUN sudo apt-get -y install r-cran-plotly
RUN sudo apt-get -y install r-cran-gridextra

RUN sudo apt-get update


RUN chmod +rwx $(pwd)/JSEQ_scRNAseq/setup/r_req.R 
RUN sudo -i Rscript $(pwd)/JSEQ_scRNAseq/setup/r_req.R 


RUN sudo apt -y install default-jdk


RUN sudo apt-get -y install samtools=1.10-3

RUN sudo apt-get update
RUN chmod +rwx $(pwd)/JSEQ_scRNAseq/setup
RUN cd JSEQ_scRNAseq/setup \
	&& git clone https://github.com/alexdobin/STAR.git --branch STAR_2.5.0a
RUN cd JSEQ_scRNAseq/setup/STAR/source \
	&& make STAR
RUN sudo apt -y install rna-star


RUN cd JSEQ_scRNAseq/setup \
	&& git clone https://github.com/OpenGene/fastp.git --branch v0.22.0
RUN cd JSEQ_scRNAseq/setup/fastp \
	&& make fastp
RUN sudo apt -y install fastp


RUN sudo apt-get install wget
RUN sudo apt-get update
RUN cd JSEQ_scRNAseq/setup \
	&& wget -O DropSeq.zip https://github.com/broadinstitute/Drop-seq/releases/download/v2.4.0/Drop-seq_tools-2.4.0.zip \
	&& unzip DropSeq \
	&& mv Drop-seq_tools-2.4.0 DropSeq \
	&& rm -r DropSeq.zip \
	&& sudo chmod +rwx DropSeq

RUN cd JSEQ_scRNAseq/setup \
	&& git clone https://github.com/broadinstitute/picard.git --branch 2.26.5
	

RUN cd JSEQ_scRNAseq/setup/picard \
	&& sudo chmod +rwx ../picard \
	&& sudo chmod +x gradlew \
	&& ./gradlew shadowJar


RUN sudo apt-get update -y

RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/analysis_mix
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/analysis_species
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/barcodes_aligment.py
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/convert_mtx_umi.py
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/converter.R
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/functions.R
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/genome_indexing
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/merge_genome.py
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/merge_reads.py
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/project_selection
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/projects
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/raport_mix.Rmd
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/raport_species.Rmd
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/rna_metrics.R
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/seurat_analysis
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/seurat_cluster_mix.R
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/seurat_cluster_species.R
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/rna_metrics.R
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/scripts/umi_extract.py

RUN mkdir $(pwd)/JSEQ_scRNAseq/projects
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/projects
RUN mkdir $(pwd)/JSEQ_scRNAseq/results
RUN sudo chmod +rwx $(pwd)/JSEQ_scRNAseq/results


WORKDIR /app/JSEQ_scRNAseq
CMD $(pwd)/scripts/docker




