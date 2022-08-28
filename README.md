## JSEQ® - single cell sequencing analysis tool

<p align="right">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/main/setup/fig/logo_jbs.PNG?raw=true" alt="drawing" width="250" />
</p>



### Authors: Jakub Kubiś & Maciej Figiel

<div align="center">
 Institute of Bioorganic Chemistry<br />
 Polish Academy of Sciences<br />
 Department of Molecular Neurobiology<br />
</div>


<br />

### Description


<div align="justify"> The JSEQ_scRNAseq® was developed for high-resolution single-cell data analysis </div>

<br />

<p align="center">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/main/setup/fig/sc.png?raw=true" alt="drawing" width="600" />
</p>

#### The single-cell method performance: A) libraries preparing  B) sequencing and analysis 
*Created in BioRender

<br />

<div align="justify"> The JSEQ_scRNAseq® bioinformatics pipeline performs all basic single-cell sequencing analysis and moreover determine higher resolution cell composition by the combinatorial CSSG algorithm involving a great number of iterations of cell cluster specific gene arrays to adjust gene combination to each cluster, where each gene occures in a part of cells inside the cluster and combination of them explain the entire cluster. Our pipeline aids to precisely discover the heterogeneity of cells and help determine cell populations particularly vulnerable to diseases. The pipeline was developed and tested on publicly available data of 812 945 cells. </div>

</br>


### Installation

#### Download:

```
git clone https://github.com/jkubis96/JSEQ_scRNAseq.git
```

#### Get JSEQ_scRNAseq directory:

```
cd JSEQ_scRNAseq
```

#### Run script:

```
./JSEQ
```

#### Choose installation destination [local / docker]:

* before choosing 'docker' make sure that you have docker on your device

```
docker
```

#### Run installation:

```
install
```

<br />

### Usage 

<div align="justify"> For all JSEQ features and how to use it in your analysis see JSEQ_manual.pdf  </div>


#### Pipeline workflow

<p align="center">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/main/setup/fig/pipeline.png?raw=true" alt="drawing" width="800" />
</p>

##### The JSEQ® pipeline was prepared and tested on AMD Ryzen Threadripper 24-Core, RAM 256GB, Ubuntu 20.04 LTS. For more information, I invite you to familiarize yourself with the manual.
