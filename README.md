## JSEQ® - single cell sequencing analysis tools

<p align="right">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/main/setup/fig/logo_jbs.PNG?raw=true" alt="drawing" width="250" />
</p>


#### NEW RELEASE: JSEQ® v2.1.1 
##### Changes:
- new versions of mouse [Release M28 (GRCm39)] and human [Release 39 (GRCh38.p13)] genomes <br />
- added configuration file <br />
- JSEQ® clustering steps optymalized <br />
- CSSG algorithm performance & resolution improved <br />
- improved charts auto-adjusting <br />
- added charts saving in vector format <br />
- adjusted JSEQ results for Seurat v4 compatibility <br />

<br />

### Authors: Jakub Kubiś & Maciej Figiel
<div align="center">
 Institute of Bioorganic Chemistry<br />
 Polish Academy of Sciences<br />
 Department of Molecular Neurobiology<br />
</div>


<br />
<br />

<div align="justify"> Mechanisms of brain development, diseases, and function are strictly related to dynamically changing cellular composition resulting from nerve cell proliferation, differentiation, and programmed death. Single-cell transcriptomes' resolution brings new and exact insights on cellular gene expressions, activated promoters, diseased populations in the brain, sequence of developmental events, regional structure, and cellular architecture. One of the approaches in single-cell is DropSeq technology (Fig.1A), based on encapsulation of single cells with single beads containing a barcode and UMI sequences. The barcode sequence is specific to every single cell (beads), and UMI is specific to the single copy of mRNA (count) in every single cell. Single-cell sequencing methods (Fig.1B)  need to acquire a hugh number of cellular transcriptomes in parallel and generate high volumes of data; however, the current bioinformatics pipelines lack of features to better investigate subtypes of brain cells in high-resolution to understand cellular heterogeneity and complexity of the brain. </div>

<br />

<p align="center">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/main/setup/fig/sc.png?raw=true" alt="drawing" width="600" />
</p>

##### Figure 1 Single-cell sequencing in DropSeq technology A) libraries preparing  B) sequencing and analysis  [Created in BioRender]

<br />


<div align="justify"> Therefore, we have constructed the JSEQ® bioinformatics pipeline to perform all basic single-cell sequencing analysis and moreover determine higher resolution cell composition by the combinatorial CSSG algorithm involving a great number of iterations of cell cluster specific gene arrays to adjust gene combination to each cluster, where each gene occures in a part of cells inside the cluster and combination of them explain the entire cluster. Our pipeline aids to precisely discover the heterogeneity of cells and help determine cell populations particularly vulnerable to brain diseases. The pipeline was developed and tested on publicly available data of 812 945 cells from several mouse brain regions, including the orbitofrontal, entorhinal, auditory and motor cortex, hypothalamus, thalamus, cerebellum, and human brain organoids. </div>

### Pipeline pathway

<p align="center">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/main/setup/fig/pipeline.png?raw=true" alt="drawing" width="900" />
</p>

##### The JSEQ® pipeline was prepared and tested on AMD Ryzen Threadripper 24-Core, RAM 256GB, Ubuntu 20.04 LTS. For more information, I invite you to familiarize yourself with the manual JSEQ.
