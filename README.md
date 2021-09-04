## JSEQ® - single cell sequencing analysis tools

<p align="right">
<img  src="https://github.com/Qubix96/JSEQ_scRNAseq/blob/main/setup/fig/logo_jbs.PNG?raw=true" alt="drawing" width="250" />
</p>

#### Authors: Jakub Kubiś & Maciej Figiel
<div align="center">
 Institute of Bioorganic Chemistry<br />
 Polish Academy of Sciences<br />
 Department of Molecular Neurobiology<br />
</div>

<br />
<br />

<div align="justify"> Mechanisms of brain development, diseases, and function are strictly related to dynamically changing cellular composition resulting from nerve cell proliferation, differentiation, and programmed death. Single-cell transcriptomes' resolution brings new and exact insights on cellular gene expressions, activated promoters, diseased populations in the brain, sequence of developmental events, regional structure, and cellular architecture. One of the approaches in single-cell is DropSeq technology, based on encapsulation (Fig.1A) single-cell with single beads including barcode and UMI sequences. The barcode sequence is specific to every single cell (beads), and UMI is specific to the single copy of mRNA (count) in every single cell. </div>

<p align="right">
<img  src="https://github.com/Qubix96/JSEQ_scRNAseq/blob/main/setup/fig/sc.png?raw=true" alt="drawing" width="250" />
</p>

##### Figure 1 Single-cell sequencing in DropSeq technology A) libraries preparing  B) sequencing and analysis  [Created in BioRender]



<div align="justify"> Single-cell sequencing methods (Fig.1B) need to acquire at least 1000+ of cellular transcriptomes in parallel and generate high volumes of data; however, the current bioinformatics pipelines lack features to investigate types of brain cells in high-resolution to understand cellular heterogeneity of the brain. Therefore, we have constructed the JSEQ® bioinformatics pipeline to perform all basic single-cell sequencing analysis and to determine higher resolution cell composition by the combinatorial CSSG algorithm involving repeated rounds of MAST-selected marker genes most accurately define the cell subtypes heterogeneity (p range 0.001 - 0.1). Our pipeline aids to precisely discover the heterogeneity of cells and help determine cell populations particularly vulnerable to brain diseases. The pipeline was developed and tested on publicly available data of 812 945 cells from several mouse brain regions, including the orbitofrontal, entorhinal, auditory and motor cortex, hypothalamus, thalamus, cerebellum, and human brain organoids. The JSEQ® fully automatic pipeline provides an easy and fast way to obtain high-resolution results for single cells using various types of input. </div>

### Pipeline pathway

<p align="center">
<img  src="https://github.com/Qubix96/JSEQ_scRNAseq/blob/main/setup/fig/pipeline.png?raw=true" alt="drawing" width="600" />
</p>

##### The JSEQ® pipeline was prepared and tested on AMD Ryzen Threadripper 24-Core, RAM 256GB, Ubuntu 20.04 LTS. For more information, I invite you to familiarize yourself with the manual JSEQ.
