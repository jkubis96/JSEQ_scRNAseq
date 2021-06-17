## JSEQ® - single cell sequencing analysis tools

<p align="center">
<img  src="https://github.com/Qubix96/Pictures/blob/main/Pictures/Logo.png?raw=true" alt="drawing" width="400" />
</p>

#### Authors: Jakub Kubiś & Maciej Figiel
<div align="center">
 Institute of Bioorganic Chemistry<br />
 Polish Academy of Sciences<br />
 Department of Molecular Neurobiology<br />
</div>


<div align="justify"> Single-cell sequencing is a modern technique of sequencing at the single-cell level. Such an approach gives a lot of possibilities to get knowledge about gene expression in different cells populations during development, carcinogenesis or after treatment with novel gen therapy. But how it works? Actually, so-called single-cell sequencing does not focus on ‘single cell’ but on whole populations of cells which we sequence. So whence obtain we information about single cell? Currently, we have some approaches for single-cell sequencing but the final effect is similar. Results of single-cell sequencing is a great number of cell transcriptome sequences with special sequence signs as ‘barcode’ and ‘UMI’ - Unique Molecular Identifiers. The barcode sequence is different for each cell during sequencing. Furthermore ‘UMI’ is independent of the barcodes sequence which allows distinguishing each copy of mRNA (count) in the cell (for each barcode) and removing duplicated values. So there are two important sequences in bioinformatics single-cell analysis, barcode (about 12 bp) specific to the single cell, and UMI (about 8 bp) specific to the single copy of mRNA (count) in every single cell. Due to this, we can talk about ‘single cell’ sequencing although during sequencing and analysis we have thousands of cells' transcriptomes. After sequencing, obtained results are conducted by multi-stage quality control, statistical analysis (eg. clustering algorithm PCA, kNN, sNN and visualization support algorithm by dimensionality reduction tSNE and UMAP). The JSEQ® fully automatic pipeline provides an easy and fast way to obtain high-resolution results for single cells using various types of input. </div>

### Pipeline pathway

<p align="center">
<img  src="https://github.com/Qubix96/Pictures/blob/main/Pictures/pipeline.png?raw=true" alt="drawing" width="600" />
</p>

##### The JSEQ® pipeline was prepared and tested on AMD Ryzen Threadripper 24-Core, RAM 256GB, Ubuntu 20.04 LTS. For more information, I invite you to familiarize yourself with the manual JSEQ.
