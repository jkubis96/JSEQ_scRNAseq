## JSEQ_scRNAseq© - single cell sequencing analysis tool




<p align="right">
<img  src="https://github.com/jkubis96/Logos/blob/main/logos/jbs_current.png?raw=true" alt="drawing" width="250" />
</p>


### Authors: Jakub Kubiś, Maciej Figiel

<div align="center">
 Institute of Bioorganic Chemistry<br />
 Polish Academy of Sciences<br />
 Department of Molecular Neurobiology<br />
</div>


<br />


## Description

<div align="justify"> The JSEQ_scRNAseq© bioinformatics pipeline performs all basic single-cell sequencing analysis and moreover determine higher resolution cell composition by the combinatorial CSSG algorithm involving a great number of iterations of cell cluster specific gene arrays to adjust gene combination to each cluster, where each gene occures in a part of cells inside the cluster and combination of them explain the entire cluster. Our pipeline aids to precisely discover the heterogeneity of cells and help determine cell populations particularly vulnerable to diseases. The pipeline was developed and tested on publicly available data of 812 945 cells. </div>

</br>

##### The process of single-cell method performance: A) libraries preparing  B) sequencing and analysis 
*Created in BioRender

<p align="center">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/main/setup/fig/sc.png?raw=true" alt="drawing" width="600" />
</p>



<br/>

#### Pipeline workflow

<p align="center">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/v2.3.1/setup/fig/JSEQ.jpg?raw=true" alt="drawing" width="1000" />
</p>

The JSEQ® pipeline was prepared and tested on AMD Ryzen Threadripper 24-Core, RAM 256GB, Ubuntu 20.04 LTS. For more information, I invite you to familiarize yourself with the manual.


## Table of contents

1. [Installation](#installation)
2. [Start application](#start)
3. [ Actions in the application](#action) \
3.1 [Genome](#genome) \
3.2 [Project](#project) \
3.2.1 [FastQ data project](#project1) \
3.2.2 [Expression matrices data project](#project2) \
3.3 [Analysis](#anal) \
3.3.1 [FastQ analysis](#anal1) \
3.3.2 [Expression matrices analysis](#anal2) \
3.3.3 [Manual analysis](#anal3) \
3.4 [Analysis parameters](#analpar) \
3.4.1 [Smart Primer](#analpar1) \
3.4.2 [Configuration file](#analpar2) \
3.4.3 [Barcodes & UMI](#analpar3) \
3.4.4 [Adapters](#analpar4) \
3.5 [Testing mode](#test) 
4. [Additional algorithms](#aa) \
4.1 [GTFtool](#aagtf) \
4.2 [Genes per cell - range estimation algorithm](#gout) \
4.3 [Component selection algorithm](#pca) 
4.4 [CSSG (Cell Subtypes Selection by Gene) algorithm](#aaCSSG)
4.5 [Cell naming algorithm](#aacn)
4.6 [Removing outlier results - algorithms](#out)
4.7 [HDMAP (High Density Manifold Approximation and Projection) - advanced visualization of CSSG subtypes](#hdmap)
5. [Used techniques](#used) 
5.1 [Ribosomal & mitochondrial gene thresholds](#used1) 
5.2 [Data normalization](#used2) 
5.3 [Variable features](#used3) 
5.4 [Dimensionality reduction (PCA)](#used4) 
5.5 [Data clustering](#used5) 
5.6 [Cluster visualization (UMAP)](#used6) 
6. [Performance testing](#perform) 
6.1 [FastQ data analysis](#perform1) 
6.2 [Expresion matrix data analysis](#perform2) 
7. [References tools](#ref)
7.1 [Tools and algorithms](#ref1) 
7.2 [Publications](#ref2) 








<br />

### 1.Installation <a id="installation"></a>

#### Download:

```
git clone https://github.com/jkubis96/JSEQ_scRNAseq.git -b v2.3.2
```

#### Get JSEQ_scRNAseq directory:

```
cd JSEQ_scRNAseq
```

#### Run script:

```
./JSEQ
```

* if you see 'Permission denied' try:

```
cd ..
chmod -R u+r+x JSEQ_scRNAseq/
cd JSEQ_scRNAseq
./JSEQ
```

#### Select 'install'


<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/1.bmp" alt="drawing" width="1000" />
</p>


```
install
```


<br />



#### Instalation options:


<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/2.bmp" alt="drawing" width="1000" />
</p>

* install - The tool will install from a Dockerfile.
* pull - The Docker container will be downloaded from DockerHub [recommended].

<br />

If Docker is not running or is not installed, an error will occur:

 <p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/2.1.bmp" alt="drawing" width="1000" />
</p>



<br />




### 2. Start application <a id="start"></a>


<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/1.bmp" alt="drawing" width="1000" />
</p>

```
start
```

<br />

 If the application was not installed previously, an error will occur:

 <p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/1.1.bmp" alt="drawing" width="1000" />
</p>

<br/>

### 3. Actions in the application <a id="action"></a>


<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/3.bmp" alt="drawing" width="1000" />
</p>


* genome - option for downloading reference genomes appropriate for the planned analysis

```
genome
```
* project - option for creating a new analytical project
```
project
```
* analysis - option to start the analysis for created projects
```
analysis
```

<br/>

### 3.1 Genomes <a id="genome"></a>


<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/4.bmp" alt="drawing" width="1000" />
</p>

* human - downloading and preparing Homo sapiens genome release GRCh37

```
human
```
* mouse - downloading and preparing Mus musculus genome release GRCm39
```
mouse
```
* custom - downloading and preparing cutom genome set in /requirements_file/genome.conf; default custom genome: Schmidtea mediterranea - PRJNA379262
```
custom
```


<br />

### genome.conf



 It is the place where the user can set a basic genome for human and mouse or add a custom genome for analysis, different from the provided options.

 ```
cd requirements_file
nano genome.conf
# or use differen text editor
 ```

 <p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/genome_config.bmp" alt="drawing" width="1000" />
</p>

Additional options:
* three_prime_utr - maximum length of 3'UTR for GTFtool
* three_prime_utr - maximum length of 5'UTR for GTFtool
* _extend - information (boolean) [T/F]; whether GTFtool should run for the selected genome

A description of GTFtool and the reasons for its development are available in the GTFtool documentation.

<br/>

### 3.2 Projects <a id="project"></a>

#### 3.2.1 FastQ data project <a id="project1"></a>


<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/5.bmp" alt="drawing" width="1000" />
</p>

* Project name – Enter your project name. If the name consists of more than one word, use '_' instead of spaces.
* Reads length – Specify the read length based on your sequencing technology. If the two read types [Read1 / Read2] have different lengths, provide the length for Read2.
* Species - Enter the name of the species for which the analysis will be run.
* Estimated number of cells – Indicate the expected number of cells.
* Marker set – Select the appropriate marker set to be used in the current analysis.

If you need more information about markers and naming processes, check the section on Naming.

Warnings: 

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/5.2.bmp" alt="drawing" width="1000" />
</p>

If any of the parameters are incorrect, the project will not be created. The user must provide all parameters again in the correct format. 

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/5.1.bmp" alt="drawing" width="1000" />
</p>

<div align="justify"> Completing the project requires adding data to the project directory. The data must be in the correct format. For a FastQ project, data should be in .fastq format for both Read1 and Read2. Multiple files can be provided and will be combined (e.g., 5x Read1_R1.fq and 5x Read2_R2.fq). Each read file should include a suffix indicating the read number: '_R1' for Read1 and '_R2' for Read2. Ensure that all Read1 files have the suffix '_R1' and all Read2 files have the suffix '_R2'. Failure to follow this format will prevent progress to the next step! </div>


<br/>

#### 3.2.2 Expression matrices data project <a id="project2"></a>

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/6.bmp" alt="drawing" width="1000" />
</p>

* Project name – Enter your project name. If the name consists of more than one word, use '_' instead of spaces.
* Species - Enter the name of the species for which the analysis will be run.
* Marker set – Select the appropriate marker set to be used in the current analysis.

If you need more information about markers and naming processes, check the section on Naming.

* Data format - Select input data format.

<br />


#### Data formats:
-count matrix - matrix of raw counts allowed in *.tsv | *.csv | *.txt formats
-normalized expression matrix - matrix of normalized data allowed in *.tsv | *.csv | *.txt formats
-sparse [genes.tsv, barcodes.tsv, matrix.mtx] - separate files for sparse matrix


Warnings: 

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/5.2.bmp" alt="drawing" width="1000" />
</p>

If any of the parameters are incorrect, the project will not be created. The user must provide all parameters again in the correct format. 

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/6.1.bmp" alt="drawing" width="1000" />
</p>

Completing the project requires adding data to the project directory. The data must be in the correct format as specified above. Failure to follow these formats will prevent progress to the next step!

<br/>




### 3.3 Analysis <a id="anal"></a>

#### 3.3.1 FastQ analysis <a id="anal1"></a>

Select the set number to be analyzed

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/7.bmp" alt="drawing" width="1000" />
</p>


Analysis progress output

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/7.1.bmp" alt="drawing" width="1000" />
</p>

Final files -> projects/'project_name'/results

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/res_1.bmp" alt="drawing" width="1000" />
</p>

Output:

-*.svg | *.png | *.pdf  - graphs of particular analysis results
-*.csv | *.xlsx - tabular data of particular analysis results 
-Completed.bam - reads aligments from STAR
-exp_matrix - directory containing different matrices with results 
-process.log.out - Log of the complete analysis
-Log.final.out - Log from the STAR analysis
-Seurat_object.Rds - input data (before analysis) in SeuratObject format
-Results.Rds - output data (after analysis) in SeuratObject format
-manual_analysis.R - scripts for adjusted manual analysis; more information can be found in the Manual Analysis section
-Report.html - full raport with graphs and description of final analysis

Raport example -> [Report.html](https://htmlpreview.github.io/?https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/example/Report_example1.html) 
or in directory 'JSEQ_scRNAseq/example/Report_example1.html'




<br/>

#### 3.3.2 Expression matrices analysis <a id="anal2"></a>

Select the set number to be analyzed

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/8.bmp" alt="drawing" width="1000" />
</p>


Analysis progress output

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/8.1.bmp" alt="drawing" width="1000" />
</p>

Final files -> projects/'project_name'/results

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/res_2.bmp" alt="drawing" width="1000" />
</p>

Output:

-*.svg | *.png | *.pdf  - graphs of particular analysis results
-*.csv | *.xlsx - tabular data of particular analysis results 
-Completed.bam - reads aligments from STAR
-exp_matrix - directory containing different matrices with results 
-process.log.out - Log of the complete analysis
-Log.final.out - Log from the STAR analysis
-Seurat_object.Rds - input data (before analysis) in SeuratObject format
-Results.Rds - output data (after analysis) in SeuratObject format
-manual_analysis.R - scripts for adjusted manual analysis; more information can be found in the Manual Analysis section
-Report.html - full raport with graphs and description of final analysis


Raport example -> [Report.html](https://htmlpreview.github.io/?https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/example/Report_example2.html) 
or in directory 'JSEQ_scRNAseq/example/Report_example2.html'


<br/>

#### 3.3.3 Manual analysis <a id="anal3"></a>

```
cd projects/'project_name'/results
```
If changes in the analysis are needed or you have to conduct more advanced analyses, run manual_analysis.R and load the data before or after the analysis to conduct your own analysis.

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/manual1.bmp" alt="drawing" width="1000" />
</p>


Users can perform different analysis steps from the beginning or only adjust the obtained results, such as cell names—all information is provided in the script as '#' tips.

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/manual2.bmp" alt="drawing" width="1000" />
</p>

<br/>

### 3.4 Analysis parameters <a id="analpar"></a>

#### 3.4.1 Smart Primer <a id="analpar1"></a>

 ```
cd requirements_file
nano smart_primer
# or use differen text editor
 ```

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/smart.bmp" alt="drawing" width="1000" />
</p>

If your single-cell technique uses different smart primers, make sure to update it accordingly.

<br/>

#### 3.4.2 Configuration file <a id="analpar2"></a>

 ```
cd requirements_file
nano config_file.conf
# or use differen text editor
 ```

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/config.bmp" alt="drawing" width="1000" />
</p>

<div align="justify"> The config_file.conf includes additional parameters necessary for the analysis, allowing users to control threshold settings, explain heterogeneity using the CSSG cell subtypes algorithm, and assess marker importance. It also allows adjustment of hardware usage depending on computational capacity. </div>

<br/> 


#### 3.4.3 Barcodes & UMI <a id="analpar3"></a>

Example of UMI/BARCODE scheme for DropSeq technology:

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/library.bmp" alt="drawing" width="1000" />
</p>


 ```
cd requirements_file
nano barcodes
# or use differen text editor
 ```

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/barcodes.bmp" alt="drawing" width="1000" />
</p>

<div align="justify"> Adjust the UMI/BARCODE layout for your analysis based on the technology used e.g., NadiaDolomite / DropSeq | 10xGenomics | etc. Select the appropriate UMI/BARCODE layout and comment out (#) the other options. The default layout is for DropSeq technology. </div>


<br/> 

#### 3.4.4 Adapters <a id="analpar4"></a>

 ```
cd requirements_file
nano Adapters.fa
# or use differen text editor
 ```

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/adapters.bmp" alt="drawing" width="1000" />
</p>

<div align="justify"> Adapters are part of the sequence necessary for proper sequencing and must be removed before analysis. New sequencing technologies typically remove adapters automatically, but additional quality checks may still be needed. If necessary, add any required adapters for your analysis or modify the existing ones. </div>

<br/>

### 3.5 Testing mode <a id="test"></a>

<div align="justify"> If the user is unsure whether the JSEQ_scRNAseq tool will work properly on their hardware or if it was installed correctly, they can run the testing mode. The tool must be installed before starting testing mode! </div>

<br/>

Testing mode includes:
* Genome downloading
* Genome annotation
* Test data downloading
* Test analysis with report creation

<br/>

In main JSEQ_scRNAseq directory run:

 ```
./testing_mode
 ```

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/t1.bmp" alt="drawing" width="1000" />
</p>


Run the tests; this can take up to several hours...


```
run
```

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/t2.bmp" alt="drawing" width="1000" />
</p>


If everything is okay, the user will see: All tests passed!


<br/>


### 4. Additional algorithms <a id="aa"></a>

#### 4.1 GTFtool <a id="aagtf"></a>

<div align="justify"> GTFtool provides an option for annotation file enrichment, which is connected to the single-cell library structure. Typically, single-cell libraries are created from the 3'UTR side, starting with the polyA tail, followed by the 3'UTR sequence, and then the rest of the transcript. Annotation files do not always include 3' UTR sequence information. To address this, GTFtool includes scripts for extending UTRs, thereby improving the mapping process to genes. Because of variations in library sequence lengths, many sequences may include UTRs, and without UTR information in the annotation file, reads may map to intergenic regions, resulting in the loss of additional count information. You can now choose to extend UTR sequences, and if you do, specify the length to extend for both (if applicable). </div>


<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/utr.bmp" alt="drawing" width="1000" />
</p>

<br/>


#### 4.2 Genes per cell - range estimation algorithm <a id="gout"></a>

<div align="justify"> An algorithm has been created that searches for the highest upper and lower density of results within the entire spectrum of gene/cell dependencies and sets the lower and upper control ranges for cells with a given expression of genes per cell. </div>

##### Graph presents example threshold on testing data

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/492bfadcad2f3d5293bc3aa3defc435b9026dd61/setup/fig/DropletQC.svg" alt="drawing" width="1000" />
</p>


<div align="justify"> The algorithm is used when the values for "down" and "up" in the config_file.conf are set to NA (default setting). Users can set their own thresholds by changing NA to numerical values. The algorithm is important because fixed values often exclude a significant portion of results. For example, if a user sets the lower threshold to 100 genes/cell and 25% of the results have 99 genes/cell, then 25% of the results would be excluded due to a difference of just one gene. </div>

<br/>

#### 4.3 Component selection algorithm <a id="pca"></a>

<div align="justify"> An algorithm has been created that checks where the variance in successive principal components does not change significantly and selects the number of principal components sufficient for the overall analysis of single-cell data clustering. </div>

##### Graph (ElbowPlot) presents example threshold for PCs

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/492bfadcad2f3d5293bc3aa3defc435b9026dd61/setup/fig/Elbow.svg" alt="drawing" width="1000" />
</p>

<br/>


<div align="justify"> In this case, the algorithm selected 10 as the appropriate number of PCs for this data set. To verify the correctness of the selected number of PCs, the JackStraw function was used to indicate the statistical significance of individual PCs (p_val = 0.05). </div>

<br />

##### Graph (JackStraPlot) presents of significient PCs

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/492bfadcad2f3d5293bc3aa3defc435b9026dd61/setup/fig/JackStrawPlot.svg" alt="drawing" width="1000" />
</p>

<br/>


#### 4.4 CSSG (Cell Subtypes Selection by Gene) algorithm <a id="aaCSSG"></a>

<div align="justify"> The CSSG algorithm is a very powerful tool based on the diversity of individual gene occurrence inside estimated cell clusters and divides them into smaller, less heterogeneous subsets. The algorithm can align combinations of single cluster-specific genes to explain a cluster's wide heterogeneity. In this approach, we exclude information about gene expression and transform it into binary data containing information about gene occurrence, where the value '1' means positive gene and '0' means negative gene when it comes to occurrence in every single cell belonging to one cluster. The algorithm uses previously selected by the MAST algorithm (or in a different way) marker genes of each cluster which allows analyzing gene combinations based on cluster-specific genes on a significant level p <= 0.05 of Wilcoxon-Test, which reduces working time and CPU and RAM usage by the algorithm. During the course of the algorithm's work are conducted iterative series of matches binary information about specific marker genes occurrence in each cluster, where two types of matrices are created. In each series of matches for successive genes, the first matrix contains a replication of one marker gene vector with dimensions the same as the second matrix containing all marker genes occurrence information minus the information about the currently studied gene. Such prepared matrices are summed to obtain each gene occurrence combination inside the cluster. The possible number of combinations is between 1 and the factorial from a max number of marker genes for the cluster at a significant level p <= 0.05 of the Wilcoxon-Test (𝑙𝑖𝑚𝐶𝑛→∞ 𝐶𝑛∈1:𝑛!). During the next steps for each row in the summed matrix are calculated statistics (loss_val [value for unexplained cell losses], ht [heterogeneity factor], adj_ht [adjusted heterogeneity factor]). If the relevant conditions are not obtained, it means loss_val reached the default level of 0.05, or if the next combination did not reduce loss_val or increased ht, the next iteration of the combination would start. Used thresholds are necessary to reduce the number of results to only important combinations for the next iterations due to a great number of possible gene combinations, which can be created during analysis that influence the longer working time of the algorithm and more CPU and RAM usage. The results matrix from the previous iteration is the second matrix for the current iteration sum of matrices; thus, we obtain a combination of more genes until the conditions for completing the analysis are met. If the conditions are reached, all gene combinations determined in the analysis, where loss_val <= Q25 [quantile], are saved to the data frame along with their statistics. Then the cluster sub-setting is created inside each cluster using the dominant expression [max(log(CPM+1)] of genes from the best genes combination based on the adj_ht statistic for a given cluster. </div>

<br />

##### CSSG scheme

<br />

<p align="center">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/v2.3.2/setup/fig/cssg.jpg?raw=true" alt="drawing" width="1000" />
</p>

<br />
<br />


##### Example of subclustering using CSSG algorithm:


<br />

<p align="center">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/v2.3.2/setup/fig/cssg_comp.jpg?raw=true" alt="drawing" width="1000" />
</p>

<br />

<br />


##### Heatmap of subtypes with the CSSG markers that were responsible for separating:

<br />

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/ab91211d96a95fd5f78035b103501495e6e2e81a/setup/fig/pheatmap_cells_populations.svg" alt="drawing" width="1000" />
</p>


<br />
<br />


##### Plot of heterogeneity of cells within subtypes:

<br />

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/ab91211d96a95fd5f78035b103501495e6e2e81a/setup/fig/cells_heterogenity.svg" alt="drawing" width="1000" />
</p>

The graph shows the x-axis percent of expressed genes for all cells included in this cell subtype and the y-axis mean expression value for all cells in this subtype. Thanks to this plot, we can check the heterogeneity/homogeneity obtained in the analysis of cell subtypes.


<br />



#### 4.5 Cell naming algorithm <a id="aacn"></a>


 <div align="justify"> A particular Cell Clusters Naming (CCN) algorithm was written for cell naming, which checked the most expressed marker for each cluster. It incluide two types of cell naming: for the user who defined a canonical set of markers for cell population naming – canonical name structure, and for the user who wants only to check the structure of the single-cell data without knowledge about cell populations lineage affiliation and attached to the cell names only automatical selected gene markers – non-canonical name structure. </div>

<br />

<p align="center">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/v2.3.2/setup/fig/naming1.jpg?raw=true" alt="drawing" width="1000" />
</p>



<div align="justify"> Canonical approach is based on three types of markers: class markers, subclass markers, and cluster-specific markers selected by MAST. In the first step, the algorithm checks the cell class based on the most expressed class marker (log(CPM+1)) and gives the name to the class. In the following steps, the algorithm gives names for cell subclass and subtypes in the same way. </div>

<div align="justify"> User markers are in the excel file in JSEQ_scRNAseq/requirements_file/markers_.xlsx. We have two types of markers: the first type is in the first sheet (cell class), and the second is in the second sheet (cell subclass). </div>


Cell class markers:
<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/markers1.bmp" alt="drawing" width="1000" />
</p>


<div align="justify"> You can change your markers depending on your data and experiments, but remember if you write your own markers for cell classes, you have to add '+' before the marker gene name. Gene markers without '+' will not be readable. It is an excellent manner to save markers without using them in analysis. </div>

Cell subtypes markers:

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/markers2.bmp" alt="drawing" width="1000" />
</p>


<div align="justify"> In this case, you need not use additional '+' for markers. If you want only class names based on canonical markers, the sheet for the subclass markers should be empty, and then the non-canonical markers for subclass will be set with algorithm. </div>

Avaiable data sets:

 ```
cd requirements_file
 ```

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/markers_sets.bmp" alt="drawing" width="1000" />
</p>

Currently available datasets are for:

* mature brain structures
* developing brain structures
* non-canonical - for non-canonical naming approaches

Users can define their own marker sets in an Excel file and use them during analysis.

Non-canonical approach is based only on cluster information and cluster-specific markers selected by MAST. 



<br />



#### 4.6 Removing outlier results - algorithms <a id="out"></a>

<div align="justify"> This pipeline contains many checkpoints that protect against lousy quality or badly clustered cells. Even though duplicates removing at the beginning and additional selection points in the pipeline were projected as another quality control step. After dividing cell populations with the CSSG algorithm, obtained cell subtypes in new clusters are checked in terms of proper names. Cell subtypes groups that were poorly marked are renamed to the correct form. Furthermore, when new cell subtypes do not express markers selected by CSSG, they drop out of the analysis. Moreover, the number of cells on the cell subtypes is controlled by the binomial test. If the number of cells is not statistically significant for a given subtype at a level of 0.1, they are excluded from further analysis. The significance level has been set to 0.1 due to potentially sparse subtypes and may be changed to a lower level in the config file. </div>

<br/>

##### Graph presents (red & gray) removed results

<p align="center">
<img src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/492bfadcad2f3d5293bc3aa3defc435b9026dd61/setup/fig/cells_type_threshold.svg">
</p>


<br />


##### Graph presents changes in the number of cells throughout the analysis


<p align="center">
<img src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/492bfadcad2f3d5293bc3aa3defc435b9026dd61/setup/fig/Cells.svg">
</p>

<br />



#### 4.7 HDMAP (High Density Manifold Approximation and Projection) - advanced visualization of CSSG subtypes <a id="hdmap"></a>

<div align="justify"> The HDMAP was developed to enhance the visualization of obtained cell subtypes by clustering subsets using the CSSG algorithm results. This approach builds on UMAP calculations with additional steps. The HDMAP uses the UMAP components 1 and 2 coordinates of the parent clusters, which are averaged and scaled to [0-1] values as centroids for these clusters. Next, the centroid information for each cell is multiplied by a distance factor of 100 for co-clusters. Subsequently, for each cluster, UMAP coordinates are computed using all gene combination results from the CSSG algorithm where loss_pval ≤ Q25 for the given cluster. Finally, the centroid information is combined with the new UMAP components 1 and 2 for each cluster, and the results are plotted. </div>



<br/>

##### UMAP plot

<p align="center">
<img src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/umap.bmp">
</p>


<br />


<br/>

##### HDMAP plot

<p align="center">
<img src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/setup/fig/hdmap.bmp">
</p>


<br />


### 5. Used techniques <a id="used"></a>

#### 5.1 Ribosomal & mitochondrial gene thresholds <a id="used1"></a>

<div align="justify"> Depending on the analysis type: single-cell or single-nuclei; the amount of mitochondrial genes per cell should differ. The default value for mitochondrial genes in JSEQ_scRNAseq is up to 20%. There is no threshold for ribosomal genes. The amounts of mitochondrial and ribosomal genes are shown in the results (graphs). Thresholds can be changed in config_file.conf or after running manual_analysis.R </div>

<br />


##### Mitochondrial & Ribosomal genes

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/327506e4901c19f155a143d89807db1ceac08544/setup/fig/Ribo%7EMito.svg
" alt="drawing" width="1000" />
</p>


<br />


##### Mitochondrial genes

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/327506e4901c19f155a143d89807db1ceac08544/setup/fig/MitoQC.svg" alt="drawing" width="1000" />
</p>


<br />

##### Ribosomal genes

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/327506e4901c19f155a143d89807db1ceac08544/setup/fig/RiboQC.svg" alt="drawing" width="1000" />
</p>



<br />


#### 5.2 Data normalization <a id="used2"></a>

The Seurata Normalize data function with "Log Normalize" normalization, and the scale factor "1e6" (CPM) was used to normalize the data.
##### Formula:
$$
CPM = \frac{\text{count of genes}}{\text{sum of counts per cell}} \times 1000000
$$

$$
\text{NormalizedData} = \log (CPM + 1)
$$

<br />


#### 5.3 Variable features <a id="used3"></a>

<div align="justify"> To calculate the most variable genes, the 'vst' (Variance Stabilizing Transformation) selection method was used with the 'equal_frequency' method based on the Seurat function FindVariableFeatures. </div>

<br />

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/327506e4901c19f155a143d89807db1ceac08544/setup/fig/variable_genes.svg" alt="drawing" width="1000" />
</p>

<br />


#### 5.4 Dimensionality reduction (PCA) <a id="used4"></a>

<div align="justify"> Principal components are a method responsible for reducing the dimensionality of 'p' numerical variables for each 'n' element, increasing interpretability without losing significant information (Jolliffe & Cadima, 2016). This method allows us to manage gene expression matrices for a large number of cells. It is very difficult to compare all genes (p) (for example, in humans and mice, it is about 30,000 genes) in all cells (n). By using the PC method, we can obtain only the important statistical information in the form of PCs, which explain the maximum amount of variance in the data set. Principal components were calculated on scaled data (using the ScaleData function) with the most variable features (genes) using the Seurat function RunPCA. </div>

<br />

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/327506e4901c19f155a143d89807db1ceac08544/setup/fig/PCA_DimPlot_class.svg" alt="drawing" width="1000" />
</p>


<br />

#### 5.5 Data clustering <a id="used5"></a>

<div align="justify"> Data clustering based on previously selected PCs uses two Seurat functions: FindNeighbors and FindClusters. These functions are based on KNN (k-nearest neighbor) with Jaccard distance and SNN (shared nearest neighbor) graph with Louvain modularity optimization. The FindClusters function is set with a resolution of 0.5, the number of starts at 10, and the number of iterations at 1000. Both algorithms are common for single-cell analysis and provide clusters connected with different cell populations (Zhu et al., 2020). In the following steps, based on these clusters, we obtain marker genes for each cluster, name cell populations using known gene markers, and divide cell populations into cell subtypes. </div>


<br />



#### 5.6 Cluster visualization (UMAP) <a id="used6"></a>

<div align="justify"> In contrast to traditional linear dimensionality reduction methods like PCA, UMAP is a non-linear method. The UMAP method, similar to t-SNE, is based on dimensionality reduction and belongs to the non-linear visualization methods. The main role of UMAP algorithms is single-cell data visualization (Narayan et al., 2020). Furthermore, UMAP seems more convenient than the t-SNE method, which can be problematic with large data sets. UMAP optimizes the embedding coordinates of individual data points using iterative algorithms and constructs a high-dimensional graph representation of the data, then optimizes a low-dimensional graph to be as structurally similar as possible. </div>

<br />

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/327506e4901c19f155a143d89807db1ceac08544/setup/fig/UMAP.svg" alt="drawing" width="1000" />
</p>


<br />



### 6. Performance testing <a id="perform"></a>


##### Inputed cells

<p align="center">
<img  src="https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/327506e4901c19f155a143d89807db1ceac08544/setup/fig/Cells.svg" alt="drawing" width="1000" />
</p>

<br/>


#### 6.1 FastQ data analysis <a id="perform1"></a>

Data for testing was used from testing_mode. The analysis report is available at example/Report_example1.html or [Report.html](https://htmlpreview.github.io/?https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/example/Report_example1.html)



#### Performance statistics:

##### Memory usage

<br />


<p align="center">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/v2.3.2/setup/fig/memory_full_analysis_real.jpeg?raw=true" alt="drawing" width="1000" />
</p>


<p align="center">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/v2.3.2/setup/fig/memory_full_analysis.jpeg?raw=true" alt="drawing" width="1000" />
</p>

<br/>

##### CPU time

<br />


<p align="center">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/v2.3.2/setup/fig/cpu_full_analysis.jpeg?raw=true" alt="drawing" width="1000" />
</p>

<br/>


#### 6.2 Expresion matrix data analysis <a id="perform2"></a>

Data for testing was used from testing_mode. The analysis report is available at example/Report_example2.html or [Report.html](https://htmlpreview.github.io/?https://raw.githubusercontent.com/jkubis96/JSEQ_scRNAseq/v2.3.2/example/Report_example2.html)



#### Performance statistics:

##### Memory usage

<br />


<p align="center">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/v2.3.2/setup/fig/memory_full_analysis_real_short.jpeg?raw=true" alt="drawing" width="1000" />
</p>


<p align="center">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/v2.3.2/setup/fig/memory_full_analysis_short.jpeg?raw=true" alt="drawing" width="1000" />
</p>


<br/>

##### CPU time

<br />

<p align="center">
<img  src="https://github.com/jkubis96/JSEQ_scRNAseq/blob/v2.3.2/setup/fig/cpu_full_analysis_short.jpeg?raw=true" alt="drawing" width="1000" />
</p>



<br/>



### 7. References <a id="ref"></a>

#### 7.1 Tools and algorithms <a id="ref1"></a>


* [DropSeqPipe](https://github.com/Hoohm/dropSeqPipe)
* [STAR](https://github.com/alexdobin/STAR)
* [Drop-seq](https://github.com/broadinstitute/Drop-seq)
* [fastp](https://github.com/OpenGene/fastp)
* [Seurat](https://github.com/satijalab/seurat)
* [UMI-tools](https://github.com/CGATOxford/UMI-tools)
* [UMAP](https://umap-learn.readthedocs.io/en/latest/)

Other tools and algorithms that are included in the aforementioned programs

<br/>



#### 7.2 Publications <a id="ref2"></a>

* Narayan, A., Berger, B., & Cho, H. (2020). Density-preserving data visualization unveils dynamic patterns of single-cell transcriptomic variability. BioRxiv, 1–50. https://doi.org/10.1101/2020.05.12.077776
* Zhu, X., Zhang, J., Xu, Y., Wang, J., Peng, X., & Li, H. D. (2020). Single-Cell Clustering Based on Shared Nearest Neighbor and Graph Partitioning. Interdisciplinary Sciences: Computational Life Sciences, 12(2), 117–130. https://doi.org/10.1007/s12539-019-00357-4
* Ian T Jolliffe, Jorge Cadima (2016) Principal component analysis: a review and recent developments. Philos Trans A Math Phys Eng Sci, 13;374(2065):20150202. https://doi.org/10.1098/rsta.2015.0202



<br/>
<br/>


### Have fun JBS©