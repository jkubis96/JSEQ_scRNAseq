#Before start you have to install required packages 
#Pipeline include Seurat version 3.1.5
#You can use other versions of seurat but some functions may not work properly
#Belowe script was adjusted for use with Seurat 4
#but due to different Seurat and R versions
#can be the difference in results in comparison to pipeline results
#More information and instruction on https://satijalab.org/seurat/articles/install.html

library(Seurat)
library(patchwork)
library(tidyverse)
library(doSNOW)
library(foreach)
library(doParallel)
library(stringr)

#README
#If you are here it means you want to improve analysis parameters or adjust obtained results
#Below is presented analytic part of the pipeline
#If you want to start the analysis from the beginning you choose UMI <- readRDS('Seurat_object.rds')
#If you want to adjust obtained results (e.g. cell names, plots) you choose UMI <- readRDS('Results.rds')

#Paths and arguments from env
#Load markers

markers <- '../../../requirements_file/markers.xlsx'

#If you use custom markers change markers file name ^

markers_class <- readxl::read_xlsx(markers, sheet = 1)
markers_subclass <- readxl::read_xlsx(markers, sheet = 2, col_names = F)

#Create new directory for manual analysis and load requirements 
{  
  
  dir.create('manual_results')
  OUTPUT <- file.path('manual_results')
  functions <- '../../../scripts/functions.R'
  cssg <- '../../../scripts/cssg.R'
  source(functions, local = T)
  source(cssg, local = T)
  species <- 'mix'
  
}

#Configuration file 

## Load configuration from config_file.conf or !!!provide other variables changing value after '<-'!!!
  
  
{
  conf_file <- read.csv(file = '../../../requirements_file/config_file.conf', header = F, sep = ':', row.names = 1)
  
  mt_per <- as.numeric(as.character(conf_file$V2[grep(pattern = 'mt_per', rownames(conf_file))]))
  
  down_tr <- as.numeric(as.character(conf_file$V2[grep(pattern = 'down', rownames(conf_file))]))
  
  up_tr <- as.numeric(as.character(conf_file$V2[grep(pattern = 'up', rownames(conf_file))]))
  
  mt_cssg <- as.character(conf_file$V2[grep(pattern = 'mt_cssg', rownames(conf_file))])
  
  s_factor <- as.numeric(as.character(conf_file$V2[grep(pattern = 's_factor', rownames(conf_file))]))
  
  m_val <- as.numeric(as.character(conf_file$V2[grep(pattern = 'm_val', rownames(conf_file))]))
  
  max_genes <- as.numeric(as.character(conf_file$V2[grep(pattern = 'max_genes', rownames(conf_file))]))
  
  max_combine <- as.numeric(as.character(conf_file$V2[grep(pattern = 'max_combine', rownames(conf_file))]))
  
  loss_pval <- as.numeric(as.character(conf_file$V2[grep(pattern = 'loss_pval', rownames(conf_file))]))
  
  p_bin <- as.numeric(as.character(conf_file$V2[grep(pattern = 'p_bin', rownames(conf_file))]))
  
}

###########################################################################################################################################################

UMI_human <- readRDS('Seurat_object_human.rds')
UMI_mice <- readRDS('Seurat_object_mice.rds')

#Data ^ - for new manual analysis 
#If you choose this option start from PART A of pipeline (code line 100-561)

UMI <- readRDS('Results.rds')

#Data ^ - for adjusting

#In some causes in your results you can see cell subtypes with wrong names obtained based on top marker
#If you want to change cell names run fallow code

#Cells renameing

Idents(UMI) <- gsub(pattern = 'Old_name', replacement = 'New_name', x = Idents((UMI)))
#				          	Write old name ^         Write new name ^

#Go to PART B and generate new plots (code line 562-684)

###########################################################################################################################################################

                                    #PART "A" OF PIPELINE FOR NEW ANALYSIS


#Graphs of counts content

UC_plot <- VlnPlot(UMI, features = c("nGenes", "nCounts"), ncol = 2)
svg(file.path(OUTPUT, "counts~genes.svg"), width=15, height=10)
UC_plot
rm(UC_plot)
dev.off()

MR_plot <- VlnPlot(UMI, features = c("RiboPercent", "MitoPercent"), ncol = 2)
svg(file.path(OUTPUT, "Ribo~Mito.svg"), width=15, height=10)
MR_plot
rm(MR_plot)
dev.off()


CG_plot <- FeatureScatter(UMI, feature1 = "nCounts", feature2 = "nGenes")
svg(file.path(OUTPUT, "counts~genes_QC.svg"), width=15, height=10)
CG_plot
rm(CG_plot)
dev.off()


rm(UMI)

###########################################################################################################################################################

#Droplet content and QC human
#Human
n_gen_human <- UMI_human@meta.data$nGenes[UMI_human@meta.data$nGenes > down_tr]

if (is.na(up_tr)) {
  n_gen_human <- as.numeric(mean(n_gen_human))*2 + 1.5*IQR(as.numeric(UMI_human@meta.data$nGenes))
} else {
  n_gen_human <- up_tr
}

QC_UMI_human <- data.frame()
QC_UMI_human <- as.data.frame(UMI_human$nGenes)
QC_UMI_human$V2 <- UMI_human$MitoPercent
QC_UMI_human$V3 <- UMI_human$RiboPercent

colnames(QC_UMI_human) <- c('nGenes','MitoPercent','RiboPercent')

QC_UMI_human$Mito_Status[QC_UMI_human$MitoPercent > mt_per] <- paste0('> ' , mt_per , '%')
QC_UMI_human$Mito_Status[QC_UMI_human$MitoPercent <= mt_per] <- 'Proper'

QC_UMI_human$nGenes_Status[UMI_human$nGenes < down_tr] <- 'Empty'
QC_UMI_human$nGenes_Status[UMI_human$nGenes > n_gen_human] <- 'Double'
QC_UMI_human$nGenes_Status[UMI_human$nGenes >= down_tr & UMI_human$nGenes <= n_gen_human] <- 'Proper'

QC_UMI_human$Ribo_Status[QC_UMI_human$RiboPercent == 0] <- '0%'
QC_UMI_human$Ribo_Status[QC_UMI_human$RiboPercent > 0] <- '> 0 %'
QC_UMI_human$species <- 'Human'

#Mouse
n_gen_mice <- UMI_mice@meta.data$nGenes[UMI_mice@meta.data$nGenes > down_tr]

if (is.na(up_tr)) {
  n_gen_mice <- as.numeric(mean(n_gen_mice))*2 + 1.5*IQR(as.numeric(UMI_mice@meta.data$nGenes))
} else {
  n_gen_mice <- up_tr
}


QC_UMI_mice <- data.frame()
QC_UMI_mice <- as.data.frame(UMI_mice$nGenes)
QC_UMI_mice$V2 <- UMI_mice$MitoPercent
QC_UMI_mice$V3 <- UMI_mice$RiboPercent

colnames(QC_UMI_mice) <- c('nGenes','MitoPercent','RiboPercent')

QC_UMI_mice$Mito_Status[QC_UMI_mice$MitoPercent > mt_per] <- paste0('> ' , mt_per , '%')
QC_UMI_mice$Mito_Status[QC_UMI_mice$MitoPercent <= mt_per] <- 'Proper'

QC_UMI_mice$nGenes_Status[QC_UMI_mice$nGenes < down_tr] <- 'Empty'
QC_UMI_mice$nGenes_Status[QC_UMI_mice$nGenes > n_gen_mice] <- 'Double'
QC_UMI_mice$nGenes_Status[QC_UMI_mice$nGenes >= down_tr & UMI_mice$nGenes <= n_gen_mice] <- 'Proper'

QC_UMI_mice$Ribo_Status[QC_UMI_mice$RiboPercent == 0] <- '0%'
QC_UMI_mice$Ribo_Status[QC_UMI_mice$RiboPercent > 0] <- '> 0 %'
QC_UMI_mice$species <- 'Mouse'


QC_UMI <- rbind(QC_UMI_mice, QC_UMI_human)

DQC <- ggplot()+
  geom_point(QC_UMI, mapping = aes(x = nGenes, y = nGenes, color = nGenes_Status))+
  ylab("Number of genes for each cells") +
  xlab("Number of genes for each cells")+
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))+
  labs(color='Droplet content') +
  facet_grid(.~QC_UMI$species)

svg(filename = file.path(OUTPUT,'DropletQC.svg'), width = 10, height = 7)
DQC
dev.off()
rm(DQC)

MQC <- ggplot()+
  geom_point(QC_UMI, mapping = aes(x = MitoPercent, y = MitoPercent , color = Mito_Status))+
  ylab("% MitoRNA") +
  xlab("% MitoRNA")+
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))+
  labs(color='% Content threshold') +
  facet_grid(.~QC_UMI$species)

svg(filename = file.path(OUTPUT,'MitoQC.svg'), width = 10, height = 7)
MQC
dev.off()
rm(MQC)

RQC <- ggplot()+
  geom_point(QC_UMI, mapping = aes(x = RiboPercent, y = RiboPercent, color = Ribo_Status))+
  ylab("% RiboRNA") +
  xlab("% RiboRNA")+
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))+
  labs(color='% Content threshold') +
  facet_grid(.~QC_UMI$species)

svg(filename = file.path(OUTPUT,'RiboQC.svg'),  width = 10, height = 7)
RQC
dev.off()
rm(RQC)
rm(QC_UMI)
rm(QC_UMI_mice)
rm(QC_UMI_human)

###########################################################################################################################################################

#Selecting right cells
#It is the place where you can change thresholds for cells for further analysis

UMI_human <- subset(UMI_human, subset = nGenes > down_tr & nGenes <= n_gen_human & MitoPercent < mt_per)
#TRESHOLDS:					                          LOWER ^               UPPER ^              MITOGENES ^

n_gen_human <- max(as.numeric(UMI_human@meta.data$nGenes))*0.95
cells_number_human <- length(Idents(UMI_human))

UMI_mice <- subset(UMI_mice, subset = nGenes > down_tr & nGenes <= n_gen_mice & MitoPercent < mt_per)
#TRESHOLDS:					                        LOWER ^             UPPER ^              MITOGENES ^

n_gen_mice <- max(as.numeric(UMI_mice@meta.data$nGenes))*0.95
cells_number_mice <- length(Idents(UMI_mice))


###########################################################################################################################################################


UMI <- merge(x = UMI_human, y = UMI_mice)
UMI <- NormalizeData(UMI, normalization.method = "LogNormalize", scale.factor = 1e6)

###########################################################################################################################################################


###########################################################################################################################################################

n_gen <- mean(n_gen_human + n_gen_mice)
UMI <- FindVariableFeatures(UMI, selection.method = "vst", nfeatures = n_gen, binning.method = 'equal_frequency')

# Identify the 10 most highly variable genes

top20 <- head(VariableFeatures(UMI), 20)

plot1 <- VariableFeaturePlot(UMI)

plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)

svg(file.path(OUTPUT, "variable_genes.svg"), width=10, height=7)
plot2
dev.off()
rm(plot2)

###########################################################################################################################################################

all.genes <- rownames(UMI)
UMI <- ScaleData(UMI, features = all.genes)

UMI <- RunPCA(UMI, features = VariableFeatures(object = UMI))

################################

Elbow <- ElbowPlot(UMI, ndims = 50)

dims <- as.data.frame(Elbow$data$stdev)

#select the most variable reduction

dim <- dim_reuction_pcs(dims)

svg(file.path(OUTPUT, "Elbow.svg"), width=10, height=7)
Elbow + geom_vline(xintercept = dim, color = 'red')
dev.off()
rm(Elbow)

###########################################################################################################################################################

UMI <- JackStraw(UMI, num.replicate = 10, dims = dim)
UMI <- ScoreJackStraw(UMI, dims = 1:dim)

#Select significient PCs
jc <- as.data.frame(UMI@reductions$pca@jackstraw@overall.p.values)
jc <- jc[jc$Score < 0.05,]
dim <- as.vector(jc$PC)


svg(file.path(OUTPUT, "JackStrawPlot.svg"), width=10, height=7)
JackStrawPlot(UMI, dims = dim)
dev.off()

UMI <- FindNeighbors(UMI, dims = dim, reduction = 'pca', nn.method="rann")
UMI <- FindClusters(UMI, resolution = 0.5, n.start = 10, n.iter = 1000)


UMI <- RunUMAP(UMI, dims = dim, umap.method = "umap-learn")

UMAP_coordinates <- as.data.frame(UMI@reductions$umap@cell.embeddings)
cluster_idents <- as.data.frame(Idents(UMI))

width <- 15 + (length(unique(Idents(UMI))))/7


species_plot <- DimPlot(UMI, reduction = "umap", group.by = 'orig.ident', raster = FALSE)
clusters_plot <- DimPlot(UMI, reduction = "umap", raster = FALSE)
mutual_plot <- species_plot + clusters_plot

svg(file.path(OUTPUT, "UMAP.svg"), width = width, height = 15)
mutual_plot
dev.off()
rm(mutual_plot)

###########################################################################################################################################################



#find markers for every cluster compared to all remaining cells, report only the positive ones

print('Searching for cluster marker genes')

UMI.markers <- FindAllMarkers(UMI, only.pos = TRUE, min.pct = 0.10, test.use = 'MAST',  logfc.threshold = 0.10, base = exp(1))

if (sum(as.numeric(levels(UMI))) != sum(unique(as.integer(UMI.markers$cluster)-1))) {
  UMI.markers <- FindAllMarkers(UMI, only.pos = TRUE, min.pct = 0.001 , logfc.threshold = 0.10, test.use = 'MAST',  base = exp(1))
}

top10 <- UMI.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

if (length(markers_subclass) != 0) {
  top10 <- top10[!toupper(top10$gene) %in% toupper(markers_subclass), ]
}

top_sig <- UMI.markers

if (mt_cssg == "exclude") {
  top_sig <- top_sig[!top_sig$gene %in% top_sig$gene[grepl('t-', top_sig$gene)],]
  top_sig <- top_sig[!top_sig$gene %in% top_sig$gene[grepl('T-', top_sig$gene)],]
  top_sig <- top_sig[!top_sig$gene %in% top_sig$gene[grepl('t.', top_sig$gene)],]
  top_sig <- top_sig[!top_sig$gene %in% top_sig$gene[grepl('T.', top_sig$gene)],]
}


subclasses_marker <- UMI.markers %>% group_by(cluster) %>% top_n(n = 1000, wt = avg_logFC)

subclasses_marker_report <- UMI.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

print('Cluster genes - DONE')

##Cells cluster naming with top genes (different between cell groups)

###########################################################################################################################################################
### Most variable genes select and cell subtypes nameing 


###########################################################################################################################################################
#Cells subtypes selection

tmp <- GetAssayData(UMI, slot = 'data')
colnames(tmp) <- UMI@active.ident

marker_df <- heterogenity_select(cells_wide_df = tmp, marker_df = top_sig, heterogenity_factor = s_factor, p_val =  m_val, max_genes =  max_genes, select_stat = 'p_val')

CSSG_df <- CSSG_markers(cells_wide_df = tmp, markers_df = marker_df$marker_df, max_combine = max_combine, loss_pval = loss_pval)
hd_factors <- hd_cluster_factors(UMI, CSSG_df)

write.table(CSSG_df, file = file.path(OUTPUT, "CSSG_marker.csv"), sep = ',')

###########################################################################################################################################################

##Cell nameing

##Create markers DF


print('Single cell naming')

cell_names <- cssg_naming(UMI, CSSG_df)

print('Naming - DONE')

###########################################################################################################################################################
#Cluster average expression for nameing

###########################################################################################################################################################

average_expression <- aggregation_num(UMI)

###########################################################################################################################################################

print('Clusters naming')

cluster_nameing(matrix_a = average_expression, markers = markers_class)

clust_names <- colnames(average_expression)


###########################################################################################################################################################
#Subclass naming

marker_list <- subcluster_naming(average_expression, markers_subclass, marker_df$heterogenity_markers_df, top10)


###########################################################################################################################################################
#Repair subclass_names

new.cluster.ids <- paste(clust_names, marker_list$sub_names)
names(new.cluster.ids) <- levels(UMI)
UMI <- RenameIdents(UMI, new.cluster.ids)

print('Naming - DONE')

###########################################################################################################################################################

colnames(average_expression) <- new.cluster.ids
dir.create(path = file.path(OUTPUT,'exp_matrix'))
write.table(average_expression, file = file.path(OUTPUT, "exp_matrix/class_average_expression.csv"), sep = ',')

###########################################################################################################################################################

#PCA plot and UMAP plot with names

#Class

Tmp_idents <- Idents(UMI)

part_name_1 <- sub(" .*", "", Tmp_idents)
part_name_2 <- sub('.*? ', "", Tmp_idents)

part_name_2 <- toupper(part_name_2)

Tmp_idents_species <- paste(part_name_1, part_name_2)
Idents(UMI) <- Tmp_idents_species


meta_data <- as.data.frame(Idents(UMI))
meta_data$idents <- rownames(meta_data)
colnames(meta_data) <- c('subclass', 'idents')
UMAP_coordinates$idents <- rownames(UMAP_coordinates)
cluster_idents$idents <- rownames(cluster_idents)
meta_data <- merge(meta_data, UMAP_coordinates, by = 'idents')
meta_data <- merge(meta_data, cluster_idents, by = 'idents')
colnames(meta_data)[5] <- 'cluster'


marker_cell_names <- meta_data[c('cluster', 'subclass')] %>% 
  distinct()


subclasses_marker <- merge(subclasses_marker, marker_cell_names, by = 'cluster')

write.table(subclasses_marker, file = file.path(OUTPUT, "markers_clusters.csv"), sep = ',')

subclasses_marker_report <- merge(subclasses_marker_report, marker_cell_names, by = 'cluster')

write.table(subclasses_marker_report, file = file.path(OUTPUT, "subclasses_marker_report.csv"), sep = ',')


width <- 20 + (length(unique(Idents(UMI))))/5


svg(file.path(OUTPUT, "PCA_DimPlot_class.svg"), width = width, height = 15)
DimPlot(UMI, reduction = "pca", raster = FALSE)
dev.off()

svg(file.path(OUTPUT, "UMAP_with_DE_gene_class.svg"), width = width, height = 15)
DimPlot(UMI, reduction = "umap", raster = FALSE) 
dev.off()

Idents(UMI) <- Tmp_idents

rm(Tmp_idents)
rm(Tmp_idents_species)


#Subtypes


new.names <- paste0(UMI@active.ident,' - ', firstup(tolower(cell_names)))


Idents(UMI) <- new.names


###########################################################################################################################################################


#Average expression matrix populations

rm(average_expression)

average_expression <- aggregation_chr(UMI)

###########################################################################################################################################################

print('Checking and renaming subtypes')


renamed_list  <- name_repairing(UMI, average_expression, markers_class, markers_subclass, species)
Idents(UMI) <- renamed_list$Renamed_idents

idents_subtypes <- as.data.frame(Idents(UMI))
idents_subtypes$idents <- rownames(idents_subtypes)
colnames(idents_subtypes) <- c('subtypes', 'idents')
meta_data <- merge(meta_data, idents_subtypes, by = 'idents')
meta_data$subtypes[grep(pattern = 'BAD!', meta_data$subtypes)] <- '-'
meta_data$subtypes[grep(pattern = 'Bad!', meta_data$subtypes)] <- '-'

write.table(meta_data, file = file.path(OUTPUT, "meta_data.csv"), sep = ',')

print('Checking - DONE')

###########################################################################################################################################################

print('QC of subtypes')

data <- bin_cell_test(p_val = p_bin, renamed_list = renamed_list)

threshold <- cell_stat_graph(data$data)

height <- 10 + (length(unique(Idents(UMI))))/5

svg(filename = file.path(OUTPUT,'cells_type_threshold.svg'), width = 15, height = height)
threshold
dev.off()
rm(threshold)

#save bad cells

bad.subnames <- c(as.character(data$bad.subnames), as.character(data$below.names))
bad.subnames <- unique(as.character(bad.subnames))
bad.subnames <- bad.subnames[bad.subnames %in% renamed_list$Renamed_idents]

UMI_unknow <- subset(UMI, idents = bad.subnames)

bad_cells <- as.data.frame(GetAssayData(object = UMI_unknow, slot = 'counts'))


write.table(bad_cells, file = file.path(OUTPUT, "exp_matrix/unknow_cells_count_matrix.csv"), sep = ',')

rm(bad_cells)

###########################################################################################################################################################

right.names <- unique(renamed_list$Renamed_idents[!as.character(renamed_list$Renamed_idents) %in% as.character(bad.subnames)])


UMI <- subset(UMI, idents = right.names)

print('DONE')

###########################################################################################################################################################
#Subtype markers selection

UMI.subtypes <- FindAllMarkers(UMI, only.pos = TRUE, min.pct = 0.20, logfc.threshold = 0.25, test.use = 'MAST',  base = exp(1))

if (length(unique(Idents(UMI))) != length(unique(UMI.subtypes$cluster))) {
  UMI.subtypes <- FindAllMarkers(UMI, only.pos = TRUE, min.pct = 0.10 , logfc.threshold = 0.25, test.use = 'MAST',  base = exp(1))
} 

if (length(unique(Idents(UMI))) != length(unique(UMI.subtypes$cluster))) {
  UMI.subtypes <- FindAllMarkers(UMI, only.pos = TRUE, min.pct = 0.001 , logfc.threshold = 0.25, test.use = 'MAST',  base = exp(1))
}

subtypes_marker <- UMI.subtypes %>% group_by(cluster) %>% top_n(n = 1000, wt = avg_logFC)
write.table(subtypes_marker, file = file.path(OUTPUT, "markers_subtypes.csv"), sep = ',')

subtypes_marker_report <- UMI.subtypes %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.table(subtypes_marker_report, file = file.path(OUTPUT, "subtypes_marker_report.csv"), sep = ',')


                                            #PART B OF PIPELINE - PLOTS SIZE ADJUSTING
#########################################################################################################################################################################################

width <- 20 + (length(unique(Idents(UMI))))/5

svg(file.path(OUTPUT, "PCA_DimPlot_subtypes.svg"), width = width, height = 15)
DimPlot(UMI, reduction = "pca", raster = FALSE)
dev.off()

svg(file.path(OUTPUT, "UMAP_with_DE_gene_subtypes.svg"), width = width, height = 15)
DimPlot(UMI, reduction = "umap", raster = FALSE) 
dev.off()

htmlwidgets::saveWidget(plotly::ggplotly(DimPlot(UMI, reduction = "umap", raster = FALSE)) , file.path(OUTPUT, "UMAP_with_DE_gene_subtypes.html"))


HDMAP <- hdmap_cordinates(UMI, hd_factors)

hd_map_plot <- plotly::ggplotly(DimPlotFactor(HDMAP))

htmlwidgets::saveWidget(hd_map_plot, file.path(OUTPUT, "HDMAP_subtypes.html"))

write.table(HDMAP, file = file.path(OUTPUT, "hdmap_cordinates.csv"), sep = ',')

#Create Expression Matrix

#Expression matrix cells

###########################################################################################################################################################

exp_stat <- heterogenity_stats(UMI)

###########################################################################################################################################################

width <- 25 + (length(unique(Idents(UMI))))/5
height <- 20 + (length(unique(Idents(UMI))))/5

cells <- ggplot(exp_stat, mapping = aes(x = mean_expression, y = positive_expression_perc, fill = names, color = names)) +
  geom_jitter() +
  facet_wrap(names~.) +
  theme(legend.position = 'none')


svg(filename = file.path(OUTPUT,'cells_heterogenity.svg'), width = width, height = height)
cells
dev.off()
rm(cells)

###########################################################################################################################################################


exp_matrix <- GetAssayData(UMI, slot = 'data')
colnames(exp_matrix) <- UMI@active.ident


write(colnames(exp_matrix), file = file.path(OUTPUT, "exp_matrix/barcodes.tsv"))
write(rownames(exp_matrix), file = file.path(OUTPUT, "exp_matrix/genes.tsv"))
Matrix::writeMM(exp_matrix, file = file.path(OUTPUT, "exp_matrix/matrix.mtx"))


###########################################################################################################################################################
#Average expression matrix populations

rm(average_expression)

average_expression <- aggregation_chr(UMI)

write.table(average_expression, file = file.path(OUTPUT, 'exp_matrix/average_expression_matrix_subtypes.csv'), sep = ',')


#save seurat output file

saveRDS(UMI, file = file.path(OUTPUT, "Results.rds"))


###########################################################################################################################################################

#Cell populations pheatmaps

if (length(markers_subclass) != 0) {
  ms<- c()
  for (m in renamed_list$subclass_marker_list) {
    ms <- c(ms, m[grepl(toupper(m), toupper(list(colnames(average_expression))))])
  }
  
  marker_list <- unique(c(marker_list$class_marker_list, marker_list$used_markers, ms))
  
}

if (length(markers_subclass) == 0) {
  
  marker_list <- unique(c(marker_list$class_marker_list, marker_list$used_markers))
  
}

width <- 30 + (length(unique(Idents(UMI))))/5
height <- 28 + (length(unique(Idents(UMI))))/5

average_expression <- average_expression[toupper(rownames(average_expression)) %in% toupper(marker_list),]
average_expression <- as.matrix(drop_na(average_expression))
average_expression <- average_expression[!rowSums(average_expression) == 0,]


pheat <- pheatmap::pheatmap(average_expression, 
                            clustering_method = 'ward.D',
                            angle_col = 270, fontsize_row = 20, fontsize_col = 20)


svg(file.path(OUTPUT, "pheatmap_cells_populations.svg"), width = width, height = height)
pheat
dev.off()
dev.off()
rm(pheat)


###########################################################################################################################################################

print('Report creating')

rmarkdown::render(input = file.path(getwd(), '../../../scripts/report_mix_manual.Rmd'), 
                  output_format = 'html_document', output_dir = OUTPUT, knit_root_dir = getwd(),
                  output_file = 'Report')

