library(Seurat)
library(patchwork)
library(tidyverse)



args <- commandArgs()

species <- FALSE


#Paths and arguments from env
{
print(args)

path <- args[6]
markers <-args[7]
species <- args[8]
seurat_umi <- file.path(path,'sc_data/')
OUTPUT <- file.path(path, 'results')
project_name <- args[9]
data <- args[10]
functions <- file.path(getwd(), 'scripts/functions.R')
source(functions)
}

#########################################################

if (species == 'human') {
  markers <- readxl::read_xlsx(markers, sheet = 1)
} else if (species == 'mice') {
  markers <- readxl::read_xlsx(markers, sheet = 2)
}

##########################################################


{ 
# Load the raw dataset by UMI
UMI_raw <- Read10X(seurat_umi, gene.column = 1)

#Create SeuratObject
UMI <- CreateSeuratObject(counts = UMI_raw, project = project_name, min.cells = 1, min.features = 1)
UMI
}

###########################ZMIENIÆ###############################
cells_number <- length(UMI@active.ident)
cells_number <- as.numeric(cells_number)
#################################################################

#Create Ribo and Mito percent stats

UMI[['MitoPercent']] <- PercentageFeatureSet(UMI, pattern = "^MT-")
UMI[['RiboPercent']] <- PercentageFeatureSet(UMI, pattern = "^Rps-") + PercentageFeatureSet(UMI, pattern = "^Rpl")

UMI@meta.data <- UMI@meta.data %>% 
  rename(nCounts = nCount_RNA) %>% 
  rename(nGenes = nFeature_RNA)

saveRDS(UMI, file = file.path(OUTPUT, "Seurat_object.rds"))

#Graphs of counts content

jpeg(file.path(OUTPUT, "counts~genes.jpeg") , units="in", width=15, height=10, res=300)
UC_plot <- VlnPlot(UMI, features = c("nGenes", "nCounts"), ncol = 2)
UC_plot
dev.off()

jpeg(file.path(OUTPUT, "Ribo~Mito.jpeg") , units="in", width=15, height=10, res=300)
MR_plot <- VlnPlot(UMI, features = c("RiboPercent", "MitoPercent"), ncol = 2)
MR_plot
dev.off()


jpeg(file.path(OUTPUT, "counts~genes_QC.jpeg") , units="in", width=15, height=10, res=300)
CG_plot <- FeatureScatter(UMI, feature1 = "nCounts", feature2 = "nGenes")
CG_plot
dev.off()

####################################################################################

#Droplet content and QC

n_gen <- as.numeric(quantile(UMI@meta.data$nGenes, 0.75))

QC_UMI <- data.frame()
QC_UMI <- as.data.frame(UMI$nGenes)
QC_UMI$V2 <- UMI$MitoPercent
QC_UMI$V3 <- UMI$RiboPercent

colnames(QC_UMI) <- c('nGenes','MitoPercent','RiboPercent')

QC_UMI$Mito_Status[QC_UMI$MitoPercent > 5] <- '> 5%'
QC_UMI$Mito_Status[QC_UMI$MitoPercent <= 5] <- 'Proper'

QC_UMI$nGenes_Status[UMI$nGenes < 200] <- 'Empty'
QC_UMI$nGenes_Status[UMI$nGenes > n_gen] <- 'Double'
QC_UMI$nGenes_Status[UMI$nGenes >= 200 & UMI$nGenes <= n_gen] <- 'Proper'

QC_UMI$Ribo_Status[QC_UMI$RiboPercent == 0] <- '0%'
QC_UMI$Ribo_Status[QC_UMI$RiboPercent > 0] <- '> 0 %'



DQC <- ggplot()+
  geom_point(QC_UMI, mapping = aes(x = nGenes, y = nGenes, color = nGenes_Status))+
  ylab("Number of genes for each cells") +
  xlab("Number of genes for each cells")+
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))+
  labs(color='Droplet content') 

ggsave(DQC, filename = file.path(OUTPUT,'DropletQC.jpeg'), width = 10, height = 7, dpi = 600)
rm(DQC)

MQC <- ggplot()+
  geom_point(QC_UMI, mapping = aes(x = MitoPercent, y = MitoPercent , color = Mito_Status))+
  ylab("% MitoRNA") +
  xlab("% MitoRNA")+
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))+
  labs(color='% Content treshold') 

ggsave(MQC, filename = file.path(OUTPUT,'MitoQC.jpeg'), width = 10, height = 7, dpi = 600)
rm(MQC)

RQC <- ggplot()+
  geom_point(QC_UMI, mapping = aes(x = RiboPercent, y = RiboPercent, color = Ribo_Status))+
  ylab("% RiboRNA") +
  xlab("% RiboRNA")+
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))+
  labs(color='% Content treshold') 

ggsave(RQC, filename = file.path(OUTPUT,'RiboQC.jpeg'),  width = 10, height = 7, dpi = 600)
rm(RQC)
rm(QC_UMI)

####################################################################################

#Selecting right cells

UMI <- subset(UMI, subset = nGenes > 200 & nGenes < n_gen & MitoPercent < 5)

#####################################################################################

UMI <- NormalizeData(UMI, normalization.method = "LogNormalize", scale.factor = 1e6)


######################################################################################

#For normalized expression matrix input

if (data == 2) {
  UMI@assays$RNA@data <- UMI@assays$RNA@counts
}

UMI <- FindVariableFeatures(UMI, selection.method = "vst", nfeatures = n_gen)

# Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(UMI), 10)

plot1 <- VariableFeaturePlot(UMI)

jpeg(file.path(OUTPUT, "variable_genes.jpeg") , units="in", width=10, height=7, res=300)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
dev.off()

#####################################################################################

all.genes <- rownames(UMI)
UMI <- ScaleData(UMI, features = all.genes)

UMI <- RunPCA(UMI, features = VariableFeatures(object = UMI))



UMI <- JackStraw(UMI, num.replicate = 100)
UMI <- ScoreJackStraw(UMI, dims = 1:20)



################################

jpeg(file.path(OUTPUT, "Elbow.jpeg") , units="in", width=10, height=7, res=300)
Elbow <- ElbowPlot(UMI)
Elbow
dev.off()

dims <- as.data.frame(Elbow$data$stdev)

#select the most variable reduction

{
  dim <- c()
  score <- c()
  element <- 0
  for (i in dims$`Elbow$data$stdev`) {
    element <- element + 1
    score <- cbind(i)
    print(dims$`Elbow$data$stdev`[element+1])
    print(i)
    if (i-0.05*i >= dims$`Elbow$data$stdev`[element+1]) {
      dim <- element
    } else break
  }
  
  
  dim <- as.numeric(element) + 1
}

#########################################################################################


jpeg(file.path(OUTPUT, "JackStrawPlot.jpeg") , units="in", width=10, height=7, res=300)
JackStrawPlot(UMI, dims = 1:(dim+1))
dev.off()

jpeg(file.path(OUTPUT, "PCA_heatmap.jpeg") , units="in", width=10, height=7, res=300)
DimHeatmap(UMI, dims = 1:dim, cells = cells_number, balanced = TRUE)
dev.off()


UMI <- FindNeighbors(UMI, dims = 1:dim)
UMI <- FindClusters(UMI, resolution = 0.5)


UMI <- RunUMAP(UMI, dims = 1:dim)


jpeg(file.path(OUTPUT, "UMAP.jpeg") , units="in", width=10, height=7, res=300)
DimPlot(UMI, reduction = "umap")
dev.off()

####################################################################################



# find markers for every cluster compared to all remaining cells, report only the positive ones
UMI.markers <- FindAllMarkers(UMI, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- UMI.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)



##Cells cluster naming with top genes (different between cell groups)

{
  top10 <- data.frame(top10$cluster, top10$gene)
  colnames(top10) <- c('cluster','gene')
  
  
  all_unique('top10','[1:2]','cluster')
  
  top10$cluster <- top10$Group.1
  top10 <- top10[,-1] 
}

##########################################################
### LogFC genes select and cell nameing ###


exp_matrix <- as.data.frame(GetAssayData(UMI, slot = 'data'))
colnames(exp_matrix) <- UMI@active.ident
exp_matrix <- as.data.frame(exp_matrix)


lvl_num <- as.numeric(max(levels(UMI))) +1 
{
  for (lvl in 1:lvl_num) {
    for (i in 1:ncol(exp_matrix)) {
      rename_df <- c()
      col <- lvl -1
      rename_df <- exp_matrix[rownames(exp_matrix) %in% top10$gene[[lvl]],]
      if (colnames(rename_df[i]) == as.numeric(col)) {
        rename_df <- rename_df[order(rename_df[,i], decreasing = T),]
        rename_df <- rownames(rename_df[1:2,])
        rename_df <- sort(rename_df)
        colnames(exp_matrix)[i] <- paste(rename_df[1], rename_df[2])
      }
    }
  } 
  
}

cell_names <- colnames(exp_matrix)
cell_nameing(matrix_a = exp_matrix, markers = markers)
colnames(exp_matrix) <- paste(colnames(exp_matrix), cell_names)



############################################
#PCA Cluster average expression for nameing

average_expression <- as.data.frame(GetAssayData(UMI, slot = 'data'))
colnames(average_expression) <- UMI@active.ident
average_expression <- as.data.frame(average_expression)
average_expression <- t(average_expression)
average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = mean, trim = 0.25)
rownames(average_expression) <- average_expression[,1]
average_expression <- average_expression[,-1]
average_expression <- t(average_expression)

cluster_nameing(matrix_a = average_expression, markers = markers)

new.cluster.ids <- colnames(average_expression)
names(new.cluster.ids) <- levels(UMI)
UMI <- RenameIdents(UMI, new.cluster.ids)

rm(average_expression)


#PCA plot and UMAP plot with names

jpeg(file.path(OUTPUT, "PCA_DimPlot.jpeg") , units="in", width=10, height=7, res=300)
DimPlot(UMI, reduction = "pca")
dev.off()

jpeg(file.path(OUTPUT, "UMAP_with_DE_gene.jpeg") , units="in", width=10, height=7, res=300)
DimPlot(UMI, reduction = "umap") 
dev.off()


#Cluster tree for show connection between groups

UMI <- BuildClusterTree(object = UMI, slot = 'data', verbose = T)

tree_cluster <- Tool(object = UMI, slot = 'BuildClusterTree')


jpeg(file.path(OUTPUT, "hierarchical_tree.jpeg") , units="in", width=10, height=7, res=300)
plot(tree_cluster)
dev.off()

#Create Expression Matrix

#Expression matrix cells


write.table(exp_matrix, file = file.path(OUTPUT, "expression_matrix_cells.csv"))

#Average expression matrix populations

average_expression <- as.data.frame(exp_matrix)
rm(exp_matrix)
average_expression <- t(average_expression)
average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = mean, trim = 0.25)
rownames(average_expression) <- average_expression[,1]
average_expression <- average_expression[,-1]
average_expression <- t(average_expression)



write.table(average_expression, file = file.path(OUTPUT, "average_expression_matrix_populations.csv"))


#save seurat output file

saveRDS(UMI, file = file.path(OUTPUT, "Results.rds"))


#######################################################################################

top_genes <- top10$gene[[1]]
index <- 1
for (i in 2:nrow(top10)){
  index <- index + 1
  top_genes <- c(top_genes, top10$gene[[index]])
}


average_expression <- as.data.frame(average_expression)
phet_exp_matrix <- average_expression[rownames(average_expression) %in% top_genes,]
jpeg(file.path(OUTPUT, "pheatmap_cells_populations.jpeg") , units="in", width=15, height=7, res=300)
pheatmap::pheatmap(phet_exp_matrix, 
                          clustering_method = 'complete', 
                          cluster_rows = F)
dev.off()


#########################################################################################

if (species %in% c('human','mice')) {
  rmarkdown::render(input = file.path(getwd(), 'scripts/raport_species.Rmd'), 
                    output_format = 'html_document', output_dir = OUTPUT, 
                    output_file = 'Raport')
} else if (species == 'mix') {
  rmarkdown::render(input = file.path(getwd(), 'scripts/raport_mix.Rmd'), 
                    output_format = 'html_document', output_dir = OUTPUT, 
                    output_file = 'Raport')
} else {
  quit()
  n
}

