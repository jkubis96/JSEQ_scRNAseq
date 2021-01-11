library(Seurat)
library(patchwork)
library(tidyverse)

options(scipen = 999)
memory.limit(30000000)
memory.size(T)

args <- commandArgs()


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

  markers_class <- readxl::read_xlsx(markers, sheet = 1)
  markers_subclass <- readxl::read_xlsx(markers, sheet = 2, col_names = F)

##########################################################


{ 
# Load the raw dataset by UMI
UMI_raw <- Read10X(seurat_umi, gene.column = 1)

#Create SeuratObject
UMI <- CreateSeuratObject(counts = UMI_raw, project = project_name, min.cells = 1, min.features = 1)
Idents(UMI) <- make.unique(as.character(head(Idents(UMI))))
}


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

n_gen <- as.numeric(quantile(sort(UMI@meta.data$nGenes, decreasing = F), 0.99))
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

UMI <- subset(UMI, subset = nGenes > 200 & nGenes <= n_gen & MitoPercent < 5)
n_gen <- as.numeric(quantile(sort(UMI@meta.data$nGenes, decreasing = F), 0.95))*0.25
cells_number <- as.numeric(length(UMI@active.ident))

#####################################################################################

UMI <- NormalizeData(UMI, normalization.method = "LogNormalize", scale.factor = 1e6)


######################################################################################

#For normalized expression matrix input
decision <- mean(UMI@assays$RNA@counts@x) < 15
if (data == 2 & decision == T) {
  UMI@assays$RNA@data@x <- UMI@assays$RNA@counts@x
} else if (data == 2 & decision == F) {
  UMI@assays$RNA@data@x <- log(UMI@assays$RNA@counts@x+1)
}

UMI <- FindVariableFeatures(UMI, selection.method = "vst", nfeatures = n_gen)

# Identify the 10 most highly variable genes

top20 <- head(VariableFeatures(UMI), 20)
top50 <- head(VariableFeatures(UMI), 50)
top100 <- head(VariableFeatures(UMI), 100)

plot1 <- VariableFeaturePlot(UMI)

jpeg(file.path(OUTPUT, "variable_genes.jpeg") , units="in", width=10, height=7, res=300)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
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
  dim <- 1
  score <- c()
  element <- 0
  for (i in dims$`Elbow$data$stdev`) {
    element <- element + 1
    if (i-i*0.01 > dims$`Elbow$data$stdev`[element+1] & element < 20 | i-i*0.03 > dims$`Elbow$data$stdev`[element+2] & element < 20) {
      dim <- dim + 1
    } else 
      break
  }
  dim <- as.numeric(dim)
}


#########################################################################################


jpeg(file.path(OUTPUT, "JackStrawPlot.jpeg") , units="in", width=10, height=7, res=300)
JackStrawPlot(UMI, dims = 1:(dim))
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
# UMI.markers <- FindAllMarkers(UMI, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# top10 <- UMI.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# 
# 
# 
# ##Cells cluster naming with top genes (different between cell groups)
# 
# {
#   top10 <- data.frame(top10$cluster, top10$gene)
#   colnames(top10) <- c('cluster','gene')
#   
#   
#   all_unique('top10','[1:2]','cluster')
#   
#   top10$cluster <- top10$Group.1
#   top10 <- top10[,-1] 
# }

##########################################################
### Most variable genes select and cell subtypes nameing ###

`%!in%` = Negate(`%in%`)
genes <- markers_subclass[1,]
exp_matrix <- as.data.frame(GetAssayData(UMI, slot = 'data'))
colnames(exp_matrix) <- UMI@active.ident
exp_matrix <- as.data.frame(exp_matrix)
exp_matrix <- exp_matrix[toupper(rownames(exp_matrix)) %!in% genes[1,],]
exp_matrix <- exp_matrix[rownames(exp_matrix) %in% top20,]


num <- as.numeric(length(genes))
index = 0
for (i in colnames(exp_matrix)) {
  index = index +1
  tmp <- as.data.frame(exp_matrix[order(exp_matrix[,index], decreasing = T), ,drop =F])
  colnames(exp_matrix)[index] <- rownames(tmp[1,])
}

cell_names <- colnames(exp_matrix)
rm(exp_matrix)


############################################
#PCA Cluster average expression for nameing

average_expression <- as.data.frame(GetAssayData(UMI, slot = 'data'))
colnames(average_expression) <- UMI@active.ident
average_expression <- as.data.frame(average_expression)
average_expression <- t(average_expression)
average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = mean)
rownames(average_expression) <- average_expression[,1]
average_expression <- average_expression[,-1]
average_expression <- t(average_expression)

cluster_nameing(matrix_a = average_expression, markers = markers_class)

clust_names <- colnames(average_expression)


##########################################


genes <- markers_subclass[1,]

num <- as.numeric(length(genes))
index = 0
for (i in colnames(average_expression)) {
  index = index +1
  rename_df <- c()
  rename_df <- average_expression[toupper(rownames(average_expression)) %in% genes[1,],]
  rename_df <- as.data.frame(rename_df[order(rename_df[,index], decreasing = T), ,drop = F])
  colnames(average_expression)[index] <- rownames(rename_df[1,])
}

new.cluster.ids <- paste(clust_names, colnames(average_expression))
names(new.cluster.ids) <- levels(UMI)
UMI <- RenameIdents(UMI, new.cluster.ids)



exp_matrix <- as.data.frame(GetAssayData(UMI, slot = 'data'))
colnames(exp_matrix) <- UMI@active.ident
exp_matrix <- as.data.frame(exp_matrix)
colnames(exp_matrix) <- paste(UMI@active.ident, cell_names)


#PCA plot and UMAP plot with names

#Class

jpeg(file.path(OUTPUT, "PCA_DimPlot_class.jpeg") , units="in", width=10, height=7, res=300)
DimPlot(UMI, reduction = "pca")
dev.off()

jpeg(file.path(OUTPUT, "UMAP_with_DE_gene_class.jpeg") , units="in", width=10, height=7, res=300)
DimPlot(UMI, reduction = "umap") 
dev.off()


#Subtypes

Idents(UMI) <- colnames(exp_matrix)

jpeg(file.path(OUTPUT, "PCA_DimPlot_subtypes.jpeg") , units="in", width=20, height=10, res=300)
DimPlot(UMI, reduction = "pca")
dev.off()

jpeg(file.path(OUTPUT, "UMAP_with_DE_gene_subtypes.jpeg") , units="in", width=20, height=10, res=300)
DimPlot(UMI, reduction = "umap") 
dev.off()


#Create Expression Matrix

#Expression matrix cells

write.table(exp_matrix, file = file.path(OUTPUT, "expression_matrix_cells.csv"), sep = ',')

#Average expression matrix populations

average_expression <- as.data.frame(exp_matrix)
average_expression <- t(average_expression)
average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = mean)
rownames(average_expression) <- average_expression[,1]
average_expression <- average_expression[,-1]
average_expression <- t(average_expression)


write.table(average_expression, file = file.path(OUTPUT, "average_expression_matrix_populations.csv"), sep = ',')


#save seurat output file

saveRDS(UMI, file = file.path(OUTPUT, "Results.rds"))


#######################################################################################

#Cell populations pheatmaps for all genes

average_expression <- as.data.frame(average_expression)
phet_exp_matrix <- average_expression[top50,]
phet_exp_matrix <- drop_na(phet_exp_matrix)
phet_exp_matrix <- scale(phet_exp_matrix)
jpeg(file.path(OUTPUT, "pheatmap_cells_populations.jpeg") , units="in", width=15, height = 10, res=300)
pheatmap::pheatmap(phet_exp_matrix, 
                   clustering_method = 'complete')
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

