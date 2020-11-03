library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)



args <- commandArgs()

species <- NULL

print(args)

path_tmp <- args[6]
path_results <- args[7]

cells_number <- args[8]
cells_number <- as.numeric(cells_number)


UMI_PATH <- file.path(path_tmp,'seurat_umi/')
OUTPUT <- file.path(path_results)

project_name_mode <- args[9]

species <- args[10]

scripts_path <- args[11]

#for mix species analysis >
path_mutual_results <- args[12]
# <

{ 
# Load the raw dataset by UMI
UMI_raw <- Read10X(UMI_PATH, gene.column = 1)


# Initialize the Seurat object for UMI
UMI <- CreateSeuratObject(counts = UMI_raw, project = project_name_mode, min.cells = 1, min.features = 1)
UMI
}


UMI[['MitoPercent']] <- PercentageFeatureSet(UMI, pattern = "^MT-")
UMI[['RiboPercent']] <- PercentageFeatureSet(UMI, pattern = "^Rps-") + PercentageFeatureSet(UMI, pattern = "^Rpl")

UMI@meta.data <- UMI@meta.data %>% 
  rename(nCounts = nCount_RNA) %>% 
  rename(nGenes = nFeature_RNA)

saveRDS(UMI, file = file.path(OUTPUT, "Seurat_object_input.rds"))

jpeg(file.path(OUTPUT, "counts~genes.jpeg") , units="in", width=15, height=10, res=300)
UC_plot <- VlnPlot(UMI, features = c("nGenes", "nCounts"), ncol = 2)
UC_plot
dev.off()

jpeg(file.path(OUTPUT, "Ribo~Mito.jpeg") , units="in", width=15, height=10, res=300)
MR_plot <- VlnPlot(UMI, features = c("RiboPercent", "MitoPercent"), ncol = 2)
MR_plot
dev.off()

###############################################################################

jpeg(file.path(OUTPUT, "counts~genes_QC.jpeg") , units="in", width=15, height=10, res=300)
CG_plot <- FeatureScatter(UMI, feature1 = "nCounts", feature2 = "nGenes")
CG_plot
dev.off()

####################################################################################

#Droplet contain QC



UMI_selected <- subset(UMI, subset = nGenes > 200 & nGenes < 2500 & MitoPercent < 5)

QC_UMI <- data.frame()
QC_UMI <- as.data.frame(UMI$nGenes)
QC_UMI$V2 <- UMI$MitoPercent
QC_UMI$V3 <- UMI$RiboPercent

colnames(QC_UMI) <- c('nGenes','MitoPercent','RiboPercent')

QC_UMI$Mito_Status[QC_UMI$MitoPercent > 5] <- '> 5%'
QC_UMI$Mito_Status[QC_UMI$MitoPercent <= 5] <- 'Proper'

QC_UMI$nGenes_Status[UMI$nGenes < 200] <- 'Empty'
QC_UMI$nGenes_Status[UMI$nGenes > 2500] <- 'Double'
QC_UMI$nGenes_Status[UMI$nGenes >= 200 & UMI$nGenes <= 2500] <- 'Proper'

QC_UMI$Ribo_Status[QC_UMI$RiboPercent == 0] <- '0%'
QC_UMI$Ribo_Status[QC_UMI$RiboPercent > 0] <- '> 0 %'



DQC <- ggplot()+
  geom_point(QC_UMI, mapping = aes(x = nGenes, y = nGenes, color = nGenes_Status))+
  ylab("Number of genes for each cells") +
  xlab("Number of genes for each cells")+
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))+
  labs(color='Droplet content') 

ggsave(DQC, filename = file.path(OUTPUT,'DropletQC.jpeg'), width = 10, height = 7, dpi = 600)

MQC <- ggplot()+
  geom_point(QC_UMI, mapping = aes(x = MitoPercent, y = MitoPercent , color = Mito_Status))+
  ylab("% MitoRNA") +
  xlab("% MitoRNA")+
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))+
  labs(color='% Content treshold') 

ggsave(MQC, filename = file.path(OUTPUT,'MitoQC.jpeg'), width = 10, height = 7, dpi = 600)

RQC <- ggplot()+
  geom_point(QC_UMI, mapping = aes(x = RiboPercent, y = RiboPercent, color = Ribo_Status))+
  ylab("% RiboRNA") +
  xlab("% RiboRNA")+
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))+
  labs(color='% Content treshold') 

ggsave(RQC, filename = file.path(OUTPUT,'RiboQC.jpeg'),  width = 10, height = 7, dpi = 600)

####################################################################################
UMI <- NormalizeData(UMI, normalization.method = "LogNormalize", scale.factor = 10000)


UMI <- FindVariableFeatures(UMI, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(UMI), 10)

plot1 <- VariableFeaturePlot(UMI)

jpeg(file.path(OUTPUT, "variable_genes.jpeg") , units="in", width=10, height=7, res=300)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
dev.off()



##########################################################################
all.genes <- rownames(UMI)
UMI <- ScaleData(UMI, features = all.genes)

UMI <- RunPCA(UMI, features = VariableFeatures(object = UMI))

print(UMI[["pca"]], dims = 1:5, nfeatures = 5)


jpeg(file.path(OUTPUT, "PCA_DimPlot.jpeg") , units="in", width=10, height=7, res=300)
DimPlot(UMI, reduction = "pca")
dev.off()

jpeg(file.path(OUTPUT, "PCA_heatmap.jpeg") , units="in", width=10, height=7, res=300)
DimHeatmap(UMI, dims = 1:15, cells = cells_number, balanced = TRUE)
dev.off()



UMI <- JackStraw(UMI, num.replicate = 100)
UMI <- ScoreJackStraw(UMI, dims = 1:20)

jpeg(file.path(OUTPUT, "JackStrawPlot.jpeg") , units="in", width=10, height=7, res=300)
JackStrawPlot(UMI, dims = 1:15)
dev.off()

jpeg(file.path(OUTPUT, "Elbow.jpeg") , units="in", width=10, height=7, res=300)
ElbowPlot(UMI)
dev.off()

UMI <- FindNeighbors(UMI, dims = 1:10)
UMI <- FindClusters(UMI, resolution = 0.5)


head(Idents(UMI), 5)


UMI <- RunUMAP(UMI, dims = 1:10)

jpeg(file.path(OUTPUT, "UMAP.jpeg") , units="in", width=10, height=7, res=300)
DimPlot(UMI, reduction = "umap")
dev.off()


# find all markers of cluster 1
cluster1.markers <- FindMarkers(UMI, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)



# find markers for every cluster compared to all remaining cells, report only the positive ones
UMI.markers <- FindAllMarkers(UMI, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
new_claster <- UMI.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)


top5 <- data.frame(new_claster$cluster, new_claster$gene)
colnames(top5) <- c('cluster','gene')


#top5$gene <- c("GAD1"  ,   "LGALS1"   ,  "S100A4"    , "ACTBP2"  ,   "AP002982.1", "C1QBP","HMGN1"  ,    "CKB"     ,   "NOL7"    ,   "GAL"     ,   "RPL26"    ,  "TXNDC17"  , "PWP1"  ,     "MT-ND2"   ,  "AURKA"    )


##Paterning via markers
{

all_unique <- function(data, range, value_check_1){
  eval(parse(text= paste0(data,'<<- aggregate(',data,range,' ,list(',
                          data,'$',value_check_1,'),FUN=list)')))

}

all_unique('top5','[1:2]','cluster')

top5$cluster <- top5$Group.1
top5 <- top5[,-1]

cell <- c()

for (i in top5$gene) {
  if ('GAD1' %in% i) {
    cell <<- c(cell, 'Inh')
  }else if ('SLC17A7' %in% i) {
    cell <<- c(cell, 'Exc')
  }else if ('SLC1A3' %in% i) {
    cell <<- c(cell, 'Non-neuronal')
  } else 
    cell <<- c(cell, 'Non')
}    

dv <- c()

for (i in top5$gene) {
  if (c('EMX1', 'TBR1', 'TBR2') %in% i) {
    dv <<- c(dv, 'Dorsal')
  }else if ('SLC17A7' %in% i) {
    dv <<- c(dv, 'Ventral')
  } else 
    dv <<- c(dv, '')
}

top5$cell <- cell

 for (i in 1:length(top5$cluster)) {
   top5$cell[i] <- paste(top5$cell[i], top5$gene[[i]][1], top5$gene[[i]][2])
}


top5$cell <- paste(top5$cell, dv)


}


new.cluster.ids <- top5$cell

names(new.cluster.ids) <- levels(UMI)
UMI <- RenameIdents(UMI, new.cluster.ids)

jpeg(file.path(OUTPUT, "UMAP_with_DE_gene.jpeg") , units="in", width=10, height=7, res=300)
DimPlot(UMI, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()


saveRDS(UMI, file = file.path(OUTPUT, "Seurat_object_output.rds"))

#Create Expression Matrix

saveRDS(UMI, file = file.path(OUTPUT, "Seurat_object_output.rds"))
AE <- AverageExpression(UMI)
AE <- as.matrix(AE[[1]])

write.csv(AE, file = file.path(OUTPUT, "expression_matrix.csv"))

#########################################################################################

if (species %in% c('human','mice')) {
  rmarkdown::render(input = file.path(scripts_path, 'raport_species.Rmd'), 
                    output_format = 'html_document', output_dir = OUTPUT, 
                    output_file = 'Raport')
} else if (species == 'mix') {
  rmarkdown::render(input = file.path(scripts_path, 'raport_mix.Rmd'), 
                    output_format = 'html_document', output_dir = path_mutual_results, 
                    output_file = 'Raport')
} else {
  quit()
  n
}

