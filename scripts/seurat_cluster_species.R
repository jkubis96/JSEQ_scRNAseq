library(Seurat)
library(patchwork)
library(tidyverse)
library(doSNOW)
library(foreach)
library(doParallel)

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
  estimated_cells <- args[11]
  functions <- file.path(getwd(), 'scripts/functions.R')
  source(functions, local = T)
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
  
  cell_input <- length(Idents(UMI))
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

jpeg(file.path(OUTPUT, "counts~genes.jpeg") , units="in", width=15, height=10, res=600)
UC_plot <- VlnPlot(UMI, features = c("nGenes", "nCounts"), ncol = 2)
UC_plot
dev.off()

jpeg(file.path(OUTPUT, "Ribo~Mito.jpeg") , units="in", width=15, height=10, res=600)
MR_plot <- VlnPlot(UMI, features = c("RiboPercent", "MitoPercent"), ncol = 2)
MR_plot
dev.off()


jpeg(file.path(OUTPUT, "counts~genes_QC.jpeg") , units="in", width=15, height=10, res=600)
CG_plot <- FeatureScatter(UMI, feature1 = "nCounts", feature2 = "nGenes")
CG_plot
dev.off()

####################################################################################

#Droplet content and QC
n_gen <- UMI@meta.data$nGenes[UMI@meta.data$nGenes > 500]
n_gen <- as.numeric(mean(n_gen))*2 + 1.5*IQR(as.numeric(UMI@meta.data$nGenes))

QC_UMI <- data.frame()
QC_UMI <- as.data.frame(UMI$nGenes)
QC_UMI$V2 <- UMI$MitoPercent
QC_UMI$V3 <- UMI$RiboPercent

colnames(QC_UMI) <- c('nGenes','MitoPercent','RiboPercent')

QC_UMI$Mito_Status[QC_UMI$MitoPercent > 5] <- '> 5%'
QC_UMI$Mito_Status[QC_UMI$MitoPercent <= 5] <- 'Proper'

QC_UMI$nGenes_Status[UMI$nGenes < 500] <- 'Empty'
QC_UMI$nGenes_Status[UMI$nGenes > n_gen] <- 'Double'
QC_UMI$nGenes_Status[UMI$nGenes >= 500 & UMI$nGenes <= n_gen] <- 'Proper'

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
  labs(color='% Content threshold') 

ggsave(MQC, filename = file.path(OUTPUT,'MitoQC.jpeg'), width = 10, height = 7, dpi = 600)
rm(MQC)

RQC <- ggplot()+
  geom_point(QC_UMI, mapping = aes(x = RiboPercent, y = RiboPercent, color = Ribo_Status))+
  ylab("% RiboRNA") +
  xlab("% RiboRNA")+
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))+
  labs(color='% Content threshold') 

ggsave(RQC, filename = file.path(OUTPUT,'RiboQC.jpeg'),  width = 10, height = 7, dpi = 600)
rm(RQC)
rm(QC_UMI)

####################################################################################

#Selecting right cells

UMI <- subset(UMI, subset = nGenes > 500 & nGenes <= n_gen & MitoPercent < 5)
n_gen <- max(as.numeric(UMI@meta.data$nGenes))*0.75
cells_number <- length(Idents(UMI))

######################################################################################
#Cells_stats

cells <- factor(c('Estimated_cells', 'Input_cells', 'Analyzed_cells'), levels = c('Estimated_cells', 'Input_cells', 'Analyzed_cells'))
cell_num <- c(as.numeric(estimated_cells), as.numeric(cell_input), as.numeric(cells_number))
df_cells <- data.frame(cells, cell_num)

cells <- ggplot(df_cells, aes(x = cells, y = cell_num, fill = cells)) +
  geom_col() +
  ylab("Number of cells") +
  xlab("Cells in analysis")+
  geom_text(aes(label = cell_num), vjust = -0.5)+
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))+
  theme_classic() +
  theme(legend.position = 'none') 


ggsave(cells, filename = file.path(OUTPUT,'Cells.jpeg'), width = 10, height = 7, dpi = 600)
rm(cells)


#####################################################################################

UMI <- NormalizeData(UMI, normalization.method = "LogNormalize", scale.factor = 1e6)


######################################################################################

#For normalized expression matrix input

decision <- mean(as.numeric(UMI@assays$RNA@counts@x)) < 15
if (data == 2 & decision == T) {
  UMI@assays$RNA@data@x <- UMI@assays$RNA@counts@x
} else if (data == 2 & decision == F) {
  UMI@assays$RNA@data@x <- log(UMI@assays$RNA@counts@x+1)
}

#######################################################################################

UMI <- FindVariableFeatures(UMI, selection.method = "vst", nfeatures = n_gen, binning.method = 'equal_frequency')

# Identify the 10 most highly variable genes

top20 <- head(VariableFeatures(UMI), 20)

plot1 <- VariableFeaturePlot(UMI)

jpeg(file.path(OUTPUT, "variable_genes.jpeg") , units="in", width=10, height=7, res=600)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2
dev.off()

#####################################################################################

all.genes <- rownames(UMI)
UMI <- ScaleData(UMI, features = all.genes)

UMI <- RunPCA(UMI, features = VariableFeatures(object = UMI))



################################

jpeg(file.path(OUTPUT, "Elbow.jpeg") , units="in", width=10, height=7, res=600)
Elbow <- ElbowPlot(UMI, ndims = 50)
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
    if (i-i*0.01 > dims$`Elbow$data$stdev`[element+1] & element < 50 | i-i*0.02 > dims$`Elbow$data$stdev`[element+2] & element < 49 | i-i*0.02 > dims$`Elbow$data$stdev`[element+3] & element < 48 | i-i*0.02 > dims$`Elbow$data$stdev`[element+4] & element < 47) {
      dim <- dim + 1
    } else 
      break
  }
  dim <- as.numeric(dim)
}


#########################################################################################

UMI <- JackStraw(UMI, num.replicate = 10, dims = dim)
UMI <- ScoreJackStraw(UMI, dims = 1:dim)


jpeg(file.path(OUTPUT, "JackStrawPlot.jpeg") , units="in", width=10, height=7, res=600)
JackStrawPlot(UMI, dims = 1:(dim))
dev.off()


UMI <- FindNeighbors(UMI, dims = 1:dim, k.param = 49, reduction = 'pca')
UMI <- FindClusters(UMI, resolution = 0.5, n.start = 10, n.iter = 1000)


UMI <- RunUMAP(UMI, dims = 1:dim, n.neighbors = 49)


jpeg(file.path(OUTPUT, "UMAP.jpeg") , units="in", width=10, height=7, res=600)
DimPlot(UMI, reduction = "umap")
dev.off()

####################################################################################



#find markers for every cluster compared to all remaining cells, report only the positive ones

print('Searching for cluster marker genes')


UMI.markers <- FindAllMarkers(UMI, only.pos = TRUE, min.pct = 0.20, logfc.threshold = 0.25, test.use = 'MAST')

if (sum(as.numeric(levels(UMI))) != sum(unique(as.integer(UMI.markers$cluster)-1))) {
  UMI.markers <- FindAllMarkers(UMI, only.pos = TRUE, min.pct = 0.20, logfc.threshold = 0.10, test.use = 'MAST')
  }

top10 <- UMI.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top100 <- UMI.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)


MAST_markers <- UMI.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

write.table(MAST_markers, file = file.path(OUTPUT, "MAST_markers_clusters.csv"), sep = ',')

print('Cluster genes - DONE')

##Cells cluster naming with top genes (different between cell groups)

##########################################################
### Most variable genes select and cell subtypes nameing 


##################################################################################
#Cells subtypes selection

print('Subsetting clusters')

iterations <- max(as.numeric(unique(UMI@active.ident)))
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


CPU <- detectCores() - 1

if (CPU > max(as.numeric(unique(UMI@active.ident)))) {
  cl <- makeCluster(max(as.numeric(unique(UMI@active.ident))))
} else if (CPU < max(as.numeric(unique(UMI@active.ident)))) {
  cl <- makeCluster(CPU)
}  

registerDoParallel(cl)
registerDoSNOW(cl)

pca_cluster_genes <- list()
clusters <- 4
tmp <- GetAssayData(UMI, slot = 'data')
colnames(tmp) <- UMI@active.ident

pca_cluster_genes <- foreach(pca_cluster = 1:max(as.numeric(unique(UMI@active.ident))), .packages =c('Seurat', 'patchwork','tidyverse'), .options.snow = opts) %dopar% {
  
  pca_cluster <- pca_cluster - 1
  
  tmp2 <- as.data.frame(tmp[, colnames(tmp) %in% pca_cluster])
  tmp2[tmp2 > 1] <- 1L
  tmp2[tmp2 < 1] <- 0L
  
  gen_cor <- unique(top100$gene[top100$cluster %in% pca_cluster])
  tmp3 <- tmp2[gen_cor, ]
  rm(tmp2)
  
  df_groups <- data.frame(1:length(colnames(tmp3)))
  df_groups <- as.data.frame(t(df_groups))
  df_groups_tmp <- data.frame(1:length(colnames(tmp3)))
  df_groups_tmp <- as.data.frame(t(df_groups_tmp))
  selec_group <- data.frame()
  df_length <- length(rownames(tmp3))
  br <- F
  for (cluster in 1:clusters) {
    if(br == T){break}
    for (i in 1:length(rownames(tmp3))) {
      if (cluster == 1) {
        tmp_elements <- tmp3[i,]
        gen_set <- rownames(tmp3)[i]
        perc_1 <- length(tmp_elements[tmp_elements == 1])/length(tmp_elements)
        perc_0 <- length(tmp_elements[tmp_elements == 0])/length(tmp_elements)
        gene_combination <- gen_set
        if (is.nan(perc_1)) {perc_1 <- 0}
        if (is.nan(perc_0)) {perc_2 <- 0}
        if (i == 1) {selec_group <- data.frame(gene_combination, perc_1, perc_0)
        } else if (i > 1) {selec_group_tmp <- data.frame(gene_combination, perc_1, perc_0)
        selec_group <- rbind(selec_group, selec_group_tmp)}
        df_groups[i,] <- tmp_elements
        rownames(df_groups)[i] <- gen_set
        if ((mean(tail(sort(selec_group$perc_0, decreasing = T), n = 10,)) < 0.01) == T) {br <- T
        break}
      }
      index_j <- 0
      if(br == T){break}
      for (j in 1:df_length) {
        if(br == T){break}
        if (cluster == 1 & grepl(rownames(tmp3)[i], rownames(tmp3)[j]) == F) {
          index_j <- index_j + 1
          tmp_elements <- colSums(rbind(tmp3[i,], tmp3[j,])) 
          gen_set <- rownames(tmp3)[i]
          gen_set <- paste(gen_set, rownames(tmp3)[j])
          perc_1 <- length(tmp_elements[tmp_elements == 1])/length(tmp_elements)
          perc_0 <- length(tmp_elements[tmp_elements == 0])/length(tmp_elements)
          gene_combination <- gen_set
          if (is.nan(perc_1)) {perc_1 <- 0}
          if (is.nan(perc_0)) {perc_2 <- 0}
          selec_group_tmp <- data.frame(gene_combination, perc_1, perc_0)
          selec_group <- rbind(selec_group, selec_group_tmp)
          df_groups[index_j,] <- tmp_elements
          rownames(df_groups)[index_j] <- gen_set
          if ((mean(tail(sort(selec_group$perc_0, decreasing = T), n = 10,)) < 0.01) == T) {br <- T
          break}
          if (j == df_length) {
            df_length <- length(rownames(df_groups))}
        } else if (cluster != 1 & grepl(rownames(tmp3)[i], rownames(df_groups)[j]) == F) {
          index_j <- index_j + 1
          colnames(df_groups) <- colnames(tmp3)
          tmp_elements <- colSums(rbind(df_groups[j,], tmp3[i,]))
          gen_set <- rownames(df_groups)[j]
          gen_set <- paste(gen_set, rownames(tmp3)[i])
          perc_1 <- length(tmp_elements[tmp_elements == 1])/length(tmp_elements)
          perc_0 <- length(tmp_elements[tmp_elements == 0])/length(tmp_elements)
          gene_combination <- gen_set
          if (is.nan(perc_1)) {perc_1 <- 0}
          if (is.nan(perc_0)) {perc_2 <- 0}
          selec_group_tmp <- data.frame(gene_combination, perc_1, perc_0)
          selec_group <- rbind(selec_group, selec_group_tmp)
          df_groups_tmp[index_j,] <- tmp_elements
          rownames(df_groups_tmp)[index_j] <- gen_set
          if (j == df_length) {
            df_groups <- df_groups_tmp
            df_length <- length(rownames(df_groups))
            if ((mean(tail(sort(selec_group$perc_0, decreasing = T), n = 10,)) < 0.01) == T ) {br <- T
            break}
          }
        }
      }
    }
  }
  
  if ((tail(selec_group$perc_0[order(selec_group$perc_0, decreasing = T)], n = 1) <= 0.05) == T) {
    selec_group <- selec_group[selec_group$perc_0 <= 0.05, ]
  } else if ((tail(selec_group$perc_0[order(selec_group$perc_0, decreasing = T)], n = 1) <= 0.1) == T) {
    selec_group <- selec_group[selec_group$perc_0 <= 0.1, ]
  } 
  
  pca_cluster_genes[paste(pca_cluster)] <- selec_group$gene_combination[order(selec_group$perc_1[1], decreasing = T)]
  
} 


close(pb)
stopCluster(cl) 

print('Single cell types marker list')
print(pca_cluster_genes)

###############################################################################

##Cell nameing

##Create markers DF

unl <- unlist(pca_cluster_genes)


check = 0
for (num in 1:length(unl)) {
  for (gen in UMI.markers$gene) {
    gen_vector <- gen[grepl(gen, unl[num])]
    cluster <- num
    if (check == 0 & cluster == 1 & length(gen_vector) == 1) {
      check = 1
      gen_vector_list <- data.frame(cluster, gen_vector)
    } else if (check == 1 & length(gen_vector) == 1){
      gen_vector <- gen[grepl(gen, unl[num])]
      cluster <- num
      gen_vector_list_tmp <- data.frame(cluster, gen_vector)
      gen_vector_list <- unique(rbind(gen_vector_list, gen_vector_list_tmp))
      
    }
  }
}

print('Single cell naming')

cell_names <- c()
for (cluster in 1:max(as.numeric(unique(UMI@active.ident)))) {
  tmp <- subset(UMI, features = gen_vector_list$gen_vector[gen_vector_list$cluster %in% cluster])
  tmp <- as.matrix(GetAssayData(tmp, slot = 'data'))
  colnames(tmp) <- UMI@active.ident
  for (i in 1:length(colnames(tmp))) {
    tmp <- tmp[order(tmp[,i], decreasing = T), ,drop =F]
    if (as.numeric(cluster-1) %in% colnames(tmp)[i] & colSums(tmp)[i] > 0) {
      cell_names[i] <- rownames(tmp)[1]
    } else if (as.numeric(cluster-1) %in% colnames(tmp)[i] & colSums(tmp)[i] == 0) {
      cell_names[i] <- 'Bad'
    }
  }
}

print('Naming - DONE')

############################################
#PCA Cluster average expression for nameing

#######################################################################################################


subset_num <- round(length(colnames(UMI))/10000)
cells_num <- round(length(colnames(UMI))/subset_num)
exp_matrix <- GetAssayData(UMI, slot = 'data')
cols <- as.data.frame(summary(UMI@active.ident, maxsum = length(unique(colnames(exp_matrix)))))
colnames(cols)[1] <- 'n'

if (round(length(colnames(UMI))) > 14999) {
  
  for (batch in 1:subset_num) {
    if (batch == 1 & round(length(colnames(UMI))) > 14999) {
      average_expression <- as.data.frame(exp_matrix[,1:cells_num])
      colnames(average_expression) <- UMI@active.ident[1:cells_num]
      average_expression <- t(average_expression)
      average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
      rownames(average_expression) <- average_expression[,1]
      average_expression <- average_expression[,-1]
      average_expression <- t(average_expression)
      
    } else if (batch == subset_num & round(length(colnames(UMI))) > 14999) {
      average_expression_tmp <- as.data.frame(exp_matrix[,(((batch - 1) * cells_num)+1):length(colnames(exp_matrix))])
      colnames(average_expression_tmp) <- UMI@active.ident[(((batch - 1) * cells_num)+1):length(colnames(exp_matrix))]
      print((((batch - 1) * cells_num)+1))
      print(length(colnames(UMI)))
      average_expression_tmp <- t(average_expression_tmp)
      average_expression_tmp <- aggregate(average_expression_tmp, by = list(rownames(average_expression_tmp)), FUN = sum)
      rownames(average_expression_tmp) <- average_expression_tmp[,1]
      average_expression_tmp <- average_expression_tmp[,-1]
      average_expression_tmp <- t(average_expression_tmp)
      average_expression <- cbind(average_expression, average_expression_tmp)
      rm(average_expression_tmp)
      average_expression <- as.data.frame(average_expression)
      average_expression <- t(average_expression)
      average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
      rownames(average_expression) <- average_expression[,1]
      average_expression <- average_expression[,-1]
      average_expression <- t(average_expression)
      average_expression <- average_expression[,order(as.numeric(colnames(average_expression)))]
      average_expression <- as.data.frame(average_expression)
      
      num_col <- 0
      for (col in average_expression) {
        num_col <- num_col + 1
        num_row <- 0
        for (mean in col) {
          num_row <- num_row + 1
          average_expression[num_row , num_col] <- as.numeric(average_expression[num_row , num_col]/cols$n[num_col])
        } 
      }
    } else if (batch > 1 & batch < subset_num & round(length(colnames(UMI))) > 14999) {
      average_expression_tmp <- as.data.frame(exp_matrix[,(((batch - 1) * cells_num)+1):(batch * cells_num)])
      colnames(average_expression_tmp) <- UMI@active.ident[(((batch - 1) * cells_num)+1):(batch * cells_num)]
      print((((batch - 1) * cells_num)+1))
      print(batch * cells_num)
      average_expression_tmp <- t(average_expression_tmp)
      average_expression_tmp <- aggregate(average_expression_tmp, by = list(rownames(average_expression_tmp)), FUN = sum)
      rownames(average_expression_tmp) <- average_expression_tmp[,1]
      average_expression_tmp <- average_expression_tmp[,-1]
      average_expression_tmp <- t(average_expression_tmp)
      average_expression <- cbind(average_expression, average_expression_tmp)
      rm(average_expression_tmp)
    } 
  }
}

if (round(length(colnames(UMI))) < 14999) {
  average_expression <- as.data.frame(exp_matrix)
  colnames(average_expression) <- UMI@active.ident
  average_expression <- t(average_expression)
  average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
  rownames(average_expression) <- average_expression[,1]
  average_expression <- average_expression[,-1]
  average_expression <- t(average_expression)
  average_expression <- average_expression[,order(as.numeric(colnames(average_expression)))]
  average_expression <- as.data.frame(average_expression)
  
  num_col <- 0
  for (col in average_expression) {
    num_col <- num_col + 1
    num_row <- 0
    for (mean in col) {
      num_row <- num_row + 1
      average_expression[num_row , num_col] <- as.numeric(average_expression[num_row , num_col]/cols$n[num_col])
      
    } 
    
  }
}

#################################################################################

print('Clusters naming')

cluster_nameing(matrix_a = average_expression, markers = markers_class)

clust_names <- colnames(average_expression)


##########################################

if (length(markers_subclass) != 0) {
  
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  
  cell_names.1 <- c()
  for (i in 1:length(colnames(average_expression))) {
    rename_df <- average_expression[rownames(average_expression) %in% markers_subclass,]
    rename_df <- as.data.frame(rename_df[order(rename_df[,i], decreasing = T), ,drop = F])
    cell_names.1[i] <- rownames(rename_df[1,])
  }
  
  
  cluster <- 0
  cell_names.2 <- c()
  for (col in 1:length(colnames(average_expression))) {
    tmp_names <- top10$gene[top10$cluster %in% cluster]
    cluster <- cluster + 1
    rename_df <- average_expression[rownames(average_expression) %in% toupper(tmp_names),]
    rename_df <- as.data.frame(rename_df[order(rename_df[,col], decreasing = T), ,drop = F])
    if (!cell_names.1[col] %in% toupper(rownames(rename_df[1,]))) {
      cell_names.2[col] <- firstup(tolower(rownames(rename_df[1,])))
    } else if (cell_names.1[col] %in% toupper(rownames(rename_df[1,]))) {
      cell_names.2[col] <- firstup(tolower(rownames(rename_df[2,])))
    }
  }
  
  colnames(average_expression) <- paste(cell_names.1, cell_names.2)
  
} else if (length(markers_subclass) == 0) {
  
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  cluster <- 0
  cell_names.1 <- c()
  cell_names.2 <- c()
  for (col in 1:length(colnames(average_expression))) {
    tmp_names <- top10$gene[top10$cluster %in% cluster]
    cluster <- cluster + 1
    rename_df <- average_expression[rownames(average_expression) %in% toupper(tmp_names),]
    rename_df <- as.data.frame(rename_df[order(rename_df[,col], decreasing = T), ,drop = F])
    cell_names.1[col] <- toupper(rownames(rename_df[1,]))
    cell_names.2[col] <- firstup(tolower(rownames(rename_df[2,])))
  }
  
  
  colnames(average_expression) <- paste(cell_names.1, cell_names.2)
  
}


subclass_marker_list_pheat <- rownames(exp_matrix)[toupper(rownames(exp_matrix)) %in% toupper(cell_names.2)]
subclass_marker_list_pheat <- unique(c(subclass_marker_list_pheat, cell_names[!grepl('Bad',cell_names)]))
#######################################################################################
#Repair subclass_names

new.cluster.ids <- paste(clust_names, colnames(average_expression))
names(new.cluster.ids) <- levels(UMI)
UMI <- RenameIdents(UMI, new.cluster.ids)

print('Naming - DONE')


#################################################################################

colnames(average_expression) <- new.cluster.ids
dir.create(path = file.path(OUTPUT,'exp_matrix'))
write.table(average_expression, file = file.path(OUTPUT, "exp_matrix/class_average_expression.csv"), sep = ',')

########################NOWE#####################################################

#PCA plot and UMAP plot with names

#Class

jpeg(file.path(OUTPUT, "PCA_DimPlot_class.jpeg") , units="in", width=10, height=7, res=600)
DimPlot(UMI, reduction = "pca")
dev.off()


jpeg(file.path(OUTPUT, "UMAP_with_DE_gene_class.jpeg") , units="in", width=10, height=7, res=600)
DimPlot(UMI, reduction = "umap") 
dev.off()


#Subtypes

new.names <- paste0(UMI@active.ident,' - ', cell_names)


Idents(UMI) <- new.names


######################################################################################################################


#Average expression matrix populations
rm(average_expression)

subset_num <- round(length(colnames(UMI))/10000)
cells_num <- round(length(colnames(UMI))/subset_num)
cols <- as.data.frame(summary(UMI@active.ident, maxsum = length(unique(colnames(exp_matrix)))))
cols <- as.data.frame(cols[order(rownames(cols)), , drop = F])
colnames(cols)[1] <- 'n'

if (round(length(colnames(UMI))) > 14999){
  
  for (batch in 1:subset_num) {
    if (batch == 1 & round(length(colnames(UMI))) > 14999) {
      average_expression <- as.data.frame(exp_matrix[,1:cells_num])
      colnames(average_expression) <- UMI@active.ident[1:cells_num]
      average_expression <- t(average_expression)
      average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
      rownames(average_expression) <- average_expression[,1]
      average_expression <- average_expression[,-1]
      average_expression <- t(average_expression)
    } else if (batch == subset_num & round(length(colnames(UMI))) > 14999) {
      average_expression_tmp <- as.data.frame(exp_matrix[,(((batch - 1) * cells_num)+1):length(colnames(exp_matrix))])
      colnames(average_expression_tmp) <- UMI@active.ident[(((batch - 1) * cells_num)+1):length(colnames(exp_matrix))]
      print((((batch - 1) * cells_num)+1))
      print(length(colnames(UMI)))
      average_expression_tmp <- t(average_expression_tmp)
      average_expression_tmp <- aggregate(average_expression_tmp, by = list(rownames(average_expression_tmp)), FUN = sum)
      rownames(average_expression_tmp) <- average_expression_tmp[,1]
      average_expression_tmp <- average_expression_tmp[,-1]
      average_expression_tmp <- t(average_expression_tmp)
      average_expression <- cbind(average_expression, average_expression_tmp)
      rm(average_expression_tmp)
      average_expression <- as.data.frame(average_expression)
      average_expression <- t(average_expression)
      average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
      rownames(average_expression) <- average_expression[,1]
      average_expression <- average_expression[,-1]
      average_expression <- t(average_expression)
      average_expression <- average_expression[,order(colnames(average_expression))]
      average_expression <- as.data.frame(average_expression)
      
      num_col <- 0
      for (col in average_expression) {
        num_col <- num_col + 1
        num_row <- 0
        for (mean in col) {
          num_row <- num_row + 1
          average_expression[num_row , num_col] <- as.numeric(average_expression[num_row , num_col]/cols$n[num_col])
        } 
      }
      
    } else if (batch > 1 & batch < subset_num & round(length(colnames(UMI))) > 14999) {
      average_expression_tmp <- as.data.frame(exp_matrix[,(((batch - 1) * cells_num)+1):(batch * cells_num)])
      colnames(average_expression_tmp) <- UMI@active.ident[(((batch - 1) * cells_num)+1):(batch * cells_num)]
      print((((batch - 1) * cells_num)+1))
      print(batch * cells_num)
      average_expression_tmp <- t(average_expression_tmp)
      average_expression_tmp <- aggregate(average_expression_tmp, by = list(rownames(average_expression_tmp)), FUN = sum)
      rownames(average_expression_tmp) <- average_expression_tmp[,1]
      average_expression_tmp <- average_expression_tmp[,-1]
      average_expression_tmp <- t(average_expression_tmp)
      average_expression <- cbind(average_expression, average_expression_tmp)
      rm(average_expression_tmp)
    } 
  }
}

if (round(length(colnames(UMI))) < 14999) {
  average_expression <- as.data.frame(exp_matrix)
  colnames(average_expression) <- UMI@active.ident
  average_expression <- t(average_expression)
  average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
  rownames(average_expression) <- average_expression[,1]
  average_expression <- average_expression[,-1]
  average_expression <- t(average_expression)
  average_expression <- average_expression[,order(colnames(average_expression))]
  average_expression <- as.data.frame(average_expression)
  
  num_col <- 0
  for (col in average_expression) {
    num_col <- num_col + 1
    num_row <- 0
    for (mean in col) {
      num_row <- num_row + 1
      average_expression[num_row , num_col] <- as.numeric(average_expression[num_row , num_col]/cols$n[num_col])
      
    } 
    
  }
}
#############################################################################################################################

print('Checking and renaming subtypes')


#remove empty cells (without markers expression)

class_marker_list <- c()
for (class in markers_class) {
  class_marker_list <- c(class_marker_list, textclean::mgsub(class, c('+'), c('')))
}


class_marker_list <- rownames(UMI)[toupper(rownames(UMI)) %in% class_marker_list]


second_matrix <- average_expression[class_marker_list,]



renamed_old.1 <- c()
renamed_new.1 <- c()
renamed_old.2 <- c()
renamed_new.2 <- c()
col_new.names <- as.character(colnames(second_matrix))

index_marker <- 0
for (marker in markers_class) {
  index_marker <- index_marker + 1
  for (col in 1:length(colnames(second_matrix))) {
    second_matrix <- as.data.frame(second_matrix[order(second_matrix[,col], decreasing = T), ,drop = F])
    if (max(second_matrix[,col]) == 0 & !grepl('Unknow', as.character(colnames(second_matrix)[col]))) {
      colnames(second_matrix)[col] <- 'Bad'
    } else if (grepl(as.character(colnames(markers_class)[index_marker]), as.character(colnames(second_matrix)[col])) & !grepl(as.character(toupper(rownames(second_matrix)[1])), as.character(list(textclean::mgsub(marker, c('+'), c('')))))) {
      renamed_old.1 <- c(renamed_old.1, colnames(second_matrix)[col])
      mark <- c()
      for (change in colnames(markers_class)) {
        mark <- c(mark, change[grepl(change, colnames(second_matrix)[col])])
      }
      n_marker = 0
      for (marker in markers_class) {
        n_marker = n_marker + 1
        if (grepl(toupper(rownames(second_matrix)[1]), list(marker))) {
          colnames(second_matrix)[col] <- gsub(pattern = mark, replacement =  colnames(markers_class)[n_marker], x = colnames(second_matrix)[col])
          renamed_new.1 <- c(renamed_new.1, colnames(second_matrix)[col])
        }
      }
    }
  }
}


Renamed_idents <- as.character(Idents(UMI))


if (length(renamed_old.1) != 0) {
  n = 0
  for (name in 1:length(renamed_new.1)) {
    n = n + 1
    Renamed_idents <- gsub(x = Renamed_idents, pattern = renamed_old.1[n], replacement = renamed_new.1[n])
  }
}




if (length(markers_subclass) != 0) {
  subclass_marker_list <- c()
  for (subclass in markers_subclass) {
    subclass_marker_list <- c(subclass_marker_list, subclass)
  }
}
  subclass_marker_list <- rownames(UMI)[toupper(rownames(UMI)) %in% subclass_marker_list]
 

  
  old.names <- colnames(second_matrix)
  
  
  second_matrix <- average_expression[subclass_marker_list,]
  colnames(second_matrix) <- old.names
  
  ################################################
  #Renamed function
  
  second_matrix <- second_matrix[!colnames(second_matrix) %in% c('Bad', 'Unknow')]
  
  for (col in 1:length(colnames(second_matrix))) {
    second_matrix <- as.data.frame(second_matrix[order(second_matrix[,col], decreasing = T), ,drop = F])
    if (second_matrix[1,col] == 0) {
      colnames(second_matrix)[col] <- 'Bad'
    } else if (!grepl(toupper(rownames(second_matrix)[1]), colnames(second_matrix)[col]) & !grepl('Bad', colnames(second_matrix)[col])) {
      renamed_old.2 <- unique(c(renamed_old.2, colnames(second_matrix)[col]))
      mark <- c()
      for (change in markers_subclass) {
        mark <- c(mark, change[grepl(change, colnames(second_matrix)[col])])
      }
      colnames(second_matrix)[col] <- gsub(pattern = mark, replacement =  toupper(rownames(second_matrix)[1]), x = colnames(second_matrix)[col])
      renamed_new.2 <- c(renamed_new.2, colnames(second_matrix)[col])
    } 
  }

  
  if (length(renamed_old.2) != 0) {
    n = 0
    for (name in 1:length(renamed_new.2)) {
      n = n + 1
      Renamed_idents <- gsub(x = Renamed_idents, pattern = renamed_old.2[n], replacement = renamed_new.2[n])
    }
  }
  
  


Idents(UMI) <- Renamed_idents

print('Checking - DONE')

###############################################

print('QC of subtypes')

subclass_names <- Renamed_idents
bad <- subclass_names[grepl('Bad', as.character(subclass_names))]
renamed.subnames <- c(as.character(renamed_new.1), as.character(renamed_new.2))
renamed.subnames <- renamed.subnames[!as.character(renamed.subnames) %in% as.character(new.names)]
renamed.subnames <- renamed.subnames[!as.character(renamed.subnames) %in% as.character(bad)]
new.subnames <- subclass_names[!as.character(subclass_names) %in% as.character(bad)]
new.subnames <- new.subnames[!as.character(new.subnames) %in% as.character(renamed.subnames)]
bad.subnames <- subclass_names[as.character(subclass_names) %in% as.character(bad)]


#############################################################################################################################

data <- as.data.frame(summary(as.factor(subclass_names), maxsum = length(unique(subclass_names))))
colnames(data)[1] <- 'n'
data$names <- rownames(data)


num <- head(sort(data$n, decreasing = T), n = 1)*0.01


below.names <- data$names[data$n < num]


data$test[data$names %in% bad.subnames] <- "Bad marked types"
data$test[data$names %in% renamed.subnames] <- "Renamed"
data$test[data$names %in% new.subnames] <- "Good marked types" 
data$test[data$names %in% below.names] <- "Non-significant < 0.01"



threshold <- ggplot(data, aes(x = n, y = reorder(names, -n), fill = test, sort = test)) +
  geom_bar(stat = 'identity') +
  ylab("Cells types") +
  xlab("Number of cells")+
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) + 
  labs(fill = "Cells threshold")

ggsave(threshold, filename = file.path(OUTPUT,'cells_type_threshold.jpeg'), units = 'in', width = 15, height = 10, dpi = 600)



#save bad cells

bad.subnames <- c(as.character(bad.subnames), as.character(below.names))
bad.subnames <- unique(as.character(bad.subnames))


UMI_unknow <- subset(UMI, idents = bad.subnames)

bad_cells <- as.data.frame(GetAssayData(object = UMI_unknow, slot = 'counts'))
colnames(bad_cells)[1:length(colnames(bad_cells))] <- 'Unknow'

write.table(bad_cells, file = file.path(OUTPUT, "exp_matrix/unknow_cells_count_matrix.csv"), sep = ',')

rm(bad_cells)

#####################################################################################################

right.names <- unique(subclass_names[!as.character(subclass_names) %in% as.character(bad.subnames)])


UMI <- subset(UMI, idents = right.names)

print('DONE')


#########################################################################################################################################################################################

jpeg(file.path(OUTPUT, "PCA_DimPlot_subtypes.jpeg") , units="in", width=10, height=7, res=600)
DimPlot(UMI, reduction = "pca")
dev.off()

jpeg(file.path(OUTPUT, "UMAP_with_DE_gene_subtypes.jpeg") , units="in", width=10, height=7, res=600)
DimPlot(UMI, reduction = "umap") 
dev.off()


#Create Expression Matrix

#Expression matrix cells

#######################3


subset_num <- round(length(colnames(UMI))/1000)
cells_num <- round(length(colnames(UMI))/subset_num)
exp_matrix_obl <- GetAssayData(UMI, slot = 'data')

if (round(length(colnames(UMI))) > 14999){
  
  for (batch in 1:subset_num) {
    if (batch == 1 & round(length(colnames(UMI))) > 14999) {
      exp_matrix_obl_tmp <- as.data.frame(exp_matrix_obl[,1:cells_num])
      colnames(exp_matrix_obl_tmp) <- UMI@active.ident[1:cells_num]
      mean_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, mean)
      exp_matrix_obl_tmp[exp_matrix_obl_tmp != 0] <- 1
      positive_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, sum)
      sum_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, length)
      positive_expression_perc <- positive_expression/sum_expression
      names <- colnames(exp_matrix_obl_tmp)
      exp_stat <- data.frame(names,mean_expression, positive_expression_perc)
    } else if (batch == subset_num & round(length(colnames(UMI))) > 14999) {
      exp_matrix_obl_tmp <- as.data.frame(exp_matrix_obl[,(((batch - 1) * cells_num)+1):length(colnames(exp_matrix_obl))])
      colnames(exp_matrix_obl_tmp) <- UMI@active.ident[(((batch - 1) * cells_num)+1):length(colnames(exp_matrix_obl))]
      mean_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, mean)
      exp_matrix_obl_tmp[exp_matrix_obl_tmp != 0] <- 1
      positive_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, sum)
      sum_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, length)
      positive_expression_perc <- positive_expression/sum_expression
      names <- colnames(exp_matrix_obl_tmp)
      exp_stat_tmp <- data.frame(names,mean_expression, positive_expression_perc)
      exp_stat <- rbind(exp_stat, exp_stat_tmp)
    } else if (batch > 1 & batch < subset_num & round(length(colnames(UMI))) > 14999) {
      exp_matrix_obl_tmp <- as.data.frame(exp_matrix_obl[,(((batch - 1) * cells_num)+1):(batch * cells_num)])
      colnames(exp_matrix_obl_tmp) <- UMI@active.ident[(((batch - 1) * cells_num)+1):(batch * cells_num)]
      mean_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, mean)
      exp_matrix_obl_tmp[exp_matrix_obl_tmp != 0] <- 1
      positive_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, sum)
      sum_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, length)
      positive_expression_perc <- positive_expression/sum_expression
      names <- colnames(exp_matrix_obl_tmp)
      exp_stat_tmp <- data.frame(names,mean_expression, positive_expression_perc)
      exp_stat <- rbind(exp_stat, exp_stat_tmp)
    } 
  }
} 

if (round(length(colnames(UMI))) < 14999){
  exp_matrix_obl_tmp <- as.data.frame(exp_matrix_obl)
  colnames(exp_matrix_obl_tmp) <- UMI@active.ident
  mean_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, mean)
  exp_matrix_obl_tmp[exp_matrix_obl_tmp != 0] <- 1
  positive_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, sum)
  sum_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, length)
  positive_expression_perc <- positive_expression/sum_expression
  names <- colnames(exp_matrix_obl_tmp)
  exp_stat <- data.frame(names,mean_expression, positive_expression_perc)
}
###############################################################################################################


cells <- ggplot(exp_stat, mapping = aes(x = mean_expression, y = positive_expression_perc, fill = names)) +
  geom_boxplot() +
  facet_wrap(names~.) +
  theme(legend.position = 'none')

ggsave(cells, filename = file.path(OUTPUT,'box_matrix.jpeg'), units = 'in', width = 20, height = 15, dpi = 600)

################################################################################################################


exp_matrix <- GetAssayData(UMI, slot = 'data')
colnames(exp_matrix) <- UMI@active.ident

write(colnames(exp_matrix), file = file.path(OUTPUT, "exp_matrix/barcodes.tsv"))
write(rownames(exp_matrix), file = file.path(OUTPUT, "exp_matrix/genes.tsv"))
Matrix::writeMM(exp_matrix, file = file.path(OUTPUT, "exp_matrix/matrix.mtx"))


################################################################################################################
#Average expression matrix populations

rm(average_expression)

subset_num <- round(length(colnames(UMI))/10000)
cells_num <- round(length(colnames(UMI))/subset_num)
cols <- as.data.frame(summary(UMI@active.ident, maxsum = length(unique(colnames(exp_matrix)))))
cols <- as.data.frame(cols[order(rownames(cols)), , drop = F])
colnames(cols)[1] <- 'n'

if (round(length(colnames(UMI))) > 14999){
  
  for (batch in 1:subset_num) {
    if (batch == 1 & round(length(colnames(UMI))) > 14999) {
      average_expression <- as.data.frame(exp_matrix[,1:cells_num])
      colnames(average_expression) <- UMI@active.ident[1:cells_num]
      average_expression <- t(average_expression)
      average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
      rownames(average_expression) <- average_expression[,1]
      average_expression <- average_expression[,-1]
      average_expression <- t(average_expression)
    } else if (batch == subset_num & round(length(colnames(UMI))) > 14999) {
      average_expression_tmp <- as.data.frame(exp_matrix[,(((batch - 1) * cells_num)+1):length(colnames(exp_matrix))])
      colnames(average_expression_tmp) <- UMI@active.ident[(((batch - 1) * cells_num)+1):length(colnames(exp_matrix))]
      print((((batch - 1) * cells_num)+1))
      print(length(colnames(UMI)))
      average_expression_tmp <- t(average_expression_tmp)
      average_expression_tmp <- aggregate(average_expression_tmp, by = list(rownames(average_expression_tmp)), FUN = sum)
      rownames(average_expression_tmp) <- average_expression_tmp[,1]
      average_expression_tmp <- average_expression_tmp[,-1]
      average_expression_tmp <- t(average_expression_tmp)
      average_expression <- cbind(average_expression, average_expression_tmp)
      rm(average_expression_tmp)
      average_expression <- as.data.frame(average_expression)
      average_expression <- t(average_expression)
      average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
      rownames(average_expression) <- average_expression[,1]
      average_expression <- average_expression[,-1]
      average_expression <- t(average_expression)
      average_expression <- average_expression[,order(colnames(average_expression))]
      average_expression <- as.data.frame(average_expression)
      
      num_col <- 0
      for (col in average_expression) {
        num_col <- num_col + 1
        num_row <- 0
        for (mean in col) {
          num_row <- num_row + 1
          average_expression[num_row , num_col] <- as.numeric(average_expression[num_row , num_col]/cols$n[num_col])
        } 
      }
      
    } else if (batch > 1 & batch < subset_num & round(length(colnames(UMI))) > 14999) {
      average_expression_tmp <- as.data.frame(exp_matrix[,(((batch - 1) * cells_num)+1):(batch * cells_num)])
      colnames(average_expression_tmp) <- UMI@active.ident[(((batch - 1) * cells_num)+1):(batch * cells_num)]
      print((((batch - 1) * cells_num)+1))
      print(batch * cells_num)
      average_expression_tmp <- t(average_expression_tmp)
      average_expression_tmp <- aggregate(average_expression_tmp, by = list(rownames(average_expression_tmp)), FUN = sum)
      rownames(average_expression_tmp) <- average_expression_tmp[,1]
      average_expression_tmp <- average_expression_tmp[,-1]
      average_expression_tmp <- t(average_expression_tmp)
      average_expression <- cbind(average_expression, average_expression_tmp)
      rm(average_expression_tmp)
    } 
  }
}

if (round(length(colnames(UMI))) < 14999) {
  average_expression <- as.data.frame(exp_matrix)
  colnames(average_expression) <- UMI@active.ident
  average_expression <- t(average_expression)
  average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
  rownames(average_expression) <- average_expression[,1]
  average_expression <- average_expression[,-1]
  average_expression <- t(average_expression)
  average_expression <- average_expression[,order(colnames(average_expression))]
  average_expression <- as.data.frame(average_expression)
  
  num_col <- 0
  for (col in average_expression) {
    num_col <- num_col + 1
    num_row <- 0
    for (mean in col) {
      num_row <- num_row + 1
      average_expression[num_row , num_col] <- as.numeric(average_expression[num_row , num_col]/cols$n[num_col])
      
    } 
    
  }
}

write.table(average_expression, file = file.path(OUTPUT, 'exp_matrix/average_expression_matrix_subtypes.csv'), sep = ',')


#save seurat output file

saveRDS(UMI, file = file.path(OUTPUT, "Results.rds"))


#######################################################################################

#Cell populations pheatmaps

ms<- c()
for (m in subclass_marker_list) {
ms <- c(ms, m[grepl(toupper(m), list(colnames(average_expression)))])
}

marker_list <- unique(c(class_marker_list, subclass_marker_list_pheat, ms))




average_expression <- average_expression[marker_list,]
average_expression <- drop_na(average_expression)
jpeg(file.path(OUTPUT, "pheatmap_cells_populations.jpeg"),units="in", width = 35, height = 30,  res=600)
pheatmap::pheatmap(average_expression, 
                   clustering_method = 'ward.D',
                   angle_col = 270, fontsize_row = 20, fontsize_col = 20)
dev.off()


#########################################################################################

print('Report creating')

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
