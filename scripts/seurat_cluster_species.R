library(Seurat)
library(patchwork)
library(tidyverse)


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
  n_gen <- as.numeric(mean(n_gen))*2 + 1.5*IQR(UMI@meta.data$nGenes)
  
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
  
  UMI <- subset(UMI, subset = nGenes > 500 & nGenes <= n_gen & MitoPercent < 5)
  n_gen <- as.numeric(max(UMI@meta.data$nGenes))*0.25
  cells_number <- length(Idents(UMI))
  
  ######################################################################################
  #Cells_stats
  
  cells <- factor(c('Estimated_cells', 'Input_cells', 'Analyzed_cells'), levels = c('Estimated_cells', 'Input_cells', 'Analyzed_cells'))
  cell_num <- c(estimated_cells, cell_input, cells_number)
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
  
  decision <- mean(UMI@assays$RNA@counts@x) < 15
  if (data == 2 & decision == T) {
    UMI@assays$RNA@data@x <- UMI@assays$RNA@counts@x
  } else if (data == 2 & decision == F) {
    UMI@assays$RNA@data@x <- log(UMI@assays$RNA@counts@x+1)
  }
  
  #######################################################################################
  
  UMI <- FindVariableFeatures(UMI, selection.method = "vst", nfeatures = n_gen)
  
  # Identify the 10 most highly variable genes
  
  top20 <- head(VariableFeatures(UMI), 20)
  top100 <- head(VariableFeatures(UMI), 100)
  
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
      if (i-i*0.01 > dims$`Elbow$data$stdev`[element+1] & element < 50 | i-i*0.02 > dims$`Elbow$data$stdev`[element+2] & element < 49 | i-i*0.02 > dims$`Elbow$data$stdev`[element+3] & element < 48) {
        dim <- dim + 1
      } else 
        break
    }
    dim <- as.numeric(dim)
  }
  
  
  #########################################################################################
  
  UMI <- JackStraw(UMI, num.replicate = 100, dims = dim)
  UMI <- ScoreJackStraw(UMI, dims = 1:dim)
  
  
  jpeg(file.path(OUTPUT, "JackStrawPlot.jpeg") , units="in", width=10, height=7, res=600)
  JackStrawPlot(UMI, dims = 1:(dim))
  dev.off()
  
  
  UMI <- FindNeighbors(UMI, dims = 1:dim)
  UMI <- FindClusters(UMI, resolution = 0.5)
  
  
  UMI <- RunUMAP(UMI, dims = 1:dim)
  
  
  jpeg(file.path(OUTPUT, "UMAP.jpeg") , units="in", width=10, height=7, res=600)
  DimPlot(UMI, reduction = "umap")
  dev.off()
  
  ####################################################################################
  
  
  
  #find markers for every cluster compared to all remaining cells, report only the positive ones
  
  
  UMI.markers <- FindAllMarkers(UMI, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'MAST')
  top1 <- UMI.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
  MAST_markers <- UMI.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
  
  write.table(MAST_markers, file = file.path(OUTPUT, "MAST_markers_clusters.csv"), sep = ',')
  
  ##Cells cluster naming with top genes (different between cell groups)
  
  ##########################################################
  ### Most variable genes select and cell subtypes nameing 
  
  
  #################ZMIANA#################################
  exp_matrix <- subset(UMI, features = top100)
  exp_matrix <- as.data.frame(GetAssayData(exp_matrix, slot = 'data'))
  colnames(exp_matrix) <- UMI@active.ident
  
  
  for (cluster in min(colnames(exp_matrix)):max(colnames(exp_matrix))) {
    tmp <- exp_matrix[colnames(exp_matrix) %in% cluster]
    clust_cells <- hclust(dist(t(tmp)), method = 'ward.D')
    clusters_cells <- as.data.frame(cutree(clust_cells, h = 20))
    tmp <- as.data.frame(t(tmp))
    tmp$clusters <- clusters_cells$`cutree(clust_cells, h = 20)`
    gen_model <- glm(formula = clusters ~ ., data = tmp)
    summary <- as.data.frame(summary.glm(gen_model)$coefficients)
    summary <- summary[-1,]
    summary <- summary[summary$`Pr(>|t|)` <= 0.001, ]
    cluster_genes <- rownames(summary)
    index <- 0
    for (i in colnames(exp_matrix)) {
      index = index +1
      if (i %in% cluster) {
        tmp <- exp_matrix[cluster_genes,]
        tmp <- as.data.frame(exp_matrix[order(exp_matrix[,index], decreasing = T), ,drop =F])
        colnames(exp_matrix)[index] <- rownames(tmp[1,])
      }
    }
  } 
  
  
  

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
  
  #####################NOWE##################################################
  
  
  
  subset_num <- round(length(colnames(UMI))/10000)
  cells_num <- round(length(colnames(UMI))/subset_num)
  exp_matrix <- GetAssayData(UMI, slot = 'data')
  cols <- as.data.frame(summary(UMI@active.ident))
  
  for (batch in 1:subset_num) {
    if (batch == 1 & round(length(colnames(UMI))) > 10000) {
      average_expression <- as.data.frame(exp_matrix[,1:cells_num])
      colnames(average_expression) <- UMI@active.ident[1:cells_num]
      average_expression <- t(average_expression)
      average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
      rownames(average_expression) <- average_expression[,1]
      average_expression <- average_expression[,-1]
      average_expression <- t(average_expression)
    } else if (batch > 1 & batch < subset_num & round(length(colnames(UMI))) > 10000) {
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
    } else if (batch == subset_num & round(length(colnames(UMI))) > 10000) {
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
          average_expression[num_row , num_col] <- as.numeric(average_expression[num_row , num_col]/cols$`summary(UMI@active.ident)`[num_col])
        } 
      }
      
    } else if (round(length(colnames(UMI))) < 10000){
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
          average_expression[num_row , num_col] <- as.numeric(average_expression[num_row , num_col]/cols$`summary(UMI@active.ident)`[num_col])
          
        } 
        
      }
    }
  }
  
  #################################################################################
  
  
  
  cluster_nameing(matrix_a = average_expression, markers = markers_class)
  
  clust_names <- colnames(average_expression)
  
  
  ##########################################
  
  
  if (length(markers_subclass) != 0) {
    
    
    index = 0
    for (i in colnames(average_expression)) {
      index = index +1
      rename_df <- average_expression[rownames(average_expression) %in% markers_subclass,]
      rename_df <- as.data.frame(rename_df[order(rename_df[,index], decreasing = T), ,drop = F])
      colnames(average_expression)[index] <- rownames(rename_df[1,])
    }
    
    
  }
  
  
  if (length(markers_subclass) == 0) {
    
    for (col in 1:length(colnames(average_expression))) {
      colnames(average_expression)[col] <- top1$gene[col]
    }
    
  }
  
  
  new.cluster.ids <- paste(clust_names, colnames(average_expression))
  names(new.cluster.ids) <- levels(UMI)
  UMI <- RenameIdents(UMI, new.cluster.ids)
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
  
  new.names <- paste(UMI@active.ident, cell_names)
  
  
  Idents(UMI) <- new.names
  
  
  select_data <- as.data.frame(summary(as.factor(new.names), maxsum = length(unique(new.names))))
  colnames(select_data)[1] <- 'n'
  select_data$names <- rownames(select_data)
  num <- (max(select_data$n)-min(select_data$n))/max(select_data$n)*mean(select_data$n, trim = 0.25)
  right.names <- select_data[select_data$n > num,]
  bad.names <- select_data[select_data$n <= num,]
  
  
  select_data$test[select_data$n > num] <- "Good marked types" 
  select_data$test[select_data$n <= num] <- "Wrong marked types"
  
  
  threshold <- ggplot(select_data, aes(x = n, y = reorder(names, -n), fill = test, sort = test)) +
    geom_bar(stat = 'identity') +
    ylab("Cells types") +
    xlab("Number of cells")+
    theme_bw() +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) + 
    labs(fill = "Cells threshold")
  
  ggsave(threshold, filename = file.path(OUTPUT,'cells_type_threshold.jpeg'), units = 'in', width = 15, height = 10, dpi = 600)
  
  
  UMI <- subset(UMI, idents = right.names$names)
  
  jpeg(file.path(OUTPUT, "PCA_DimPlot_subtypes.jpeg") , units="in", width=30, height=15, res=600)
  DimPlot(UMI, reduction = "pca")
  dev.off()
  
  jpeg(file.path(OUTPUT, "UMAP_with_DE_gene_subtypes.jpeg") , units="in", width=30, height=15, res=600)
  DimPlot(UMI, reduction = "umap") 
  dev.off()
  
  
  #Create Expression Matrix
  
  #Expression matrix cells
  
  #######################3
  
  
  subset_num <- round(length(colnames(UMI))/1000)
  cells_num <- round(length(colnames(UMI))/subset_num)
  exp_matrix_obl <- GetAssayData(UMI, slot = 'data')
  cols <- as.data.frame(summary(UMI@active.ident))
  
  for (batch in 1:subset_num) {
    if (batch == 1 & round(length(colnames(UMI))) > 10000) {
      exp_matrix_obl_tmp <- as.data.frame(exp_matrix_obl[,1:cells_num])
      colnames(exp_matrix_obl_tmp) <- UMI@active.ident[1:cells_num]
      mean_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, mean)
      exp_matrix_obl_tmp[exp_matrix_obl_tmp != 0] <- 1
      positive_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, sum)
      sum_expression <- apply(exp_matrix_obl_tmp, MARGIN = 2, length)
      positive_expression_perc <- positive_expression/sum_expression
      names <- colnames(exp_matrix_obl_tmp)
      exp_stat <- data.frame(names,mean_expression, positive_expression_perc)
    } else if (batch > 1 & batch < subset_num & round(length(colnames(UMI))) > 10000) {
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
    } else if (batch == subset_num & round(length(colnames(UMI))) > 10000) {
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
    } else if (round(length(colnames(UMI))) < 10000){
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
  }
  
  ###############################################################################################################
  
  
  cells <- ggplot(exp_stat, mapping = aes(x = mean_expression, y = positive_expression_perc, fill = names)) +
    geom_boxplot() +
    facet_wrap(names~.) +
    theme(legend.position = 'none')
  
  ggsave(cells, filename = file.path(OUTPUT,'box_matrix.jpeg'), units = 'in', width = 15, height = 10, dpi = 600)
  
  ################################################################################################################
  
  
  exp_matrix <- GetAssayData(UMI, slot = 'data')
  
  write(colnames(exp_matrix), file = file.path(OUTPUT, "exp_matrix/barcodes.tsv"))
  write(rownames(exp_matrix), file = file.path(OUTPUT, "exp_matrix/genes.tsv"))
  Matrix::writeMM(exp_matrix, file = file.path(OUTPUT, "exp_matrix/matrix.mtx"))
  
  
  ################################################################################################################
  #Average expression matrix populations
  
  
  subset_num <- round(length(colnames(UMI))/10000)
  cells_num <- round(length(colnames(UMI))/subset_num)
  cols <- as.data.frame(summary(UMI@active.ident))
  cols <- as.data.frame(cols[order(rownames(cols)), , drop = F])
  
  for (batch in 1:subset_num) {
    if (batch == 1 & round(length(colnames(UMI))) > 10000) {
      average_expression <- as.data.frame(exp_matrix[,1:cells_num])
      colnames(average_expression) <- UMI@active.ident[1:cells_num]
      average_expression <- t(average_expression)
      average_expression <- aggregate(average_expression, by = list(rownames(average_expression)), FUN = sum)
      rownames(average_expression) <- average_expression[,1]
      average_expression <- average_expression[,-1]
      average_expression <- t(average_expression)
    } else if (batch > 1 & batch < subset_num & round(length(colnames(UMI))) > 10000) {
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
    } else if (batch == subset_num & round(length(colnames(UMI))) > 10000) {
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
          average_expression[num_row , num_col] <- as.numeric(average_expression[num_row , num_col]/cols$`summary(UMI@active.ident)`[num_col])
        } 
      }
      
    } else if (round(length(colnames(UMI))) < 10000){
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
          average_expression[num_row , num_col] <- as.numeric(average_expression[num_row , num_col]/cols$`summary(UMI@active.ident)`[num_col])
          
        } 
        
      }
    }
  }
  
  
  write.table(average_expression, file = file.path(OUTPUT, "exp_matrix/average_expression_matrix_subtypes.csv"), sep = ',')
  
  
  #save seurat output file
  
  saveRDS(UMI, file = file.path(OUTPUT, "Results.rds"))
  
  
  #######################################################################################
  
  #Cell populations pheatmaps
  
  class_marker_list <- c()
  for (class in markers_class) {
    print(class)
    class_marker_list <- c(class_marker_list, textclean::mgsub(class, c('+'), c('')))
  }
  class_marker_list <- rownames(UMI)[toupper(rownames(UMI)) %in% class_marker_list]
  
  subclass_marker_list <- c()
  for (subclass in markers_subclass) {
    subclass_marker_list <- c(subclass_marker_list, subclass)
  }
  subclass_marker_list <- rownames(UMI)[toupper(rownames(UMI)) %in% subclass_marker_list]
  
  
  
  if (length(markers_subclass) != 0) {
    
    
    marker_list <- unique(c(class_marker_list, subclass_marker_list))
    
  } else if (length(markers_subclass) == 0) {
    
    
    marker_list <- unique(c(class_marker_list, top1$gene))
    
  }
  
  
  
  average_expression <- average_expression[marker_list,]
  jpeg(file.path(OUTPUT, "pheatmap_cells_populations.jpeg"),units="in", width = 40, height = 25,  res=600)
  pheatmap::pheatmap(average_expression, 
                     clustering_method = 'ward.D',
                     angle_col = 270, fontsize_row = 20, fontsize_col = 20)
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
  
  