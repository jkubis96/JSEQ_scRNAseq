library(Seurat)
library(tidyverse)
library(doSNOW)
library(foreach)
library(doParallel)
library(stringr)

args <- commandArgs()

#Paths and arguments from env
{
  
  print(args)
  
  path <- args[6]
  markers <-args[7]
  species <- args[8]
  if (tolower(species) == 'human') {
    species = 'Homo sapiens'
  }
  seurat_umi <- file.path(path,'sc_data/')
  OUTPUT <- file.path(path, 'results')
  project_name <- args[9]
  data <- args[10]
  estimated_cells <- args[11]
  dir.create(path = file.path(OUTPUT,'matrices'))
  dir.create(path = file.path(OUTPUT,'matrices/sparse'))
  dir.create(path = file.path(OUTPUT,'figures'))
  dir.create(path = file.path(OUTPUT,'markers'))
  dir.create(path = file.path(OUTPUT,'rds'))
  dir.create(path = file.path(OUTPUT,'metadata'))
  
  
  
  
  
}

#Configuration file



{
  
  
  
  conf_file <- read.csv(
    file = file.path(getwd(), 'requirements_file/config_file.conf'),
    header = FALSE,
    sep = ':',
    row.names = 1,
    comment.char = "#"
  )
  
  
  mt_per <- as.numeric(as.character(conf_file$V2[grep(pattern = 'mt_per', rownames(conf_file))]))
  
  down_tr <- as.numeric(as.character(conf_file$V2[grep(pattern = 'down', rownames(conf_file))]))
  
  up_tr <- as.numeric(as.character(conf_file$V2[grep(pattern = 'up', rownames(conf_file))]))
  
  mt_cssg <- as.logical(conf_file$V2[grep(pattern = 'mt_cssg', rownames(conf_file))])
  
  s_factor <- as.numeric(as.character(conf_file$V2[grep(pattern = 's_factor', rownames(conf_file))]))
  
  m_val <- as.numeric(as.character(conf_file$V2[grep(pattern = 'm_val', rownames(conf_file))]))
  
  max_genes <- as.numeric(as.character(conf_file$V2[grep(pattern = 'max_genes', rownames(conf_file))]))
  
  max_combine <- as.numeric(as.character(conf_file$V2[grep(pattern = 'max_combine', rownames(conf_file))]))
  
  loss_val <- as.numeric(as.character(conf_file$V2[grep(pattern = 'loss_val', rownames(conf_file))]))
  
  p_bin <- as.numeric(as.character(conf_file$V2[grep(pattern = 'p_bin', rownames(conf_file))]))
  
  scale_factor <- as.numeric(as.character(conf_file$V2[grep(pattern = 'scale_factor', rownames(conf_file))]))
  
  n_features <- as.numeric(as.character(conf_file$V2[grep(pattern = 'n_features', rownames(conf_file))]))
  
  heterogeneity <- as.character(as.character(conf_file$V2[grep(pattern = 'heterogeneity', rownames(conf_file))]))
  
  drop_sub <- as.logical(as.character(conf_file$V2[grep(pattern = 'drop', rownames(conf_file))]))
  
  min_c <- as.numeric(as.character(conf_file$V2[grep(pattern = 'min_c', rownames(conf_file))]))
  
  c_res <- as.numeric(as.character(conf_file$V2[grep(pattern = 'c_res', rownames(conf_file))]))
  
  top_m <- as.numeric(as.character(conf_file$V2[grep(pattern = 'top_m', rownames(conf_file))]))
  
  
  if(length(data) == 0) {data = 3}
  
}

###########################################################################################################################################################

markers_class <- readxl::read_xlsx(markers, sheet = 1)
markers_subclass <- readxl::read_xlsx(markers, sheet = 2, col_names = F)

###########################################################################################################################################################

{
  # Load the raw dataset by UMI
  UMI_raw <- Read10X(seurat_umi, gene.column = 1)
  
  #Create SeuratObject
  UMI <- CreateSeuratObject(counts = UMI_raw, project = project_name, min.cells = 1, min.features = 1)
  
  cell_input <- ncol(UMI)
}


UMI@meta.data$orig.ident  <- make.unique(as.character(names(Idents(UMI))))

###########################################################################################################################################################


UMI[["MitoPercent"]] <- PercentageFeatureSet(UMI, pattern = "^(MT-|MT\\.)")
UMI[["RiboPercent"]] <- PercentageFeatureSet(UMI, pattern = "^(RPS|RPL|MRPL|MRPS|RS-)")


UMI@meta.data <- UMI@meta.data %>%
  rename(nCounts = nCount_RNA) %>%
  rename(nGenes = nFeature_RNA)


#Graphs of counts content

UC_plot <- VlnPlot(UMI, features = c("nGenes", "nCounts"), ncol = 2)
svg(file.path(OUTPUT, "figures/counts~genes.svg"), width=15, height=10)
UC_plot
rm(UC_plot)
dev.off()


CG_plot <- FeatureScatter(UMI, feature1 = "nCounts", feature2 = "nGenes")
svg(file.path(OUTPUT, "figures/counts~genes_QC.svg"), width=15, height=10)
CG_plot
rm(CG_plot)
dev.off()


if (!any(is.na(UMI@meta.data$MitoPercent)) && !any(is.na(UMI@meta.data$RiboPercent))) {
  
  MR_plot <- VlnPlot(UMI, features = c("RiboPercent", "MitoPercent"), ncol = 2)
  svg(file.path(OUTPUT, "figures/Ribo~Mito.svg"), width=15, height=10)
  print(MR_plot)
  rm(MR_plot)
  dev.off()
  
}



###########################################################################################################################################################

#Droplet content and QC

thresholds <- CSSG.toolkit::outlires(UMI@meta.data$nGenes)


if (is.na(up_tr)) {
  n_gen <- thresholds$thresholds[length(thresholds$thresholds)]
} else {
  n_gen <- up_tr
}


if (is.na(down_tr)) {
  down_tr <- thresholds$thresholds[1]
} else {
  down_tr <- down_tr
}


QC_UMI <- data.frame()
QC_UMI <- as.data.frame(UMI$nGenes)


colnames(QC_UMI) <- c('nGenes')

QC_UMI$nGenes_Status[UMI$nGenes < down_tr] <- 'Empty'
QC_UMI$nGenes_Status[UMI$nGenes > n_gen] <- 'Double'
QC_UMI$nGenes_Status[UMI$nGenes >= down_tr & UMI$nGenes <= n_gen] <- 'Proper'




DQC <- ggplot()+
  geom_point(QC_UMI, mapping = aes(x = nGenes, y = nGenes, color = nGenes_Status))+
  ylab("Number of genes for each cell") +
  xlab("Number of genes for each cell")+
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))+
  geom_vline(xintercept = down_tr, color = 'green') + annotate("text", x=down_tr - 40 , y=as.integer(round(max(QC_UMI$nGenes) / 1.3)), label= as.character(paste0('> ', as.character(round(down_tr)))) , angle=90) +
  geom_vline(xintercept = n_gen, color = 'red') + annotate("text", x=n_gen - 40  , y=as.integer(round(max(QC_UMI$nGenes) / 1.3)), label= as.character(paste0('< ', as.character(round(n_gen)))), angle=90) +
  labs(color='Droplet content') +
  theme_classic()


svg(filename = file.path(OUTPUT,'figures/DropletQC.svg'), width = 10, height = 7)
DQC
dev.off()
rm(DQC)


svg(filename = file.path(OUTPUT,'figures/DropletQC_hist.svg'), width = 10, height = 7)
thresholds$plot
dev.off()


rm(QC_UMI)
rm(thresholds)


###########################################################################################################################################################

#Selecting right cells

if (!any(is.na(UMI@meta.data$MitoPercent))) {
  
  UMI <- subset(UMI, subset = nGenes > down_tr & nGenes <= n_gen & MitoPercent < mt_per)
  
} else {
  
  UMI <- subset(UMI, subset = nGenes > down_tr & nGenes <= n_gen)
  
  
}


qc_cells <- ncol(UMI)




###########################################################################################################################################################

if (data == 2) {
  UMI@assays$RNA@data <- UMI@assays$RNA@counts
} else {
  UMI <- NormalizeData(UMI, normalization.method = "LogNormalize", scale.factor = scale_factor)
}


###########################################################################################################################################################

UMI <- FindVariableFeatures(UMI, selection.method = "vst", nfeatures = n_features, binning.method = 'equal_frequency')

# Identify the 10 most highly variable genes

top20 <- head(VariableFeatures(UMI), 20)

plot1 <- VariableFeaturePlot(UMI)

plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)

svg(file.path(OUTPUT, "figures/variable_genes.svg"), width=10, height=7)
plot2
dev.off()
rm(plot2)

###########################################################################################################################################################

all.genes <- rownames(UMI)
UMI <- ScaleData(UMI, features = all.genes)


var_features = VariableFeatures(object = UMI)

if (mt_cssg == F) {
  var_features = var_features[!grepl("^(MT-|MT\\.)", toupper(var_features))]
}

UMI <- RunPCA(UMI, features = VariableFeatures(object = UMI))

###########################################################################################################################################################


Elbow <- ElbowPlot(UMI, ndims = 50)

dims <- as.data.frame(Elbow$data$stdev)

#select the most variable reduction

dim <- CSSG.toolkit::dim_reuction_pcs(dims)


svg(file.path(OUTPUT, "figures/Elbow.svg"), width=10, height=7)
Elbow <- Elbow + geom_vline(xintercept = dim, color = 'red') +   
  geom_text(aes(x = dim + 3, y = round(max(Elbow$data$stdev)/1.5,0), label = paste("Dim =", dim)), color = 'red', vjust = -1)
Elbow
dev.off()
rm(Elbow)

###########################################################################################################################################################

UMI <- JackStraw(UMI, num.replicate = 10, dims = dim)
UMI <- ScoreJackStraw(UMI, dims = 1:dim)

#Select significient PCs
jc <- as.data.frame(UMI@reductions$pca@jackstraw@overall.p.values)
jc <- jc[jc$Score < 0.05,]
dim <- as.vector(jc$PC)


svg(file.path(OUTPUT, "figures/JackStrawPlot.svg"), width=10, height=7)
JackStrawPlot(UMI, dims = dim)
dev.off()

UMI <- FindNeighbors(UMI, dims = dim, reduction = 'pca')


UMI <- FindClusters(UMI, resolution = c_res, n.start = 10, n.iter = 1000)


UMI <- RunUMAP(UMI, dims = dim, umap.method = "umap-learn")

#CHANGE


width <- 10 + (length(unique(Idents(UMI))))/5

svg(file.path(OUTPUT, "figures/UMAP_clusters.svg"), width = width, height = 10)
DimPlot(UMI, reduction = "umap", raster = FALSE)
dev.off()

###########################################################################################################################################################


########## METADATA STEP

meta_data <- as.data.frame(Idents(UMI))
colnames(meta_data)[1] <- 'clusters'
meta_data$barcodes <- as.character(rownames(meta_data))
UMAP_coordinates <- as.data.frame(Embeddings(UMI, reduction = "umap"))
meta_data$UMAP1 <- as.numeric(UMAP_coordinates[,1])
meta_data$UMAP2 <- as.numeric(UMAP_coordinates[,2])

#########


#find markers for every cluster compared to all remaining cells, report only the positive ones

print('Searching for clusters / subclasses marker genes')

# clusters /  subclasses





sc_project <- CSSG.toolkit::create_project_from_seurat(UMI)

# get markers  for subclass naming

print('Clusters stats')

sc_project <- CSSG.toolkit::get_cluster_stats(sc_project = sc_project, type = 'primary', only_pos = TRUE, min_pct = 0.05)


# elect markers for subclass naming

print('Clusters naming - subclasses')

sc_project <- CSSG.toolkit::namign_genes_selection(sc_project, type = 'primary', top_n = top_m,
                                                   p_val = m_val, select_stat = "p_val",
                                                   mito_content = FALSE, ribo_content = FALSE)



sc_project <- CSSG.toolkit::subclass_naming(sc_project = sc_project, 
                                            class_markers = markers_class, 
                                            subclass_markers = markers_subclass, 
                                            species = species, 
                                            chunk_size = 5000)

# unique(sc_project@names$subclass)

print('Composition binomial - subclasses')

# print(sc_project@names$subclass) # debug

data <- CSSG.toolkit::bin_cell_test(p_val = p_bin, names = sc_project@names$subclass, min_cells = min_c)

# print(data) # debug

print('Composition - thresholds')

# print(data$data) # debug

threshold <- CSSG.toolkit::cell_stat_graph(data$data, include_ns = TRUE)

height <- 4 + (length(unique(threshold$data$names)))/5
svg(file.path(OUTPUT, "figures/subclasses_composition.svg"), width = 7, height = height)
threshold
dev.off()


data_avg <- CSSG.toolkit::get_avg_data(sc_project = sc_project, type = 'subclasses', data = 'norm')

# print(data_avg) # debug

print('Save - subclasses avg')

write.table(data_avg, file = file.path(OUTPUT, "matrices/subclasses_average_expression.csv"), sep = ',')

rm(data_avg)
rm(data)
rm(threshold)

################################################################################

print('Get names markers - subclasses')


# print(sc_project@names$subclass) # debug


markers <- CSSG.toolkit::get_names_markers(sc_project, type = 'subclasses')


print('Subclasses heatmaps')

height <- 10 + (length(unique(markers)))/5
width <- 15 + (length(unique(Idents(UMI))))/5



svg(file.path(OUTPUT, "figures/heatmap_cells_subclasses.svg"), width = width, height = height)
CSSG.toolkit::marker_heatmap(sc_project, type = 'subclasses', markers = markers,
                             angle_col = 270,
                             fontsize_row = 15,
                             fontsize_col = 15,
                             font_labels = 15,
                             clustering_method = 'complete',
                             x_axis = 'Cells',
                             y_axis = 'Genes [log(CPM +1)]',
                             scale = FALSE)

dev.off()




svg(file.path(OUTPUT, "figures/heatmap_cells_subclasses_scaled.svg"), width = width, height = height)
CSSG.toolkit::marker_heatmap(sc_project, type = 'subclasses', markers = markers,
                             angle_col = 270,
                             fontsize_row = 15,
                             fontsize_col = 15,
                             font_labels = 15,
                             clustering_method = 'complete',
                             x_axis = 'Cells',
                             y_axis = 'Genes scaled([log(CPM +1)])',
                             scale = TRUE)
dev.off()





################################################################################


print('Subclasses renaming')

Idents(UMI) <- sc_project@names$subclass


width <- 15 + (length(unique(Idents(UMI))))/5


svg(file.path(OUTPUT, "figures/PCA_DimPlot_subclasses.svg"), width = width, height = 10)
DimPlot(UMI, reduction = "pca", raster = FALSE)
dev.off()


svg(file.path(OUTPUT, "figures/UMAP_DimPlot_subclasses.svg"), width = width, height = 10)
DimPlot(UMI, reduction = "umap", raster = FALSE)
dev.off()




print('Searching for clusters / subclasses marker genes - DONE')


################################################################################

print('Searching for CSSG markers')

if (tolower(heterogeneity) == 'var') {
  
  sc_project <- CSSG.toolkit::heterogeneity_select_variance(sc_project = sc_project,
                                                            heterogeneity_factor = s_factor,
                                                            max_genes = max_genes,
                                                            min_occ = 5,
                                                            min_exp = 0.01,
                                                            rep_factor = 0.2,
                                                            mito_content = mt_cssg)
  
} else {
  
  sc_project <- CSSG.toolkit::heterogeneity_select_specificity(sc_project = sc_project,
                                                               type = 'primary',
                                                               heterogeneity_factor = s_factor,
                                                               p_val =  m_val,
                                                               max_genes =  max_genes,
                                                               select_stat = 'p_val',
                                                               min_occ = 5,
                                                               mito_content = mt_cssg)
}

sc_project <- CSSG.toolkit::CSSG_markers(sc_project = sc_project,
                                         max_combine = max_combine,
                                         loss_val = loss_val)



marker_to_write <- sc_project@metadata$primary_markers

mapper <- data.frame(
  primary  = sc_project@names$primary,
  subclass = sc_project@names$subclass,
  stringsAsFactors = FALSE
)

mapper <- distinct(mapper)

marker_to_write$subclass <- mapper$subclass[ match(marker_to_write$cluster, mapper$primary) ]

write.table(marker_to_write, file = file.path(OUTPUT, "markers/markers_subclasses.csv"), sep = ',')


rm(marker_to_write)
rm(mapper)



print('Searching for CSSG markers - DONE')


print('Subtypes naming')

# subtypes
sc_project <- CSSG.toolkit::subtypes_naming(sc_project = sc_project, markers_class = markers_class, markers_subclass = markers_subclass,  species = species)


thr_data <- CSSG.toolkit::bin_cell_test(p_val = p_bin,
                                        names = sc_project@names$repaired,
                                        min_cells = min_c)

threshold <- CSSG.toolkit::cell_stat_graph(thr_data$data, include_ns = TRUE)


height <- 7 + (length(unique(threshold$data$names)))/5
svg(file.path(OUTPUT, "figures/subtypes_composition.svg"), width = 10, height = height)
threshold
dev.off()

# metadata fill

meta_data$subclass <- as.character(sc_project@names$subclass)
meta_data$subtypes <- as.character(sc_project@names$subtypes)
meta_data$subtypes[grepl('BAD!', toupper(meta_data$subtypes))] <- 'Undefined'


write.table(sc_project@metadata$cssg_markers, file = file.path(OUTPUT, "markers/CSSG_marker.csv"), sep = ',')

###########################################################################################################################################################

# subsetns analysis

print('Renames subtypes')

Idents(UMI) <- sc_project@names$subtypes


if (drop_sub) {
  
  select_list <- thr_data$data$names[thr_data$data$test %in% c("Good marked types", "Renamed")]
  
} else {
  
  select_list <- thr_data$data$names[thr_data$data$test %in% c("Good marked types", "Renamed", "Non-significant")]
  
}

print('Select list set')

meta_data$reduced <- TRUE
meta_data$reduced[meta_data$subtypes %in% select_list] <- FALSE

centroids <- meta_data %>%
  group_by(subtypes) %>%
  summarise(
    umap1_centroid = mean(UMAP1),
    umap2_centroid = mean(UMAP2)
  )


meta_data <- meta_data %>%
  left_join(centroids, by = "subtypes")


factor <- 0.6

meta_data <- meta_data %>%
  mutate(
    fUMAP1 = UMAP1 + factor * (umap1_centroid - UMAP1),
    fUMAP2 = UMAP2 + factor * (umap2_centroid - UMAP2)
  )



meta_data <- meta_data %>%
  group_by(subtypes) %>%
  mutate(n_subtype = n()) %>%
  ungroup()

meta_data <- meta_data %>%
  arrange(desc(n_subtype))


meta_data <- meta_data[,!colnames(meta_data) %in% c('umap1_centroid', 'umap2_centroid')]


print('Metadata writing')

write.table(meta_data, file = file.path(OUTPUT, "metadata/metadata.csv"), sep = ',')



meta_data_plot <- meta_data[meta_data$subtypes %in% select_list,]



print('Subtypes UMAP')

plot <- ggplot(meta_data_plot, aes(x = fUMAP1, y = fUMAP2, color = subtypes)) +
  geom_point(size = 0.8, alpha = 0.8) +
  xlab('fUMAP_1') +
  ylab('fUMAP_2') +
  labs(colour = "Cell subtypes") +
  theme_classic()



htmlwidgets::saveWidget(plotly::ggplotly(plot) , file.path(OUTPUT, "figures/UMAP_subtypes.html"))



width <- 15 + (length(unique(Idents(UMI))))/5


svg(file.path(OUTPUT, "figures/UMAP_subtypes.svg"), width = width, height = 15)
plot
dev.off()


################################ Reducing

print('Projects subsetting')

UMI <- subset(UMI, idents = select_list)
sc_project <- CSSG.toolkit::subset_project(sc_project = sc_project, type = 'subtypes', select_list = select_list)




markers <- CSSG.toolkit::get_names_markers(sc_project, type = 'subtypes')

height <- 10 + (length(unique(markers)))/5
width <- 15 + (length(unique(Idents(UMI))))/5


print('Subtypes heatmaps')


svg(file.path(OUTPUT, "figures/heatmap_cells_subtypes.svg"), width = width, height = height)
CSSG.toolkit::marker_heatmap(sc_project, type = 'subtypes', markers = markers,
                             angle_col = 270,
                             fontsize_row = 15,
                             fontsize_col = 15,
                             font_labels = 15,
                             clustering_method = 'complete',
                             x_axis = 'Cells',
                             y_axis = 'Genes [log(CPM +1)]',
                             scale = FALSE)

dev.off()




svg(file.path(OUTPUT, "figures/heatmap_cells_subtypes_scaled.svg"), width = width, height = height)
CSSG.toolkit::marker_heatmap(sc_project, type = 'subtypes', markers = markers,
                             angle_col = 270,
                             fontsize_row = 15,
                             fontsize_col = 15,
                             font_labels = 15,
                             clustering_method = 'complete',
                             x_axis = 'Cells',
                             y_axis = 'Genes scaled([log(CPM +1)])',
                             scale = TRUE)
dev.off()




data_avg <- CSSG.toolkit::get_avg_data(sc_project = sc_project, type = 'subtypes', data = 'norm')

write.table(data_avg, file = file.path(OUTPUT, "matrices/sybtypes_average_expression.csv"), sep = ',')




print('Subtypes naming - DONE')



###########################################################################################################################################################
#Subtype markers selection


sc_project <- CSSG.toolkit::get_cluster_stats(sc_project = sc_project, type = 'subtypes', only_pos = TRUE)



marker_to_write <- sc_project@metadata$subtypes_markers
colnames(marker_to_write)[colnames(marker_to_write) == 'cluster'] <- 'subtype'

mapper <- data.frame(
  primary  = sc_project@names$primary,
  subtype = sc_project@names$subtypes,
  stringsAsFactors = FALSE
)

mapper <- distinct(mapper)

marker_to_write$cluster <- mapper$primary[match(marker_to_write$subtype, mapper$subtype) ]
write.table(marker_to_write, file = file.path(OUTPUT, "markers/markers_subtypes.csv"), sep = ',')

rm(marker_to_write)
rm(mapper)

############################################################################################################################################################

#Cells_stats

cells <- factor(c('Estimated_cells', 'Input_cells', 'Analyzed_cells', 'Output_cells'), levels = c('Estimated_cells', 'Input_cells', 'Analyzed_cells', 'Output_cells'))
cell_num <- c(as.numeric(estimated_cells), as.numeric(cell_input), as.numeric(qc_cells), ncol(UMI))


df_cells <- data.frame(cells, cell_num)

cells <- ggplot(df_cells, aes(x = cells, y = cell_num, fill = cells)) +
  geom_col() +
  ylab("Number of cells") +
  xlab("Cells in analysis")+
  geom_text(aes(label = cell_num), vjust = -0.5)+
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))+
  theme_classic() +
  theme(legend.position = 'none')


svg(filename = file.path(OUTPUT,'figures/Cells.svg'), width = 10, height = 7)
cells
dev.off()
rm(cells)



###########################################################################################################################################################

print('Matrix saving ')


exp_matrix <- GetAssayData(UMI, slot = 'data')
colnames(exp_matrix) <- UMI@active.ident

write(colnames(exp_matrix), file = file.path(OUTPUT, "matrices/sparse/barcodes.tsv"))
write(rownames(exp_matrix), file = file.path(OUTPUT, "matrices/sparse/genes.tsv"))
Matrix::writeMM(exp_matrix, file = file.path(OUTPUT, "matrices/sparse/matrix.mtx"))

print('Matrix saving - DONE')



print('Project saving ')


#save seurat output file

saveRDS(UMI, file = file.path(OUTPUT, "rds/Results.rds"))

print('Project saving - DONE')



###########################################################################################################################################################

print('Report creating')


rmarkdown::render(input = file.path(getwd(), 'scripts/report_species.Rmd'),
                  output_format = 'html_document', output_dir = OUTPUT,
                  output_file = 'Report')
