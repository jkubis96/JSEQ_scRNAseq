#' ---
#' title:  plot_violine.R
#' author: Sebastian Mueller (sebm_at_posteo.de)
#' date:   2018-04-10
#' ---
### for debug
# If you wish to access the snakefile object first invoke snakemake and save the session automatically
# Since there are no debug flags to my knowledge, just uncomment the line below and run snakemake which
# creates an R object that can be loaded into a custom R session
# save.image(file="R_workspace_debug.rdata")
# load("R_workspace_debug.rdata")
#### /debug



options(warn = -1)
library(plyr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE) # Dataframe manipulation
library(Matrix, quietly = TRUE, warn.conflicts = FALSE) # Sparse matrices
library(stringr, quietly = TRUE, warn.conflicts = FALSE)
library(RColorBrewer, quietly = TRUE, warn.conflicts = FALSE)
library(devtools, quietly = TRUE, warn.conflicts = FALSE)
library(Seurat, quietly = TRUE, warn.conflicts = FALSE)
library(plotly, quietly = TRUE, warn.conflicts = FALSE)


count_matrix <- Read10X(data.dir = file.path('../MappingResultsHuman/read/'), gene.column = 1)
umi_matrix <- Read10X(file.path('../MappingResultsHuman/umi/'),gene.column = 1)

input <- c('sample1', '270', '74', 'batch1')
design <- data.frame(input, stringsAsFactors = T)
t <- t(as.matrix(design))
design <- as.data.frame(t)
colnames(design) <- c('samples','expected_cells','read_length','batch')
rownames(design) <- NULL


metaData <- data.frame(cellNames = colnames(umi_matrix)) %>% 
  mutate(samples = factor(str_replace(cellNames, "_[^_]*$", ""))) %>%
  mutate(barcode = factor(str_replace(cellNames, ".+_", "")))

metaData$samples <- factor(paste(design$batch,metaData$samples))
metaData$cellNames <- factor(paste(design$samples,metaData$cellNames))

rownames(metaData) <- metaData$cellNames

seuratobj <- CreateSeuratObject(counts = umi_matrix, meta.data = metaData)
seuratobj <- SetAllIdent(object = seuratobj, id = "samples")
# relabel cell idenity (https://github.com/satijalab/seurat/issues/380)
seuratobj@meta.data$orig.ident <- seuratobj@meta.data$samples

mycount <- CreateSeuratObject(counts = count_matrix, meta.data = metaData)
mycount@meta.data$orig.ident <- mycount@meta.data$samples
meta.data <- seuratobj@meta.data
meta.data$nCounts <- mycount@meta.data$nUMI
seuratobj@meta.data <- meta.data
rm(mycount)


meta.data <- meta.data %>% 
  rename(nCounts = nCount_RNA) %>% 
  rename(nUMI = nFeature_RNA)
  
  

# mytheme <- theme_bw(base_size = 9) +
mytheme <- theme_bw() +
  theme(
    legend.position = "right",
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 300, hjust = 0)
  )
theme_set(mytheme)

# predefined ggplot layers for subsequent plots
gglayers <- list(
  geom_smooth(method = "loess"),
  geom_point(size = .5),
  scale_y_continuous(
    labels = scales::unit_format(unit = "", scale = 1e-3, digits = 2),
    breaks = scales::pretty_breaks(n = 8)
  ),
  scale_x_continuous(
    labels = scales::unit_format(unit = "", scale = 1e-3, digits = 2),
    breaks = scales::pretty_breaks(n = 8)
  )
)

gg <- ggplot(meta.data, aes(x = nUMI, y = nCounts, color = orig.ident)) +
  #   coord_trans(y="log10",x = "log10") +
  gglayers +
  geom_abline(intercept = 0, slope = 1) +
  labs(
    title = "UMI counts vs raw Counts",
    subtitle = "Number of UMIs and raw Counts for each Bead",
    x = "Number of UMIs per Bead [k]",
    y = "Number of Counts per Bead [k]"
  )

# dev.new()
# htmlwidgets::saveWidget(ggplotly(gg), file.path(getwd(),snakemake@output$html_umivscounts))
ggsave(gg, file = file.path(getwd(), snakemake@output$pdf_umivscounts), width = 12, height = 7)
################################TUTAJ SKONCZYLEM##############################





# how about unaligned reads/UMI?
# Note(Seb): raw.data is actually filtered data i.e. nr of genes likely to be smaller than input data!
mito.gene.names <- grep("^mt-", rownames(seuratobj@meta.data), value = TRUE, ignore.case = TRUE)
sribo.gene.names <- grep("^Rps", rownames(seuratobj@meta.data), value = TRUE, ignore.case = TRUE)
lribo.gene.names <- grep("^Rpl", rownames(seuratobj@meta.data), value = TRUE, ignore.case = TRUE)

col.total <- Matrix::colSums(seuratobj@meta.data)
meta.data$col.total <- col.total

seuratobj.top_50 <- apply(seuratobj@raw.data, 2, function(x) sum(x[order(x, decreasing = TRUE)][1:50]) / sum(x))
# mycount.top_50 <- apply(mycount@raw.data, 2, function(x) sum(x[order(x, decreasing = TRUE)][1:50])/sum(x))

seuratobj <- AddMetaData(seuratobj, Matrix::colSums(seuratobj@raw.data[sribo.gene.names, ]) / col.total, "pct.sribo")
seuratobj <- AddMetaData(seuratobj, Matrix::colSums(seuratobj@raw.data[lribo.gene.names, ]) / col.total, "pct.lribo")
seuratobj <- AddMetaData(seuratobj, Matrix::colSums(seuratobj@raw.data[unique(c(sribo.gene.names, lribo.gene.names)), ]) / col.total, "pct.Ribo")
seuratobj <- AddMetaData(seuratobj, Matrix::colSums(seuratobj@raw.data[mito.gene.names, ]) / col.total, "pct.mito")
seuratobj <- AddMetaData(seuratobj, seuratobj.top_50, "top50")
tmp <- seuratobj@meta.data$nUMI / seuratobj@meta.data$nGene
names(tmp) <- rownames(seuratobj@meta.data)
seuratobj <- AddMetaData(seuratobj, tmp, "umi.per.gene")


gg <- VlnPlot(seuratobj,
  c("nUMI", "nGene", "top50", "umi.per.gene", "pct.Ribo", "pct.mito"),
  x.lab.rot = TRUE, do.return = TRUE
)
# ggsave(gg,file=file.path("violinplots_comparison_UMI.pdf"),width=18,height=18)
ggsave(gg, file = snakemake@output$pdf_violine, width = 18, height = 18)
# gg <- VlnPlot(mycount,c("nUMI", "nGene", "top50", "count.per.gene","pct.Ribo", "pct.mito"), x.lab.rot = TRUE, do.return = TRUE)
# ggsave(gg,file=file.path("violinplots_comparison_count.pdf"),width=18,height=18)

# gg <- GenePlot(object = seuratobj, gene1 = "nUMI", gene2 = "nGene")
# ggsave(gg,file=file.path("violinplots_comparison.pdf"),width=18,height=18)


gg <- ggplot(meta.data, aes(x = nUMI, y = nGene, color = orig.ident)) +
  gglayers +
  labs(
    title = "Genes (pooled mouse and human set) vs UMIs for each bead",
    x = "Number of UMIs per Bead [k]",
    y = "Number of Genes per Bead [k]"
  )

# dev.new()
# htmlwidgets::saveWidget(ggplotly(gg),
# file.path(getwd(), snakemake@output$html_umi_vs_gene))
ggsave(gg, file = snakemake@output$pdf_umi_vs_gene, width = 12, height = 7)



################################################################################
## same for Counts instead UMIs (using mycount object)
gg <- ggplot(meta.data, aes(x = nCounts, y = nGene, color = orig.ident)) +
  gglayers +
  labs(
    title = "Genes (pooled mouse and human set) vs Counts for each bead",
    x = "Number of Counts per Bead [k]",
    y = "Number of Genes per Bead [k]"
  )

# dev.new()
# htmlwidgets::saveWidget(ggplotly(gg),
#                       file.path(getwd(), snakemake@output$html_count_vs_gene))

ggsave(gg, file = snakemake@output$pdf_count_vs_gene, width = 12, height = 7)


# head(meta.data,2)
#                              nGene nUMI                    cellNames         samples      barcode expected_cells read_length  batch      orig.ident pct.sribo  pct.lribo  pct.Ribo  pct.mito     top50 umi.per.gene
# sample1_GAGTCTGAGGCG     6    6 sample1_GAGTCTGAGGCG sample1 GAGTCTGAGGCG            100         100 batch1 sample1 0.0000000 0.00000000 0.0000000 0.0000000 1.0000000     1.000000
# sample1_CAGCCCTCAGTA   264  437 sample1_CAGCCCTCAGTA sample1 CAGCCCTCAGTA            100         100 batch1 sample1 0.0389016 0.07551487 0.1144165 0.0228833 0.5102975     1.655303

# saving snakemake meta information into misc slot so all can be exported as one object
seuratobj@misc <- snakemake
# exporting R Seurat objects into summary/R_Seurat_objects.rdata
saveRDS(seuratobj, file = file.path(snakemake@output$R_objects))

if (debug_flag) {
  save.image(file = file.path(path_debug, "plot_violin_workspace.rdata"))
}
