library(Matrix)

#convert gene count matric to seurat obj

args <- commandArgs()

format <- args[6]
path <- args[7]
input_name <- args[8]


if (format == 'txt') {
  m_counts <-  read.csv(file = file.path(path, input_name), header = T, row.names = 1, sep = '\t')
  }else if (format == 'tsv') {
    m_counts <-  read.csv(file = file.path(path, input_name), header = T, row.names = 1, sep = '\t')
   }else if (format == 'csv') {
      m_counts <-  read.csv(file = file.path(path, input_name), header = T, row.names = 1)
}

m_counts <- as.matrix(m_counts)
sparse_m <- Matrix(m_counts , sparse = T )
writeMM(obj = sparse_m, file=file.path(path,"matrix.mtx"))

genes <- rownames(m_counts)
write.table(x = as.data.frame(genes), file = file.path(path, 'genes.tsv'), row.names = F, col.names = F)
barcodes <- colnames(m_counts)
write.table(x = as.data.frame(barcodes), file = file.path(path, 'barcodes.tsv'), row.names = F, col.names = F)
  
