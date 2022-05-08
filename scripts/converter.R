library(Matrix) 

#convert gene count matrix to seurat obj

args <- commandArgs()

format <- args[6]
path <- args[7]
input_name <- args[8]
tmp <- args[9]
con <- args[10]

if(con == 'raw') {

m_counts <-  read.csv(file = file.path(tmp, 'umi_expression.tsv'), header = T, row.names = 1, sep = '\t')

m_counts <- as.matrix(m_counts)
sparse_m <- Matrix(m_counts , sparse = T )
writeMM(obj = sparse_m, file=file.path(path,"matrix.mtx"))

genes <- rownames(m_counts)
genes <- gsub(' ', '', genes)
genes <- make.unique(genes)
write.table(x = as.data.frame(genes), file = file.path(path, 'genes.tsv'), row.names = F, col.names = F)
barcodes <- colnames(m_counts)
barcodes <- gsub(' ', '', barcodes)
write.table(x = as.data.frame(barcodes), file = file.path(path, 'barcodes.tsv'), row.names = F, col.names = F)

} else {


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
genes <- gsub(' ', '', genes)
genes <- make.unique(genes)
write.table(x = as.data.frame(genes), file = file.path(path, 'genes.tsv'), row.names = F, col.names = F)
barcodes <- colnames(m_counts)
barcodes <- gsub(' ', '', barcodes)
write.table(x = as.data.frame(barcodes), file = file.path(path, 'barcodes.tsv'), row.names = F, col.names = F)
  
}
