library(Matrix)
library(data.table)
set.seed(123)


args <- commandArgs()

format <- args[6]
path <- args[7]
input_name <- args[8]
tmp <- args[9]
con <- args[10]

read_and_process <- function(filepath, sep = "\t") {
  m_counts <- fread(filepath, header = TRUE, sep = sep, data.table = FALSE)
  rownames(m_counts) <- m_counts[, 1]
  m_counts <- m_counts[, -1]
  m_counts <- as.matrix(m_counts)
  sparse_m <- Matrix(m_counts, sparse = TRUE)
}

write_output_files <- function(sparse_m, path) {
  writeMM(obj = sparse_m, file = file.path(path, "matrix.mtx"))

  genes <- rownames(sparse_m)
  genes <- gsub(" ", "", genes)
  genes <- make.unique(genes)
  fwrite(as.data.table(genes), file = file.path(path, "genes.tsv"), col.names = FALSE, sep = "\t")

  barcodes <- colnames(sparse_m)
  barcodes <- gsub(" ", "", barcodes)
  fwrite(as.data.table(barcodes), file = file.path(path, "barcodes.tsv"), col.names = FALSE, sep = "\t")
}

if (con == "raw") {
  sparse_m <- read_and_process(file.path(tmp, "umi_expression.tsv"), sep = "\t")
} else {
  if (format == "txt" || format == "tsv") {
    sparse_m <- read_and_process(file.path(path, input_name), sep = "\t")
  } else if (format == "csv") {
    sparse_m <- read_and_process(file.path(path, input_name), sep = ",")
  }
}

write_output_files(sparse_m, path)
