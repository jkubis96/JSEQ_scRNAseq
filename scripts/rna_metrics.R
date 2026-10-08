args <- commandArgs()
path <- args[6]
results <- args[7]

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
library(viridis)
debug_flag <- FALSE



#### /debug

metrics <- file.path(path, "/scRNAmetrics.txt")

mydata <- read.csv(
  file = metrics, header = T,
  stringsAsFactors = F, skip = 6, sep = "\t"
)
mydata <- mydata[order(mydata$PF_ALIGNED_BASES, decreasing = T), ]
mydata_pct <- mydata[, c(
  "READ_GROUP",
  "PCT_INTERGENIC_BASES",
  "PCT_UTR_BASES",
  "PCT_RIBOSOMAL_BASES",
  "PCT_INTRONIC_BASES",
  "PCT_CODING_BASES"
)]
colnames(mydata_pct) <- c("Cell Barcode", "Intergenic", "UTR", "Ribosomial", "Intronic", "Coding")


mydata_long_pct <- mydata_pct %>% gather("Read Overlap", fraction, -"Cell Barcode")
mydata_long_pct$`Cell Barcode` <- factor(mydata_long_pct$`Cell Barcode`,
  levels = factor(unique(mydata_long_pct$`Cell Barcode`))
)
mydata_long_pct$`Read Overlap` <- factor(mydata_long_pct$`Read Overlap`,
  levels = unique(mydata_long_pct$`Read Overlap`)
)

p <- ggplot(mydata_long_pct, aes(x = `Cell Barcode`, y = fraction, fill = `Read Overlap`)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, size = 8, vjust = 0.05), legend.position = "bottom") +
  labs(x = "Barcodes", y = "%Bases") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis(discrete = TRUE, option = "viridis")


ggsave(filename = file.path(results, "scRNAmetrics.jpeg"), plot = p, dpi = 600)
