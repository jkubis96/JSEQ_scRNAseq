#install.packages('ggplot2')

#install.packages('tidyr')

#install.packages('gridExtra')
 
#install.packages('grid')

#install.packages('viridis')

args <- commandArgs()
project_name_mode <- args[6]

library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(viridis)
debug_flag <- FALSE



#### /debug

metrics <- file.path(getwd(),'/projects/',project_name_mode,'/tmp/scRNAmetrics.txt')

mydata <- read.csv(file = metrics, header = T,
                   stringsAsFactors = F, skip = 6, sep = "\t")
mydata <- mydata[order(mydata$PF_ALIGNED_BASES, decreasing = T), ]
mydata_pct <- mydata[, c("READ_GROUP",
                         "PCT_INTERGENIC_BASES",
                         "PCT_UTR_BASES",
                         "PCT_RIBOSOMAL_BASES",
                         "PCT_INTRONIC_BASES",
                         "PCT_CODING_BASES")
                     ]
colnames(mydata_pct) = c('Cell Barcode', 'Intergenic', 'UTR', 'Ribosomial', 'Intronic', 'Coding')

mydata <- mydata[, c("READ_GROUP",
                     "INTERGENIC_BASES",
                     "UTR_BASES",
                     "RIBOSOMAL_BASES",
                     "INTRONIC_BASES",
                     "CODING_BASES")
                 ]
colnames(mydata) = c('Cell Barcode', 'Intergenic', 'UTR', 'Ribosomial', 'Intronic', 'Coding')

# converting into long format for ploting
mydata_long <- mydata %>% gather("Read Overlap", count, -"Cell Barcode")

# Keep the original order of the barcodes using factor and levels.
mydata_long$`Cell Barcode` <- factor(mydata_long$`Cell Barcode`,
                                 levels = factor(unique(mydata_long$`Cell Barcode`)))
mydata_long$`Read Overlap` <- factor(mydata_long$`Read Overlap`,
                                   levels = unique(mydata_long$`Read Overlap`))

p1 <- ggplot(mydata_long, aes(x = `Cell Barcode`, y = count, fill = `Read Overlap`)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0), legend.position = "none")
p1 <- p1 + labs(title = paste(nrow(mydata),
                              "selected barcodes for"
                              ),
                x = "Barcodes", y = "Bases")
p1 <- p1 + theme(axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank())
p1 <- p1 + scale_y_continuous(labels = scales::scientific)
p1 <- p1 + scale_fill_viridis(discrete = TRUE, option = "viridis")


mydata_long_pct <- mydata_pct %>% gather("Read Overlap", fraction, -"Cell Barcode")
# Keep the original order of the barcodes using factor and levels.
mydata_long_pct$`Cell Barcode` <- factor(mydata_long_pct$`Cell Barcode`,
                                     levels = factor(unique(mydata_long_pct$`Cell Barcode`)))
mydata_long_pct$`Read Overlap` <- factor(mydata_long_pct$`Read Overlap`,
                                       levels = unique(mydata_long_pct$`Read Overlap`))

p2 <- ggplot(mydata_long_pct, aes(x = `Cell Barcode`, y = fraction, fill = `Read Overlap`)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, size=8, vjust = 0.05), legend.position = "bottom") +
  labs(x = "Barcodes", y = "%Bases") +
  scale_y_continuous(labels = scales::percent) + scale_fill_viridis(discrete = TRUE, option = "viridis")
# This allows to align the main plots so that we can relate both directly with the label from the bottom one.
gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)
pdf(file = file.path(getwd(),'/projects/',project_name_mode,'/results/scRNAmetrics.pdf'), width = 16, height = 13)
grid::grid.newpage()
grid::grid.draw(rbind(gp1, gp2, size = "last"))
dev.off()

ggsave(filename = file.path(getwd(),'/projects/',project_name_mode,'/results/scRNAmetrics.jpeg'), dpi = 600)

