 
Sys.setenv(R_INSTALL_STAGED = FALSE)

install.packages('remotes', dependencies = TRUE, ask = FALSE)
install.packages('patchwork', version = '1.1.1', dependencies = TRUE, ask = FALSE)
# remotes::install_version("Seurat", version = "3.1.5", dependencies = TRUE, ask = FALSE)
install.packages('viridis', version = "0.5.1", dependencies = TRUE, ask = FALSE)
install.packages('textclean', version = "0.9.3", dependencies = TRUE, ask = FALSE)
install.packages('ape', version = '5.4.1', dependencies = TRUE, ask = FALSE)
install.packages('umap', dependencies = TRUE, ask = FALSE)

BiocManager::install("MAST", ask = FALSE)


