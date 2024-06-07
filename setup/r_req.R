 
Sys.setenv(R_INSTALL_STAGED = FALSE)

install.packages('patchwork', version = '1.1.1', dependencies = FALSE, ask = FALSE)
# install.packages('viridis', version = "0.5.1", dependencies = FALSE, ask = FALSE)
install.packages('textclean', version = "0.9.3", dependencies = FALSE, ask = FALSE)
# install.packages('ape', version = '5.4.1', dependencies = FALSE, ask = FALSE)
install.packages('umap', dependencies = FALSE, ask = FALSE)
BiocManager::install("MAST", ask = FALSE)


