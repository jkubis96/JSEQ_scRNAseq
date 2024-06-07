Sys.setenv(R_INSTALL_STAGED = FALSE)

install.packages('patchwork', version = '1.1.1')
# install.packages('viridis', version = "0.5.1", dependencies = FALSE, ask = FALSE)
install.packages('textclean', version = "0.9.3")
# install.packages('ape', version = '5.4.1', dependencies = FALSE, ask = FALSE)
install.packages('umap')
BiocManager::install("MAST")


