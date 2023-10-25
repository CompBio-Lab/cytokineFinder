# CRAN
install.packages("tidyverse")
install.packages("here")
install.packages("ggplotify")
install.packages("ggpubr")
install.packages("pheatmap")

# Bioconductor
install.packages("BiocManager", repos = "http://cran.us.r-project.org")
data <- c("GEOquery")
gsva <- c("GSVA")
eda <- c("mixOmics")
limma <- c("limma")
fgsea <- c("fgsea")
bioconductor_packages <- c(data, gsva, eda, limma, fgsea)
BiocManager::install(bioconductor_packages)