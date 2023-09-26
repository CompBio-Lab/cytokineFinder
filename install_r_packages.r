install.packages("devtools")
require("devtools")

install.packages("BiocManager", repos = "http://cran.us.r-project.org")

data <- c("GEOquery")
gsva <- c("GSVA")
eda <- c("mixOmics")
limma <- c("limma")
fgsea <- c("fgsea")
bioconductor_packages <- c(data, gsva, eda, limma, fgsea)
BiocManager::install(bioconductor_packages, version = "3.17")

install.packages("parallel")