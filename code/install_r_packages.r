install.packages("devtools")
require("devtools")

install.packages("BiocManager", repos = "http://cran.us.r-project.org")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GSVA")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mixOmics")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("fgsea")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")