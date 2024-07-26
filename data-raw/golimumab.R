# Load libraries
library(cytokineFindeR)
library(tidyverse)

## code to prepare `golimumab` dataset goes here
dbs_all <- dbs_all
combined_data <- retrieve_geo("GSE92415")
dim(combined_data$metadata); dim(combined_data$annotations); dim(combined_data$eset);
all(rownames(combined_data$metadata) == colnames(combined_data$eset))

# in this specific dataset, I subset to just golimumab treated samples
golimumab <- subset(combined_data$metadata, `treatment:ch1` == "golimumab")
eset <- combined_data$eset[rownames(combined_data$annotations), rownames(golimumab)]
all(rownames(eset) == rownames(combined_data$annotations))

# subset the gene list 
gensym <- sapply(strsplit(combined_data$annotations$`Gene Symbol`, "///"), trimws) 

id_gensym <- tibble(probeids = rep(rownames(combined_data$annotations), 
                                   sapply(gensym, length)),
                    gensym = unlist(gensym)) %>%
  filter(gensym %in% unique(unlist(dbs_all)))

golimumab_eset <- clean_eset(combined_data$eset, 
                   gene_list_df = id_gensym)

# subset Expression Set data on the treatment 
golimumab_eset <- golimumab_eset[,colnames(combined_data$eset) %in% rownames(golimumab)]

# Combine metadata and expression matrix into a list
golimumab <- list(
  metadata = golimumab,
  eset = golimumab_eset
)

usethis::use_data(golimumab, overwrite = TRUE)
