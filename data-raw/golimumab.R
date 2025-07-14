# Load libraries
library(cytokineFindeR)
library(tidyverse)

# ULCERATIVE COLITIS anti-TNF DRUG paired data
## retrieve from GEO and clean metadata and annotations
golimumab <- retrieve_geo("GSE92415")

metadata <- golimumab$GSE92415_series_matrix.txt.gz$metadata %>%
  filter(`treatment:ch1` %in% c("golimumab")) 

# Get metadata for paired samples 
# this will be used to filter the columns for the eset
md_paired <- metadata |> 
  group_by(`subject:ch1`) |>
  filter(n()>1) |>
  select(title, geo_accession, `subject:ch1`,`visit:ch1`)

# get obs_id and condition
cond <- factor(md_paired$`visit:ch1`, levels = c("Week 0", "Week 6"))
obs_id <- md_paired$`subject:ch1`

# clean up eset:
gensym <- sapply(strsplit(golimumab$GSE92415_series_matrix.txt.gz$annotations$`Gene Symbol`, "///"), trimws) 

id_gensym <- tibble(probeids = rep(rownames(golimumab$GSE92415_series_matrix.txt.gz$annotations), 
                                   sapply(gensym, length)),
                    gensym = unlist(gensym)) 


## apply annotations to eset and QC
eset <- clean_eset(golimumab$GSE92415_series_matrix.txt.gz$eset, 
                   id_gensym)

match <- match(md_paired$geo_accession, colnames(eset))
eset <- eset[,match,drop=FALSE] 


## Filter zero variance genes 
# filter non-zero variance genes
eset <- eset[apply(eset,1,sd) > 0,]
#hist(eset)

# No obvious spike so did not apply a low gexprs filter

golimumab <- list(md = md_paired,
                  qc_eset = eset,
                  cond = cond,
                  obs_id = obs_id)

usethis::use_data(golimumab, overwrite = TRUE)
