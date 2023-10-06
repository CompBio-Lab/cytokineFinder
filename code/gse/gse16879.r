library(GEOquery)
gds <- getGEO("GSE16879")
saveRDS(gds, "gse16879.rds")
