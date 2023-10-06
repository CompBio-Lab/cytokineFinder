library(GEOquery)
gds <- getGEO("GSE226244")
saveRDS(gds, "gse226244.rds")
