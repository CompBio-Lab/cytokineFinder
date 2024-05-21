library(GEOquery)
gds <- getGEO("GSE220652")
saveRDS(gds, "gse220652.rds")
