library(tidyverse)
library(GEOquery)

load(here::here("data/ligand_receptor_db.RData"))
dbs_all = list(baderlab=baderlab, nichenet=nichenet, 
               fantom5=fantom5, citedb=citedb, all_dbs=all_dbs)

source("code/funcs.R")

# retrieve GEO data set and clean data
geo_data <- "GSE92415"
series_matrix <- paste0(geo_data,"_series_matrix.txt.gz")
# geo <- getGEO(geo_data, GSEMatrix=TRUE)
geo <- readRDS("code/gse/gse92415.rds") ## sockeye

e1 <- geo[[series_matrix]]
phenoData <- pData(e1)
ann <- e1@featureData@data
ann <- ann[ann$`Gene Symbol` != "", ]
exp <- exprs(e1)
dim(phenoData); dim(ann); dim(exp);
all(rownames(ann) == rownames(exp)); all(rownames(phenoData) == colnames(exp))

golimumab <- subset(phenoData, `treatment:ch1` == "golimumab")
eset <- exp[rownames(ann), rownames(golimumab)]
all(rownames(eset) == rownames(ann))

gensym <- sapply(strsplit(ann$`Gene Symbol`, "///"), trimws)
id_gensym <- data.frame(gensym = unlist(gensym),
                        probeids = rep(rownames(ann), sapply(gensym, length)))

X = eset[id_gensym$probeids, ] %>% 
  as.data.frame() %>% 
  mutate(genesym = id_gensym$gensym) %>% 
  group_by(genesym) %>% 
  summarise(across(everything(), mean, na.rm=TRUE))

eset <- as.matrix(X[,-1])
rownames(eset) <- X$genesym
y = golimumab$`visit:ch1`
obs_id = golimumab$`subject:ch1`

## rm genes from db if not in eset
dbs <- lapply(dbs_all, function(db){
  db <- db[names(db) %in% rownames(eset)]
  db <- lapply(db, function(i){
    intersect(i, rownames(eset))
  })
  db[sapply(db, length) > 0]
})

# Fast GSEA, run GSEA for all datasets and rank
cores <- 4

funs <- list(cfgsea_p = cfgsea_p,
             cgsva_p = cgsva_p,
             cpca_p = cpca_p,
             cpcr_p = cpcr_p,
             cgsvar_p = cgsvar_p)

run_all = function(eset, y, obs_id, dbs, cores, funs){
  lapply(funs, function(fun) fun(eset, y, obs_id, dbs, cores))
}

result <- run_all(eset, y, obs_id, dbs, cores, funs)

xranks <- lapply(result, function(i){
  100 - 100*round(sapply(i, function(j){ which(names(j) == "TNF")})/sapply(i, length), 2)
})

saveRDS(result, paste0("../results/",geo_data,".rds"))
