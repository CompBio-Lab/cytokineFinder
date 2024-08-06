preprocess <- function(eset, dbs = dbs_all) {
  # check NA in data
  ifelse(eset)
  # check 0 variable ligands from dbs and 0 variable 
  # Run PCA to get the first PC
  pc <- sapply(db, function(ligand){
    tryCatch({
      genexp <- t(eset[intersect(rownames(eset), ligand), , drop=FALSE])
      prcomp(genexp, center = TRUE, scale. = TRUE, rank. = 1)$x[, "PC1"]  
    }, error = function(e) NA)
  })
  # receptors (specific for PCA)
  # Check overlapping genes within the eset and dbs
}