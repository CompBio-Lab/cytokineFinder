#' Title
#'
#' @param eset 
#' @param dbs 
#'
#' @return
#' 
#'
#' @examples

preprocess <- function(eset, dbs = dbs_all) {
  # create a list of unique receptor genes across dbs
  unique_receptors <- unique(unlist(dbs_all))
  
  # 1) assess overlap:
  # clean up database and check intersect with eset
  # recursively check empty element remove ligand from database 
  dbs <- extract_from_db(eset = eset, dbs = dbs_all)
  
  # 2) check 0 variable ligands from dbs and intersect with eset 
  # Run PCA to get the first PC
  pc <- sapply(dbs, function(ligand){
    tryCatch({
      genexp <- t(eset[intersect(rownames(eset), ligand), , drop=FALSE])
      prcomp(genexp, center = TRUE, scale. = TRUE, rank. = 1)$x[, "PC1"]  
    }, error = function(e) NA)
  })
  # receptors (specific for PCA)
  # Check overlapping genes within the eset and dbs
}