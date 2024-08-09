#' Preprocess the eset
#'
#' @param eset 
#' @param dbs 
#'
#' @return The filtered expression set that contains 0 variability ligands
#' 
#'
#' @examples

preprocess_eset <- function(eset, db) {
  # create a list of unique receptor genes across dbs
  unique_receptors <- unique(unlist(dbs_all))
  
   
  dbs <- filter_db_against_eset(eset = eset, dbs = dbs_all)

  # 2) check 0 variable ligands from dbs and intersect with eset 
  # Run PCA to get the first PC
  pc <- sapply(db, function(ligand){
    tryCatch({
      genexp <- t(eset[intersect(rownames(eset), ligand), , drop=FALSE])
      prcomp(genexp, center = TRUE, scale. = TRUE, rank. = 1)$x[, "PC1"]  
    }, error = function(e) NA)
  })
  # remove 0 variability ligands
  pc <- pc[sapply(pc, length) > 1] %>% 
    do.call(rbind, .)
  
  eset_filt <- eset[rownames(eset) %in% rownames(pc),] 
  return(eset_filt)
 }