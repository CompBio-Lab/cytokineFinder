#' Preprocess the eset
#'
#' @param eset 
#' @param dbs 
#'
#' @return The filtered expression set that contains 0 variability ligands
#' @export
#'
#' @examples

preprocess_eset <- function(eset, dbs) {
  # 1) filter the dbs against the eset 
  filtered_dbs <- filter_db_against_eset(eset, dbs)
  
  # 2) check 0 variable ligands from dbs and intersect with eset 
  # Run PCA to get the first PC
  filtered_dbs_zeroVariance <- remove_zero_variance_ligands(eset, dbs)
  
  unique_receptors <- unique(unlist(filtered_dbs_zeroVariance))
  filtered_eset <- eset[intersect(rownames(eset), unique_receptors),]

  return(list(filtered_eset, 
              filtered_dbs_zeroVariance))
}

#' Filter database and extract ligands-receptors pairs that match the eset data
#'
#' @description
#' Provide a list of databases and match against receptors of interest. 
#' Remove ligands that do not have 
#'
#' @details
#' Iterate through the expression matrix input and match a vector of genes
#' of interest.
#'
#' @param eset A numeric matrix that represents the expression set
#' @param dbs A nested list of lists containing the databases, the ligand genes, 
#' and the receptor genes for each ligand gene
#' @export
#'
#' @return A list of lists showing a set of ligand genes where each ligand
#' contains a list of genes that were subset defined as "receptor" genes
#' @examples
#' 

filter_db_against_eset <- function(eset, dbs) {
  # For all databases, search for ligand that matched cytokine and eset
  genes <- rownames(eset)
  
  # Iterate through each database and subset based on receptor genes available
  # in the eset
  output_db <- lapply(dbs, function(db) {
    # Update list of receptors against the genes
    db <- lapply(db, function(receptors) { 
      intersect(receptors, genes) 
    })
    db[sapply(db, length) > 0]
  })
  return(output_db)
}


#' Title
#'
#' @param eset 
#' @param dbs 
#'
#' @return
#' @export
#'
#' @examples

remove_zero_variance_ligands <- function(eset, dbs) {
  dbs_list_removedZeroVariance <- 
    lapply(dbs, function(db) {
      pc <- lapply(db, function(ligand){
        tryCatch({
          genexp <- t(eset[intersect(rownames(eset), ligand), , drop=FALSE])
          prcomp(genexp, center = TRUE, scale. = TRUE, rank. = 1)$x[, "PC1"]
        }, error = function(e) NA)
      })
      db[sapply(pc, length) > 1]
    })
  # remove 0 variability ligands
  return(dbs_list_removedZeroVariance)
}
