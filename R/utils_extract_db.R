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
    filter_db <- db[sapply(db, length) > 0]
    return(filter_db)
  })
  
  return(output_db)
  }

