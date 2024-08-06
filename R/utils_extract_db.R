#' Helper function to filter database and extract ligands-receptors pairs match
#'
#' @description
#' Provide a list of databases and match ligands of interest.
#'
#' @details
#' Iterate through the expression matrix input and match a vector of genes
#' of interest.
#'
#' @param cytokine An optional vector of cytokines that can be added to address
#' some specifically named receptors that may have been omitted due to character
#' changes such as TNFA vs TNF.
#'
#' @return A list of lists showing a set of ligand genes where each ligand
#' contains a list of genes that were subset defined as "receptor" genes

extract_from_db <- function(cytokine = NULL, eset, dbs) {
  # For all databases, search for ligand that matched cytokine and eset
  genes <- rownames(eset)
  receptors_interest <- c(cytokine, genes)
  
  # Now we could get databases based on the rownames of eset and cytokine
  output_db <- lapply(dbs, function(db) {
    # Then filter on db as well
    db <- db[names(db) %in% receptors_interest] |>
      lapply(function(i) { intersect(i, genes) })
    filter_db <- db[sapply(db, length) > 0]
  })
  return(output_db)
}

