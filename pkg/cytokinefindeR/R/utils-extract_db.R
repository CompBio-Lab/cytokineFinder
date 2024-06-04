# Helper fun to filter db and extract ligands-receptors pairs match
extract_db <- function(cytokine = NULL, eset, dbs) {
  # For all databases, search for ligand that matched cytokine and eset
  genes <- rownames(eset)
  receptors_interest <- c(cytokine, genes)
  # Now we could databases based on the rownames of eset and cytokine
  output_db <- lapply(dbs, function(db) {
    # Then filter on db as well
    db <- db[names(db) %in% receptors_interest] |>
      lapply(function(i) { intersect(i, genes) })
    filter_db <- db[sapply(db, length) > 0]
    return(filter_db)
  })
  return(output_db)
}
