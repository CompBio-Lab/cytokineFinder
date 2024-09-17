#' Title
#'
#' @param eset Expression Set object containing gene expression data.
#' @param treatment Treatment response variable
#' @param db ligand-receptor database
#'
#' @return 
#' @export
#'
#' @examples
#' # Load the database
#' # dbs_all <- dbs_all
#' # cpca_plsda(eset, treatment, dbs_all$baderlab)

pca_plsda <- function(eset, treatment, db){
  pc <- lapply(db, function(ligand){
    tryCatch({
      genexp <- t(eset[intersect(rownames(eset), ligand), , drop=FALSE])
      prcomp(genexp, center = TRUE, scale. = TRUE, rank. = 1)$x[, "PC1"]  
    }, error = function(e) NA)
  })
  # remove 0 variability ligands
  pc <- pc[sapply(pc, length) > 1] %>% 
    do.call(rbind, .)
  # fit to plsda
  fit <- mixOmics::plsda(t(pc), treatment)
  coef <- abs(mixOmics::selectVar(fit, comp=1)$value$value.var)
  names(coef) <- rownames(mixOmics::selectVar(fit, comp=1)$value)
  return(enframe(coef[order(coef, decreasing = TRUE)],
                 name = "ligand",
                 value = "coef")
         )
}
