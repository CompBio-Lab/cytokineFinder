#' Calculate top ligands using a PCA approach and use for weights to 
#'
#' @param eset Expression Set object containing gene expression data.
#' @param design Design matrix generated from create_design()
#' @param db Ligand-receptor database
#' 
#' @return List of differentially expressed ligands ordered by p-values
#' @export
#'
#' @examples
#' 
#' @importFrom stats prcomp
#' @importFrom limma eBayes
#' @importFrom limma lmFit
#' @importFrom limma topTable
#' @importFrom tibble enframe

cpca <- function(eset, design, db){
  # Check if the design matrix is a data frame or matrix
  if (!is.data.frame(design) && !is.matrix(design)) {
    stop("The design argument must be a data frame or matrix.")
  }
  
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
  
  # run DEA using the design matrix integrated from previous create_design()
  fit <- eBayes(lmFit(pc, design))
  top <- topTable(fit, coef = 2, n = nrow(fit))
  pval <- top$P.Value
  names(pval) <- rownames(top)
  return(enframe(pval[order(pval)], name = "ligand", value = "pval"))
}