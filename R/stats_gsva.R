#' Gene Set Variation Analysis for cytokines (cGSVA)
#'
#' @param eset Expression Set object containing gene expression data.
#' @param design Design matrix generated from create_design()
#' @param db Ligand-receptor database
#'
#' @return
#' @export
#'
#' @examples

cgsva <- function(eset, design, db){
  gsvapar <- GSVA::gsvaParam(eset, db, maxDiff = TRUE)
  gsva_eset <- GSVA::gsva(gsvapar)
  
  # Run DEA
  fit <- limma::eBayes(limma::lmFit(gsva_eset, design))
  top <- limma::topTable(fit, coef = 2, n = nrow(fit))
  pval <- top$P.Value
  names(pval) <- rownames(top)
  return(pval[order(pval)])
}