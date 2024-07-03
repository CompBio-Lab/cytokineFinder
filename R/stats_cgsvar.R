#' Calculate coefficients from variable selection using GSVA and PLS-DA
#' 
#' Uses a multiple regression approach when computing receptor weights with
#' Partial Least Squares - Discriminant Analysis and pools these receptors
#' together to compute a value for the ligand or cytokine given a gene set.
#'
#' @param eset Expression Set object containing gene expression data.
#' @param y Response variable
#' @param obs_id Observation ID or sample if looking there are biological replicates
#' @param db ligand-receptor database
#'
#' @return A named vector of ligands indicating importance 
#' @export
#'
#' @examples
#' 
#' @importFrom GSVA gsvaParam
#' @importFrom GSVA gsva
#' @importFrom mixOmics plsda
#' @importFrom mixOmics selectVar

cgsvar <- function(eset, y, obs_id, db){
  gsvapar <- gsvaParam(eset, db, maxDiff = TRUE)
  gsva_eset <- gsva(gsvapar)
  fit <- plsda(gsva_eset, y)
  coef <- abs(selectVar(fit, comp=1)$value$value.var)
  names(coef) <- rownames(selectVar(fit, comp=1)$value)
  return(coef[order(coef, decreasing = TRUE)])
}