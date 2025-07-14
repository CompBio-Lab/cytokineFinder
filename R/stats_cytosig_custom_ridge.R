#' Run CytoSig Custom Ridge Regression function
#'
#' @param eset Expression matrix of microarray or RNAseq experiments
#' @param design design matrix list that indicates conditions to compare and if paired design
#' @param obs_id if paired design, indicate observation ID to map biological replicates to sample of origin if paired
#' @param correlation if paired design, indicate correlation blocks for paired data
#' @param beta_coef beta matrix for ridge regression from published CytoSig
#' 
#' @return CytoSig ridge regression results table
#' @export
#'
#' @examples
#' 
#' 

cytosig_custom_ridge <- function(eset, design, 
                                 obs_id = NULL, 
                                 correlation = NULL,
                                 beta_coef) {
  
  # Force evaluation of all parameters

  logfc <- run_limma(eset, design, obs_id, correlation) %>%
    dplyr::select(genes,logFC) %>%
    tibble::column_to_rownames("genes")
  com_genes <- intersect(rownames(beta_coef), rownames(logfc))
  bulk <- as.matrix(logfc[com_genes, ])
  sig <- as.matrix(beta_coef[com_genes, ])
  
  # create adjacency matrix for cytoSig
  beta1 <- solve(crossprod(sig, sig) + diag(ncol(sig)))
  beta2 <- crossprod(sig, bulk)
  # get ridge results
  ridge_res <- data.frame(coef = beta1 %*% beta2) %>%
    rownames_to_column("ligand")
  
  return(ridge_res)
}  