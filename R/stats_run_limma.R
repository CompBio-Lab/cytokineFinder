#' Run limma in unpaired or paired mode
#'
#' @param eset expression matrix
#' @param design design matrix using the create_design function which indicates blocks for paired samples
#' @param obs_id required to indicate the blocks that map biological replicates to the same sample
#' @param correlation the average estimated inter-duplicate correlation. The average is the trimmed mean of the individual correlations on the atanh-transformed scale. 
#'
#' @returns top table that runs differential gene expression analysis. Main purpose for this is to get logFC of all genes for CytoSig input.
#' @export
#'
#' @examples
run_limma <- function (eset, design, obs_id = NULL, correlation = NULL) 
{
  if (!is.null(obs_id)) {
    fit <- limma::lmFit(eset, design, block = obs_id, correlation = correlation)
    message("fitting model with paired samples.")
  }
  else {
    fit <- limma::lmFit(eset, design)
    message("fitting model without paired sample consideration.")
  }
  efit <- limma::eBayes(fit)
  top <- limma::topTable(efit, coef = 2, number = nrow(efit)) %>%
    rownames_to_column(var = "genes")
  return(top)
}
