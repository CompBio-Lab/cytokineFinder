#' Gene Set Enrichment Analysis
#'
#' @param eset Expression Set object containing gene expression data.
#' @param design Design matrix generated from create_design()
#' @param db Ligand-receptor database
#' @param obs_id A vector of observation IDs to indicate if it's paired data
#' @param correlation Add a correlation block based on the dupcor package for paired analysis
#'
#' @return a table with GSEA results. Each row corresponds to a ligand
#' @export
#'
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @importFrom fgsea fgsea
#' 
#' @examples
#' # This is part of a series of enrichment analysis methods
#' # Basic usage:
#' \dontrun{
#' gsea_res <- cfgsea(eset, design, dbs)
#' }

cfgsea <- function(eset, design, db, 
                   obs_id = NULL, correlation = NULL) {
  # generate linear model from limma for DEA
  # First check if experiment samples are paired
  if (!is.null(obs_id)) {
    fit <- lmFit(eset, design, block = obs_id, correlation = correlation)
    message("fitting model with paired samples.")
  } else {
    fit <- lmFit(eset, design) 
    message("fitting model without paired sample consideration.")
  }
  efit <- eBayes(fit)
  # get topTable
  top <- topTable(efit, coef = 2, number = nrow(efit))
  
  # create named vector of t-stats for each gene
  stats <- top$t
  names(stats) <- rownames(top)
  
  # Get the list of lengths of receptors (gene sets) 
  # representing each ligand
  length_receptors <- sapply(db, length)
  
  # Run fgsea
  fgsea_results <- fgsea(pathways = db,
                         stats = stats,
                         minSize = min(length_receptors),
                         maxSize = max(length_receptors)
                         )
  fgsea_results <- dplyr::rename(fgsea_results, ligand = pathway)
  return(fgsea_results)
}
