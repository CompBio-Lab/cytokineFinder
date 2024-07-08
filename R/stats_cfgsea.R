#' Title
#'
#' @param eset Expression Set object containing gene expression data.
#' @param design Design matrix generated from create_design()
#' @param db Ligand-receptor database
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

cfgsea <- function(eset,design, db) {
  fit <- lmFit(eset, design)
  ebayes <- eBayes(fit)
  top <- topTable(ebayes, coef = 2)
  
  # create named vector of t-stats for each gene
  stats <- top$t
  names(stats) <- rownames(top)
  
  # Get the list of lengths of receptors (gene sets) 
  # representing each ligand
  length_receptors <- sapply(db, length)
  
  fgsea_results <- fgsea(pathways = db,
                         stats = stats,
                         minSize = min(length_receptors),
                         maxSize = max(length_receptors)
                         )
  return(fgsea_results)
}
