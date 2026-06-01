#' Run the complete cytokineFindeR workflow
#'
#' Preprocesses the expression matrix, creates the design matrix (with optional
#' paired-sample handling via \code{duplicateCorrelation}), and runs all
#' requested LRI methods across all databases.
#'
#' @param study_data Named list with required elements:
#'   \itemize{
#'     \item \code{qc_eset}: genes x samples numeric matrix (log-normalised counts)
#'     \item \code{cond}: factor or character vector of condition labels
#'     \item \code{obs_id}: (optional) donor ID vector for paired analysis
#'   }
#' @param databases Named list of LRI databases (e.g. \code{dbs_cytosig}).
#' @param methods Character vector of methods to run. Supported:
#'   \code{"cfgsea"}, \code{"pca_limma"}, \code{"gsva_limma"}.
#'
#' @return \code{study_data} with added elements:
#'   \itemize{
#'     \item \code{design}: output of \code{create_design()}
#'     \item \code{benchmarks}: \code{BenchmarkResults} object
#'     \item \code{preprocessing}: output of \code{preprocess_eset()}
#'   }
#' @export
#'
#' @importFrom rlang %||%
run_cytokinefinder <- function(study_data, databases, methods) {

  required <- c("qc_eset", "cond")
  missing  <- setdiff(required, names(study_data))
  if (length(missing) > 0)
    stop("Missing required elements: ", paste(missing, collapse = ", "))

  eset   <- study_data$qc_eset
  obs_id <- study_data$obs_id %||% NULL

  # Validate that the condition vector matches the number of samples up front,
  # so the dimension mismatch is reported clearly rather than surfacing deep
  # inside limma (where run_lri_methods would swallow it via tryCatch).
  if (length(study_data$cond) != ncol(eset))
    stop("row dimension of design doesn't match column dimension")

  # Preprocess: filter databases against expressed genes, remove zero-variance ligands
  preprocess        <- preprocess_eset(eset, databases)
  eset_f            <- preprocess$eset_f
  dbs_f             <- preprocess$dbs_f
  study_data$preprocessing <- preprocess

  # Design matrix + duplicate correlation (paired) or plain OLS (unpaired)
  if (!is.null(obs_id)) {
    design_results <- create_design(study_data$cond, obs_id = obs_id, eset = eset)
    benchmark_results <- run_lri_methods(
      eset        = eset_f,
      design      = design_results$design,
      dbs         = dbs_f,
      methods     = methods,
      treatment   = study_data$cond,
      obs_id      = obs_id,
      correlation = design_results$dupcor$consensus.correlation  # explicit field
    )
  } else {
    design_results <- create_design(study_data$cond, eset = eset)
    benchmark_results <- run_lri_methods(
      eset      = eset_f,
      design    = design_results$design,
      dbs       = dbs_f,
      methods   = methods,
      treatment = study_data$cond
    )
  }

  study_data$design     <- design_results
  study_data$benchmarks <- benchmark_results
  return(study_data)
}
