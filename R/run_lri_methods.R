#' Run LRI benchmarking methods across databases
#'
#' Dispatches \code{cfgsea}, \code{pca_limma}, and \code{gsva_limma} across all
#' databases in \code{dbs} using parallel execution via \code{future}.
#' PLSDA methods have been removed.
#'
#' @param eset Numeric matrix (genes x samples), pre-filtered by
#'   \code{preprocess_eset()}.
#' @param design Design matrix from \code{create_design()}.
#' @param dbs Named list of databases (each a named list of receptor gene vectors).
#' @param methods Character vector of method names. Must be one or more of:
#'   \code{"cfgsea"}, \code{"pca_limma"}, \code{"gsva_limma"}.
#' @param treatment Condition vector (used for logging only).
#' @param obs_id Optional donor ID vector for paired analysis.
#' @param correlation Optional within-block correlation.
#' @param verbose Logical; print progress messages. Default \code{FALSE}.
#'
#' @return A \code{BenchmarkResults} object (named list: method -> database -> data.frame).
#' @export
#'
#' @importFrom future plan multicore multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom stats setNames
run_lri_methods <- function(eset, design, dbs, methods,
                            treatment   = NULL,
                            obs_id      = NULL,
                            correlation = NULL,
                            verbose     = FALSE) {

  # Refactor R3: flat registry — no PLSDA, no grepl() dispatch
  method_registry <- list(
    cfgsea     = cfgsea,
    pca_limma  = pca_limma,
    gsva_limma = gsva_limma
  )

  methods <- match.arg(methods, choices = names(method_registry), several.ok = TRUE)

  tryCatch(
    future::plan(if (.Platform$OS.type == "unix") future::multicore
                 else future::multisession),
    error = function(e) {
      message("Parallel unavailable, using sequential: ", e$message)
      future::plan(future::sequential)
    }
  )

  results <- list()

  for (method_name in methods) {
    method <- method_registry[[method_name]]

    method_results <- future_lapply(names(dbs), function(database) {
      if (verbose)
        message("Method: ", method_name, "  DB: ", database)

      result <- tryCatch(
        method(eset, design, dbs[[database]],
               obs_id = obs_id, correlation = correlation),
        error = function(e) {
          message("Error in ", method_name, "/", database, ": ", e$message)
          NULL
        }
      )
      list(database = database, result = result)
    }, future.seed = TRUE)

    results[[method_name]] <- setNames(
      lapply(method_results, `[[`, "result"),
      sapply(method_results, `[[`, "database")
    )
  }

  structure(results, class = "BenchmarkResults")
}
