#' Convert nested benchmark results to a flat tibble
#'
#' Handles both the LRI benchmark structure
#' (\code{cytokine -> benchmarks -> method -> database -> data.frame}) and the
#' CytoSig structure (\code{cytokine -> method -> data.frame}).
#'
#' @param results_list Nested list of benchmark results.
#' @param study_type Character string identifying the study (e.g. \code{"treatment"}).
#' @param has_benchmarks_layer Logical. \code{TRUE} for LRI results (default);
#'   \code{FALSE} for CytoSig results.
#' @param method_class Character string for the \code{class} column.
#'   Defaults to \code{"LRI"} when \code{has_benchmarks_layer = TRUE} and
#'   \code{"CytoSig_ridge"} otherwise. Refactor R5: no longer hardcoded.
#'
#' @return A tibble with columns: study_type, cytokine, method, database, class,
#'   ligand_tables.
#' @export
#'
#' @importFrom purrr imap_dfr
#' @importFrom tibble tibble
benchlist_to_tbl <- function(results_list,
                             study_type,
                             has_benchmarks_layer = TRUE,
                             method_class         = NULL) {

  # Refactor R5: explicit class label, not hardcoded
  if (is.null(method_class))
    method_class <- if (has_benchmarks_layer) "LRI" else "CytoSig_ridge"

  imap_dfr(results_list, function(cytokine_val, cytokine_name) {
    data_to_process <- if (has_benchmarks_layer) {
      if (!"benchmarks" %in% names(cytokine_val)) return(NULL)
      cytokine_val$benchmarks
    } else {
      cytokine_val
    }

    imap_dfr(data_to_process, function(method_val, method_name) {
      imap_dfr(method_val, function(db_val, db_name) {
        if (!is.data.frame(db_val)) return(NULL)
        tibble(
          study_type    = study_type,
          cytokine      = cytokine_name,
          method        = method_name,
          database      = db_name,
          class         = method_class,
          ligand_tables = list(db_val)
        )
      })
    })
  })
}
