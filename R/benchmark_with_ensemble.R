#' Create ensemble results combining LRI and CytoSig rankings
#'
#' Joins LRI-based results with CytoSig results and computes a mean-rank
#' ensemble score. Uses \code{ties.method = "min"} consistently (refactor R1).
#'
#' @param master_tbl Tibble from \code{benchlist_to_tbl()}, with columns:
#'   study_type, cytokine, method, database, class, ligand_tables.
#' @param ensemble_method Currently only \code{"mean_rank"}.
#' @param pval_col_lri P-value column name in LRI results. Default \code{"pval"}.
#' @param pval_col_cytosig P-value column name in CytoSig results. Default \code{"pval"}.
#'
#' @return \code{master_tbl} (LRI rows only) with added columns:
#'   overlap_count, ensemble_table, ensemble_rank, lri_rank, cytosig_rank.
#' @export
#'
#' @importFrom dplyr filter left_join select mutate inner_join
#' @importFrom purrr map2_int pmap map map_dbl
create_ensemble_results <- function(master_tbl,
                                    ensemble_method  = "mean_rank",
                                    pval_col_lri     = "pval",
                                    pval_col_cytosig = "pval") {

  master_tbl %>%
    filter(class == "LRI") %>%
    left_join(
      master_tbl %>%
        filter(class == "CytoSig_ridge") %>%
        select(study_type, cytokine, cytosig_table = ligand_tables),
      by = c("study_type", "cytokine")
    ) %>%
    mutate(
      overlap_count = map2_int(ligand_tables, cytosig_table, function(lri, cyto) {
        if (is.null(lri) || is.null(cyto)) return(0L)
        length(intersect(lri$ligand, cyto$ligand))
      }),

      ensemble_data = pmap(
        list(ligand_tables, cytosig_table, cytokine),
        function(lri, cyto, target_cytokine) {

          no_data <- list(ensemble_table = NULL, ensemble_rank = NA_real_,
                          lri_rank = NA_real_, cytosig_rank = NA_real_,
                          source = "no_data")

          if (is.null(cyto)) return(no_data)

          # Refactor R1: assign_pct_rank() with ties.method="min" everywhere
          if (is.null(lri) || !target_cytokine %in% lri$ligand) {
            cyto_ranked <- cyto %>%
              mutate(cytosig_rank = assign_pct_rank(.data[[pval_col_cytosig]]))
            target_row <- cyto_ranked %>% filter(ligand == target_cytokine)
            if (nrow(target_row) == 0) return(no_data)
            return(list(
              ensemble_table = cyto_ranked,
              ensemble_rank  = target_row$cytosig_rank,
              lri_rank       = NA_real_,
              cytosig_rank   = target_row$cytosig_rank,
              source         = "cytosig_only"
            ))
          }

          ensemble_table <- inner_join(lri, cyto, by = "ligand",
                                       suffix = c("_lri", "_cytosig")) %>%
            mutate(
              lri_rank      = assign_pct_rank(.data[[paste0(pval_col_lri, "_lri")]]),
              cytosig_rank  = assign_pct_rank(.data[[paste0(pval_col_cytosig, "_cytosig")]]),
              ensemble_rank = (lri_rank + cytosig_rank) / 2
            )

          target_row <- ensemble_table %>% filter(ligand == target_cytokine)
          if (nrow(target_row) == 0)
            return(c(list(ensemble_table = ensemble_table), no_data[-1]))

          list(
            ensemble_table = ensemble_table,
            ensemble_rank  = target_row$ensemble_rank,
            lri_rank       = target_row$lri_rank,
            cytosig_rank   = target_row$cytosig_rank,
            source         = "ensemble"
          )
        }
      ),

      ensemble_table = map(ensemble_data, "ensemble_table"),
      ensemble_rank  = map_dbl(ensemble_data, "ensemble_rank"),
      lri_rank       = map_dbl(ensemble_data, "lri_rank"),
      cytosig_rank   = map_dbl(ensemble_data, "cytosig_rank")
    ) %>%
    select(-cytosig_table, -ensemble_data)
}


# ---------------------------------------------------------------------------
# Refactor R1: shared percentile-rank helper used by ensemble and ranking scripts
# ---------------------------------------------------------------------------

#' Convert a numeric vector to percentile ranks (100 = best)
#'
#' Lower values are ranked better (suitable for p-values). NA values are
#' excluded from ranking and returned as NA in the output.
#'
#' @param values Numeric vector.
#' @param decreasing Logical. If \code{TRUE}, higher values rank better
#'   (suitable for coefficients such as \code{abs(coef)} from CytoSig).
#'   Default \code{FALSE}.
#' @param ties.method Passed to \code{rank()}. Default \code{"min"}.
#'
#' @return Numeric vector of the same length as \code{values}, with values in
#'   (0, 100]. NA inputs produce NA outputs.
#' @export
assign_pct_rank <- function(values, decreasing = FALSE, ties.method = "min") {
  out <- rep(NA_real_, length(values))
  ok  <- !is.na(values)
  if (!any(ok)) return(out)
  v   <- if (decreasing) -values[ok] else values[ok]
  r   <- rank(v, ties.method = ties.method)
  n   <- sum(ok)
  out[ok] <- 100 * (1 - (r - 1) / n)
  out
}
