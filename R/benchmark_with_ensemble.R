#' Create ensemble results combining LRI-based and CytoSig methods
#'
#' This function combines results from LRI-based methods with CytoSig results 
#' using a rank-based ensemble approach. It performs an inner join between 
#' LRI and CytoSig results and computes ensemble rankings based on p-values.
#'
#' @param master_tbl A tibble containing benchmark results from benchlist_to_tbl(), 
#'   with columns: study_type, cytokine, method, database, class, ligand_tables
#' @param ensemble_method Character string specifying ensemble method. 
#'   Currently supports "mean_rank" (default)
#' @param pval_col_lri Character string specifying the p-value column name 
#'   in LRI method results. Default is "pval"
#' @param pval_col_cytosig Character string specifying the p-value column name 
#'   in CytoSig method results. Default is "pval"
#'
#' @return A tibble with additional ensemble columns:
#'   \itemize{
#'     \item overlap_count: Number of overlapping ligands between LRI and CytoSig
#'     \item ensemble_table: Combined LRI and CytoSig results with rankings
#'     \item ensemble_rank: Mean rank score combining LRI and CytoSig rankings
#'     \item lri_rank: Percentile rank from LRI method (100 = best, 0 = worst)  
#'     \item cytosig_rank: Percentile rank from CytoSig method (100 = best, 0 = worst)
#'   }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Filters for LRI-based methods (class == "LRI")
#'   \item Joins with CytoSig results (class == "CytoSig_Web") by study_type and cytokine
#'   \item Computes overlap between LRI and CytoSig ligand sets
#'   \item Creates ensemble rankings using percentile ranks (lower p-value = higher rank)
#'   \item Falls back to CytoSig-only results when no LRI data is available
#' }
#'
#' Ranking system: Converts p-values to percentile ranks where 100 represents 
#' the best (lowest p-value) and 0 represents the worst (highest p-value).
#'
#' @export
#'
#' @importFrom dplyr filter left_join select mutate inner_join
#' @importFrom purrr map2_int pmap map map_dbl
#'
#' @examples
#' \dontrun{
#' # Assuming you have benchmark results
#' results <- run_cytokinefinder_workflow(study_data, databases, methods)
#' 
#' # Convert to tibble format  
#' master_tbl <- benchlist_to_tbl(results$benchmarks, "my_study", FALSE)
#' 
#' # Create ensemble results
#' ensemble_results <- create_ensemble_results(
#'   master_tbl, 
#'   ensemble_method = "mean_rank",
#'   pval_col_lri = "pval",
#'   pval_col_cytosig = "pval"
#' )
#' 
#' # View ensemble rankings
#' ensemble_results %>% 
#'   select(cytokine, method, database, ensemble_rank, lri_rank, cytosig_rank)
#' }
create_ensemble_results <- function(master_tbl, 
                                    ensemble_method = "mean_rank",
                                    pval_col_lri = "pval",
                                    pval_col_cytosig = "pval") {
  
  master_tbl %>%
    filter(class == "LRI") %>%
    left_join(
      master_tbl %>% 
        filter(class == "CytoSig_Web") %>%
        select(study_type, cytokine, cytosig_table = ligand_tables),
      by = c("study_type", "cytokine")
    ) %>%
    mutate(
      # Count overlap
      overlap_count = map2_int(ligand_tables, cytosig_table, function(lri, cyto) {
        if(is.null(lri) || is.null(cyto)) return(0L)
        length(intersect(lri$ligand, cyto$ligand))
      }),
      
      # Create ensemble table with rank-based approach
      ensemble_data = pmap(list(ligand_tables, cytosig_table, cytokine), function(lri, cyto, target_cytokine) {
        if(is.null(cyto)) return(list(ensemble_table = NULL, ensemble_rank = NA, lri_rank = NA, cytosig_rank = NA, source = "no_data"))
        
        # Fallback to cytosig only if no LRI data
        if(is.null(lri) || !target_cytokine %in% lri$ligand) {
          cyto_ranked <- cyto %>%
            mutate(cytosig_rank = 100 * (1 - (rank(.data[[pval_col_cytosig]], ties.method = "min") - 1) / n()))
          
          target_row <- cyto_ranked %>% filter(ligand == target_cytokine)
          if(nrow(target_row) == 0) return(list(ensemble_table = NULL, ensemble_rank = NA, lri_rank = NA, cytosig_rank = NA, source = "no_data"))
          
          return(list(
            ensemble_table = cyto_ranked,
            ensemble_rank = target_row$cytosig_rank,
            lri_rank = NA,
            cytosig_rank = target_row$cytosig_rank,
            source = "cytosig_only"
          ))
        }
        
        # Create ensemble table (inner join) and compute all ranks at once
        ensemble_table <- inner_join(lri, cyto, by = "ligand", suffix = c("_lri", "_cytosig")) %>%
          mutate(
            # Compute ranks: lower p-value = higher rank (better)
            lri_rank = 100 * (1 - (rank(.data[[paste0(pval_col_lri, "_lri")]], ties.method = "min") - 1) / n()),
            cytosig_rank = 100 * (1 - (rank(.data[[paste0(pval_col_cytosig, "_cytosig")]], ties.method = "min") - 1) / n()),
            ensemble_rank = (lri_rank + cytosig_rank) / 2)
        
        # Extract target cytokine results
        target_row <- ensemble_table %>% filter(ligand == target_cytokine)
        if(nrow(target_row) == 0) return(list(ensemble_table = ensemble_table, ensemble_rank = NA, lri_rank = NA, cytosig_rank = NA, source = "no_data"))
        
        return(list(
          ensemble_table = ensemble_table,
          ensemble_rank = target_row$ensemble_rank,
          lri_rank = target_row$lri_rank,
          cytosig_rank = target_row$cytosig_rank,
          source = "ensemble"
        ))
      }),
      
      # Extract individual components
      ensemble_table = map(ensemble_data, "ensemble_table"),
      ensemble_rank = map_dbl(ensemble_data, "ensemble_rank"),
      lri_rank = map_dbl(ensemble_data, "lri_rank"),
      cytosig_rank = map_dbl(ensemble_data, "cytosig_rank"),
    ) %>%
    select(-cytosig_table, -ensemble_data)
}