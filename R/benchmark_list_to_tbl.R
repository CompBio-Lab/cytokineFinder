#' Convert nested benchmark results to tibble format
#' 
#' @param results_list Nested list of benchmark results (cytokine -> benchmarks -> method -> database -> data) or CytoSig Results (cytokine -> method -> data)
#' @param study_type Character string identifying the study type
#' @param has_benchmarks_layer boolean to handle benchmark and cytosig structures
#' 
#' @return A tibble with flattened benchmark results
#' Convert nested results to tibble format (handles both benchmark and cytosig structures)
#' @importFrom purrr imap_dfr
#' @importFrom tibble tibble
#' 
#' @export

benchlist_to_tbl <- function(results_list,
                             study_type, 
                             has_benchmarks_layer = TRUE) {
  imap_dfr(results_list, function(cytokine_val, cytokine_name) {
    # Determine the data to process based on structure
    if(has_benchmarks_layer) {
      # Benchmark structure: cytokine -> benchmarks -> method -> database -> data
      if(!"benchmarks" %in% names(cytokine_val)) return(NULL)
      data_to_process <- cytokine_val$benchmarks
    } else {
      # CytoSig structure: cytokine -> method -> database -> data
      data_to_process <- cytokine_val
    }
    imap_dfr(data_to_process, function(method_val, method_name) {
      imap_dfr(method_val, function(db_val, db_name) {
        if(!is.data.frame(db_val)) return(NULL)
        
        tibble(study_type = study_type,
               cytokine = cytokine_name,
               method = method_name,
               database = db_name,
               class = ifelse(has_benchmarks_layer == TRUE, "LRI", "CytoSig_Web"),
               ligand_tables = list(db_val))
      })
    })
  })
}