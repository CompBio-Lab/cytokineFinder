#' Extract specific ligands and merge results based on a chosen metric
#'
#' This function extracts specified ligands from a `BenchmarkResults` object 
#' and merges the results into a single data frame using the specified metric.
#'
#' @param benchmark_results A BenchmarkResults object containing nested results
#' @param ligands A vector of ligands to filter BenchmarkResults against
#' @param metric The metric column to use for sorting and ranking ("padj", "pval", "coef")
#'
#' @return A data frame containing extracted results for specified ligands
#' @export
#' @import dplyr
#' @import purrr
#' @examples
#' \dontrun{
#' results_df <- extract_ligands(
#'     benchmark_results = results,
#'     ligands = c("LigandA", "LigandB"),
#'     metric = "coef")
#' }
extract_ligands <- function(benchmark_results, ligands, metric = "padj") {
  # Provide a vector of cytokines to find
  ligands <- as.vector(ligands)
  
  results_df <- map_dfr(
    names(benchmark_results),
    function(method_name) {
      map_dfr(
        names(benchmark_results[[method_name]]),
        function(db_name) {
          process_method_db(
            df = benchmark_results[[method_name]][[db_name]],
            method_name = method_name,
            db_name = db_name,
            metric = metric
          )
        }
      )
    }
  )
  
  summary_df <- summarize_df(results_df, metric) %>%
    filter(ligand %in% ligands)
  
  return(summary_df)
}

#' @rdname extract_ligands
#' @param df Input data frame
#' @param method_name Name of the method
#' @param db_name Name of the database
#' @param metric Metric column to use for sorting
#'
#' @return Processed data frame with ranking
#' @export

process_method_db <- function(df, method_name, db_name, metric) {
  if (!metric %in% colnames(df)) {
    stop("Metric column '", metric, "' not found in data frame")
  }
  
  sorted_df <- if (metric == "coef") {
    df %>% arrange(desc(.data[[metric]]))
  } else {
    df %>% arrange(.data[[metric]])
  }
  
  sorted_df %>%
    mutate(
      database = db_name,
      method = method_name,
      # rank is based on top table list n of all ligands and the ligand's relative position 
      rank = 100 * (1 - (row_number() / n()))
    )
}

#' @rdname extract_ligands
#' @param results_df Combined results data frame
#' @param metric Metric column to summarize
#'
#' @return Summary data frame in long format
#' @export

summarize_df <- function(results_df, metric) {
  reshaped <- reshape_metric(results_df, metric, metric)
  reshaped %>% arrange(method, type, rank)
}

#' @rdname extract_ligands
#' @param df Input data frame
#' @param metric Column name to reshape
#' @param type Type label for the metric
#'
#' @return Long format data frame with value and type columns
#' @export
reshape_metric <- function(df, metric, type) {
  df %>%
    select(ligand, all_of(metric), method, database, rank) %>%
    rename(value = !!metric) %>%
    mutate(type = type) %>%
    filter(!is.na(value))
}