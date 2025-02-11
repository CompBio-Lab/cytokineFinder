#' Extract specific ligands and merge results based on a chosen metric
#'
#' This function extracts specified ligands from a `BenchmarkResults` object 
#' and merges the results into a single data frame using the specified metric.
#'
#' @param benchmark_results A BenchmarkResults object containing nested results
#' @param ligands A vector of ligands to filter BenchmarkResults against
#' @param metric A vector of metrics to select from ("padj", "pval", "coef")
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

extract_ligands <- function(benchmark_results, ligands, metrics = c("padj", "coef")) {
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
            metrics = metrics
          )
        }
      )
    }
  )
  
  summary_df <- results_df %>%
    filter(ligand %in% ligands, metric_type %in% metrics)
  
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

process_method_db <- function(df, method_name, db_name, metrics) {
  # Validate available metrics
  available_metrics <- intersect(metrics, colnames(df))
  if (length(available_metrics) == 0) {
    warning("No specified metrics found in ", method_name, "/", db_name)
    return(NULL)
  }
  
  # Process each metric independently
  map_dfr(available_metrics, function(metric) {
    # Determine sort direction
    sorted_df <- if (metric == "coef") {
      df %>% arrange(desc(.data[[metric]]))
    } else {
      df %>% arrange(.data[[metric]])
    }
    
    # Calculate metric-specific ranks
    sorted_df %>%
      mutate(
        database = db_name,
        method = method_name,
        metric_type = metric,
        rank = 100 * (1 - (row_number() / n())),
        value = .data[[metric]]
      ) %>%
      select(ligand, value, method, database, rank, metric_type)
  })
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