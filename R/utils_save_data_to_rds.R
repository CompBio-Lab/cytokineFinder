#' Save data to RDS or RDA format
#'
#' This function fetches GEO data using the geoQuery package 
#' and saves it as an RDS
#'
#' @param geo_id A character string representing the GEO ID 
#' to be downloaded
#' @param dir_path A character string representing the directory path 
#' where the RDS files should be saved. DEFAULT: data/
#'
#' @return The GEO data if successful; otherwise, return NULL
#' @export
#'
#' @examples
#' 
#' @importFrom GEOquery getGEO
#' @importFrom Biobase exprs pData featureData

save_to_rds_from_geo <- function(geo_id, dir_path = "data/"){
  tryCatch({
    saveRDS(geo_data, paste0(dir_path, geo_id, ".rds"))
    return(geo_data)
  }, error = function(e) {
    message("Error in fetching GEO data for ID: ", geo_id, " - ", e$message)
    return(NULL)
  })
}

#' Save multiple GEO datasets in batch
#' 
#' This function saves multiple GEO datasets to separate RDS files by
#' calling `save_to_rds` on each GEO ID.
#' 
#' @param geo_ids A vector of character strings representing the GEO IDs
#' @param dir_path A character string representing the directory path 
#' where the RDS files should be saved. DEFAULT: data/
#'
#' @return A list of GEO with GEO IDs as names
#' @export
#'
#' @examples

save_rds_in_batch <- function(geo_ids, dir_path = "data/"){
  geo_datasets <- lapply(geo_ids, function(geo_id) {
    save_to_rds_from_geo(geo_id, dir_path)
  })
  names(geo_datasets) <- geo_ids
  return(geo_datasets)
}