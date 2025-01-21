#' Retrieve GEO data using the Bioconductor package GEOquery, clean it up, and store in a df list
#'
#' @param geo_ids A character string representing the GEO ID
#'
#' @return
#' @export
#' @name retrieve_geo
#'
#' @examples
#'
#' @importFrom GEOquery getGEO
#' @importFrom Biobase pData
#' @importFrom Biobase exprs

retrieve_geo <- function(geo_id){
  geo_data <- tryCatch({
    getGEO(geo_id, GSEMatrix = TRUE)
  }, error = function(e) {
    message("Error in fetching GEO data for ID: ", geo_id, 
            " - ", e$message)
    return(NULL)
  })
  e1 <- geo_data
  # handle retrieval if e1 returns a list
  combined_data <- lapply(e1, function(i){
    metadata <- pData(i)
    eset <- exprs(i)
    annotations <- i@featureData@data
    annotations <- annotations[annotations$`Gene Symbol` != "", ]
    list(metadata = metadata, eset = eset, annotations = annotations)
  })
  if (any(sum(rapply(combined_data, nrow)) == 0)) {
    message("Warning: One or more data frames           
            in the combined data list are empty.")
  }
  return(combined_data)
}

# retrieve_geo <- function(geo_id) {
#   # Try to fetch the GEO data
#   geo_data <- tryCatch({
#     getGEO(geo_id, GSEMatrix = TRUE)
#     }, error = function(e) {
#       message("Error in fetching GEO data for ID: ", geo_id, " - ", e$message)
#       return(NULL)
#       })
#     
#   e1 <- geo_data[[paste0(geo_id, "_series_matrix.txt.gz")]]
#   
#   # Extract phenotype data
#   metadata <- pData(e1)
#   # Extract expression data
#   eset <- exprs(e1)
#   
#   # Extract annotation data
#   annotations <- e1@featureData@data
#   annotations <- annotations[annotations$`Gene Symbol` != "", ]
#   
#   # Create a list to store all data
#   combined_data <- list(
#     metadata = metadata,
#     eset = eset,
#     annotations = annotations
#     )
#   
#   # Check if any of the data frames are empty
#   if (any(sapply(combined_data, nrow) == 0)) {
#     message("Warning: One or more data frames 
#             in the combined data list are empty.")
#   }
#   
#   return(combined_data)
#   
# }