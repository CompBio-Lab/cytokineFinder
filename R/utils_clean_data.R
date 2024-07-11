#' Title
#'
#' @param geo_ids A character string representing the GEO ID
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom GEOquery getGEO
#' @importFrom Biobase pData
#' @importFrom Biobase exprs

clean_geo_dataset <- function(geo_id) {
  tryCatch({
    geo_data <- getGEO(geo_id, GSEMatrix = TRUE)
    e1 <- geo_data[[paste0(geo_id, "_series_matrix.txt.gz")]]
    
    # Extract phenotype data
    metadata <- pData(e1)
    
    # Extract expression data
    eset <- exprs(e1)
    
    # Extract annotation data
    annotations <- e1@featureData@data
    annotations <- annotations[annotations$`Gene Symbol` != "", ]
    
    # Create a list to store all data
    combined_data <- list(
      metadata = metadata,
      eset = eset,
      annotations = annotations
    )
    
    return(combined_data)
  }, error = function(e) {
    message("Error in fetching GEO data for ID: ", geo_id, " - ", e$message)
    return(NULL)
  })
}
