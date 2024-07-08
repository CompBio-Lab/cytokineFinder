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

clean_geo_dataset <- function(geo_id) {
  tryCatch({
    geo_data <- getGEO(geo_id, GSEMatrix = TRUE)
    e1 <- geo_data[[1]]
    
    # Extract phenotype data
    pheno_data <- pData(e1)
    
    # Extract expression data
    exp_data <- exprs(e1)
    
    # Extract annotation data
    ann_data <- e1@featureData@data
    ann_data <- ann_data[ann_data$`Gene Symbol` != "", ]
    
    # Create a list to store all data
    combined_data <- list(
      phenotype = pheno_data,
      expression = exp_data,
      annotation = ann_data
    )
    
    return(combined_data)
  }, error = function(e) {
    message("Error in fetching GEO data for ID: ", geo_id, " - ", e$message)
    return(NULL)
  })
}
