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

retrieve_geo = function(geo_id){
  geo_data <- tryCatch({
    getGEO(geo_id, GSEMatrix = TRUE)
  }, error = function(e) {
    message("Error in fetching GEO data for ID: ", geo_id, 
            " - ", e$message)
    return(NULL)
  })
  e1 <- geo_data
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


#' Given an expression set, make sure that data is 
#' cleaned by gene name using a provided gene list data frame
#'
#' @param eset Expression Set as a numeric matrix
#' @param gene_list_df 
#'
#' @return
#' @export
#' @name clean_eset
#'
#' @examples
#' 

clean_eset <- function(eset, gene_list_df){
  # clean eset against list of probe genes
  # combine probes that bind to multiple genes
  
  X <- eset[gene_list_df$probeids, ] %>% 
    # convert to a data frame to transform probe IDs to genes
    # Take a mean of all probes that match to the Gene ID
    as.data.frame() %>% 
    mutate(genesym = gene_list_df$gensym) %>% 
    group_by(genesym) %>% 
    summarise(across(everything(), ~ mean(.x, na.rm=TRUE)))
  
  eset <- as.matrix(X[,-1])
  rownames(eset) <- X$genesym
  
  return(eset)
  
}
