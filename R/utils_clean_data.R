#' Given an expression set, make sure that data is 
#' cleaned by gene name using a provided gene list data frame
#'
#' @param eset Expression Set as a numeric matrix
#' @param gene_list_df 
#'
#' @return cleaned up expression matrix
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
