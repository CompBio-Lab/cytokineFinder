#' Create design matrix for the data based on the outcome and the samples
#'
#' @param y 
#' @param obs_id 
#'
#' @return
#' @export
#'
#' @examples

create_design <- function(y, obs_id){
  if(is.null(obs_id)) {
    #for unpaired datasets
    design <- stats::model.matrix(~y)
  } else{
    design <- stats::model.matrix(~y+obs_id)
  }
  return(design)
}