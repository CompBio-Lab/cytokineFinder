#' Helper fun to create the design matrix
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
    design <- model.matrix(~y)
  } else{
    design <- model.matrix(~y+obs_id)
  }
  return(design)
}
