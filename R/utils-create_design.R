#' Helper fun to create the design matrix
#'
#' @description
#' Provide a design matrix for fitting the linear model
#'
#'


create_design <- function(y, obs_id){
  if(is.null(obs_id)) {
    #for unpaired datasets
    design <- model.matrix(~y)
  } else{
    design <- model.matrix(~y+obs_id)
  }
  return(design)
}
