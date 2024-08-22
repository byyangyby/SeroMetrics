
#' Calculate the gini coefficient of the titers
#'
#' @param titer a vector of titer values
#' @param input.log.trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param adjust how to adjust titers, default is the minimum of titers - 1
#'
#' @return gini
#'
#' @export
gini = function(titer, input.log.trans, base = 2, adj = 10, adjust = NULL){

  if(input.log.trans){
    log_titer = titer
  } else {
    log_titer = log(titer/adj, base)
  }

  if (is.null(adjust)) {
    adjust <- min(min(log_titer), 0)
  }

  log_titer = log_titer - adjust

  rt_val = Gini(log_titer)

  return(rt_val)

}

