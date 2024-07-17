#' Calculate geometric mean titer against a panel of viruses
#'
#'
#' @param titer titer value
#' @param input.log.trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param output.log.trans whether the output titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @return GMT
#'

gmt = function(titer, input.log.trans = TRUE, base = 2, adj = 10, output.log.trans = TRUE){

  if(input.log.trans){
    log_titer = titer
  } else {
    log_titer = log(titer/adj, base)
  }

  rt_val = mean(log_titer, na.rm=TRUE)

  if(!output.log.trans){
    rt_val = (base^rt_val)*adj
  }

  return(rt_val)

}
