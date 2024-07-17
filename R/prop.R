#' Calculate proportion of titers over a threshold
#'
#' @param titer titer value#'
#' @param threshold the titers above which are counted, defalut is 2
#' @param threshold.log.trans whether the threshold value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param input.log.trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10

#' @return PROP
#' @export
prop = function(titer, threshold = 2, threshold.log.trans = TRUE, input.log.trans = TRUE, base = 2, adj = 10){

  if(input.log.trans){
    log_titer = titer
  } else {
    log_titer = log(titer/adj, base)
  }

  if(threshold.log.trans){
    log_threshold = threshold
  } else {
    log_threshold = log(threshold/adj, base)
  }

  len <- sum(!is.na(log_titer))
  x = log_titer >= log_threshold

  rt_val = sum(x, na.rm=TRUE)/(len)

  return(rt_val)

}

