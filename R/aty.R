
#' Calculate average distance weighted by titers
#'
#' @param titer a vector of titer values
#' @param distance.weight a vector of the distance used to sorting titers. Eg. year of isolation or genetic distance
#' @param input.log.trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param adjust how to adjust titers, default is the minimum of titers

#' @return ATY
#'


aty = function(titer, distance.weight,
               input.log.trans, base = 2, adj = 10, adjust = NULL){

  if(input.log.trans){
    log_titer = titer
  } else {
    log_titer = log(titer/adj, base)
  }

  if(is.null(adjust)){
    adjust = min(log_titer)
  }

  len <- length(distance.weight)

  log_titer = log_titer - adjust

  y1 <- log_titer[1:(len-1)]
  y2 <- log_titer[2:len]
  t1 <- distance.weight[1:(len-1)]
  t2 <- distance.weight[2:len]

  auc = sum((y1 + y2)/2*(t1 - t2), na.rm=TRUE)
  rt_val = sum((y1 + y2)*(t1 + t2)/4*(t1 - t2), na.rm=TRUE)/auc
  if(is.infinite(rt_val)){rt_val <- NA}

  return(rt_val)

}

