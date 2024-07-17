#' Calculate the width of the titers
#'
#' @param titer titer value
#' @param distance.weight a vector the distance used to sorting titers. Eg. year of isolation or amino acid distance
#' @param threshold the titers above which are counted, default is 2
#' @param threshold_trans whether the threshold value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param input.log.trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' #' @param adjust how to adjust titers, default is the minimum of titers
#'
#' @return width
#' @export
width = function(titer, distance.weight, threshold = 2, threshold_trans = TRUE, input.log.trans = TRUE, base = 2, adj = 10, adjust = NULL){

  if(input.log.trans){
    log_titer = titer
  } else {
    log_titer = log(titer/adj, base)
  }

  if(threshold_trans){
    z = threshold
  } else {
    z = log(threshold/adj, base)
  }

  if(is.null(adjust)){
    adjust = min(log_titer)
  }

  len <- length(distance.weight)

  log_titer = log_titer - adjust

  z = z + 1

  y.1 <- log_titer[1:(len-1)]
  y.2 <- log_titer[2:len]
  t.1 <- distance.weight[1:(len-1)]
  t.2 <- distance.weight[2:len]

  t.1GreaterThanZ <- y.1>=z
  t.2GreaterThanZ <- y.2>=z

  width <- rep(0, len-1)
  n = sum(t.2-t.1)

  # condition 1
  con <- t.1GreaterThanZ & t.2GreaterThanZ # both greater than z
  width[con] <- t.2[con] - t.1[con]

  # condition 2
  con <- !(t.1GreaterThanZ) & t.2GreaterThanZ # yi+1 greater than z & yi less than z
  width[con] <- t.2[con] - (((z-y.2)*t.1 + (y.1-z)*t.2)/(y.1-y.2))[con]

  # condition 3
  con <- t.1GreaterThanZ & !(t.2GreaterThanZ) # yi+1 less than z & yi greater than z
  width[con] <- (((z-y.2)*t.1 + (y.1-z)*t.2)/(y.1-y.2))[con] - t.1[con]

  rt_val = sum(width)/n

  return(rt_val)

}

