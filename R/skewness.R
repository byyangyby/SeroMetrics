
#' Calculate the skewness of the titers
#'
#' @param titer a vector of titer values
#' @param distance.weight a vector the distance used to sorting titers. Eg. year of isolation or amino acid distance
#' @param input.log.trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param adjust how to adjust titers, default is min(min(log_titer), 0)
#'
#' @return skewness
#'
#' @examples
#' titer <- c(3, 5, 6, 4, 2)
#' year  <- c(1968, 1972, 1977, 1987, 1997)
#' skew(titer, year, input.log.trans = TRUE)
#'
#' @export
skew = function(titer, distance.weight, input.log.trans, base = 2, adj = 10, adjust = NULL){

  if(input.log.trans){
    log_titer = titer
  } else {
    log_titer = log(titer/adj, base)
  }

  if (is.null(adjust)) {
    adjust <- min(min(log_titer), 0)
  }

  log_titer_transf = log_titer - adjust

  dist_df <- unlist(lapply(seq_along(log_titer), function(i) {
    rep(distance.weight[i], log_titer_transf[i])
  }))

  rt_val = skewness(dist_df, na.rm = TRUE)

  return(rt_val)

}

