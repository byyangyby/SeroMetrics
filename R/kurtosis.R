
#' Calculate the kurtosis of the titers
#'
#' @param titer a vector of titer values
#' @param distance.weight a vector the distance used to sorting titers. Eg. year of isolation or amino acid distance
#' @param input.log.trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param adjust how to adjust titers, default is the minimum of titers
#'
#' @return kurtosis
#'
#' @export
kurt = function(titer, distance.weight, input.log.trans, base = 2, adj = 10, adjust = NULL){

  if(input.log.trans){
    log_titer = titer
  } else {
    log_titer = log(titer/adj, base)
  }

  if(is.null(adjust)){
    adjust = min(log_titer)
  }

  log_titer_transf = log_titer - adjust

  dist_df = lapply(seq_along(log_titer), function(i){

    rep(distance.weight[i], log_titer_transf[i])

  }) %>%
    do.call(
      "c", .
    ) %>%
    unlist

  rt_val = kurtosis(dist_df, na.rm = TRUE)

  return(rt_val)

}

