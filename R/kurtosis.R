
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
#' @examples
#' titer <- c(3, 5, 6, 4, 2)
#' year  <- c(1968, 1972, 1977, 1987, 1997)
#' kurt(titer, year, input.log.trans = TRUE)
#'
#' @export
kurt = function(titer, distance.weight, input.log.trans, base = 2, adj = 10, adjust = NULL){

  if(input.log.trans){
    log_titer = titer
  } else {
    log_titer = log(titer/adj, base)
  }

  if(is.null(adjust)){
    adjust = min(log_titer, na.rm = TRUE)
  }

  log_titer_transf = log_titer - adjust

  dist_df <- lapply(seq_along(log_titer), function(i) {
    t <- round(log_titer_transf[i])
    if (is.na(t) || t <= 0L) return(numeric(0))
    rep(distance.weight[i], t)
  }) |>
    (\(x) do.call("c", x))() |>
    unlist()

  rt_val = kurtosis(dist_df, na.rm = TRUE)

  return(rt_val)

}

