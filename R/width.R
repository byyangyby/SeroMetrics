#' Calculate the width of the titers
#'
#' @param titer a vector of titer values
#' @param distance.weight a vector of the distance used to sorting titers. Eg. year of isolation or amino acid distance
#' @param threshold titer cutoff above which a `distance.weight` value is considered "covered". Default is 2.
#'        Uses the same convention as `prop()`: with the defaults (`threshold_trans = TRUE`, `adj = 10`, `base = 2`),
#'        `threshold = 2` means a linear titer of `10 * 2^2 = 40`. Set `threshold_trans = FALSE` to pass a raw titer.
#' @param threshold_trans whether the threshold value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param input.log.trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param adjust how to adjust titers, default is min(min(log_titer), 0)
#'
#' @return width
#'
#' @details
#' Width measures the breadth of the antibody response: the fraction of the
#' `distance.weight` axis (e.g. year-of-isolation range) over which the titer
#' exceeds `threshold`. For consecutive points where the titer crosses the
#' threshold, the crossing location is found by linear interpolation between
#' the two `distance.weight` values.
#'
#' Threshold semantics mirror `prop()`: with the defaults, `threshold = 2`
#' is the same as "titer >= 40". Setting `threshold_trans = FALSE` lets you
#' pass a linear titer directly (e.g. `threshold = 40`).
#'
#' Note: internally `log_titer` is shifted by `adjust` before the comparison,
#' so passing a non-default `adjust` shifts the effective threshold by
#' `adjust` log units relative to the documented value above.
#'
#' @examples
#' titer <- c(3, 5, 6, 4, 2)
#' year  <- c(1968, 1972, 1977, 1987, 1997)
#' width(titer, year)
#'
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

  if (is.null(adjust)) {
    adjust <- min(min(log_titer), 0)
  }

  len <- length(distance.weight)

  log_titer = log_titer - adjust

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

