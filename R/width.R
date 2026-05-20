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
#' @param adjust how to adjust titers, default is the minimum of titers
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
#' **Effective threshold.** Internally `log_titer` is shifted by `adjust`
#' (`log_titer - adjust`) and the threshold is incremented by 1
#' (`z = threshold + 1`) before the comparison. With `threshold_trans = TRUE`,
#' the cutoff that gets counted as "covered" is therefore:
#'
#' \deqn{linear\ titer \geq adj \times base^{(threshold + 1 + adjust)}}
#'
#' For the package defaults (`adj = 10`, `base = 2`):
#'
#' | `threshold` | `adjust = 0` | `adjust = -1` | `adjust = -2` |
#' |---|---|---|---|
#' | 0 | titer >= 20  | titer >= 10 | titer >= 5  |
#' | 2 | titer >= 80  | titer >= 40 | titer >= 20 |
#' | 4 | titer >= 320 | titer >= 160 | titer >= 80 |
#'
#' The shift makes the cutoff *relative* to the participant's minimum titer
#' when `adjust` is left at its default. Pass an explicit `adjust` to align
#' the cutoff with a fixed linear titer across participants
#' (e.g. `adjust = -1` makes `width(threshold = N)` and `prop(threshold = N)`
#' match numerically).
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

