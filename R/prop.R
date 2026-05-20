#' Calculate proportion of titers over a threshold
#'
#' @param titer a vector of titer values
#' @param threshold titer cutoff. Default is 2. With the defaults
#'        (`threshold.log.trans = TRUE`, `adj = 10`, `base = 2`),
#'        `threshold = 2` means a linear titer of `10 * 2^2 = 40`.
#'        Set `threshold.log.trans = FALSE` to pass a raw titer
#'        (e.g. `threshold = 40, threshold.log.trans = FALSE`).
#' @param threshold.log.trans whether the threshold value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param input.log.trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#'
#' @return proportion of titers above threshold
#'
#' @details
#' Returns the fraction of (non-missing) titers that are greater than or
#' equal to `threshold`. The threshold convention is shared with `width()`:
#' both functions interpret `threshold = N` as "linear titer >= adj * base^N"
#' by default, so the same numeric value of `threshold` selects the same
#' cutoff in both metrics.
#'
#' @examples
#' titer <- c(3, 5, 6, 4, 2)
#' prop(titer, threshold = 4)
#'
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

