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
#' equal to `threshold`. Unlike `width()`, `prop()` does **not** apply
#' any `adjust`-based shift, so the effective cutoff is fixed:
#'
#' \deqn{linear\ titer \geq adj \times base^{threshold}}
#'
#' For the package defaults (`adj = 10`, `base = 2`):
#'
#' | `threshold` | linear cutoff |
#' |---|---|
#' | 0 | titer >= 10 |
#' | 2 | titer >= 40 |
#' | 3 | titer >= 80 |
#'
#' To pass a linear titer directly, set `threshold.log.trans = FALSE`
#' (e.g. `prop(titer, threshold = 40, threshold.log.trans = FALSE)`).
#'
#' Note: `width()` uses a related-but-different convention. For
#' `width(threshold = N)` to count the same titers as `prop(threshold = N)`,
#' pass `adjust = -1` to `width()`. See `?width` for the full table.
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

