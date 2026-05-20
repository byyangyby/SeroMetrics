
#' Calculate the gini coefficient of the titers
#'
#' @param titer a vector of titer values
#' @param input.log.trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param adjust how to adjust titers, default is min(min(log_titer), 0)
#'
#' @return gini coefficient
#'
#' @examples
#' titer <- c(3, 5, 6, 4, 2)
#' gini(titer, input.log.trans = TRUE)
#'
#' @export
gini = function(titer, input.log.trans, base = 2, adj = 10, adjust = NULL){

  if(input.log.trans){
    log_titer = titer
  } else {
    log_titer = log(titer/adj, base)
  }

  if (is.null(adjust)) {
    adjust <- min(min(log_titer), 0)
  }

  log_titer = log_titer - adjust

  x <- sort(log_titer)
  n <- length(x)
  rt_val <- (2 * sum(x * seq_len(n)) - (n + 1) * sum(x)) / (n * sum(x))

  return(rt_val)

}

