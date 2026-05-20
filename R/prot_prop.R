
#' Calculate the proportion of titers weighted by protection levels
#'
#' @param titer a vector of titer values
#' @param min_titer the smallest titer value(s), default is the minimum value in the titer vector
#' @param weight a vector of protection weights indexed by integer log-titer levels (lowest to highest).
#'        Default is `c(0.00, 0.05, 0.10, 0.30, 0.60, 0.75, 0.85, 0.90, 0.95, 0.99)`,
#'        taken from Hoa et al. (2022) — see Details.
#' @param input.log.trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param input.pct.trans whether the input weight value is transformed to percentages, default is TRUE
#'        TRUE = percentage transformed, FALSE = not percentage transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @return protected proportion
#'
#' @details
#' Each titer is mapped to a protection weight by its integer log-titer level
#' (counting up from `min_titer`). Titers at or below `min_titer` receive
#' `weight[1]`; titers at or above `min_titer + length(weight) - 1` receive
#' `weight[length(weight)]`; values in between use the corresponding entry.
#' The function returns the mean weight across all titers — interpretable as
#' the average probability of protection for that individual.
#'
#' The default `weight` vector
#' `c(0.00, 0.05, 0.10, 0.30, 0.60, 0.75, 0.85, 0.90, 0.95, 0.99)`
#' is the probability of protection for hemagglutination-inhibition (HI)
#' titer levels used in Hoa et al. (2022). Provide a different `weight`
#' vector to use a different titer-to-protection mapping.
#'
#' @references
#' Hoa LNM, Sullivan SG, Mai LQ, et al. Influenza A(H1N1)pdm09 But Not
#' A(H3N2) Virus Infection Induces Durable Seroprotection: Results From
#' the Ha Nam Cohort. *J Infect Dis.* 2022;226(1):59-69.
#' \doi{10.1093/infdis/jiaa293}
#'
#' @examples
#' titer <- c(3, 5, 6, 4, 2)
#' prot_prop(titer)
#'
#' @export
prot_prop = function(titer, min_titer = NULL, weight = c(0.00,0.05,0.10,0.30,0.60,0.75,0.85,0.90,0.95,0.99),
                     input.log.trans = TRUE, input.pct.trans = TRUE, base = 2, adj = 10){

  if(input.log.trans){
    log_titer = titer
  } else {
    log_titer = log(titer/adj, base)
  }

  if(input.pct.trans){
    pct_weight = weight
  } else {
    pct_weight = weight/100
  }

  if (is.null(min_titer)) {
    min_titer <- floor(min(log_titer, na.rm = TRUE))
  }

  len = length(weight)
  sum_prot = 0

  for (t in log_titer) {
    if (t <= min_titer) {
      sum_prot = sum_prot + pct_weight[1]
    } else if (t > min_titer + len - 1) {
      sum_prot = sum_prot + pct_weight[len]
    } else {
      sum_prot = sum_prot + pct_weight[t - min_titer + 1]
    }
  }

  rt_val = sum_prot/(length(titer))
  return(rt_val)
}
