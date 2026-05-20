#' @importFrom dplyr group_by summarise ungroup rename across first
#'   left_join select distinct all_of
#' @importFrom rlang .data :=
#' @importFrom moments skewness kurtosis
"_PACKAGE"

utils::globalVariables(c("titer", "id"))
