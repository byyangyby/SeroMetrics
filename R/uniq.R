#' Select the samples according to the requirements
#'
#'
#' @param data data frame
#' @param columns columns to select
#' @return uniq
#'

uniq = function(data, columns = NULL){
  
  library(dplyr)
  
  if (is.null(columns)) {
    columns <- names(data)
  }
  
  result <- data %>%
    select(all_of(columns)) %>% 
    distinct()
  
  return(result)
  
}
