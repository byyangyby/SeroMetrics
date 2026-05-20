#' Select the samples according to the requirements
#'
#'
#' @param data data frame
#' @param columns columns to select
#' @return a data frame of unique rows for the selected columns
#'
#' @examples
#' df <- data.frame(id = c(1, 1, 2), year = c(2000, 2000, 2001),
#'                  titer = c(5, 5, 6))
#' uniq(df)
#' uniq(df, columns = c("id", "year"))
#'
#' @export
uniq = function(data, columns = NULL){
  
  if (is.null(columns)) {
    columns <- names(data)
  }
  
  result <- data |>
    select(all_of(columns)) |> 
    distinct()
  
  return(result)
  
}
