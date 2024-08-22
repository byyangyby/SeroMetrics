
#' Function to transform data frame
#'
#' @param df a data frame
#' @param panel the panel of data transformation, panel = 1 means transpose some columns, panel = 2 means combining group columns
#' @param col_start the column number from which the transformation starts, default is NULL
#' @param col_end the column number at which the transformation ends, default is NULL
#' @param strain_col the column name of the strains in the output, default is "strain"
#' @param val_col the column name of the titer values in the output, default is "titer"
#' @param group_cols the group columns intended to be combined, default is NULL
#' @param group_col the new group column name, default is NULL
#' @param def the default linking between group column names to form the new column name, default is "_"
#' @param drop whether to drop the original group columns, default is TRUE
#'
#' @return a transformed data frame, with the column in the range [col_start, col_end] shifted to two columns
#' @export
#'
dataTransform <- function(df, panel, col_start = NULL, col_end = NULL, strain_col = "strain", val_col = "titer",
                          group_cols = NULL, group_col = NULL, def = "_", drop = 1){

  if(panel == 1){
    df <- df %>%
      pivot_longer(
        cols = col_start:col_end,
        names_to = strain_col,
        values_to = val_col
      )
  }

  if(panel == 2){

    if(is.null(group_col)){
      group_col <- paste(group_cols, collapse = def)
    }

    df[[group_col]] <- apply(df[, group_cols], 1, paste, collapse = def)

    if(drop){
      df <- df %>% select(-all_of(group_cols))
    }

  }

  return(df)

}
