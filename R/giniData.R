#' Function to format data and calculate gini
#'
#' @param data a data frame containing the necessary columns (participant_id, titer)
#' @param part_col a vector identifying participant id
#' @param weight_col a vector of the distance used to sorting titers. Eg. year of isolation or genetic distance
#' @param val_col a vector of target values. Eg. titers
#' @param group_col a vector identifying the sample group, default is NULL
#' @param var_trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param adjust how to adjust titers, default is the minimum of titers - 1
#' @param gini_col the column name of the gini coefficient column, default is "gini_coefficient"
#'
#' @return a data frame containing the individual gini values
#'
giniData <- function(data, part_col,val_col,group_col = NULL, var_trans = 1, base = 2, adj = 10, adjust = NULL, gini_col = "gini_coefficient") {

  source('R/gini.R')

  if(is.null(group_col)){
    stat <- unique(data[, part_col, drop = FALSE])

    gini_coefficient <- sapply(seq_len(nrow(stat)), function(i) {
      giniCal(stat[i, part_col], NULL, data, part_col, val_col, group_col, var_trans, base, adj, adjust)
    })
  }else{
    stat <- unique(data[, c(part_col, group_col), drop = FALSE])

    gini_coefficient <- sapply(seq_len(nrow(stat)), function(i) {
      giniCal(stat[i, part_col], stat[i, group_col], data, part_col, val_col, group_col, var_trans, base, adj, adjust)
    })
  }
  stat <- cbind(stat, gini_coefficient)

  colnames(stat)[colnames(stat) == "gini_coefficient"] <- gini_col

  return(stat)
}

#' Calculate the gini for a given participant
#'
#' @param i a participant ID
#' @param j a sample group number
#' @param data the input data frame
#' @param part_col a vector identifying participant id
#' @param weight_col a vector of the distance used to sorting titers. Eg. year of isolation or genetic distance
#' @param val_col a vector of target values. Eg. titers
#' @param group_col a vector identifying the sample group
#' @param var_trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base a parameter passed to the `gini` function
#' @param adj a parameter passed to the `gini` function
#' @param adjust how to adjust titers, default is the minimum of titers - 1
#'
#' @return the gini value for the given participant
giniCal <- function(i, j, data, part_col, val_col, group_col, var_trans, base, adj, adjust = NULL) {

  if(is.null(j)){
    df <- data[data[[part_col]] == as.numeric(i), ]
  }else{
    if(is.numeric(data[[part_col]])){
      if(is.character(data[[group_col]])){
        df <- data[data[[part_col]] == as.numeric(i) & data[[group_col]] == as.character(j), ]
      }else if(is.numeric(data[[group_col]])){
        df <- data[data[[part_col]] == as.numeric(i) & data[[group_col]] == as.numeric(j), ]
      }
    }else if(is.character(data[[part_col]])){
      if(is.character(data[[group_col]])){
        df <- data[data[[part_col]] == as.character(i) & data[[group_col]] == as.character(j), ]
      }else if(is.numeric(data[[group_col]])){
        df <- data[data[[part_col]] == as.character(i) & data[[group_col]] == as.numeric(j), ]
      }
    }
  }


  y <- df[[val_col]]

  gini <- gini(y, var_trans, base, adj, adjust)

  return(gini)
}
