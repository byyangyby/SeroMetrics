#' Function to format data and calculate max_titer
#'
#' @param data a data frame containing the necessary columns (participant_id, titer)
#' @param part_col a vector identifying participant id
#' @param val_col a vector of target values. Eg. titers
#' @param group_col a vector identifying the sample group, default is NULL
#' @param var_trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param output_trans whether the output titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param max_titer_col the column name of the max titer column, default is "Max_titer"
#'
#' @return a data frame containing the individual max_titer values
#'
max_titerData <- function(data, part_col, val_col, group_col = NULL, var_trans = 1, output_trans = 1, base = 2, adj = 10, max_titer_col = "Max_titer") {

  source('R/max_titer.R')

  if(is.null(group_col)){
    stat <- unique(data[, part_col, drop = FALSE])
    Max_titer <- sapply(seq_len(nrow(stat)), function(i) max_titerCal(stat[i, part_col], NULL, data, part_col, val_col, group_col, var_trans, output_trans, base, adj))
  }else{
    stat <- unique(data[, c(part_col, group_col), drop = FALSE])
    Max_titer <- sapply(seq_len(nrow(stat)), function(i) {
      max_titerCal(stat[i, part_col], stat[i, group_col], data, part_col, val_col, group_col, var_trans, output_trans, base, adj)
    })
  }

  stat <- cbind(stat, Max_titer)

  colnames(stat)[colnames(stat) == "Max_titer"] <- max_titer_col

  return(stat)
}

#' Calculate the max_titer for a given participant
#'
#' @param i a participant ID
#' @param j a sample group number
#' @param data the input data frame
#' @param part_col a vector identifying participant id
#' @param val_col a vector of target values. Eg. titers
#' @param group_col a vector identifying the sample group
#' @param var_trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param output.log.trans whether the output titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base a parameter passed to the `max_titer` function
#' @param adj a parameter passed to the `max_titer` function
#'
#' @return the max_titer value for the given participant
max_titerCal <- function(i,j, data, part_col, val_col, group_col = NULL, var_trans, output_trans, base, adj) {

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

  max_titer <- max_titer(y, var_trans, base, adj, output_trans)

  return(max_titer)
}
