#' Function to format data and calculate the weighted prop
#'
#' @param data a data frame containing the necessary columns (participant_id titer)
#' @param part_col a vector identifying participant id
#' @param val_col a vector of target values. Eg. titers
#' @param group_col a vector identifying the sample group, default is NULL
#' @param var_trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param threshold_trans whether the threshold value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param threshold the titers above which are counted, default is 2
#' @param prop_col the column name of the proportion column, default is "Proportion"
#'
#' @return a data frame containing the individual weighted prop
#'
propData <- function(data, part_col, val_col, group_col = NULL, var_trans = 1, threshold_trans = 1, base = 2, adj = 10, threshold = 2, prop_col = "Proportion") {

  source('serometric/R/prop.R')

  if(is.null(group_col)){
    stat <- unique(data[, part_col, drop = FALSE])

    Proportion <- sapply(seq_len(nrow(stat)), function(i) {
      propCal(stat[i, part_col], NULL, data, part_col, val_col, group_col = NULL, var_trans, threshold_trans, base, adj, threshold)
    })
  }else{
    stat <- unique(data[, c(part_col, group_col), drop = FALSE])

    Proportion <- sapply(seq_len(nrow(stat)), function(i) {
      propCal(stat[i, part_col], stat[i, group_col], data, part_col, val_col, group_col, var_trans, threshold_trans, base, adj, threshold)
    })
  }

  stat <- cbind(stat, Proportion)

  colnames(stat)[colnames(stat) == "Proportion"] <- prop_col

  return(stat)
}

#' Calculate the weighted prop for a given participant
#'
#' @param i a participant ID
#' @param j a sample group number
#' @param data the input data frame
#' @param part_col a vector identifying participant id
#' @param val_col a vector of target values. Eg. titers
#' @param group_col a vector identifying the sample group, default is NULL
#' @param var_trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param threshold_trans whether the threshold value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param threshold the titers above which are counted, defalut is 2
#'
#' @return the weighted prop for the given participant
propCal <- function(i, j, data, part_col, val_col, group_col, var_trans, threshold_trans, base, adj, threshold) {

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

  prop <- prop(y, threshold, threshold_trans, var_trans, base, adj)

  return(prop)
}
