#' Function to format data and calculate GMT
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
#'
#' @return a data frame containing the individual GMT values
#'
gmtData <- function(data, part_col, val_col, group_col = NULL, var_trans = 1, output_trans = 1, base = 2, adj = 10) {

  source('R/gmt.R')

  if(is.null(group_col)){
    stat <- unique(data[, part_col, drop = FALSE])

    GMT <- sapply(seq_len(nrow(stat)), function(i) {
      gmtCal(stat[i, part_col], NULL, data, part_col, val_col, group_col, var_trans, output_trans, base, adj)
    })
  }else{
    stat <- unique(data[, c(part_col, group_col), drop = FALSE])

    GMT <- sapply(seq_len(nrow(stat)), function(i) {
      gmtCal(stat[i, part_col], stat[i, group_col], data, part_col, val_col, group_col, var_trans, output_trans, base, adj)
    })
  }

  stat <- cbind(stat, GMT)

  return(stat)
}

#' Calculate the GMT for a given participant
#'
#' @param i a participant ID
#' @param j a sample group number
#' @param data the input data frame
#' @param part_col a vector identifying participant id
#' @param val_col a vector of target values. Eg. titers
#' @param group_col a vector identifying the sample group
#' @param var_trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param output_trans whether the output titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base a parameter passed to the `gmt` function
#' @param adj a parameter passed to the `gmt` function
#'
#' @return the GMT value for the given participant
gmtCal <- function(i, j, data, part_col, val_col, group_col, var_trans, output_trans, base, adj) {

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

  gmt <- gmt(y, var_trans, base, adj, output_trans)

  return(gmt)
}
