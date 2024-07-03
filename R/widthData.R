#' Function to format data and calculate the weighted width
#'
#' @param data a data frame containing the necessary columns (participant_id, isolation_year, titer)
#' @param part_col a vector identifying participant id
#' @param weight_col a vector of the distance used to sorting titers. Eg. year of isolation or genetic distance
#' @param val_col a vector of target values. Eg. titers
#' @param group_col a vector identifying the sample group, default is NULL
#' @param var_trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param threshold_trans whether the threshold value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param threshold the titers above which are counted, default is 2
#' @param mode how to calculate values if weight_col is repeated, default is NULL
#' @param adjust how to adjust titers, default is the minimum of titers
#'
#' @return a data frame containing the individual weighted width
#'
widthData <- function(data, part_col, weight_col, val_col, group_col = NULL, var_trans = 1, threshold_trans = 1, base = 2, adj = 10, threshold = 2, mode = NULL, adjust = NULL) {

  source('sermetric/R/width.R')

  if(is.null(group_col)){
    stat <- unique(data[, part_col, drop = FALSE])

    Width <- sapply(seq_len(nrow(stat)), function(i) {
      widthCal(stat[i, part_col], NULL, data, part_col, weight_col, val_col, group_col, var_trans, threshold_trans, base, adj, threshold, mode, adjust)
    })
  }else{
    stat <- unique(data[, c(part_col, group_col), drop = FALSE])

    Width <- sapply(seq_len(nrow(stat)), function(i) {
      widthCal(stat[i, part_col], stat[i, group_col], data, part_col, weight_col, val_col, group_col, var_trans, threshold_trans, base, adj, threshold, mode, adjust)
    })
  }

  stat <- cbind(stat, Width)

  return(stat)
}

#' Calculate the weighted width for a given participant
#'
#' @param i a participant ID
#' @param j a sample group number
#' @param data the input data frame
#' @param part_col a vector identifying participant id
#' @param weight_col a vector of the distance used to sorting titers. Eg. year of isolation or genetic distance
#' @param val_col a vector of target values. Eg. titers
#' @param group_col a vector identifying the sample group, default is NULL
#' @param var_trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param threshold_trans whether the threshold value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param threshold the titers above which are counted, default is 2
#' @param mode how to calculate values if weight_col is repeated, default is NULL
#' @param adjust how to adjust titers, default is the minimum of titers
#'
#' @return the weighted width for the given participant
widthCal <- function(i, j, data, part_col, weight_col, val_col, group_col, var_trans, threshold_trans, base, adj, threshold, mode, adjust) {

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

  if(!is.null(mode)){
    if(mode == "mean"){
      df <- df %>%
        group_by(across(all_of(weight_col))) %>%
        summarise(
          titer = mean(.data[[val_col]], na.rm = TRUE),
          id = first(.data[[part_col]])
        )%>%
        ungroup() %>%
        rename(!!val_col := titer, !!part_col := id)
    }else if(mode == "max"){
      df <- df %>%
        group_by(across(all_of(weight_col))) %>%
        summarise(
          titer = max(.data[[val_col]], na.rm = TRUE),
          id = first(.data[[part_col]])
        )%>%
        ungroup() %>%
        rename(!!val_col := titer, !!part_col := id)
    }else if(mode == "min"){
      df <- df %>%
        group_by(across(all_of(weight_col))) %>%
        summarise(
          titer = min(.data[[val_col]], na.rm = TRUE),
          id = first(.data[[part_col]])
        )%>%
        ungroup() %>%
        rename(!!val_col := titer, !!part_col := id)
    }
  }

  df <- df[order(df[[weight_col]]), ]

  x <- df[[weight_col]]
  y <- df[[val_col]]

  width <- width(y, x, threshold, threshold_trans, var_trans, base, adj, adjust)

  return(width)
}
