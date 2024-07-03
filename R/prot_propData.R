#' Function to format data and calculate the weighted prot_prop
#'
#' @param data a data frame containing the necessary columns (participant_id titer)
#' @param part_col a vector identifying participant id
#' @param val_col a vector of target values. Eg. titers
#' @param group_col a vector identifying the sample group, default is NULL
#' @param min.y the smallest titer value(s), default is the minimum value in the titer vector(calculated later, temporarily being NULL in the formula)
#' @param weight a vector of weight values according to titer values, default is c(0.00,0.05,0.10,0.30,0.60,0.75,0.85,0.90,0.95,0.99)
#' @param var_trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param weight_trans whether the input weight value is transformed to percentages, default is TRUE
#'        TRUE = percentage transformed, FALSE = not percentage transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param weight a parameter passed to the `prot_prop` function
#'
#' @return a data frame containing the individual weighted prot_prop
#'
prot_propData <- function(data, part_col, val_col, group_col = NULL, min.y = NULL, weight = c(0.00,0.05,0.10,0.30,0.60,0.75,0.85,0.90,0.95,0.99), var_trans = 1, weight_trans = 1, base = 2, adj = 10) {

  source('sermetric/R/prot_prop.R')

  if(is.null(group_col)){
    stat <- unique(data[, part_col, drop = FALSE])

    Prot_prop <- sapply(seq_len(nrow(stat)), function(i) {
      prot_propCal(stat[i, part_col], NULL, data, part_col, val_col, group_col, min.y, weight, var_trans, weight_trans, base, adj)
    })
  }else{
    stat <- unique(data[, c(part_col, group_col), drop = FALSE])

    Prot_prop <- sapply(seq_len(nrow(stat)), function(i) {
      prot_propCal(stat[i, part_col], stat[i, group_col], data, part_col, val_col, group_col, min.y, weight, var_trans, weight_trans, base, adj)
    })
  }

  stat <- cbind(stat, Prot_prop)

  return(stat)
}

#' Calculate the weighted prot_prop for a given participant
#'
#' @param i a participant ID
#' @param j a sample group number
#' @param data the input data frame
#' @param part_col a vector identifying participant id
#' @param val_col a vector of target values. Eg. titers
#' @param group_col a vector identifying the sample group, default is NULL
#' @param min.y the smallest titer value(s), default is the minimum value in the titer vector
#' @param weight a vector of weight values according to titer values, default is c(0.00,0.05,0.10,0.30,0.60,0.75,0.85,0.90,0.95,0.99)
#' @param var_trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param weight_trans whether the input weight value is transformed to percentages, default is TRUE
#'        TRUE = percentage transformed, FALSE = not percentage transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#'
#' @return the weighted prot_prop for the given participant
prot_propCal <- function(i, j,data, part_col, val_col, group_col, min.y, weight, var_trans, weight_trans, base, adj) {

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

  prot_prop <- prot_prop(y, min.y, weight, var_trans, weight_trans, base, adj)

  return(prot_prop)
}
