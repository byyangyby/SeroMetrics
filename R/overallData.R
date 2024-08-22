#' Function to format all required metrics of the data
#'
#' @param data a data frame
#' @param part_col a vector identifying participant id
#' @param weight_col a vector of the distance used to sorting titers. Eg. year of isolation or genetic distance
#' @param val_col a vector of target values. Eg. titers
#' @param group_col a vector identifying the sample group, default is NULL
#' @param min.y the smallest titer value(s), default is the minimum value in the titer vector(calculated later, temporarily being NULL in the formula)
#' @param weight a vector of weight values according to titer values, default is c(0.00,0.05,0.10,0.30,0.60,0.75,0.85,0.90,0.95,0.99)
#' @param var_trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param max_titer_output_trans whether the output max titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param gmt_output_trans whether the output gmt value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param weight_trans whether the input weight value is transformed to percentages, default is TRUE
#'        TRUE = percentage transformed, FALSE = not percentage transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param mode how to calculate values if weight_col is repeated, default is NULL
#' @param adjust how to adjust titers, default is the minimum of titers
#' @param required_metrics a list of metrics that required by the user. Note that when inputting proportion and width, it should be in the format of "Prop"/"Width" followed by a log_transformed threshold, eg. "Prop2"/"Width2". The other metrics should be exactly "ATY","AUC","gini_coefficient","GMT","Kurtosis","Max_titer","Prot_prop", or "Skewness", i.e. inputting "aty","gini", etc. is considered invalid.
#' @param col_names the column names of the output metrics, default is required_metrics
#'
#' @return a data frame containing all metrics
#' @export
overallData <- function(data, part_col, weight_col, val_col, group_col = NULL,
                        min.y = NULL, weight = c(0.00,0.05,0.10,0.30,0.60,0.75,0.85,0.90,0.95,0.99),
                        var_trans = 1, max_titer_output_trans = 1, gmt_output_trans = 1, weight_trans = 1,
                        base = 2, adj = 10, mode = NULL, adjust = NULL,
                        required_metrics = c("ATY","AUC","gini_coefficient","GMT","Kurtosis","Max_titer","Prop2","Prot_prop","Skewness","Width2"),
                        col_names = NULL) {


  if(is.null(col_names)){
    col_names = required_metrics
  }

  for (i in seq_along(required_metrics)) {
    metric <- required_metrics[i]
    col_name <- col_names[i]
    if (metric == "ATY") {
      Result <- atyData(data = data,part_col = part_col,weight_col = weight_col, val_col = val_col, group_col = group_col, var_trans = var_trans, base = base, adj = adj, mode = mode, adjust = adjust, aty_col = col_name)
    }else if (metric == "AUC"){
      Result <- aucData(data = data, part_col = part_col,weight_col = weight_col, val_col = val_col, group_col = group_col, var_trans = var_trans, base = base, adj = adj, mode = mode, adjust = adjust, auc_col = col_name)
    }else if (metric == "gini_coefficient"){
      Result <- giniData(data = data, part_col = part_col, val_col = val_col,group_col = group_col, var_trans = var_trans, base = base, adj = adj, adjust = adjust, gini_col = col_name)
    }else if (metric == "GMT"){
      Result <- gmtData(data = data, part_col = part_col, val_col = val_col, group_col = group_col, var_trans = var_trans, output_trans = gmt_output_trans, base = base, adj = adj, gmt_col = col_name)
    }else if (metric == "Kurtosis"){
      Result <- kurtosisData(data = data, part_col = part_col,weight_col = weight_col, val_col = val_col,group_col = group_col, var_trans = var_trans, base = base, adj = adj, mode = mode, adjust = adjust, kurt_col = col_name)
    }else if (metric == "Max_titer"){
      Result <- max_titerData(data = data, part_col = part_col, val_col = val_col, group_col = group_col, var_trans = var_trans, output_trans = max_titer_output_trans, base = base, adj = adj, max_titer_col = col_name)
    }else if (substr(metric, 1, 4) == "Prop"){
      threshold <- as.numeric(substr(metric, 5, nchar(metric)))
      Result <- propData(data = data, part_col = part_col, val_col = val_col, group_col = group_col, var_trans = var_trans, threshold_trans = 1, base = base, adj = adj, threshold = threshold, prop_col = col_name)
    }else if(metric == "Prot_prop"){
      Result <- prot_propData(data = data, part_col = part_col, val_col = val_col, group_col = group_col, min.y = min.y, weight = weight, var_trans = var_trans, weight_trans = weight_trans, base = base, adj = adj, prot_prop_col = col_name)
    }else if(metric == "Skewness"){
      Result <- skewnessData(data = data, part_col = part_col,weight_col = weight_col, val_col = val_col, group_col = group_col, var_trans = var_trans, base = base, adj = adj , mode = mode, adjust = adjust, skew_col = col_name)
    }else if (substr(metric, 1, 5) == "Width"){
      threshold <- as.numeric(substr(metric, 6, nchar(metric)))
      Result <- widthData(data = data, part_col = part_col,weight_col = weight_col, val_col = val_col, group_col = group_col, var_trans = var_trans, threshold_trans = 1, base = base, adj = adj, threshold = threshold, mode = mode, adjust = adjust, width_col = col_name)
    }

    if (!exists("Overall") || is.null(Overall)) {
      Overall <- Result
    } else {
      if (is.null(group_col)){
        Overall <- left_join(Overall, Result, by = part_col)
      }else{
        Overall <- left_join(Overall, Result, by = c(part_col, group_col))
      }
    }
  }
  return(Overall)

}
