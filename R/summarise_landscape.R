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
#'
#' @examples
#' \donttest{
#' data(Fonville)
#' summarise_landscape(Fonville, part_col = "Subject Number", weight_col = "isolation_year",
#'             val_col = "titer", group_col = "Sample Year",
#'             required_metrics = c("GMT", "AUC"))
#' }
#'
#' @export
summarise_landscape <- function(data, part_col, weight_col, val_col, group_col = NULL,
                        min.y = NULL, weight = c(0.00,0.05,0.10,0.30,0.60,0.75,0.85,0.90,0.95,0.99),
                        var_trans = TRUE, max_titer_output_trans = TRUE, gmt_output_trans = TRUE, weight_trans = TRUE,
                        base = 2, adj = 10, mode = NULL, adjust = NULL,
                        required_metrics = c("ATY","AUC","gini_coefficient","GMT","Kurtosis","Max_titer","Prop2","Prot_prop","Skewness","Width2"),
                        col_names = NULL) {


  if(is.null(col_names)){
    col_names = required_metrics
  }

  Overall <- NULL

  valid_metrics <- c("ATY", "AUC", "gini_coefficient", "GMT", "Kurtosis",
                     "Max_titer", "Prot_prop", "Skewness")

  for (i in seq_along(required_metrics)) {
    metric <- required_metrics[i]
    col_name <- col_names[i]

    is_prop  <- substr(metric, 1, 4) == "Prop"  && nchar(metric) > 4
    is_width <- substr(metric, 1, 5) == "Width" && nchar(metric) > 5
    if (!(metric %in% valid_metrics || is_prop || is_width)) {
      stop(paste0(
        "Unknown metric: \"", metric, "\". ",
        "Valid values are: ", paste(valid_metrics, collapse = ", "),
        ', or "Prop<n>" / "Width<n>" (e.g. "Prop2", "Width2").'
      ))
    }

    if (metric == "ATY") {
      Result <- aty_data(data = data,part_col = part_col,weight_col = weight_col, val_col = val_col, group_col = group_col, var_trans = var_trans, base = base, adj = adj, mode = mode, adjust = adjust, aty_col = col_name)
    }else if (metric == "AUC"){
      Result <- auc_data(data = data, part_col = part_col,weight_col = weight_col, val_col = val_col, group_col = group_col, var_trans = var_trans, base = base, adj = adj, mode = mode, adjust = adjust, auc_col = col_name)
    }else if (metric == "gini_coefficient"){
      Result <- gini_data(data = data, part_col = part_col, val_col = val_col,group_col = group_col, var_trans = var_trans, base = base, adj = adj, adjust = adjust, gini_col = col_name)
    }else if (metric == "GMT"){
      Result <- gmt_data(data = data, part_col = part_col, val_col = val_col, group_col = group_col, var_trans = var_trans, output_trans = gmt_output_trans, base = base, adj = adj, gmt_col = col_name)
    }else if (metric == "Kurtosis"){
      Result <- kurtosis_data(data = data, part_col = part_col,weight_col = weight_col, val_col = val_col,group_col = group_col, var_trans = var_trans, base = base, adj = adj, mode = mode, adjust = adjust, kurt_col = col_name)
    }else if (metric == "Max_titer"){
      Result <- max_titer_data(data = data, part_col = part_col, val_col = val_col, group_col = group_col, var_trans = var_trans, output_trans = max_titer_output_trans, base = base, adj = adj, max_titer_col = col_name)
    }else if (substr(metric, 1, 4) == "Prop"){
      threshold <- as.numeric(substr(metric, 5, nchar(metric)))
      Result <- prop_data(data = data, part_col = part_col, val_col = val_col, group_col = group_col, var_trans = var_trans, threshold_trans = TRUE, base = base, adj = adj, threshold = threshold, prop_col = col_name)
    }else if(metric == "Prot_prop"){
      Result <- prot_prop_data(data = data, part_col = part_col, val_col = val_col, group_col = group_col, min.y = min.y, weight = weight, var_trans = var_trans, weight_trans = weight_trans, base = base, adj = adj, prot_prop_col = col_name)
    }else if(metric == "Skewness"){
      Result <- skewness_data(data = data, part_col = part_col,weight_col = weight_col, val_col = val_col, group_col = group_col, var_trans = var_trans, base = base, adj = adj , mode = mode, adjust = adjust, skew_col = col_name)
    }else if (substr(metric, 1, 5) == "Width"){
      threshold <- as.numeric(substr(metric, 6, nchar(metric)))
      Result <- width_data(data = data, part_col = part_col, weight_col = weight_col, val_col = val_col, group_col = group_col, var_trans = var_trans, threshold_trans = TRUE, base = base, adj = adj, threshold = threshold, mode = mode, adjust = adjust, width_col = col_name)
    }

    if (is.null(Overall)) {
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
