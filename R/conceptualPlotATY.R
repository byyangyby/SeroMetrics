#'
#' Conceptual plot of an antibody profile regarding ATY.
#'
#' @param i, panel i
#' @param ids, a vector of randomly selected participant id
#' @param df, dataframe
#' @param part_col a vector identifying participant id
#' @param weight_col a vector of the distance used to sorting titers. Eg. year of isolation or genetic distance
#' @param val_col a vector of target values. Eg. titers
#' @param group_col a vector identifying the sample group, default is NULL
#' @param groups a list of group numbers, default is NULL
#' @param age_col a vector of the age of participants, used for plotting the title of the plot, default is NULL
#' @param var_trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param mode how to calculate values if weight_col is repeated, default is "mean"
#' @param adjust how to adjust titers when they are too small to plot, default is NULL
#' @param x.min the start point of the x-axis, default is 1965
#' @param x.max the end point of the x-axis, default is 2015
#' @param x.diff.1 the larger interval of the segments on the x-axis, default is 5
#' @param x.diff.2 the smaller interval of the segments on the x-axis, default is 1
#' @param y.min the start point of the x-axis, default is 0
#' @param y.max the end point of the x-axis, default if 10
#' @param y.diff the interval of the segments on the y-axis, default is 1
#' @param legends the legends of the plot, default is NULL
#' @param xlab the label of the x-axis, default is "Year of isolation"
#' @param ylab the label of the y-axis, default is "Titers"
#' @param color.1 the color of the vertical unweighted average line, default is "#5bd0c8"
#' @param color.2 the color of the vertical ATY line, default is "#5bd0c8"
#' @param color.line the color of the curve, default is "#999999"
#' @param line.type.1 the type of the vertical unweighted average line, default is 2
#' @param line.type.2 the type of the vertical ATY line, default is 1
#' @param cex the font size of the legend, default is 1.5
#' @param cex.lab the font size of the lab of x- and y-axis, default is 2.0
#' @param cex.main the font size of the title of the figures, default is 2.0
#' @param cex.axis the font size of the segments of x- and y-axis, default is 1.5
#' @param mar the margin of the figure, default is c(6,6,5,2)
#' @param mgp the layout of the figure, default is c(3,1,4)
#' @param width.dataline the width of the curve, default is 2.0
#' @param width.1 the width of the vertical unweighted average line, default is 3
#' @param width.2 the width of the vertical ATY line, default is 3
#' @param legend.font the font size of the legend, default is 2
#'
#' @return the conceptual plot regarding ATY
#' @export

conceptualPlotATY <- function(i, ids, df, part_col, weight_col, val_col, group_col = NULL, groups = NULL,
                              age_col = NULL, var_trans = 1, base = 2, adj = 10, mode = "mean", adjust = NULL,
                              x.min = 1965,x.max = 2015,x.diff.1 = 5, x.diff.2 = 1,y.min = 0,y.max = 8, y.diff = 1,
                              legends = c("C"), xlab = "Year of isolation", ylab = "Titers",
                              color.1 = "#5bd0c8", color.2 = "#5bd0c8", color.line = "#999999", line.type.1 = 2, line.type.2 = 1,
                              cex = 1.5, cex.lab = 2.0, cex.main = 2.0, cex.axis = 1.5, mar = c(6,6,5,2), mgp = c(3,1,4),
                              width.dataline = 2, width.1 = 3, width.2 = 3, legend.font = 2){

  id = ids[i]

  legend = legends[i]

  if(!var_trans){
    df[[val_col]] = log(df[[val_col]]/adj,base)
  }

  if(!is.null(adjust)){
    df[[val_col]] = df[[val_col]] - adjust
  }

  if(!is.null(group_col)){

    if (mode == "mean") {
      df <- df %>%
        select(all_of(part_col), all_of(group_col), all_of(weight_col), all_of(val_col), all_of(age_col)) %>%
        group_by(across(all_of(c(part_col, group_col, weight_col, age_col)))) %>%
        summarise(across(all_of(val_col), mean, na.rm = TRUE)) %>%
        ungroup()
    } else if (mode == "max") {
      df <- df %>%
        select(all_of(part_col), all_of(group_col), all_of(weight_col), all_of(val_col), all_of(age_col)) %>%
        group_by(across(all_of(c(part_col, weight_col, age_col)))) %>%
        summarise(across(all_of(val_col), max, na.rm = TRUE)) %>%
        ungroup()
    } else if (mode == "min") {
      df <- df %>%
        select(all_of(part_col), all_of(group_col), all_of(weight_col), all_of(val_col), all_of(age_col)) %>%
        group_by(across(all_of(c(part_col, group_col, weight_col, age_col)))) %>%
        summarise(across(all_of(val_col), min, na.rm = TRUE)) %>%
        ungroup()
    }

    group = groups[i]
    con <- df[[part_col]] == id & df[[group_col]] == group

  }else{

    if (mode == "mean") {
      df <- df %>%
        select(all_of(part_col), all_of(weight_col), all_of(val_col), all_of(age_col)) %>%
        group_by(across(all_of(c(part_col, weight_col, age_col)))) %>%
        summarise(across(all_of(val_col), mean, na.rm = TRUE)) %>%
        ungroup()
    } else if (mode == "max") {
      df <- df %>%
        select(all_of(part_col), all_of(weight_col), all_of(val_col), all_of(age_col)) %>%
        group_by(across(all_of(c(part_col, weight_col, age_col)))) %>%
        summarise(across(all_of(val_col), max, na.rm = TRUE)) %>%
        ungroup()
    } else if (mode == "min") {
      df <- df %>%
        select(all_of(part_col), all_of(weight_col), all_of(val_col), all_of(age_col)) %>%
        group_by(across(all_of(c(part_col, weight_col, age_col)))) %>%
        summarise(across(all_of(val_col), min, na.rm = TRUE)) %>%
        ungroup()
    }

    con <- df[[part_col]] == id

  }

  used <- df[con, c(weight_col, val_col)]
  used <- used[order(used[[weight_col]]),]

  # par(mar = c(4,6,2,0))
  par(mar = mar, mgp = mgp)

  if(!is.null(age_col)){

    age <- unique(df[[age_col]][con])

    plot(NULL,
         xlim = c(x.min, x.max),
         ylim = c(y.min, y.max),
         las = 1,
         xlab = xlab,
         ylab = ylab,
         main = paste(age, " years",
                      sep = ""),
         axes = FALSE,
         cex.lab = cex.lab,
         cex.main= cex.main)

  }else{
    plot(NULL,
         xlim = c(x.min, x.max),
         ylim = c(y.min, y.max),
         las = 1,
         xlab = xlab,
         ylab = ylab,
         axes = FALSE,
         cex.lab = cex.lab)
  }

  # obserbed data
  x <- used[[weight_col]]
  y <- used[[val_col]]

  lines(x, y, lwd = width.dataline, col = color.line, type = "b")

  selected_df <- df[con,]

  stat <- atyData(selected_df,part_col = part_col, weight_col = weight_col, val_col = val_col, var_trans = 1, mode = "mean")

  aty <- as.numeric(stat$ATY)

  aty.mean <- mean(used[[weight_col]])

  segments(aty.mean, y.min, aty.mean, y.max, col = color.1, lwd = width.1, lty = line.type.1)
  segments(aty, y.min, aty, y.max, col = color.2, lwd = width.2, lty = line.type.2)

  axis(1, seq(x.min, x.max, x.diff.1),
       pos = y.min - 0.5,
       cex.axis = cex.axis)
  axis(1, seq(x.min, x.max, x.diff.2),
       pos = y.min - 0.5, labels = NA,
       tcl = -0.2, col.ticks = 1, col=NA,
       cex.axis = cex.axis)

  axis(2, seq(y.min, y.max, y.diff),
       pos = x.min - 0.5,
       las=1,
       labels = seq(y.min, y.max, y.diff),
       cex.axis = cex.axis)

  mtext(legend,
    side = 2,
    line = 4.5,
    at = y.max+(y.max-y.min)*0.1,
    outer = FALSE,
    las = 1,
    cex = cex,
    font = legend.font)

}
