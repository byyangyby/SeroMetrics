#'
#' Conceptual plot of an antibody profile regarding AUC.
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
#' @param color.1 the color of the area above the 0-axis and below the curve, default is "#7b896a"
#' @param color.2 the color of the area below the 0-axis and above the curve, default is "#c3bc1a"
#' @param color.line the color of the curve and the vertical lines from the data points to the 0-axis, default is "#999999"
#' @param line.type the type of the lines, default is 3
#' @param cex the font size of the legend, default is 1.5
#' @param cex.lab the font size of the lab of x- and y-axis, default is 2.0
#' @param cex.main the font size of the title of the figures, default is 2.0
#' @param cex.axis the font size of the segments of x- and y-axis, default is 1.5
#' @param mar the margin of the figure, default is c(6,6,5,2)
#' @param mgp the layout of the figure, default is c(3,1,4)
#' @param width.dataline the width of the curve, default is 2
#' @param legend.font the font size of the legend, default is 2
#'
#' @return the conceptual plot regarding AUC
#' @export

conceptualPlotAUC <- function(i, ids, df, part_col, weight_col, val_col, group_col = NULL, groups = NULL,
                              age_col = NULL, var_trans = 1,base = 2, adj = 10, mode = "mean", adjust = NULL,
                              x.min = 1965,x.max = 2015,x.diff.1 = 5, x.diff.2 = 1,
                              y.min = 0,y.max = 8, y.diff = 1, legends = c("A"), xlab = "Years of isolation", ylab = "Titers",
                              color.1 = "#7b896a", color.2 = "#c3bc1a", color.line = "#999999", line.type = 3,
                              cex = 1.5, cex.lab = 2.0, cex.main = 2.0, cex.axis = 1.5, mar = c(6,6,5,2), mgp = c(3,1,4),
                              width.dataline = 2, legend.font = 2){

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
  used <- used[order(used[[weight_col]]), ]

  cols <- c(color.1, # green
            color.2) # dark green


  # observed data
  x <- used[[weight_col]]
  y <- used[[val_col]]

  par(mar = mar,
      mgp = mgp)

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
         cex.lab = cex.lab,
         cex.main= cex.main)
  }

  for (i in 1:(length(x) - 1)) {

    if (y[i] >= 0 & y[i + 1] >= 0){
      polygon(x = c(x[i],x[i + 1],x[i + 1],x[i]),
              y = c(0,0,y[i + 1],y[i]),
              border = NA,
              col = adjustcolor(cols[1], alpha.f = 0.5))
    }else if (y[i] < 0 & y[i + 1] < 0){
      polygon(x = c(x[i],x[i + 1],x[i + 1],x[i]),
              y = c(0,0,y[i + 1],y[i]),
              border = NA,
              col = adjustcolor(cols[2], alpha.f = 0.5))
    }else if (y[i] >= 0 & y[i + 1] < 0){
      x_intersect <- x[i] + (0 - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
      polygon(x = c(x[i], x_intersect, x[i]),
              y = c(0, 0, y[i]),
              border = NA,
              col = adjustcolor(cols[1], alpha.f = 0.5))
      polygon(x = c(x[i + 1], x[i + 1], x_intersect),
              y = c(y[i + 1], 0, 0),
              border = NA,
              col = adjustcolor(cols[2], alpha.f = 0.5))
    }else if (y[i] < 0 & y[i + 1] >= 0){
      x_intersect <- x[i] + (0 - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
      polygon(x = c(x[i], x_intersect, x[i]),
              y = c(0, 0, y[i]),
              border = NA,
              col = adjustcolor(cols[2], alpha.f = 0.5))
      polygon(x = c(x[i + 1], x[i + 1], x_intersect),
              y = c(y[i + 1], 0, 0),
              border = NA,
              col = adjustcolor(cols[1], alpha.f = 0.5))
    }
  }

  lines(x, y, lwd = width.dataline, col = color.line, type = "b")

  len <- length(x)

  sapply(1:(len-1), function(i){

    polygon(x = c(x[i], x[i], x[i+1], x[i+1]),
            y = c(y[i],0,  0, y[i + 1]),
            col = NULL,
            border = color.line,
            lty = line.type)

  })

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
