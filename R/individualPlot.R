#' Function to plot an individual antibody profile.
#'
#' @param i, panel i
#' @param ids, a vector of randomly selected participant id
#' @param df, dataframe
#' @param part_col a vector identifying participant id, default is "id"
#' @param weight_col a vector of the distance used to sorting titers. Eg. year of isolation or genetic distance
#' @param val_col a vector of target values. Eg. titers
#' @param group_col a vector identifying the sample group, default is NULL
#' @param groups a list of group numbers, default is NULL
#' @param age_col a vector of the age of participants, used for plotting the title of the plot, default is NULL
#' @param birth_col a vector identifying the birth year of the participants, used for plotting the birth year of the participant, default is NULL
#' @param sample_years a list of sample years for plotting vertical bars, default is NULL
#' @param mode how to calculate values if weight_col is repeated, default is "mean"
#' @param var_trans whether the input titer value is log transformed, default is TRUE
#'        TRUE = log transformed, FALSE = not log transformed
#' @param base log base used to calculate the log transformed value, default is 2
#' @param adj dividing factor before taking log transformation, default is 10
#' @param adjust how to adjust titers when they are too small to plot, default is NULL
#' @param x.min the start point of the x-axis, default is 1965
#' @param x.max the end point of the x-axis, default is 2015
#' @param x.diff.1 the larger interval of the segments on the x-axis, default is 5
#' @param x.diff.2 the smaller interval of the segments on the x-axis, default is 1
#' @param y.min the start point of the x-axis, default is 0
#' @param y.max the end point of the x-axis, default if 10
#' @param y.diff the interval of the segments on the y-axis, default is 1
#' @param legends the legends of the plot, default is "A"
#' @param xlab the label of the x-axis, default is "Distance"
#' @param ylab the label of the y-axis, default is "Titers"
#' @param color the color of the curve, the data points, the vertical bar of sample year, and the vertical line of birth year, default is "#3f93aa"
#' @param color.main the color of the area under the curve, default is "#867b7b"
#' @param line.type the line type of the curve, default is 1
#' @param point.type the point type of the data points, default is 1
#' @param line.type.birth the line type of the birth line, default is 2
#' @param cex the font size of the legend, default is 1.5
#' @param cex.lab the font size of the lab of x- and y-axis, default is 2.0
#' @param cex.main the font size of the title of the figures, default is 2.0
#' @param cex.axis the font size of the segments of x- and y-axis, default is 1.5
#' @param mar the margin of the figure, default is c(6,6,5,2)
#' @param mgp the layout of the figure, default is c(3,1,4)
#' @param width.dataline the width of the curve, default is 2.0
#' @param width.points the width of the data points, default is 1.5
#' @param legend.font the font size of the legend, default is 2
#'
#' @return a plot of individual antibody profile
#' @export
#'

individualPlot <- function(i, ids, df, part_col, weight_col, val_col, group_col = NULL, groups = NULL,
                           age_col = NULL, birth_col = NULL, sample_years = NULL,
                           mode = "mean", var_trans = 1, base = 2, adj = 10, adjust = NULL,
                           x.min = 1965, x.max = 2015, x.diff.1 = 5, x.diff.2 = 1,
                           y.min = 0, y.max = 10, y.diff = 1, legends = c("A"), xlab = "Distance", ylab = "Titers",
                           color = "#3f93aa", color.main = "#867b7b", line.type = 1, point.type = 1, line.type.birth = 2,
                           cex = 1.5, cex.lab = 2.0, cex.main = 2.0, cex.axis = 1.5, mar = c(6,6,5,2), mgp = c(3,1,4),
                           width.dataline = 2.0, width.points = 1.5, legend.font = 2){

  legend = legends[i]

  if(!var_trans){
    df[[val_col]] = log(df[[val_col]]/adj, base)
  }

  if(!is.null(adjust)){
    df[[val_col]] = df[[val_col]] - adjust
  }

  if(!is.null(sample_years)){
    sample_year = sample_years[i]
  }

  id <- ids[i]

  if(!is.null(group_col)){
    group <- groups[i]
    con <- df[[part_col]] == id & df[[group_col]] == group
  }else{
    con <- df[[part_col]] == id
  }

  data <- NULL

  data$year <- unique(df[[weight_col]])

  data$titers <- sapply(data$year, function(x) {

    if(mode == "mean"){
      a <- mean(df[[val_col]][con & df[[weight_col]] == x], na.rm = TRUE)
    }else if(mode == "min"){
      a <- min(df[[val_col]][con & df[[weight_col]] == x])
    }else if(mode == "max"){
      a <- max(df[[val_col]][con & df[[weight_col]] == x])
    }

    return(a)

  })

  data <- data.frame(data)
  data <- data[!is.na(data$titers), ]
  data <- data[order(data$year), ]


  # splined data
  sm <- spline(data$year, data$titers, method = "natural")

  data.spline <- data.frame(year = sm$x,
                            titers = sm$y)

  # replace negative with 0
  data.spline$titers[data.spline$titers<0] <- 0

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
         cex.lab = cex.lab,
         cex.main= cex.main)

  }

  if(!is.null(sample_years)){
    # time period of visit
    rect(sample_year - 0.5,
         y.min,
         sample_year + 0.5,
         y.max,
         col = adjustcolor(color,
                           alpha.f = 0.2),
         border = FALSE)

  }

  # under curve area
  rowMins <- sapply(1:nrow(data.spline), function(i){
    max(min(data.spline[i, 2]), 0)
  })

  polygon(x = c(data.spline$year, rev(data.spline$year)),
          y = c(rep(0, length(data.spline$year)),
                rev(rowMins)),
          border = NA,
          col = adjustcolor(color.main,
                            alpha.f = 0.5))

  x <- data.spline$year
  y <- data.spline$titers

  # Draw points
  points(data$year, data$titers,
         pch = point.type,
         col = color,
         cex = cex,
         lwd = width.points)

  # Draw lines
  lines(x, data.spline$titers,
        lty = line.type,
        col = color,
        lwd = width.dataline)

  # time period before birth
  if(!is.null(birth_col)){

    birth <- min(unique(df[[birth_col]][df[[part_col]] == id]))

    segments(birth, y.min, birth, y.max, lty = line.type.birth, col = color)


  }

  # add axis
  axis(1, seq(x.min, x.max, x.diff.1),
       pos = y.min - 0.5, cex.axis = cex.axis)
  axis(1, seq(x.min, x.max, x.diff.2),
       pos = y.min - 0.5, labels = NA,
       tcl = -0.2, col.ticks = 1,
       col=NA, cex.axis = cex.axis)

  axis(2, seq(y.min, y.max, y.diff),
       pos = x.min - 0.5,
       las=1,
       labels = c(0, seq(y.min+1, y.max, y.diff)),
       cex.axis = cex.axis)

  mtext(legend,
        side = 2,
        line = 3.5,
        at = y.max + (y.max - y.min)*0.05,
        outer = FALSE,
        las = 1,
        cex = cex,
        font = 2)

}
