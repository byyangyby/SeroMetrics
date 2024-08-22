#' Function to plot comparisons between two antibody profiles.
#'
#' @param i, panel i
#' @param ids, a vector of randomly selected participant id
#' @param df, dataframe
#' @param part_col a vector identifying participant id, default is "id"
#' @param weight_col a vector of the distance used to sorting titers. Eg. year of isolation or genetic distance
#' @param val_col a vector of target values. Eg. titers
#' @param group_col a vector identifying the sample group, default is NULL
#' @param groups.1 a list of group numbers of the first group of profiles, default is NULL
#' @param groups.2 a list of group numbers of the second group of profiles, default is NULL
#' @param age_col a vector of the age of participants, used for plotting the title of the plot, default is NULL
#' @param birth_col a vector identifying the birth year of the participants, used for plotting the birth year of the participant, default is NULL
#' @param sample_years.1 a list of sample years for plotting vertical bars of the first profiles, default is NULL
#' @param sample_years.2 a list of sample years for plotting vertical bars of the second profiles, default is NULL
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
#' @param color.1 the color of the curve, the data points, the vertical bar of sample year, and the vertical line of birth year, for the first group. default is "#3f93aa"
#' @param color.2 the color of the curve, the data points, the vertical bar of sample year, and the vertical line of birth year, for the second group. default is "#fa8072"
#' @param color.main the color of the area under the curves, default is "#867b7b"
#' @param color.change.1 the color of the area above the curve of the first group and below the curve of the second group, default is "#590059"
#' @param color.change.2 the color of the area above the curve of the second group and below the curve of the first group, default is "#8BAD23"
#' @param line.type.1 the line type of the curve of the first group, default is 1
#' @param line.type.2 the line type of the curve of the second group, default is 3
#' @param point.type.1 the point type of the data points of the first group, default is 1
#' @param point.type.2 the point type of the data points of the second group, default is 2
#' @param line.type.birth.1 the line type of the birth line of the first group, default is 2
#' @param line.type.birth.2 the line type of the birth line of the second group, default is 4
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

dualPlot <- function(i, ids.1, ids.2, df, part_col, weight_col, val_col, group_col = NULL, groups.1 = NULL, groups.2 = NULL,
                     age_col = NULL, birth_col = NULL, sample_years.1 = NULL, sample_years.2 = NULL,
                     mode = "mean", var_trans = 1, base = 2, adj = 10, adjust = NULL,
                     x.min = 1965, x.max = 2015, x.diff.1 = 5, x.diff.2 = 1,
                     y.min = 0, y.max = 10, y.diff = 1, legends = c("A"), xlab = "Distance", ylab = "Titers",
                     color.1 = "#3f93aa", color.2 = "#fa8072", color.main = "#867b7b", color.change.1 = "#590059", color.change.2 = "#8BAD23",
                     line.type.1 = 1, line.type.2 = 3, point.type.1 = 1, point.type.2 = 2, line.type.birth.1 = 2, line.type.birth.2 = 4,
                     cex = 1.5, cex.lab = 2.0, cex.main = 2.0, cex.axis = 1.5, mar = c(6,6,5,2), mgp = c(3,1,4),
                     width.dataline = 2.0, width.points = 1.5, legend.font = 2){

  legend = legends[i]

  if(!var_trans){
    df[[val_col]] = log(df[[val_col]]/adj, base)
  }

  if(!is.null(adjust)){
    df[[val_col]] = df[[val_col]] - adjust
  }

  if(!is.null(sample_years.1)){
    sample_year.1 = sample_years.1[i]
    sample_year.2 = sample_years.2[i]
  }

  id.1 <- ids.1[i]
  id.2 <- ids.2[i]

  if(!is.null(group_col)){

    group.1 <- groups.1[i]
    group.2 <- groups.2[i]

    con.1 <- df[[part_col]] == id.1 & df[[group_col]] == group.1
    con.2 <- df[[part_col]] == id.2 & df[[group_col]] == group.2

  }else{

    con.1 <- df[[part_col]] == id.1
    con.2 <- df[[part_col]] == id.2

  }

  data <- NULL

  data$year <- unique(df[[weight_col]])
  data$titers.1 <- sapply(data$year, function(x) {

    if(mode == "mean"){
      a <- mean(df[[val_col]][con.1 & df[[weight_col]] == x], na.rm = TRUE)
    }else if(mode == "min"){
      a <- min(df[[val_col]][con.1 & df[[weight_col]] == x])
    }else if(mode == "max"){
      a <- max(df[[val_col]][con.1 & df[[weight_col]] == x])
    }

  })

  data$titers.2 <- sapply(data$year, function(x) {

    if(mode == "mean"){
      a <- mean(df[[val_col]][con.2 & df[[weight_col]] == x], na.rm = TRUE)
    }else if(mode == "min"){
      a <- min(df[[val_col]][con.2 & df[[weight_col]] == x])
    }else if(mode == "max"){
      a <- max(df[[val_col]][con.2 & df[[weight_col]] == x])
    }

  })

  data$delta_titers <- sapply(data$year, function(x) {
    if(mode == "mean"){
      a <- mean(df[[val_col]][con.1 & df[[weight_col]] == x], na.rm = TRUE)
    }else if(mode == "min"){
      a <- min(df[[val_col]][con.1 & df[[weight_col]] == x])
    }else if(mode == "max"){
      a <- max(df[[val_col]][con.1 & df[[weight_col]] == x])
    }
    if(mode == "mean"){
      b <- mean(df[[val_col]][con.2 & df[[weight_col]] == x], na.rm = TRUE)
    }else if(mode == "min"){
      b <- min(df[[val_col]][con.2 & df[[weight_col]] == x])
    }else if(mode == "max"){
      b <- max(df[[val_col]][con.2 & df[[weight_col]] == x])
    }
    return(b-a)
  })

  data <- data.frame(data)

  rownames(data) <- NULL

  data <- data[!is.na(data$titers.1), ]
  data <- data[order(data$year), ]

  # splined data
  sm.1 <- spline(data$year, data$titers.1, method = "natural")
  sm.2 <- spline(data$year, data$titers.2, method = "natural")

  data.spline <- data.frame(year = sm.1$x,
                            titers.1 = sm.1$y,
                            titers.2 = sm.2$y)

  # replace negative with 0
  data.spline$titers.1[data.spline$titers.1<0] <- 0
  data.spline$titers.2[data.spline$titers.2<0] <- 0

  data.spline$increase <- data.spline$titers.2 - data.spline$titers.1
  data.spline$increase[data.spline$increase<=0] <- 0
  data.spline$decrease <- data.spline$titers.2 - data.spline$titers.1
  data.spline$decrease[data.spline$decrease>=0] <- 0

  cols <- c(color.1, color.2)

  cols.change <- c(color.change.1, color.change.2)

  par(mar = mar, mgp = mgp)

  if(!is.null(age_col)){

    age.1 <- unique(df[[age_col]][con.1])
    age.2 <- unique(df[[age_col]][con.2])

    plot(NULL,
         xlim = c(x.min, x.max),
         ylim = c(y.min, y.max),
         las = 1,
         xlab = c("Distance"),
         ylab = c("Titers"),
         main = paste(age.1, " years vs. ", age.2, " years",
                      sep = ""),
         axes = FALSE,
         cex.lab = cex.lab,
         cex.main= cex.main)

  }else{

    plot(NULL,
         xlim = c(x.min, x.max),
         ylim = c(y.min, y.max),
         las = 1,
         xlab = c("Distance"),
         ylab = c("Titers"),
         axes = FALSE,
         cex.lab = cex.lab,
         cex.main= cex.main)

  }

  # time period of visit
  if(!is.null(sample_years.1)){
    rect(sample_year.1 - 0.5,
         y.min,
         sample_year.1 + 0.5,
         y.max,
         col = adjustcolor(cols[1],
                           alpha.f = 0.2),
         border = FALSE)

    rect(sample_year.2 - 0.5,
         y.min,
         sample_year.2 + 0.5,
         y.max,
         col = adjustcolor(cols[2],
                           alpha.f = 0.2),
         border = FALSE)
  }

  # # titer reference line
  # segments(x.min, 2,
  #          x.max, 2,
  #          lty = 2)
  # segments(x.min, 4,
  #          x.max, 4,
  #          lty = 3)

  # brown area
  rowMins <- sapply(1:nrow(data.spline), function(i){
    max(min(data.spline[i, 2:3]), 0)
  })

  polygon(x = c(data.spline$year, rev(data.spline$year)),
          y = c(rep(0, length(data.spline$year)),
                rev(rowMins)),
          border = NA,
          col = adjustcolor(color.main,
                            alpha.f = 0.5))
  for(j in 1:2){

    x <- data.spline$year
    y <- data.spline$titers.1 + data.spline[ ,c("increase", "decrease")[j]]

    polygon(x = c(x, rev(x)),
            y = c(data.spline$titers.1, rev(y)),
            border = NA,
            col = adjustcolor(cols.change[j],
                              alpha.f = 0.7))

    if(j == 1){

      points(data$year, data$titers.1,
             pch = point.type.1,
             col = cols[j],
             cex = cex,
             lwd = width.points)

      lines(x, data.spline$titers.1,
            lty = line.type.1,
            col = cols[j],
            lwd = width.dataline)

    }else if(j == 2){

      points(data$year, data$titers.2,
             pch = point.type.2,
             col = cols[j],
             cex = cex,
             lwd = width.points)

      lines(x, data.spline$titers.2,
            lty = line.type.2,
            col = cols[j],
            lwd = width.dataline)

    }

  }

  # time period before birth
  if(!is.null(birth_col)){
    birth.1 <- min(unique(df[[birth_col]][df[[part_col]] == id.1]))
    birth.2 <- min(unique(df[[birth_col]][df[[part_col]] == id.2]))

    segments(birth.1, y.min, birth.1, y.max, lty = line.type.birth.1, col = color.1)
    segments(birth.2, y.min, birth.2, y.max, lty = line.type.birth.2, col = color.2)

  }

  # add axis
  axis(1, seq(x.min, x.max, x.diff.1),
       pos = y.min - 0.5, cex.axis = cex.axis)
  axis(1, seq(x.min, x.max, x.diff.2),
       pos = y.min - 0.5, labels = NA,
       tcl = -0.2, col.ticks = 1,
       col=NA, cex.axis = cex.axis)

  axis(2, seq(y.min, y.max, 1),
       pos = x.min - 0.5,
       las=1,
       labels = c(0, seq(y.min+1, y.max, 1)),
       cex.axis = cex.axis)

  mtext(legend,
        side = 2,
        line = 3.5,
        at = y.max + (y.max - y.min)*0.05,
        outer = FALSE,
        las = 1,
        cex = cex,
        font = legend.font)

}
