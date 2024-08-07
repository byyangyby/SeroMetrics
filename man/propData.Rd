% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/propData.R
\name{propData}
\alias{propData}
\title{Function to format data and calculate the weighted prop}
\usage{
propData(
  data,
  part_col,
  val_col,
  group_col = NULL,
  var_trans = 1,
  threshold_trans = 1,
  base = 2,
  adj = 10,
  threshold = 2,
  prop_col = "Proportion"
)
}
\arguments{
\item{data}{a data frame containing the necessary columns (participant_id titer)}

\item{part_col}{a vector identifying participant id}

\item{val_col}{a vector of target values. Eg. titers}

\item{group_col}{a vector identifying the sample group, default is NULL}

\item{var_trans}{whether the input titer value is log transformed, default is TRUE
TRUE = log transformed, FALSE = not log transformed}

\item{threshold_trans}{whether the threshold value is log transformed, default is TRUE
TRUE = log transformed, FALSE = not log transformed}

\item{base}{log base used to calculate the log transformed value, default is 2}

\item{adj}{dividing factor before taking log transformation, default is 10}

\item{threshold}{the titers above which are counted, default is 2}

\item{prop_col}{the column name of the proportion column, default is "Proportion"}
}
\value{
a data frame containing the individual weighted prop
}
\description{
Function to format data and calculate the weighted prop
}
