% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skewnessData.R
\name{skewCal}
\alias{skewCal}
\title{Calculate the skewness for a given participant}
\usage{
skewCal(
  i,
  j,
  data,
  part_col,
  weight_col,
  val_col,
  group_col,
  var_trans,
  base,
  adj,
  mode,
  adjust
)
}
\arguments{
\item{i}{a participant ID}

\item{j}{a sample group number}

\item{data}{the input data frame}

\item{part_col}{a vector identifying participant id}

\item{weight_col}{a vector of the distance used to sorting titers. Eg. year of isolation or genetic distance}

\item{val_col}{a vector of target values. Eg. titers}

\item{group_col}{a vector identifying the sample group, default is NULL}

\item{var_trans}{whether the input titer value is log transformed, default is TRUE
TRUE = log transformed, FALSE = not log transformed}

\item{base}{log base used to calculate the log transformed value, default is 2}

\item{adj}{dividing factor before taking log transformation, default is 10}

\item{mode}{how to calculate values if weight_col is repeated, default is NULL}

\item{adjust}{how to adjust titers, default is the minimum of titers}
}
\value{
the skewness value for the given participant
}
\description{
Calculate the skewness for a given participant
}
