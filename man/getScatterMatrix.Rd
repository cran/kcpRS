% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{getScatterMatrix}
\alias{getScatterMatrix}
\title{Get the matrix of optimized scatters used in locating the change points.}
\usage{
getScatterMatrix(II_, X_, H_)
}
\arguments{
\item{II_}{A D x N matrix where D is the maximum no. of segments (Kmax+1) and N is the no. of windows}

\item{X_}{An N x r dataframe where N is the no. of windows and r the no. of running statistics monitored}

\item{H_}{A D x N matrix where D is the maximum no. of segments (Kmax+1) and N is the no. of windows}
}
\value{
\item{II}{A matrix of optimized scatters}

\item{H}{A matrix of candidate changes point locations}

\item{medianK}{Median of the pairwise Euclidean distances}
}
\description{
Get the matrix of optimized scatters used in locating the change points.
}
