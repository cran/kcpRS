% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runMean.R
\name{runMean}
\alias{runMean}
\title{Running Means}
\usage{
runMean(data, wsize = 25)
}
\arguments{
\item{data}{\emph{N} x \emph{v} dataframe where \emph{N} is the no. of time points and \emph{v} the no. of variables}

\item{wsize}{Window size}
}
\value{
Running means time series
}
\description{
Extracts the running means by sliding a window comprised of \code{wsize} time points, and in each window, the mean for each variable is computed.
Each time the window is slid, the oldest time point is discarded and the latest time point is added.
}
\examples{
phase1=cbind(rnorm(50,0,1),rnorm(50,0,1)) #phase1: Means=0
phase2=cbind(rnorm(50,1,1),rnorm(50,1,1)) #phase2: Means=1
X=rbind(phase1,phase2)
RS=runMean(data=X,wsize=25)
ts.plot(RS, gpars=list(xlab="Window", ylab="Means", col=1:2,lwd=2))

}
