\name{kcpRS-package}
\alias{kcpRS-package}
\docType{package}
\title{KCP on the running statistics}
\details{
This package contains the function \code{\link{kcpRS}} that can accept a user-defined function, \code{RS_fun}, which should derive the running statistics of interest. For examples, see \code{\link{runMean}}, \code{\link{runVar}}, \code{\link{runAR}} and \code{\link{runCorr}}. \code{kcpRS} performs a full change point analysis on the running statistics starting from locating the optimal change points given \emph{k}, significance testing if \emph{k}>0, and finally, determining the most optimal \emph{k}. This function calls the function \code{\link{kcpa}} to find the most optimal change points given \emph{k} and then the \code{\link{permTest}} function to carry out the permutation test. The model selection step is embedded in the \code{kcpRS} function.

This package also contains the function \code{\link{kcpRS_workflow}} which carries out a stepwise change point analysis to flag changes in 4 basic time series statistics: mean, variance, autocorrelation (lag 1) and correlations.

Two illustrative data sets are included: \code{\link{MentalLoad}} and \code{\link{CO2Inhalation}}
}
\description{Flagging change points on a user-specified running statistics through KCP (Kernel Change Point) detection. A KCP permutation test is first implemented to confirm whether there is at least one change point (\emph{k}>0) in the running statistics. If this permutation test is significant, a model selection procedure is implemented to choose the most optimal number of change points.}

\author{Jedelyn Cabrieto (\email{jed.cabrieto@kuleuven.be}) and Kristof Meers

For the core KCP analysis, the authors built upon the codes from the Supplementary Material available in doi:10.1080/01621459.2013.849605 by Matteson and James (2012).

}


\references{
Arlot, S., Celisse, A., & Harchaoui, Z. (2019). A kernel multiple change-point algorithm via model selection. \emph{Journal of Machine Learning Research}, 20(162), 1-56.

Cabrieto, J., Tuerlinckx, F., Kuppens, P., Grassmann, M., & Ceulemans, E. (2017). Detecting correlation changes in multivariate time series: A comparison of four non-parametric change point detection methods. \emph{Behavior Research Methods}, 49, 988-1005. doi:10.3758/s13428-016-0754-9

Cabrieto, J., Tuerlinckx, F., Kuppens, P., Hunyadi, B., & Ceulemans, E. (2018). Testing for the presence of correlation changes
in a multivariate time series: A permutation based approach. \emph{Scientific Reports}, 8, 769, 1-20. doi:10.1038/s41598-017-19067-2

Cabrieto, J., Tuerlinckx, F., Kuppens, P., Wilhelm, F., Liedlgruber, M., & Ceulemans, E. (2018). Capturing correlation changes by applying kernel change point
detection on the running correlations. \emph{Information Sciences}, 447, 117-139. doi:10.1016/j.ins.2018.03.010

Cabrieto, J., Adolf, J., Tuerlinckx, F., Kuppens, P., & Ceulemans, E. (2018). Detecting long-lived autodependency changes in a multivariate system via change point detection
and regime switching models. \emph{Scientific Reports}, 8, 15637, 1-15. doi:10.1038/s41598-018-33819-8
}

\seealso{
\code{\link{kcpRS}}

\code{\link{kcpRS_workflow}}

\code{\link{MentalLoad}}

\code{\link{CO2Inhalation}}

}
