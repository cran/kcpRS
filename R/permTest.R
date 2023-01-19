
#' KCP Permutation Test
#'
#' The KCP permutation test implements the variance test and the variance drop test to determine if there is at least one change point in the running statistics
#'
#' @param data data \emph{N} x \emph{v} dataframe where \emph{N} is the number of time points and \emph{v} the number of variables
#' @param RS_fun Running statistics function: Should require the time series and \code{wsize} as input and return a dataframe of running statistics
#' as output. This output dataframe should have rows that correspond to the time windows and columns that correspond to the variable(s) on which the running statistics were computed.
#' @param wsize Window size
#' @param nperm Number of permutations to be used in the permutation test
#' @param Kmax Maximum number of change points desired
#' @param alpha Significance level of the permutation test
#' @param varTest If FALSE, only the variance DROP test is implemented, and if TRUE, both the variance and the variance DROP tests are implemented.

#' @return \item{sig}{Significance of having at least one change point. 0 - Not significant, 1- Significant}
#' \item{p_var_test}{P-value of the variance test.}
#' \item{p_varDrop_test}{P-value of the variance drop test.}
#' \item{perm_rmin}{A matrix of minimized variance criterion for the permuted data.}
#' \item{perm_rmin_without_NA}{A matrix of minimized variance criterion for the permuted data without NA values.}
#' @importFrom stats var
#' @importFrom stats na.omit
#' @importFrom foreach foreach %dopar%
#' @import roll
#' @import Rcpp
#' @references Cabrieto, J., Tuerlinckx, F., Kuppens, P., Hunyadi, B., & Ceulemans, E. (2018). Testing for the presence of correlation changes
#' in a multivariate time series: A permutation based approach. \emph{Scientific Reports}, 8, 769, 1-20. doi:10.1038/s41598-017-19067-2


permTest <-
  function(data,
           RS_fun,
           wsize = 25,
           nperm = 1000,
           Kmax = 10,
           alpha = .05,
          varTest = FALSE ) {
    myfun <- function(data, wsize = 25,FUN) {FUN(data, wsize)} #function that allows the user to input RS_fun

    #variance criterion: original data
    RS <- myfun(data, wsize, FUN = RS_fun)
    output <- kcpa(RunStat = RS,Kmax,wsize)
    kcp_RS_result <- output$kcpSoln
    rmin <- kcp_RS_result$Rmin
    alpha_per_test <- ifelse(isTRUE(varTest), alpha/2 , alpha)
    
    #variance criterion:permuted data
    N <- nrow(data)
    v <- ncol(data)
    #perm_rmin<-NULL
    perm_rmin <- matrix(0, nrow = nperm, ncol = Kmax + 1)
    p=NULL


    perm_rmin <- foreach(p = 1:nperm,.packages = c('roll', 'Rcpp', 'kcpRS'),.combine = rbind) %dopar% {
        set.seed(444+p)
        Xperm <- as.data.frame(matrix(0, nrow = N, ncol = v))
        RS_temp <- NULL
        kcp_RS_result_temp <- NULL

        indeces <- sample(N, N, replace = FALSE)#permuted index
        Xperm <- data[indeces, ]#permuted data
        RS_temp <- myfun(data = Xperm,wsize,FUN = RS_fun)#running stat of permuted data
        if (anyNA(RS_temp) || !isFinite(RS_temp)) { 
          rep(NA, Kmax+1)
        } else {
          output <- kcpa(RunStat = RS_temp,Kmax,wsize)  #kcp solution for the RS of permuted data
          if (output$medianK != 0) {
            kcp_RS_result_temp <- output$kcpSoln
            Rmin_temp_data <- kcp_RS_result_temp$Rmin #variance criterion of permuted data
            as.vector(Rmin_temp_data)
          } else {
            rep(NA, Kmax+1)
          }
        }
      }

    perm_rmin <- as.data.frame(perm_rmin, row.names = FALSE)
    colnames(perm_rmin) <- paste0("K=", 0:Kmax)

    
    perm_rmin_without_NA = na.omit(perm_rmin)
    
    if (nrow(perm_rmin_without_NA) != 0) {
      if (isTRUE(varTest)) {
        #variance test
        var <- rmin[1] #variance of the original RS, Rmin at K<-0
        var_perm <- perm_rmin_without_NA[, 1] #variances of the RS from permuted data: distribution of the variances
        p_var_test <- length(var_perm[var_perm > var]) / dim(perm_rmin_without_NA)[1]
      }
  
      #variance drop test
      var_drop <- abs(apply(perm_rmin_without_NA, 1, diff)) #computes slopes (drop in Rmin) for each permuted data
      max_vdrop_perm <- apply(var_drop, 2, max) #computes max slope for each permuted data: distribution of max var drop
      max_vdrop <- max(abs(diff(rmin))) #max slope of original data
      p_varDrop_test <-
        length(max_vdrop_perm[max_vdrop_perm > max_vdrop]) / dim(perm_rmin_without_NA)[1] #p-value, variance drop test
    } else {
      p_var_test = 0
      p_varDrop_test = 0
    }
    #significance
    if (isTRUE(varTest)) {
      sig = ifelse(p_var_test < alpha_per_test |
                     p_varDrop_test < alpha_per_test,1,0)

       output <- list(
        sig = sig,
        p_var_test = p_var_test,
        p_varDrop_test = p_varDrop_test,
        perm_rmin = perm_rmin,
        perm_rmin_without_NA = perm_rmin_without_NA
       )
    }

    if (isFALSE(varTest)) {
      sig = ifelse(p_varDrop_test < alpha_per_test, 1, 0)

      output <- list(sig = sig,
                     p_varDrop_test = p_varDrop_test,
                     perm_rmin = perm_rmin,
                     perm_rmin_without_NA = perm_rmin_without_NA
      )
    }
    return(output)

  }