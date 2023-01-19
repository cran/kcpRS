#' @rdname kcpRS
#' @param kcp_details If \emph{TRUE}, then the matrix of optimal change points solutions given \emph{k} is displayed.
#' If FALSE, then this output is suppressed.
#' @export

print.kcpRS <- function(x, kcp_details = TRUE, ...) {

  if (x$medianK == 0) {
    cat("\n")
    cat("    The KCP-RS analysis cannot be performed, as the median Euclidean distance between all pairs of running statistics is 0.")
    cat("\n")
  } else {
    if (x$nperm > 0) {
      cat("\n")
      cat("    Number of change points detected based on grid search:", x$BestK, "\n")
      if (x$BestK > 0) {
        cat("    Change point location(s):", x$changePoints, "\n")
      }
      cat("\n")
      cat("    Number of change points detected based on scree test:", x$changePoints_scree_test, "\n")
      cat("\n")
      
      if (isTRUE(x$varTest)) {
        if (x$BestK > 0) {
          cat(
            "    KCP permutation test is significant: At least one of the two subtests is significant.",
            "\n"
          )
        }
        if (x$BestK == 0) {
          cat(
            "    The KCP permutation test is NOT significant: Neither of the two subtests is significant.",
            "\n"
          )
        }
        cat("    Overall significance level:", x$alpha , "\n")
        cat("    Significance level of each subtest:", x$subTest_alpha , "\n")
        cat("         P-value of the variance test:", x$p_var_test, "\n")
        cat("         P-value of the variance drop test:", x$p_varDrop_test, "\n")
        cat("\n")
      }
  
      if (isFALSE(x$varTest)) {
        if (x$BestK > 0) {
          cat("    KCP permutation test is significant:", "\n")
        }
        if (x$BestK == 0) {
          cat("    The KCP permutation test is NOT significant:", "\n")
        }
        cat("        Significance level:", x$subTest_alpha , "\n")
        cat("        P-value of the variance drop test:", x$p_varDrop_test, "\n")
        cat("\n")
      }
  
      if (x$nperm_success < x$nperm) {
        cat("    Warning: Only", x$nperm_success, "out of", x$nperm, "permutation tests succeeded,", x$nperm - x$nperm_success, "failed. Possible reasons are:\n")
        cat("                  - NA values in the running statistics of some of the permutations\n");
        cat("                  - Inf values in the running statistics of some of the permutations\n");
        cat("                  - The median Euclidean distance between all pairs of running statistics of some of the permutations is 0\n\n");
      }
  
    }
  
    else {
      cat("    KCP permutation test was not implemented.", "\n")
      cat("\n")
    }
  
  
    if (kcp_details == TRUE) {
      cat("    Optimal change points given k:", "\n")
      space = rep("", nrow(x$CPs_given_K))              #space column
      y = cbind(space, x$CPs_given_K)                    #add the column to the dataframe to be indented
      colnames(y) = c("    ", colnames(x$CPs_given_K))   #the width of the col name is the width of the indention
      print(y, row.names = F)
    }
  }
}
