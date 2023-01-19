isFinite <- function(data) {
  
  d = as.matrix(data)
  l = length(d)
  
  sum(is.finite(d)) == l
}