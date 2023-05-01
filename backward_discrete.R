backward_discrete <- function(y, mu, p, a) {
  n <- nrow(p)
  t <- length(y)
  beta <- matrix(0, nrow = t, ncol = n)
  beta[t, ] <- rep(1, n)
  
  for (k in rev(seq_len(t-1))) {
    beta[k, ] <- t(p %*% (beta[k + 1,] * a[, y[k + 1] + 1])) }
  return(list(proba = sum(beta[1, ] * mu * a[, y[1] + 1]), beta = beta))
}
