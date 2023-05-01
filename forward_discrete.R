forward_discrete <- function(y, mu, p, a) {
  n <- nrow(p)
  t <- length(y)
  alpha <- matrix(0, nrow = t, ncol = n)
  alpha[1,] <- a[, y[1] + 1] * mu
  
  for (k in seq_len(t)[-1]) {
    alpha[k, ] <- a[, y[k] + 1] * (alpha[k - 1, ] %*% p)  }
  return(list(proba = sum(alpha[t, ]), alpha = alpha))
}
