naif_discrete <- function(y, mu, p, a, states) {
  t <- length(y)
  arrangements <- expand.grid(rep(list(states), t))
  proba <- numeric(nrow(arrangements))
  if (t == 1) {
    return(sum(mu * a[, y[1] + 1]))
  }
  for (i in 1 : nrow(arrangements)) {
    proba[i] <- mu[arrangements[i, 1] + 1]
    for (j in 1 : (t - 1)) {
      proba[i] <- proba[i] * p[arrangements[i, j] + 1, arrangements[i, j + 1] + 1] * a[arrangements[i, j] + 1, y[j] + 1]
    }
  }
  return(sum(proba * a[arrangements[, t] + 1, y[t] + 1]))
}
