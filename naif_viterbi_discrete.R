proba <- function(y, x, mu, p, a) {
  
  result <- mu[x[1] + 1]*a[x[1] + 1, y[1] + 1]
  for (i in 2:length(y)) {
    result <- result * p[x[i - 1] + 1, x[i] + 1] * a[x[i] + 1, y[i] + 1]
  }
  return(result)
}


naif_viterbi_discrete <- function(y, mu, p, a, states) {
  
  t <- length(y)
  chainpossible <- as.matrix(expand.grid(rep(list(states), t)))
  result <- rep(0, nrow(chainpossible))
  
  for (i in 1:nrow(chainpossible)) {
    result[i] <- proba(y, chainpossible[i, ], mu, a, p)
  }
  probmax <- max(result)
  
  return(list(prob = probmax, path = chainpossible[which(result == probmax), ]))
}
