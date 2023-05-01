viterbi_discrete <- function(y, mu, p, a) {
  
  # y: vecteur des observations
  # states: vecteur des etats caches possibles
  # mu: vecteur de loi initiale
  # p: matrice de transition
  # a: matrice d'emission
  
  t <- length(y)
  n <- ncol(p)
  states <- 0:(n-1)
  
  # Initialisation
  delta <- matrix(0, nrow = t, ncol = n)
  psi <- matrix(NA, nrow = t, ncol = n)  # initialisation avec NA
  delta[1,] <- mu * a[,y[1] + 1]
  
  # Recurrence
  for (k in 2:t) {
    for (j in 1:n) {
      tmp <- delta[k-1,] * p[,j] * a[j, y[k] + 1]
      delta[k, j] <- max(tmp)
      psi[k, j] <- which.max(tmp)
    }
  }
  
  # Fin
  p_star <- max(delta[t,])
  s_star <- which.max(delta[t,])
  
  # Back-Tracking du chemin optimal
  path <- integer(t)
  path[t] <- s_star
  for (k in (t-1):1) {
    path[k] <- psi[k+1,path[k+1]]
  }
  
  # Renvoi des resultats
  list(prob = p_star, path = states[path])
}
