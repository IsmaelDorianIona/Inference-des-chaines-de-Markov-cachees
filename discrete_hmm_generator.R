discrete_hmm_generator <- function(mu, p, a, t) {
  
  n <- length(mu) # nombre d'etats caches possibles
  k <- ncol(a) # nombre d'etats observes possibles
  
  x <- rep(0, t) #  sequence cachee
  y <- rep(0, t) # sequence observee
  
  # initialisation du premier etat cache
  x[1] <- sample(0:(n-1), size = 1, prob = mu)
  # initialisation du premier etat observe
  y[1] <- sample(0:(k-1), size = 1, prob = a[x[1] + 1,])
  
  # iteration sur le reste de la sequence
  for(i in 2:t) {
    # tirage de l'etat cache courant selon l'etat precedent
    x[i] <- sample(0:(n-1), size = 1, prob = p[x[i-1] + 1,])
    # tirage de l'etat observe courant selon l'etat cache courant
    y[i] <- sample(0:(k-1), size = 1, prob = a[x[i] + 1,])
  }
  
  return(list(hidden = x, observed = y))
}