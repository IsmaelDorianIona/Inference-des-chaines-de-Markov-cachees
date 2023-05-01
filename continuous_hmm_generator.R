continuous_hmm_generator <- function(mu, p, m, sigma, t){
  
  n <- length(mu) # nombre d'etats caches possibles
  x <- rep(0, t) #  sequence cachee
  y <- rep(0, t) # sequence observee
  
  # initialisation du premier etat cache
  x[1] <- sample(0:(n-1), size = 1, prob = mu)
  # initialisation du premier etat observe
  y[1] <- rnorm(1, m[x[1] + 1], sigma)
  
  # iteration sur le reste de la sequence
  for(i in 2:t) {
    # tirage de l'etat cache courant selon l'etat precedent
    x[i] <- sample(0:(n-1), size = 1, prob = p[x[i-1] + 1,])
    # tirage de l'etat observe courant selon l'etat cache courant
    y[i] <- rnorm(1, m[x[i] + 1], sigma)
  }
  
  return(data.frame(hidden = x, observed = y))
}
