source(file = "forward_continuous.R")

backward_continuous <- function(y, mu, p, m, sigma) {
  n <- length(mu) # nombre d'etats caches
  t <- length(y) #longueur de la sequence
  
  #Recuperation des constantes de normalisation
  c <- forward_continuous(y, mu, p, m, sigma)$c_value
  
  beta <- matrix(0, nrow = t, ncol = n)
  beta[t, ] <- rep(1, n) # Initialisation
  for (k in rev(seq_len(t-1))) {
    beta[k, ] <- t(p %*% (beta[k + 1,] * dnorm(y[k + 1], m, sigma)))
    beta[k, ] <- beta[k,] / c[k + 1]
  }
  # return(beta)
  return(list(proba = sum(log(c)), beta = beta))
}

#Le resultat renvoye pour la probabilite est en echelle logarithmique pour
#eviter les arrondis de R.
#Ainsi, la proabilite d'obtenir la sequence d'observations
#est egale à l'exponentielle du resultat renvoye par 
#backward_continuous()$proba