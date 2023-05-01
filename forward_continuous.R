forward_continuous <- function(y, mu, p, m, sigma) {
  n <- nrow(p)
  t <- length(y)
  alpha <- matrix(0, nrow = t, ncol = n)
  c <- rep(0, t) #constante de normalisation
  
  alpha[1,] <- dnorm(y[1], m, sigma) * mu
  c[1] <- sum(alpha[1,])
  alpha[1,] <- alpha[1,] / c[1]
  
  for (k in seq_len(t)[-1]) {
    alpha[k,] <- dnorm(y[k], m, sigma) * (alpha[k - 1, ] %*% p)
    c[k] <- sum(alpha[k,])
    alpha[k,] <- alpha[k,] / c[k]
  }
  return(list(proba = sum(log(c)), alpha = alpha, c_value = c))
}

#Le resultat renvoye pour la probabilite est en echelle logarithmique pour
#eviter les arrondis de R.
#Ainsi, la proabilite d'obtenir la sequence d'observations
#est egale à l'exponentielle du resultat renvoye par 
#forward_continuous()$proba