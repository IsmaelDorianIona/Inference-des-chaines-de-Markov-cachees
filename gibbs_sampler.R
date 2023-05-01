source(file = "backward_continuous.R")

count_transitions <- function(vec, i, j){
  temp <- diff(vec)
  val <- sum(temp == j - i & vec[2:length(vec)] == j)
  return(val)
}

# Fonction qui simule les paramètres m_i selon (3)
simulate_m <- function(x, y, sigma2, n, kappa, xi) {
  m <- numeric(n)
  for(i in 0:(n-1)){
    s_i <- sum(y[which(x == i)])
    n_i <- sum(x == i)
    m[i+1] <- rnorm(1,  (s_i+kappa*xi*sigma2)/(n_i+kappa*sigma2), sqrt(sigma2/(n_i+kappa*sigma2)))
  }
  return(m)
}


# Fonction qui simule sigma^2 selon (4)
simulate_sigma2 <- function(x, y, m, alpha, beta) {
  
  t <- length(y)
  a <- alpha + t/2
  b <- beta + sum((y - m[x+1])**2)/2
  sigma2 <- 1/rgamma(1, a, b)
  return(sigma2)
  
}


# Fonction qui simule les paramètres p_{i,*} selon (2)
simulate_p <- function(x, n) {
  
  p <- matrix(0, nrow=n, ncol=n) 
  for (k in 1:n) { 
    n_k <- sapply(0:(n-1), count_transitions, vec = x, i = k-1) 
    p[k,] <- rdirichlet(1, n_k + rep(1, n)) 
  } 
  return(p)
  
}


# Fonction qui simule les paramètres mu_i selon (1)
simulate_mu <- function(x, n) {
  
  alpha_mu <- rep(1, n) + as.integer(x[1] == 0:(n-1))
  mu <- rdirichlet(1, alpha_mu)
  return(mu)
  
}


# Fonction qui simule la chaîne de Markov cachée X selon (6) et (7)
simulate_X <- function(y, m, p, mu, n, sigma, b) {
  
  t <- length(y)
  x <- rep(0, t)
  
  #Simuler X_1
  x[1] <- sample(0:(n-1), 1, prob = mu * dnorm(y[1], m, sigma) * b[1,])
  
  #Simuler X_2:T
  for (k in 2:t) {
    x[k] <- sample(0:(n-1), size = 1, prob = b[k-1,] * dnorm(y[k], m, sigma) * p[x[k-1]+1,])
  }
  return(x)
}


gibbs_sampler <- function(y, n, alpha = 2, g = .2, iter = 10000, burnin = 1000) {
  #Initialisation des parameters
  t <- length(y)
  r <- max(y) - min(y)
  h <- 10 / r**2
  kappa <- 1 / r**2
  xi <- (min(y)+max(y)) / 2
  m <- min(y) + (1:n - 1)*r/n + r/(2*n)
  x <- rep(0,t)
  for(k in 1:t){
    x[k] <- which.min((y[k] - m)**2) - 1
  }
  sigma2 <- mean((y - m[x+1])**2)
  p <- matrix(0, nrow=n, ncol=n)
  mu <- rep(1/n, n)
  beta <- r**2 * g / 10
  # Initialisation de la matrice de transition
  for (k in 1:n) { 
    n_k <- sapply(0:(n-1), count_transitions, vec = x, i = k-1) 
    p[k,] <- n_k / sum(n_k)
  } 
  
  #Creation d'une liste pour stocker les valeurs des parametres
  m_list <- vector(mode = "list", length = iter)
  mu_list <- vector(mode = "list", length = iter)
  p_list <- vector(mode = "list", length = iter)
  sigma2_list <- vector(mode = "numeric", length = iter)
  x_list <- vector(mode = "list", length = iter)
  
  #Gibbs step
  for (k in 1:(iter + burnin)) {
    #Recuperation des probabilites backward
    b <- backward_continuous(y, mu, p, m, sqrt(sigma2))$beta
    # Step 1: mise a jour de m
    m <- simulate_m(x, y, sigma2, n, kappa, xi)
    m_list[[k]] <- m
    # Step 2: mise a jour de sigma2
    sigma2 <- simulate_sigma2(x, y, m, alpha, beta)
    sigma2_list[k] <- sigma2
    # Step 3: mise a jour de  beta
    beta <- rgamma(1, g+alpha, h+1/sigma2)
    # Step 4: mise a jour de  P:
    p <- simulate_p(x, n)
    p_list[[k]] <- p
    # Step 6: mise a jour de mu:
    mu <- simulate_mu(x, n)
    mu_list[[k]] <- mu
    # Step 7: update X:
    x <- simulate_X(y, m, p, mu, n, sqrt(sigma2), b)
    x_list[[k]] <- x
    
  }
  
  #Renvoie des valeurs pour m, mu, et p apres la periode de burn-in
  return(list(m = m_list[(burnin+1):(burnin+iter)],
              mu = mu_list[(burnin+1):(burnin+iter)],
              p = p_list[(burnin+1):(burnin+iter)],
              sigma2 = sigma2_list[(burnin+1):(burnin+iter)],
              x = x_list[(burnin+1):(burnin+iter)]))
}
