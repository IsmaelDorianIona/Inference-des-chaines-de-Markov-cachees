source(file = "forward_continuous.R")
source(file = "backward_continuous.R")

count_transitions <- function(vec, i, j){
  temp <- diff(vec)
  val <- sum(temp == j - i & vec[2:length(vec)] == j)
  return(val)
}

em_continuous <- function(y, mu, p, m, sigma2, n_iter=100, tol=1e-6) {
  
  t <- length(y)
  n <- nrow(p)
  mu_list <- list(mu)
  m_list <- list(m)
  p_list <- list(p)
  sigma2_list <- list(sigma2)
  
  # Iteration de l algorithme
  for (iter in 1:n_iter) {
    
    # E-step
    
    gamma <- matrix(0, nrow = n, ncol = t)
    xi <- array(0, dim = c(n, n, t-1))
    
    # Forward-backward
    alpha <- forward_continuous(y, mu, p, m, sqrt(sigma2))$alpha
    beta <- backward_continuous(y, mu, p, m, sqrt(sigma2))$beta
    
    
    # Calcul des quantites utiles pour les updates
    for (k in 1:t) {
      
      gamma[,k] <- alpha[k,] * beta[k,] / sum(alpha[k,] * beta[k,])
      
    }
    
    for (k in 1:(t-1)) {
      for (i in 1:n) {
        for (j in 1:n) {
          xi[i, j, k] <- alpha[k, i] * p[i, j] * dnorm(y[k+1], m[j], sqrt(sigma2)) * beta[k + 1, j]
        }
      }
      xi[, , k] <- xi[, , k] / sum(sum(xi[, , k]))
    }
    
    # M-step: updates
    mu <- gamma[,1]
    for (i in 1:n) {
      for (j in 1:n) {
        p[i,j] <- sum(xi[i,j,])
      }
      p[i,] <- p[i,] / sum(p[i,])
      m[i] <- sum(gamma[i,] * y) / sum(gamma[i,])
    }
    sigma2 <- 0
    for(i in 1:n){
      sigma2 <- sigma2 + sum(gamma[i,]*(y - m[i])**2)
    }
    sigma2 <- sigma2 / t
    
    # Stockage des valeurs itermediaires
    mu_list[[iter+1]] <- mu
    m_list[[iter+1]] <- m
    p_list[[iter+1]] <- p
    sigma2_list[[iter+1]] <- sigma2
    
    # verification de la convergence
    ll <- sum(log(apply(alpha, 1, sum)))
    ll_prev <- ll
  }
  
  return(list(mu = mu_list, m = m_list, p = p_list, sigma2 = sigma2_list))
  
}