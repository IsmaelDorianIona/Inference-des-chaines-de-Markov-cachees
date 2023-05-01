source(file = "forward_discrete.R")
source(file = "backward_discrete.R")

em_iteration <- function(y, mu, p, a) {
  n <- ncol(p)
  t <- length(y)

  #gamma  
  frwd <- forward_discrete(y, mu, p, a)
  bcwd <- backward_discrete(y, mu, p, a)
  
  gamma<- frwd[[2]] * bcwd[[2]]/frwd[[1]]
  
  
  #xi
  ksi = array(0, dim = c(n, n, t))
  for (i in 1:n) {
    for (u in 1:(t-1)) {
      ksi[i, , u] = frwd[[2]][u, i] * bcwd[[2]][u + 1, ] * p[i, ] * a[, y[u + 1] + 1]/frwd[[1]]
    }
  }
  
  #mu
  mu_star <- gamma[1, ]
  
  #P
  p_star <- matrix(0,nrow = n, ncol = n)
  for (i in 1 : n) {
    p_star[i, ] <- rowSums(ksi[i, , ])/sum(gamma[1 : t - 1, i])
  }
  
  #A
  m <- length(a[1,])
  a_star <- matrix(0, nrow = n, ncol = m)
  for (j in 1 : m) {
    a_star[, j] <- colSums(gamma * (y + 1 == j))/colSums(gamma)
  }
  
  return(list(mu_star, p_star, a_star))
}



em_discrete <- function(y, mu, p, a, n_iter = 100, tol = 1e-6) {
  for (iter in 1 : n_iter) {
    upd <- em_iteration(y, mu, a, p)
    if (norm(as.matrix(upd[[1]]) - as.matrix(mu)) < tol & norm(upd[[2]] - p) < tol & norm(upd[[3]] - a) < tol) {
      return(upd)
    }
    else {
      mu <- upd[[1]]
      p <- upd[[2]]
      a <- upd[[3]]
    }
  }
  return(upd)
}