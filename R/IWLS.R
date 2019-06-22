
# write a test for w_func that checks it has appropriate dimensions
#
w_func <- function(eta, sigma_t) {
  w <- (dh(eta)) ^ 2 / sigma_t
  return(diag(c(w)))
}

fisher_func <- function(X, W_t, M) {
  out <- t(X) %*% W_t %*% X + solve(M)
  return(out)
}

y_wgl_func <- function(eta_t, y) {
  out <- eta_t + ((y - h(eta_t)) / dh(eta_t))
  return(out)
}

mu_func <- function(X, Ft, Wt, yt_wgl, M, m) {
  out <- solve(Ft) %*% t(X) %*% Wt %*% yt_wgl + solve(M) %*% m
  return(out)
}

proposalfunction <- function(mu, sigma) {
  out <- rmvnorm(1, mu, sigma)
  return(as.vector(out))
}


cond_proposaldensity <- function(beta, mu, sigma) {
  out <- dmvnorm(beta, mu, sigma, log = T)
  return(out)
}
