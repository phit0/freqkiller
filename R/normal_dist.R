#response function
h <- function(eta) {
  out <- eta
  return(out)
}

#first derivative of h
dh <- function(eta) {
  out <- rep(1, dim(eta)[1])
  return(out)
}

#beta prior
prior_func <- function(beta_t, m, M) {
  out <- mvtnorm::dmvnorm(beta_t,
                          mean = m,
                          sigma = M,
                          log = T)
  return(out)
}

#likelihood for beta (without assuming independence)
loglik_func <- function(eta_t, sigma_t, y) {
  out <- mvtnorm::dmvnorm(y, eta_t, diag(y)*sigma_t, log = T)
  return(out)
}

