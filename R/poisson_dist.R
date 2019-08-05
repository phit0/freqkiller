
#response function
h2 <- function(eta) {
  out <- exp(eta)
  return(out)
}

#first derivative of h
dh2 <- function(eta) {
  out <- exp(e)
  return(out)
}

#likelihood for beta (without assuming independence)
loglik_func <- function(eta_t, sigma_t, y) {
  out <- dpois(y, lambda = h(eta_t), log = T)
  return(sum(out))
}

