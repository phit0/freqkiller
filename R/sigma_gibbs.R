sigma_gibbs <- function(a_t, b_t) {
  out <- 1/rgamma(1, shape = a_t, rate = b_t)
  return(out)
}

a_func <- function(y, a_t) {
  n <- length(y)
  a_new <- (n / 2) + a_t
  return(a_new)
}

b_func <- function(y, eta, b_t) {
  b_new <- b_t + 0.5 * t(y - eta) %*% (y - eta)
  return(b_new)
}
