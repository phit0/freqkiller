sigma_gibbs <- function(y, eta, a_t, b_t) {
  a <- a_func(y, a_t)
  b <- b_func(y, eta, b_t)
  sigma_new <- 1/rgamma(1, shape = a, rate = b)
  return(sigma_new)
}

a_func <- function(y, a_t) {
  n <- length(y)
  a_new <- (n / 2) + a_t
  return(a_new)
}

b_func <- function(y, eta, b_t) {
  b_new <- b_t + t(y - eta) %*% (y - eta)
  return(b_new)
}
