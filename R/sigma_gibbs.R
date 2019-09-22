# samples a value for sigma2 from an inverse gamma distribution
sigma_gibbs <- function(a_t, b_t) {
  out <- 1 / rgamma(1, shape = a_t, rate = b_t)
  return(out)
}

# updates tha parameter b from the inverse gamma distribution.
# a is a fixed value specified in the code "metroNorm"
b_func <- function(beta_t, y, X, b0) {
  eta <- X %*% beta_t
  b_new <- b0 + 0.5 * t(y - eta) %*% (y - eta)
  return(b_new)
}
