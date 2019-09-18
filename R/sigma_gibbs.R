#' Title
#'
#' @noRd
#'
sigma_gibbs <- function(a_t, b_t) {
  out <- 1 / rgamma(1, shape = a_t, rate = b_t)
  return(out)
}

b_func <- function(beta_t, y, X, b0) {
  eta <- X %*% beta_t
  b_new <- b0 + 0.5 * t(y - eta) %*% (y - eta)
  return(b_new)
}
