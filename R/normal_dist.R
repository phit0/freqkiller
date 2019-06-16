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

#sigma function
sigma_vec <- function(y) {
  out <- rep(1, length(y))
  return(out)
}

#beta prior
prior <- function(beta_t){
  out <- mvtnorm::dmvnorm(beta_t, mean = m(beta_t), sigma = M(beta_t), log = T)
  return(out)
}

#likelihood for beta
like_beta <- function(beta_t, X, y) {
  mean = X%*%beta_t
  sigma = diag(length(y))
  out <- dmvnorm(y, mean, sigma, log=T)
  return(out)
}

#posterior for beta
post_beta <- function(beta_t, X, y){
  return(like_beta(beta_t, X, y) + prior(beta_t))
}
