#response function
h <- function(eta) {
  out <- eta
  return(out)
}

#first derivative of h
dh <- function(eta) {
  out <- 1
  return(out)
}

#sigma function
sigma_vec <- function(beta_t) {
  out <- rep(1, length(y))
  return(out)
}

#beta prior
prior <- function(beta_t){
  out <- dmvnorm(beta_t, mean = m(beta_t), sigma = M(beta_t), log = T)
  return(out)
}

#likelihood for beta
like_beta <- function(beta_t) {
  out <- dmvnorm(y, mean = X%*%beta_t , sigma = diag(length(y)),log=T)
  return(out)
}

#posterior for beta
post_beta <- function(beta_t){
  return(like_beta(beta_t) + prior(beta_t))
}
