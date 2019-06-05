#response function h
h <- function(eta) {
  return(exp(eta) / (1 + exp(eta)))
}

#first derivative of h
dh <- function(eta) {
  return(exp(eta) / (1 + exp(eta))^2)
}

#parameter pi (execution of h)
pi = function(beta_t) {
  output <- h(X%*%beta_t)
  for (i in 1:length(output)) {
    if (is.nan(output[i])) {
      output[i] <- 1.0
    }
  }
  return(output)
}

#loglikelihood L(y|beta)
loglike_beta = function(beta_t) {
  result = sum(log(pi(beta_t)^y * (1- pi(beta_t))^(1-y)))
  return(result)
}

#posterior for priors m = 0 and M = I only!!!
post_beta = function(beta) {
  return(exp(-0.5*t(beta)%*%beta + loglike_beta(beta)))
}
