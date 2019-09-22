# Response function for Bernoulli dependent variables
h <- function(x){
  return(1 / (1 + exp(-x)))
}

# First derivative of the response function
dh <- function(x){
  out = exp(x) / (1 + exp(x))^2
  # using 1-time L`hospital
  result <- ifelse(is.nan(out), 0, out)
  return(result)
}

# Returns the diagonal vector of the Weights matrix W
w_func <- function(sigma2_t, beta_t, y, X, dist) {
  n = length(y)
  out <- switch(
    dist,
    "normal" = diag(n) / sigma2_t,
    "poisson" = as.vector(exp(X %*% beta_t)),
    "bernoulli" = as.vector(dh(X %*% beta_t))
    )
  return(out)
}

# Returns the expected Fisher Information
fisher_func <- function(sigma2_t, beta_t, y, X, M_1, dist) {
  out <- switch(
    dist,
    "normal" = (1/sigma2_t) * crossprod(X) + M_1,
    "poisson" = t(X * w_func(sigma2_t, beta_t, y, X, dist)) %*% X + M_1,
    "bernoulli" = t(X * w_func(sigma2_t, beta_t, y, X, dist)) %*% X + M_1)
  return(out)
}

# Returns the working observations
y_wgl_func <- function(beta_t, y, X, dist) {
  out <- switch(
    dist,
    "normal" = y,
    "poisson" = X %*% beta_t + (y / exp(X %*% beta_t) - 1),
    "bernoulli" = X %*% beta_t + (y - h(X %*% beta_t)) / dh(X %*% beta_t)
    )
  return(out)
}

# Returns the expectation of the proposal density
mu_func <- function(sigma2_t, beta_t, y, X, M_1, m, dist) {
  out <- switch(
    dist,
    "normal" = solve(fisher_func(sigma2_t, beta_t, y, X, M_1, dist)) %*%
      (t(X) %*% y_wgl_func(beta_t, y, X, dist))/sigma2_t + M_1 %*% m,

    "poisson" = chol2inv(chol(fisher_func(sigma2_t, beta_t, y, X, M_1, dist))) %*%
      (t(X * w_func(sigma2_t, beta_t, y, X, dist)) %*%
         y_wgl_func(beta_t, y, X, dist) + M_1 %*% m),

    "bernoulli" = chol2inv(chol(fisher_func(sigma2_t, beta_t, y, X, M_1, dist))) %*%
      (t(X * w_func(sigma2_t, beta_t, y, X, dist)) %*%
         y_wgl_func(beta_t, y, X, dist) + M_1 %*% m)
    )
  return(out)
}

# Prior density for the parameter beta
prior_func <- function(beta_t, m, M_1, M_det) {
  n <- length(beta_t)
  out <- -0.5 * n * log(2 * pi) - 0.5 * log(M_det) -
    0.5 * t(beta_t - m) %*% M_1 %*% (beta_t - m)

  return(out)
}

# Returns a random vector from the proposal density
proposalfunction <- function(mu, sigma2) {
  out <- mvrnorm(1, mu, sigma2)
  return(as.vector(out))
}

# Proposal density
cond_proposaldensity <- function(beta, mu, Fisher) {
  n <- length(beta)
  out <- -0.5 * n * log(2*pi) - 0.5 * log(1 / det(Fisher)) -
    0.5 * t(beta - mu) %*% Fisher %*% (beta - mu)

  return(out)
}


#likelihood function for beta (assuming independence)
loglik_func <- function(beta_t, sigma2_t, y, X, dist) {
  out <- switch(
    dist,
    "normal" = sum(dnorm(y, X%*%beta_t, sqrt(sigma2_t), log = T)),
    "poisson" = sum(dpois(y, lambda = exp(X %*% beta_t), log = T)),
    "bernoulli" = sum(dbinom(y, size = 1, prob  = h(X %*% beta_t), log = T))
    )
  return(out)
}

# Sets initial values for beta if desired by the user
beta_init <- function(mf, dist){
  out <- switch(
    dist,
    "normal" = lm(mf)$coefficients,
    "poisson" = glm(mf, family = poisson(link = log))$coefficients,
    "bernoulli" = glm(mf, family = binomial(link = logit))$coefficients
    )
}
