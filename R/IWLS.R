w_func <- function(sigma2_t,beta_t) {
  out <- switch(dist,
                "normal" = diag(length(y))/ sigma2_t,
                "poisson" = diag(exp(X%*%beta_t)))
  return(out)
}

fisher_func <- function(sigma2_t,beta_t) {
  out <- switch(dist,
                "normal" = (1/sigma2_t) * t(X) %*% X + M_1,
                "poisson" = t(X)%*%w_func(sigma2_t,beta_t))
  return(out)
}

y_wgl_func <- function(beta_t) {
  out <- switch(dist,
                "normal" = y,
                "poisson" = X%*%beta_t + ((y - exp(X%*%beta_t)) / exp(X%*%beta_t)))
  return(out)
}

mu_func <- function(sigma2_t,beta_t) {
  out <- switch(dist,
                "normal" = solve(fisher_func(sigma2_t,beta_t)) %*% (t(X) %*% w_func(sigma2_t,beta_t) %*% y_wgl_func(beta_t) + M_1 %*% m),
                "poisson" = solve(fisher_func(sigma2_t,beta_t))%*% (t(X) %*% w_func(sigma2_t,beta_t)%*% y_wgl_func(beta_t) + M_1 %*% m))
  return(out)
}

#beta prior
prior_func <- function(beta_t) {
  #out <- dmvnorm(beta_t, mean = m, sigma = M, log = T)
  n <- length(beta_t)
  out <- -0.5 * n * log(2*pi) -0.5 * log(M_det) - 0.5 * t(beta_t - m) %*% M_1 %*% (beta_t - m)
  return(out)
}


proposalfunction <- function(mu, sigma2) {
  out <- rmvnorm(1, mu, sigma2)
  return(as.vector(out))
}

cond_proposaldensity <- function(beta, mu, Fisher) {
  #out <- dmvnorm(beta, mu, solve(Fisher), log = T)
  n <- length(beta)
  out <- -0.5 * n * log(2*pi) - 0.5 * log(1 / det(Fisher)) - 0.5 * t(beta - mu)%*%Fisher%*%(beta - mu)
  return(out)
}


#likelihood for beta (without assuming independence)
loglik_func <- function(beta_t, sigma2_t) {
  out <- switch(dist,
                "normal" = sum(dnorm(y, X%*%beta_t, sigma2_t, log = T)),
                "poisson" = sum(dpois(y, lambda = exp(eta_t, dist), log = T)))

  return(out)
}


