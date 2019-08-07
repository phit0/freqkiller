
# write a test for w_func that checks it has appropriate dimensions

w_func <- function(sigma_t,beta_t) {
  out <- switch(dist,
                "normal" = diag(length(y))/ sigma_t,
                "poisson" = "CODE FEHLT HIER!")
  return(out)
}

fisher_func <- function(sigma_t, beta_t) {
  out <- switch(dist,
                "normal" = (1/sigma_t)*t(X) %*% X + M_1,
                "poisson" = "CODE FEHLT HIER Hola papi <3")
  return(out)
}

y_wgl_func <- function(sigma_t, beta_t) {
  out <- switch(dist,
         "normal" = y,
         "poisson" = eta_t + ((y - exp(eta_t, dist)) / exp(eta_t, dist)))
  return(out)
}

mu_func <- function(sigma_t, beta_t) {
  out <- switch(dist,
                "normal" = solve(fisher_func(sigma_t, beta_t)) %*% t(X) %*% w_func(sigma_t,beta_t) %*% y_wgl_func(sigma_t,beta_t) + M_1 %*% m,
                "poisson" = "CODE FEHLT HIER")
  return(out)
}

#beta prior
prior_func <- function(beta_t) {
  #out <- dmvnorm(beta_t, mean = m, sigma = M, log = T)
  n <- length(beta_t)
  out <- -0.5 * n * log(2*pi) -0.5 * log(M_det) - 0.5 * t(beta_t - m) %*% M_1 %*% (beta_t - m)
  return(out)
}

proposalfunction <- function(mu, sigma) {
  out <- rmvnorm(1, mu, sigma)
  return(as.vector(out))
}

cond_proposaldensity <- function(beta, mu, Fisher) {
  #out <- dmvnorm(beta, mu, solve(Fisher), log = T)
  n <- length(beta)
  out <- -0.5 * n * log(2*pi) - 0.5 * log(1 / det(Fisher)) - 0.5 * t(beta - mu)%*%Fisher%*%(beta - mu)
  return(out)
}


#likelihood for beta ( assuming independence)
loglik_pois <- function(beta_t, sigma_t) {
  out <- sum(dpois(y, lambda = exp(eta_t), log = T))
  return(out)
}

loglik_norm <- function(beta_t, sigma_t) {
 out <- dnorm(y, X%*%beta_t, sigma_t, log = T)
 return(sum(out))
}





