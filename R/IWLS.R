
# write a test for w_func that checks it has appropriate dimensions

w_func <- function(sigma_t,beta_t) {
  out <- switch(dist,
                "normal" = diag(length(y))/ sigma_t,
                "poisson" = "CODE FEHLT HIER!")
  return(out)
}

fisher_func <- function(sigma_t,beta_t) {
  out <- switch(dist,
                "normal" = (1/sigma_t)*t(X) %*% X + M_1,
                "poisson" = "CODE FEHLT HIER")
  return(out)
}

y_wgl_func <- function(sigma_t, beta_t) {
  out <- switch(dist,
         "normal" = y,
         "poisson" = eta_t + ((y - exp(eta_t, dist)) / exp(eta_t, dist)))
  return(out)
}

mu_func <- function(sigma_t,beta_t) {
  out <- switch(dist,
                "normal" = solve(fisher_func(sigma_t,beta_t)) %*% t(X) %*% w_func(sigma_t,beta_t) %*% y_wgl_func(sigma_t,beta_t) + M_1 %*% m,
                "poisson" = "CODE FEHLT HIER")
  return(out)
}

#beta prior
prior_func <- function(beta_t) {
  out <- dmvnorm(beta_t, mean = m, sigma = M, log = T)
  return(out)
}

proposalfunction <- function(mu, sigma) {
  out <- rmvnorm(1, mu, sigma)
  return(as.vector(out))
}

cond_proposaldensity <- function(beta, mu, sigma) {
  out <- dmvnorm(beta, mu, sigma, log = T)
  return(out)
}


#likelihood for beta (without assuming independence)
loglik_func <- function(beta_t, sigma_t) {
  out <- switch(dist,
                "normal" = dmvnorm(y, X%*%beta_t, diag(length(y))*sigma_t, log = T),
                "poisson" = sum(dpois(y, lambda = exp(eta_t, dist), log = T)))

  return(out)
}

