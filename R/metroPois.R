
metroPois <- function(formula, beta_start, anzahl_sim, m, M){

  X <- model.matrix(formula)
  y <- as.matrix(model.frame(formula)[paste(formula[2])])[,1]

  chain <- array(dim = c(anzahl_sim + 1, length(beta_start)))
  chain[1,] <- beta_start

  for (i in 1:anzahl_sim) {

    eta_t <- X%*%chain[i, ]
    lambda_t <- h(eta_t)

    # IWLS
    W_t <- w_func(eta_t, lambda_t)
    F_t <- fisher_func(X, W_t, M)
    y_wgl_t <- y_wgl_func(eta_t, y)
    mu_t <- mu_func(X, F_t, W_t, y_wgl_t, M, m)

    # Pick proposal
    proposal <- proposalfunction(mu_t, sigma = solve(F_t))

    #Update eta
    eta_star <- X%*%proposal
    #lambda_star <- h(eta_star)

    # IWLS
    W_star <- w_func(eta_star, lambda_t)
    F_star <- fisher_func(X, W_star, M)
    y_wgl_star <- y_wgl_func(eta_star, y)
    mu_star <- mu_func(X, F_star, W_star, y_wgl_star, M, m)

    q_cond_star <- cond_proposaldensity(chain[i,], mu_star, solve(F_star))
    q_cond_t <- cond_proposaldensity(proposal, mu_t, solve(F_t))

    #Posterior
    prior_t <- prior_func(chain[i,], m, M)
    prior_star <- prior_func(proposal, m, M)

    loglik_t <- loglik_func(eta_t, lambda_t, y)
    loglik_star <- loglik_func(eta_star, lambda_t, y)

    alpha <- min(c((prior_star + loglik_star + q_cond_star)
                   / (prior_t + loglik_t + q_cond_t), 1))

    if (runif(1) < alpha) {
      chain[i+1,] <- proposal
    }else{
      chain[i+1,] <- chain[i,]
    }

  }
  return(chain)
}
