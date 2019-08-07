
metroNorm <- function(formula, sigma_start, beta_start, a0, b0, anzahl_sim, m, M, dist){

  X <- model.matrix(formula)
  y <- as.matrix(model.frame(formula)[paste(formula[2])])[,1]

  chain <- array(dim = c(anzahl_sim + 1, length(beta_start)))
  chain[1,] <- beta_start
  eta_t <- X%*%beta_start

  s_chain <- array(dim = anzahl_sim + 1)
  s_chain[1] <- sigma_start
  sigma_t <- sigma_start

  a_t <- a0
  b_t <- b0

  for (i in 1:anzahl_sim) {

    # IWLS
    W_t <- w_func(eta_t, sigma_t, dist)
    F_t <- fisher_func(X, W_t, M)
    y_wgl_t <- y_wgl_func(eta_t, y, dist)
    mu_t <- mu_func(X, F_t, W_t, y_wgl_t, M, m)

    # Pick proposal
    proposal <- proposalfunction(mu_t, sigma = solve(F_t))

    #Update eta
    eta_star <- X%*%proposal

    # IWLS
    W_star <- w_func(eta_star, sigma_t, dist)
    F_star <- fisher_func(X, W_star, M)
    y_wgl_star <- y_wgl_func(eta_star, y, dist)
    mu_star <- mu_func(X, F_star, W_star, y_wgl_star, M, m)

    q_cond_star <- cond_proposaldensity(chain[i,], mu_star, solve(F_star))
    q_cond_t <- cond_proposaldensity(proposal, mu_t, solve(F_t))

    #Posterior
    prior_t <- prior_func(chain[i,], m, M)
    prior_star <- prior_func(proposal, m, M)

    loglik_t <- loglik_func(eta_t, sigma_t, y, dist)
    loglik_star <- loglik_func(eta_star, sigma_t, y, dist)

    alpha <- min(c((prior_star + loglik_star + q_cond_star)
                   / (prior_t + loglik_t + q_cond_t), 1))

    if (runif(1) < alpha) {
      chain[i+1,] <- proposal
    }else{
      chain[i+1,] <- chain[i,]
    }
    eta_t <- X%*%chain[i + 1,]

    # update sigma
    a_t <- a_func(y, a_t)
    b_t <- b_func(y, eta_t, b_t)
    sigma_t <- sigma_gibbs(a_t, b_t)
    s_chain[i+1] <- sigma_t
  }
return(cbind(chain, s_chain))
}
