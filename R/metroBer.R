metroBer <- function(formula, beta_start, anzahl_sim, m, M){

  X <- model.matrix(formula)
  y <- as.matrix(model.frame(formula)[paste(formula[2])])[,1]

  chain <- array(dim = c(anzahl_sim + 1, length(beta_start)))
  chain[1,] <- beta_start

  for (i in 1:anzahl_sim) {

    beta_t <- chain[i,]
    sigma2_t <- dh(X%*%beta_t) #this holds due to var(y) = h(eta)*(1-h(eta)) = dh(eta)

    # IWLS
    F_t <- fisher_func(sigma2_t, beta_t)
    mu_t <- mu_func(sigma2_t, beta_t)

    # Pick proposal
    proposal <- proposalfunction(mu_func(sigma2_t,  beta_t), sigma = solve(fisher_func(sigma2_t,beta_t)))

    # IWLS
    F_star <- fisher_func(sigma2_t, proposal)
    mu_star <- mu_func(sigma2_t, proposal)

    q_cond_star <- cond_proposaldensity(chain[i,], mu_star, F_star)
    q_cond_t <- cond_proposaldensity(proposal, mu_t, F_t)  #invert in function to avoid re-inverting

    #Posterior
    prior_t <- prior_func(chain[i,])
    prior_star <- prior_func(proposal)

    loglik_t <- loglik_func(chain[i,], sigma2_t)
    loglik_star <- loglik_func(proposal, sigma2_t)

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
