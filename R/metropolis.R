
proposalfunction = function(x){
  mu <- mu_func(x)
  sigma <- solve(fisher_func(x))
  output <- rmvnorm(1, mu, sigma)
  return(as.vector(output))
}

#Funktion q(beta_t gegeben beta_stern) oder umgekehrt (wie man will))
post_proposal2 <- function(x,y) {
  output <- dmvnorm(x, mean = mu_func(y), sigma = solve(fisher_func(y)), log=T)
  return(output)
}

metrohas <- function(startvalue, anzahl_sim){
  chain <- array(dim = c(anzahl_sim + 1, length(startvalue)))
  chain[1,] <- startvalue

  for (i in 1:anzahl_sim) {
    proposal <- proposalfunction(chain[i,])

    enumerator <- (post_beta(proposal) + post_proposal2(chain[i,], proposal))
    nominator <- (post_beta(chain[i, ]) + post_proposal2(proposal, chain[i, ]))
    alpha <- min(c(enumerator / nominator, 1))

    if (runif(1) < alpha) {
      chain[i+1,] <- proposal
    }else{
      chain[i+1,] <- chain[i,]
    }
  }
  print("DONE")
  return(chain)
}
