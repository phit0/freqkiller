summary.metrohas <- function(result) {
  chain <- result$chain
  res <- list()
  res$coefficients <- matrix(NA, nrow = ncol(chain), ncol = 5)
  colnames(res$coefficients) <- c("mean", "median", "sd", "2.5% quantile", "97.5% quantile")
  rownames(res$coefficients) <- colnames(chain)
  # calculation of the chain statistics
  for (i in 1:ncol(chain)) {
    res$coefficients[i, "mean"] <- mean(chain[, i])
    res$coefficients[i, "median"] <- median(chain[, i])
    res$coefficients[i, "sd"] <- sd(chain[, i])
    res$coefficients[i, "2.5% quantile"] <- quantile(chain[, i], probs = 0.025)
    res$coefficients[i, "97.5% quantile"] <- quantile(chain[, i], probs = 0.975)
  }
  # output
  print(paste("number of iterations: ", result$number_it))
  print(paste(c("starting values for beta: ", result$beta_start), collapse = " ", sep = ","))
  class(res) <- "summary.metrohas"
  print(res)
  if (result$thinning == 0) {
    print('No thinning was  performed (!)')
  }else{
    print(paste("Thinning was performed at lags ", result$thinning))
  }
}


