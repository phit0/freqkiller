summary.frequentistkiller <- function(result) {
  chain <- result$chain
  res <- list()

  ## call
  res$call <- paste("Call: \n  frequentistkiller(formula = ", deparse(result$formula),
                    ",", "number_it =", result$number_it, ",", "dist =", result$dist,")\n\n")

   ## start values
  res$beta_start <- paste("Starting values for beta: ",
                          deparse(round(result$beta_start, digits = 4)),"\n\n")

  ## startvalues of the gibbs sampler (normal only)
  if (result$dist == "normal") {
    res$gibbs <- paste("Prior parameters for variance: ", result$a0, result$b0, "\n\n")
  }
  ## Summary statistics for the chain
  res$t1 <- "Summary statistics for the sample:\n"
  if (is.null(ncol(chain))) { # if univariate
    res$coefficients <- matrix(c(mean(chain), median(chain), sd(chain),
                         quantile(chain, probs = 0.025), quantile(chain, probs = 0.975)),
                         nrow = 1)
    colnames(res$coefficients) <- c("mean", "median", "sd", "2.5% quantile", "97.5% quantile")
    rownames(res$coefficients) <- ifelse(is.null(names(chain)), "x1", names(chain))
  }else{ # if multivariate
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
  }

  ## print prior parameters
  res$br1 <- rep("~", 30)
  res$t2 <- "\nPrior asumptions for estimated parameters:\n"
  res$t3 <- "covariance matrix \"M\":\n"
  res$M <- result$M
  res$t4 <- "expected value \"m\":\n"
   res$m <- result$m
  res$br2 <-  rep("~", 30)

 ## Burnin
  res$burnin <- paste("\n\nBurn in iterations: ", result$burnin, "\n")
  ## thinning
  if (result$thinning == 1) {
    res$thin <- "No thinning was  performed (!)\n"
  }else{
    res$thin <- paste("Thinning: Every k =", result$thinning_lag, " sample was retained\n")
  }
  res$br3 <- rep("-", 30)
  ## Number of samples drawn
  res$samples <- paste("\nNumber of remaining samples: ", nrow(result$chain))

  #printing
  lapply(res[1:4], cat)
  print(res$coefficients, quote = FALSE)
  lapply(res[6:8], cat)
  print(res$M, quote = FALSE)
  cat(res$t4)
  print(res$m, quote = FALSE)
  lapply(res[12:16], cat)
  invisible(res)
}


