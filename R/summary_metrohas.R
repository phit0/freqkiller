summary.metrohas <- function(result) {
  chain <- result$chain
  res <- list()
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

  ###########################################################
  ####                 output             ###################
  ###########################################################

  ## call
  cat("\nCall: \nmetrohas(formula = ", deparse(result$formula),
      ",", "number_it =", result$number_it, ",", "dist =", result$dist, ")\n")

  ## start values
  cat("\nStartvalues for beta: ", round(result$beta_start, digits = 4),
              collapse = " ", "\n")
  ## startvalue for sigma of the gibbs sampler (normal only)
  if (result$dist == "normal") {
    cat("\nStartvalue for variance: ", result$sigma2_start, "\n")
  }

  ## print the coefficients
  cat("\nCoefficients:\n")
  print(res$coefficients, digits = 4, print.gap = 2)

  ## print prior parameters
  cat("\n", rep("~", 30), "\n")
  cat("Prior asumptions for estimated parameters: \nCovariance \'M\': \n")
  print(result$M)
  cat("\nExpected value \'m\': \n")
  print(result$m)
  cat(rep("~", 30), "\n")

  ## Number of iterations
  cat("\nNumber of Metropolis-Hastings iterations: ",
      (result$number_it + result$burnin) * result$thinning_lag, "\n")
  ## Burnin
  cat("Burn in iterations: ", result$burnin, "\n")
  ## thinning
  if (result$thinning == 1) {
    cat("No thinning was  performed (!)\n")
  }else{
    cat("Thinning: Every k =", result$thinning_lag, " sample was retained\n")
  }
  cat(rep("-", 30), "\n")
  ## Number of samples drawn
  cat("Number of remaining samples: ", result$number_it, "\n")

  invisible(res)
}


