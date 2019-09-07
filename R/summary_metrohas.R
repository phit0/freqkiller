summary.metrohas <- function(result) {
  chain <- result$chain
  res <- list()
  if (is.null(ncol(chain))) { # if univariate
    res$coefficients <- matrix(c(mean(chain), median(chain), sd(chain),
                         quantile(chain, probs = 0.025), quantile(chain, probs = 0.975)),
                         nrow = 1)
    colnames(res$coefficients) <- c("mean", "median", "sd", "2.5% quantile", "97.5% quantile")
    rownames(res$coefficients) <- ifelse(is.null(names(x1)), "x1", names(x1))
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
  if (result$dist == "normal") { # startvalue for sigma of the gibbs sampler (normal only)
    cat("\nStartvalue for variance: ", result$sigma2_start, "\n")
  }
  ## Burnin
  cat("\nBurn in iterations: ", result$burnin, "\n")

  ## print the coefficients
  cat("\nCoefficients:\n")
  print(res$coefficients, digits = 4, print.gap = 2)

  ## print prior parameters
  cat(rep("~", 30), "\n")
  cat("Prior asumptions for estimated parameters: \nCovariance Matrix \'M\': \n")
  print(result$M)
  cat("\nExpected value \'m\': \n")
  print(result$m)
  cat(rep("~", 30), "\n")

  if (result$thinning == 0) { # thinning
    cat("\n No thinning was  performed (!)")
  }else{
    cat("\n", "Thinning was performed with k = ", result$thinning)
  }

  cat("\n", "Number of Metropolis-Hastings iterations: ", result$number_it, "\n")

  # define class for summary object
  class(res) <- "summary.metrohas"
}


