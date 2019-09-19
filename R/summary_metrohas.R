#' Title
#'
#' @param result Object of class frequentistkiller
#'
#' @return list with all the elements of the summary output to the console
#' @export
#'
#' @examples
summary.frequentistkiller <- function(result) {
  chain <- result$chain
  res <- list()

  ## call
  res$call <- paste("Call: frequentistkiller(", deparse(result$formula), ",",
                    "number_it =", result$number_it,",", "dist =", result$dist,")")

   ## start values
  res$beta_start <- paste("Starting values for beta: ",
                          deparse(round(result$beta_start, digits = 4)),"")

  ## startvalues of the gibbs sampler (normal only)
  if (result$dist == "normal") {
    res$gibbs <- paste("Prior parameters for variance: ", result$a0, result$b0, "")
  }
  ## Summary statistics for the chain
  res$t1 <- "Summary statistics for the sample:"
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
  res$t2 <- "Prior asumptions for estimated parameters:"
  res$t3 <- "covariance matrix \"M\":"
  res$M <- result$M
  res$t4 <- "expected value \"m\":"
   res$m <- result$m
  res$br2 <-  rep("~", 30)

 ## Burnin
  res$burnin <- paste("Burn in iterations: ", result$burnin, "")
  ## thinning
  if (result$thinning == 1) {
    res$thin <- "No thinning was  performed (!)"
  }else{
    res$thin <- paste("Thinning: Every k =", result$thinning_lag, " sample was retained")
  }
  res$br3 <- rep("-", 30)
  ## Number of samples drawn
  res$samples <- paste("Number of remaining samples: ", nrow(result$chain))

  lapply(res, FUN = print, quote = FALSE)
  invisible(res)
}


