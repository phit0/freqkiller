test_that("Algorithm runs", {
  x <- runif(10)
  beta <- c(0.2, 3.8)
  y <- beta[1] + beta[2] * x + rnorm(10)
  test_chain <- metrohas(y ~ x, "normal",
                         beta_start = c(0, 0),
                         number_it = 1000)
  plot(test_chain$chain[, 1], type ="l")
})


