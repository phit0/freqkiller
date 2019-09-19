test_that("Algorithm runs", {
  set.seed(42)
  x <- runif(10)
  beta <- c(0.2, 3.8)
  y <- beta[1] + beta[2] * x + rnorm(10)
  test_chain <- frequentistkiller(y ~ x, "normal", number_it = 2000)
 expect_equal(median(test_chain$chain[, 1]), -0.991213)
 expect_error(frequentistkiller(y ~ x, "normal",
                                beta_start = c(2, 4),
                                number_it = 2000), NA)
})


