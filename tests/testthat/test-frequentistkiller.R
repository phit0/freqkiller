context("Algorithm run")

test_that("Algorithm runs", {
  set.seed(42)
  n <- 10
  a <- runif(n)
  eta <- 0.2 + 2 * a
  df1 <- data.frame(y = 0.2 + 3.8 * a + rnorm(n), x = a)
  df2 <- data.frame(y = rpois(n, exp(eta)), x = a)
  df3 <- data.frame(y = rbinom(10, 1, exp(eta)/(1 + exp(eta))), x = a)

  expect_error(frequentistkiller(y ~ x, df1, dist = "normal",
                                number_it = 1000,
                                burnin = 0), NA)
  expect_error(frequentistkiller(y ~ x, df2, dist = "poisson",
                                 number_it = 1000,
                                 burnin = 0), NA)
  expect_error(frequentistkiller(y ~ x, df3, dist = "bernoulli",
                                 number_it = 1000,
                                 burnin = 0), NA)
})


