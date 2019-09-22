context("Algorithm runs")

test_that("Algorithm and summary method work", {
  set.seed(42)
  n <- 10
  a <- runif(n)
  eta <- 0.2 + 2 * a
  df1 <- data.frame(y = 0.2 + 3.8 * a + rnorm(n), x = a)
  df2 <- data.frame(y = rpois(n, exp(eta)), x = a)
  df3 <- data.frame(y = rbinom(10, 1, exp(eta)/(1 + exp(eta))), x = a)

  expect_error((t1 <- frequentistkiller(y ~ x - 1, df1, dist = "normal",
                                number_it = 1000,
                                burnin = 0)), NA)
  expect_error((t2 <- frequentistkiller(y ~ x, df2, dist = "poisson",
                                 number_it = 1000,
                                 burnin = 0)), NA)
  expect_error((t3 <- frequentistkiller(y ~ x, df3, dist = "bernoulli",
                                 number_it = 1000,
                                 burnin = 0)), NA)
  expect_error(summary(t1, notify = FALSE), NA)
  expect_error(summary(t2, notify = FALSE), NA)
  expect_error(summary(t3, notify = FALSE), NA)
})
