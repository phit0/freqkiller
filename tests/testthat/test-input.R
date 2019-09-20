context("check input")

set.seed(42)
a <- runif(10)
df <- data.frame(y = 0.2 + 3.8 * a + rnorm(10), x = a)

test_that("Data type can be used", {
  y1 <- c(rep("a", 4), rep("b", 6))
  x1 <- runif(10)
  expect_error(
    frequentistkiller(y1 ~ x1, dist = "normal", number_it = 1000),
    "response variable has to be numeric")
})


test_that("correct distribution name", {
  expect_error(
    frequentistkiller(y ~ x, df, dist = "typo", beta_start = c(1,1), number_it = 1000),
  "Wrong distribution name. Choose one of the implemented distributions:
         \"normal\", \"poisson\" or \"bernoulli\".")
  })

test_that("correct beta", {
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", beta_start = "typo",
  number_it = 100),
  "beta_start can be either a numeric
               vector of appropriate length or default \"ml_estimate\"")

  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", beta_start = 1, number_it = 1000),
    "length of \"beta_start\" must equal the number of covariables in the model")
})

test_that("m and M have correct dimensions", {
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", m = c(1,1,1), number_it = 1000),
    "\"m\" must be of the same lenght as \"beta_start\"!")

  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", M = diag(3), number_it = 1000),
    paste("\"M\" must be a square matrix of ", 2, "x", 2))
})

test_that("a0 and b0 are either numeric or integer scalar values", {
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", a0 = "d", number_it = 1000),
    "a0 and b0 have to be either numeric or integer scalar values")
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", b0 = c(1,1), number_it = 1000),
    "a0 and b0 have to be either numeric or integer scalar values")
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", b0 = "r", number_it = 1000),
    "a0 and b0 have to be either numeric or integer scalar values")
})

test_that("iteration arguments are positive numbers", {
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", number_it = "d"),
    "\"number_it\" \"burnin\" and \"thinning_lag\"
         must be positive integer numbers")
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", number_it = 1000, burnin = "d"),
    "\"number_it\" \"burnin\" and \"thinning_lag\"
         must be positive integer numbers")
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", number_it = 1000,
                      thinning_lag = "d"),
    "\"number_it\" \"burnin\" and \"thinning_lag\"
         must be positive integer numbers")
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", number_it = 1000,
                      thinning_lag = -1),
    "\"number_it\" \"burnin\" and \"thinning_lag\"
         must be positive integer numbers")
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", number_it = 1000,
                       burnin = -10),
    "\"number_it\" \"burnin\" and \"thinning_lag\"
         must be positive integer numbers")
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", number_it = 0),
    "\"number_it\" \"burnin\" and \"thinning_lag\"
         must be positive integer numbers")
})

test_that("iteration arguments are scalar integers", {
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", number_it = 1000,
                      thinning_lag = 10.9),
    "\"number_it\" \"burnin\" and \"thinning_lag\"
         have to be scalar integers")
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", number_it = 0.9),
    "\"number_it\" \"burnin\" and \"thinning_lag\"
         have to be scalar integers")
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", number_it = 1000, burnin = 10.1),
    "\"number_it\" \"burnin\" and \"thinning_lag\"
         have to be scalar integers")
})

test_that( "number iterations are larger than burnin and thinning lag",{
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", number_it = 1000, burnin = 1000),
    "\"number_it\" must be larger than \"burnin\" and \"thinning_lag\"")
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", number_it = 1000, burnin = 1002),
    "\"number_it\" must be larger than \"burnin\" and \"thinning_lag\"")
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", number_it = 1000,
                      thinning_lag = 1000),
    "\"number_it\" must be larger than \"burnin\" and \"thinning_lag\"")
  expect_error(
    frequentistkiller(y ~ x, df, dist = "normal", number_it = 1000,
                      thinning_lag = 1002),
    "\"number_it\" must be larger than \"burnin\" and \"thinning_lag\"")
})









