context("")


# test_that("sigma_vec is positive",{
#   expect_true(sigma_vec(c(0, 0)) >= 0)
#   expect_true(sigma_vec(0) >= 0)
# })
rnorm(4)

test_that("w_func works", {
t1 <- sum(w_func(sigma2_t = 1, beta_t =  c(1,2), y = c(-0.9, 0.7, -0.35, -0.85),
                 X = cbind(1, c(1, 2, 3, 4)), dist = "normal"))
expect_equal(t1, 4)
t2 <- sum(w_func(sigma2_t = 1, beta_t =  c(1,2), y = c(1, 0, 0, 1),
           X = cbind(1, c(1, 2, 3, 4)), dist = "bernoulli"))
expect_equal(t2, 0.05285832)
t3 <- sum(diag(w_func(sigma2_t = 1, beta_t =  c(1,2), y = c(0, 1, 0, 0),
                 X = cbind(1, c(1, 2, 3, 1.5)), dist = "poisson")))
expect_equal(t3, 1319.73)
})

test_that("w_func returns matrix of correct dimensions", {
  d1 <- w_func(sigma2_t = 1, beta_t =  c(1,2), y = c(-0.9, 0.7, -0.35, -0.85),
         X = cbind(1, c(1, 2, 3, 4)), dist = "normal")
  d2 <- w_func(sigma2_t = 1, beta_t =  c(1,2), y = c(1, 0, 0, 1),
                   X = cbind(1, c(1, 2, 3, 4)), dist = "bernoulli")
  d3 <- w_func(sigma2_t = 1, beta_t =  c(1,2), y = c(0, 1, 0, 0),
                        X = cbind(1, c(1, 2, 3, 1.5)), dist = "poisson")
  d <- list(d1, d2, d3)
  for (i in 1:3) {
    expect_equal(dim(d[[i]])[1], dim(d[[i]])[2])
    expect_equal(dim(d[[i]])[1], 4)
  }
})

