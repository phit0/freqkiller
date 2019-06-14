context("")

test_that("Covariance matrix M of prior for beta a squared Matrix", {
  b <- c(2, 3)
  expect_equal(dim(M(b))[1], dim(M(b))[2])
})

test_that("Covariance matrix M of prior for beta is an invertible Matrix", {
  b <- c(2, 3)
  D <- M(b)
  expect_equal((D%*%solve(D)),  diag(ncol(D))) #inverse*D is identity
  expect_false(round(det(D), 10e12) == 0) #nonzero determinant
  expect_equal(qr(D)$rank, ncol(D)) #full rank matrix
})


# test_that("sigma_vec is positive",{
#   expect_true(sigma_vec(c(0, 0)) >= 0)
#   expect_true(sigma_vec(0) >= 0)
# })



