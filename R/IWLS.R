
# write a test for w_func that checks it has appropriate dimensions
#
w_func <- function(beta_t, X, y) {
  s <- sigma_gibbs(y, eta, a_t, b_t)
  w = (dh(X%*%beta_t))^2 / s
  return(diag(c(w)))
}

fisher_func <- function(beta_t, X, y) {
  Wt <- w_func(beta_t, X, y)
  Ft <- t(X)%*%Wt%*%X + solve(M(beta_t))
  return(Ft)
}

y_wgl_func <- function(beta_t, X, y) {
  eta_t <- X%*%beta_t
  y_wgl <- eta_t + ((y-h(eta_t)) / dh(eta_t))
  return(y_wgl)
}

mu_func <- function(beta_t, X, y) {
  Wt <- w_func(beta_t, X, y)
  Ft <- fisher_func(beta_t, X, y)
  yt_wgl <- y_wgl_func(beta_t, X, y)
  mu <- solve(Ft)%*%t(X)%*%Wt%*%yt_wgl + solve(M(beta_t))%*%m(beta_t)
  return(mu)
}

proposalfunction <- function(beta_t, X, y){
  mu <- mu_func(beta_t, X, y)
  sigma <- solve(fisher_func(beta_t, X, y))
  output <- rmvnorm(1, mu, sigma)
  return(as.vector(output))
}


#' conditional proposal density q
#'
#' @param x beta_t or
#' @param y beta_star
#'
#' @return q(beta_t | beta_star) or q(beta_star| beta_t)
#' @export
#'
#' @examples cond_proposaldensity(beta_t, beta_star)
cond_proposaldensity <- function(beta1, beta2, X, y ) {
  output <- dmvnorm(beta1,
                    mean = mu_func(beta2, X, y),
                    sigma = solve(fisher_func(beta2, X, y)), log=T)
  return(output)
}
