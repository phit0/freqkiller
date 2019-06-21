m <- function(beta_t) {
  return(rep(0,length(beta_t)))
}

M <- function(beta_t) {
  return(1*diag(length(beta_t)))
}

# write a test for w_func that checks it has appropriate dimensions
#
w_func <- function(eta, sigma_t) {
  w <- (dh(eta))^2 / sigma_t
  return(diag(c(w)))
}

fisher_func <- function(beta_t, X, W_t) {
  Ft <- t(X)%*%Wt%*%X + solve(M(beta_t))
  return(Ft)
}

y_wgl_func <- function(beta_t, X, y) {
  eta_t <- X%*%beta_t
  y_wgl <- eta_t + ((y-h(eta_t)) / dh(eta_t))
  return(y_wgl)
}

mu_func <- function(beta_t, X, y) {
  Wt <- w_func(eta, a_t, b_t, y)
  Ft <- fisher_func(beta_t, X,a_t, b_t, y)
  yt_wgl <- y_wgl_func(beta_t, X, y)
  mu <- solve(Ft)%*%t(X)%*%Wt%*%yt_wgl + solve(M(beta_t))%*%m(beta_t)
  return(mu)
}

proposalfunction <- function(beta_t, X, y){
  mu <- mu_func(beta_t, X, y)
  sigma <- solve(fisher_func(beta_t, X,a_t, b_t, y))
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
