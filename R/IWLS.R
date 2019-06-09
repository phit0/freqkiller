m <- function(beta_t) {
  return(rep(0,length(beta_t)))
}

M <- function(beta_t) {
  return(1*diag(length(beta_t)))
}

w_func <- function(beta_t) {
  s <- sigma_vec(beta_t)
  w = (dh(X%*%beta_t))^2 / s
  return(diag(c(w)))
}

fisher_func <- function(beta_t) {
  Wt <- w_func(beta_t)
  Ft <- t(X)%*%Wt%*%X + solve(M(beta_t))
  return(Ft)
}

y_wgl_func <- function(beta_t) {
  eta_t <- X%*%beta_t
  y_wgl <- eta_t + (y-h(X%*%beta_t)) / dh(X%*%beta_t)
  return(y_wgl)
}

mu_func <- function(beta_t) {
  Wt <- w_func(beta_t)
  Ft <- fisher_func(beta_t)
  yt_wgl <- y_wgl_func(beta_t)
  mu <- solve(Ft)%*%t(X)%*%Wt%*%yt_wgl + solve(M(beta_t))%*%m(beta_t)
  return(mu)
}

#' sample from IWLS proposal density
#'
#' @param beta_t
#'
#' @return beta_star
#' @export
#'
#' @examples proposalfuntion(c(1, 2, 3))
proposalfunction <- function(beta_t){
  mu <- mu_func(beta_t)
  sigma <- solve(fisher_func(beta_t))
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
cond_proposaldensity <- function(x,y) {
  output <- dmvnorm(x, mean = mu_func(y), sigma = solve(fisher_func(y)), log=T)
  return(output)
}
