#Pior Parameters
m <- function(beta_t) {
  return(rep(0,length(beta_t)))
}

M <- function(beta_t) {
  return(diag(length(beta_t)))
}

sigma_vec <- function(beta_t) {
  pi(beta_t)*(1-pi(beta_t))
}

w_func <- function(beta_t) {
  w <- (dh(X%*%beta_t))^2 / sigma_vec(beta_t)
  return(diag(c(w)))
}

fisher_func <- function(beta_t) {
  Wt <- w_func(beta_t)
  Ft <- t(X)%*%Wt%*%X + solve(M(beta_t))
  return(Ft)
}

y_wgl_func <- function(beta_t) {
  eta_t <- X%*%beta_t
  y_wgl <- eta_t + (y-h(X%*%beta_t))/dh(X%*%beta_t)
  return(y_wgl)
}

mu_func <- function(beta_t) {
  Wt <- w_func(beta_t)
  Ft <- fisher_func(beta_t)
  yt_wgl <- y_wgl_func(beta_t)
  mu <- solve(Ft)%*%t(X)%*%Wt%*%yt_wgl + solve(M(beta_t))%*%m(beta_t)
  return(mu)
}
