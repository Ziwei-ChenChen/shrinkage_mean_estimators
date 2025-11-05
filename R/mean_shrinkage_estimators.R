St_MSh_estimator <- function(x) {
  Tn <- nrow(x)
  N_T <- ncol(x)
  x_bar <- colMeans(x)
  a1T_hat <- sum(x_bar^2) / N_T
  trace_Sn <- sum(diag(cov(x)))
  a2T_hat <- a1T_hat + trace_Sn/(Tn * N_T)
  c_star_hat <- a1T_hat / a2T_hat
  theta_St_MSh <- c_star_hat * x_bar
  return(as.vector(theta_St_MSh))
}

D_MSh_estimator <- function(x) {
  Tn <- nrow(x)
  N_T <- ncol(x)
  x_bar <- colMeans(x)
  Sn <- cov(x)
  c_star_hat <- (x_bar^2) / (x_bar^2 + (1/Tn) * diag(Sn))
  theta_D_MSh <- c_star_hat * x_bar
  return(as.vector(theta_D_MSh))
}

O_LSh_estimator <- function(x) {
  Tn <- nrow(x)
  N_T <- ncol(x)
  x_bar <- colMeans(x)
  a1T_hat <- sum(x_bar^2) / N_T
  trace_Sn <- sum(diag(cov(x)))
  a2T_hat <- a1T_hat + trace_Sn/(Tn * N_T)
  s_i <- rowSums(x)
  term1 <- sum(s_i^2)
  term2 <- ((sum(s_i))^2 - term1) / (Tn - 1)
  d2_hat <- (term1 - term2) / (Tn^2 * N_T^2)
  grand_mean <- mean(x_bar)
  d3_hat <- sum((x_bar - grand_mean)^2) / N_T
  tilde_a1_hat <- a2T_hat - a1T_hat - d2_hat
  tilde_a2_hat <- tilde_a1_hat + d3_hat
  delta_star_hat <- tilde_a1_hat / tilde_a2_hat
  theta_O_LSh <- (1 - delta_star_hat) * x_bar +
    delta_star_hat * rep(grand_mean, N_T)
  return(as.vector(theta_O_LSh))
}

T_LSh_estimator <- function(x) {
  Tn <- nrow(x)
  N_T <- ncol(x)
  x_bar <- colMeans(x)
  a1T_hat <- sum(x_bar^2) / N_T
  trace_Sn <- sum(diag(cov(x)))
  a2T_hat <- a1T_hat + trace_Sn/(Tn * N_T)
  s_i <- rowSums(x)
  term1 <- sum(s_i^2)
  term2 <- ((sum(s_i))^2 - term1) / (Tn - 1)
  d2_hat <- (term1 - term2) / (Tn^2 * N_T^2)
  grand_mean <- mean(x_bar)
  d3_hat <- sum((x_bar - grand_mean)^2) / N_T
  d4_hat <- grand_mean^2
  tilde_a1_hat <- a2T_hat - a1T_hat - d2_hat
  tilde_a2_hat <- tilde_a1_hat + d3_hat
  delta_star_hat <- tilde_a1_hat / tilde_a2_hat
  xi_star_hat <- 1 - (d2_hat/(d2_hat + d4_hat))/delta_star_hat
  theta_T_LSh <- (1 - delta_star_hat)*x_bar +
    delta_star_hat * xi_star_hat * rep(grand_mean, N_T)
  return(as.vector(theta_T_LSh))
}

