##### Test for Assumption (4) #####

source("library.R")

### PARAMETER INITIALIZATION ###
S <- 200
p <- 3
mu_0 <- matrix(c(0, 0, 0), ncol = 1)
Sigma <- matrix(c(1, 0.5, 0.5,
                  0.5, 1, 0.5,
                  0.5, 0.5, 1), nrow = 3)
mu_seq <- seq(-2, 2, 0.2)
n_seq <- c(10, 100, 1000)

xi <- c(0, 0, 0)
Omega <- matrix(c(1, 1.5, 0.5,
                  1.5, 1, 0.5,
                  0.5, 0.5, 1), nrow = 3)
alpha <- matrix(c(0,.5,0), ncol = 1)
nu <- 10

# SN params:
dp_sn_initial <- list(xi = c(xi), Omega = Omega, alpha = c(alpha))
cp_sn_list <- dp2cp(dp_sn_initial, "SN")
# ST params:
dp_st_initial <- list(xi = c(xi), Omega = Omega, alpha = c(alpha), nu = 10)
cp_st_list <- dp2cp(dp_st_initial, "ST")

power_table_1 <- tibble(mu = rep(NA_real_, S * length(mu_seq) * length(n_seq)),
                        s = NA_real_, 
                        reject_mv_mv = NA_real_,
                        reject_sn_sn = NA_real_,
                        reject_st_st = NA_real_,
                        reject_mv_sn = NA_real_,
                        reject_mv_st = NA_real_,
                        reject_sn_st = NA_real_,
                        n = NA_integer_)


### POWER STUDY ###
set.seed(697)
i <- 1
for (mu in mu_seq) {
  mu_a <- matrix(c(0, 0, mu), ncol = 1)
  for (n in n_seq) {
    for (s in 1:S) {
      # multivariate normal 1
      x <-  mvtnorm::rmvnorm(n, c(mu_0), Sigma) 
      x_bar <- matrix(colMeans(x), ncol = 1)
      S_x <- cov(x)
      
      # multivariate normal 2
      y <-  mvtnorm::rmvnorm(n, c(mu_a), Sigma) 
      y_bar <- matrix(colMeans(y), ncol = 1)
      S_y <- cov(y)
      
      # multivariate skewed normal 1
      cp_sn_list$mean <- c(mu_0)
      cp_sn_list$var.cov <- Sigma
      dp_sn_list <- cp2dp(cp_sn_list, "SN")
      u <- sn::rmsn(n = n, dp = dp_sn_list)
      u_bar <- matrix(colMeans(u), ncol = 1)
      S_u <- cov(u)
      
      # multivariate skewed t normal 1
      cp_st_list$mean <- c(mu_0)
      cp_st_list$var.cov <- Sigma
      dp_st_list <- cp2dp(cp_st_list, "ST")
      v <-sn::rmst(n = n, dp = dp_st_list)
      v_bar <- matrix(colMeans(v), ncol = 1)
      S_v <- cov(v)
      
      # multivariate skewed normal 2
      cp_sn_list$mean <- c(mu_a)
      cp_sn_list$var.cov <- Sigma
      dp_sn_list <- cp2dp(cp_sn_list, "SN")
      z <- sn::rmsn(n = n, dp = dp_sn_list)
      z_bar <- matrix(colMeans(z), ncol = 1)
      S_z <- cov(z)
      
      # multivariate skewed t normal 2
      cp_st_list$mean <- c(mu_a)
      cp_st_list$var.cov <- Sigma
      dp_st_list <- cp2dp(cp_st_list, "ST")
      w <-sn::rmst(n = n, dp = dp_st_list)
      w_bar <- matrix(colMeans(w), ncol = 1)
      S_w <- cov(w)
      
      
      # multivariate normal 1 - multivariate normal 2
      power_table_1$reject_mv_mv[i] <- 
        (t(x_bar-y_bar) %*% solve((1/n) * (S_x + S_y)) %*% (x_bar-y_bar)) > qchisq(0.99, p)
      # multivariate skewed normal 1 - multivariate skewed normal 2
      power_table_1$reject_sn_sn[i] <- 
        (t(u_bar-z_bar) %*% solve((1/n) * (S_u + S_z)) %*% (u_bar-z_bar)) > qchisq(0.99, p)
      # multivariate skewed t normal 1 -  multivariate skewed t normal 2
      power_table_1$reject_st_st[i] <- 
        (t(v_bar-w_bar) %*% solve((1/n) * (S_v + S_w)) %*% (v_bar-w_bar)) > qchisq(0.99, p)
      # multivariate normal 1 - multivariate skewed normal 2
      power_table_1$reject_mv_sn[i] <- 
        (t(x_bar-z_bar) %*% solve((1/n) * (S_x +S_z )) %*% (x_bar-z_bar)) > qchisq(0.99, p)
      # multivariate normal 1 - # multivariate skewed t normal 2
      power_table_1$reject_mv_st[i] <- 
        (t(x_bar-w_bar) %*% solve((1/n) * (S_x + S_w)) %*% (x_bar-w_bar)) > qchisq(0.99, p)
      # multivariate skewed normal 1 -  # multivariate skewed t normal 2
      power_table_1$reject_sn_st[i] <- 
        (t(u_bar-w_bar) %*% solve((1/n) * (S_u + S_w)) %*% (u_bar-w_bar)) > qchisq(0.99, p)
      power_table_1$mu[i] <- mu
      power_table_1$s[i] <- s
      power_table_1$n[i] <- n
      i <- i + 1
    }
  }
}

save(power_table_1, file = "data/power_table_1.Rdata")
