##### Test for Assumption (4) and Assumption (5) #####

source("library.R")

### PARAMETER INITIALIZATION ###
S <- 200
p <- 3
mu_0 <- matrix(c(0, 0, 0), ncol = 1)
Sigma_0 <-  matrix(c(1, 0, 0,
                     0, 1, 0,
                     0, 0, 1), nrow = 3)
mu_seq <- seq(-2, 2, 0.2)
n_seq <- c(10, 100, 1000)
sigma_seq <- seq(1, 10, 1)

xi <- c(0, 0, 0)
Omega <-  matrix(c(1, 0.5, 0.5,
                   0.5, 1, 0.5,
                   0.5, 0.5, 1), nrow = 3)
alpha <- matrix(c(0,.5,0), ncol = 1)
nu <- 10

# SN params:
dp_sn_initial <- list(xi = c(xi), Omega = Omega, alpha = c(alpha))
cp_sn_list <- dp2cp(dp_sn_initial, "SN")
# ST params:
dp_st_initial <- list(xi = c(xi), Omega = Omega, alpha = c(alpha), nu = 10)
cp_st_list <- dp2cp(dp_st_initial, "ST")

power_table_3 <- tibble(mu = rep(NA_real_, S * length(mu_seq) * length(sigma_seq) * length(n_seq)),
                        s = NA_real_, 
                        sigma = NA_real_, 
                        reject_mv_mv = NA_real_,
                        reject_mv_sn = NA_real_,
                        reject_mv_st = NA_real_,
                        reject_sn_st = NA_real_,
                        n = NA_integer_)

### POWER STUDY ###

set.seed(697)
i <- 1
for (mu in mu_seq) {
  mu_a <- matrix(c(0, 0, mu), ncol = 1)
  for (sigma in sigma_seq){
    Sigma_a <- sigma*Sigma_0
    for (n in n_seq) {
      for (s in 1:S) {
        # multivariate normal 1
        x <-  mvtnorm::rmvnorm(n, c(mu_0), Sigma_0) 
        x_bar <- matrix(colMeans(x), ncol = 1)
        S_x <- cov(x)
        
        # multivariate normal 2
        y <-  mvtnorm::rmvnorm(n, c(mu_a), Sigma_a) 
        y_bar <- matrix(colMeans(y), ncol = 1)
        S_y <- cov(y)
        
        # multivariate skewed normal 1
        cp_sn_list$mean <- c(mu_a)
        cp_sn_list$var.cov <- Sigma_a
        dp_sn_list <- cp2dp(cp_sn_list, "SN")
        z <- sn::rmsn(n = n, dp = dp_sn_list)
        z_bar <- matrix(colMeans(z), ncol = 1)
        S_z <- cov(z)
        
        # multivariate skewed t normal 1
        cp_st_list$mean <- c(mu_a)
        cp_st_list$var.cov <- Sigma_a
        dp_st_list <- cp2dp(cp_st_list, "ST")
        w <-sn::rmst(n = n, dp = dp_st_list)
        w_bar <- matrix(colMeans(w), ncol = 1)
        S_w <- cov(w)
        
        # multivariate skewed normal 2
        cp_sn_list$mean <- c(mu_0)
        cp_sn_list$var.cov <- Sigma_0
        dp_sn_list <- cp2dp(cp_sn_list, "SN")
        u <- sn::rmsn(n = n, dp = dp_sn_list)
        u_bar <- matrix(colMeans(u), ncol = 1)
        S_u <- cov(u)
        
        power_table_3$reject_mv_mv[i] <- 
          (t(x_bar-y_bar) %*% solve((1/n) * (S_x + S_y)) %*% (x_bar-y_bar)) > qchisq(0.99, p)
        power_table_3$reject_mv_sn[i] <- 
          (t(x_bar-z_bar) %*% solve((1/n) * (S_x + S_z)) %*% (x_bar-z_bar)) > qchisq(0.99, p)
        power_table_3$reject_mv_st[i] <- 
          (t(x_bar-w_bar) %*% solve((1/n) * (S_x + S_w)) %*% (x_bar-w_bar)) > qchisq(0.99, p)
        power_table_3$reject_sn_st[i] <- 
          (t(u_bar-w_bar) %*% solve((1/n) * (S_u + S_w)) %*% (u_bar-w_bar)) > qchisq(0.99, p)
        
        power_table_3$mu[i] <- mu
        power_table_3$sigma[i] <- sigma
        power_table_3$s[i] <- s
        power_table_3$n[i] <- n
        i <- i + 1
      }
    }
  }
}

save(power_table_3, file = "data/power_table_3.Rdata")
