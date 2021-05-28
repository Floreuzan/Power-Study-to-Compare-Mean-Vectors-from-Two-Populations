##### Test for Assumption (5) #####

source("library.R")

### PARAMETER INITIALIZATION ###
S <- 2
p <- 3
mu_0 <- matrix(c(0, 0, 0), ncol = 1)
Sigma_0 <- matrix(c(1, 0.5, 0.5,
                    0.5, 1, 0.5,
                    0.5, 0.5, 1), nrow = 3)
n_seq <- c(10, 100, 1000)
mu_seq <- seq(-2, 2, 0.2)
sigma_seq <- seq(1, 10, 1)

power_table_2 <- tibble(mu = rep(NA_real_, S * length(mu_seq) * length(sigma_seq) * length(n_seq)),
                        sigma = NA_real_, 
                        s = NA_real_, 
                        reject = NA_real_,
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
        
        x <-  mvtnorm::rmvnorm(n, c(mu_0), Sigma_0) 
        x_bar <- matrix(colMeans(x), ncol = 1)
        S_x <- cov(x)
        
        y <- mvtnorm::rmvnorm(n, c(mu_a), Sigma_a) 
        y_bar <- matrix(colMeans(y), ncol = 1)
        S_y <- cov(y)
        
        power_table_2$reject[i] <- 
          (t(x_bar-y_bar) %*% solve((1/n) * (S_x + S_y)) %*% (x_bar-y_bar)) > qchisq(0.99, p)
        
        power_table_2$mu[i] <- mu
        power_table_2$sigma[i] <- sigma
        power_table_2$s[i] <- s
        power_table_2$n[i] <- n
        i <- i + 1
        
      }
    }
  }
}

save(power_table_2, file = "data/power_table_2.Rdata")