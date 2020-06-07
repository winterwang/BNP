crp_gibbs <- function(data, alpha = 0.01, mu0, sigma0, sigma_y, c_init, maxIters = 1000) {
  require(mvtnorm)
  
  data_dim <- ncol(data)
  
  N <- nrow(data)
  
##%######################################################%##
#                                                          #
####                       Priors                       ####
#                                                          #
##%######################################################%##
  tau0 <- solve(sigma0)
  tau_y <- solve(sigma_y)
  
  z <- c_init
  n_k <- as.vector(table(z))
  Nclust <- length(n_k)
  
##%######################################################%##
#                                                          #
####  Chinese Restaurant Process Gibbs sampler begins   ####
#                                                          #
##%######################################################%##

  
  res <- matrix(NA, nrow = N, ncol = maxIters) #cluster membership storage
  pb <- txtProgressBar(min = 0, max = maxIters, style = 3)
  
  for(iter in 1:maxIters) {
    for(n in 1:N) {
      c_i <- z[n] 
      n_k[c_i] <- n_k[c_i] - 1
      if (n_k[c_i] == 0) 
      {
        n_k[c_i] <- n_k[Nclust]
        loc_z <- (z == Nclust)
        z[loc_z] <- c_i
        n_k <- n_k[ -Nclust ]
        Nclust <- Nclust - 1
      }
      
      z[n] <- -1
      
      logp <- rep(NA, Nclust + 1)
      
      for (c_i in 1:Nclust) {
        tau_p <- tau0 + n_k[c_i] * tau_y
        sig_p <- solve(tau_p)
        
        loc_z <- which(z == c_i)
        if(length(loc_z) > 1) {
          sum_data <- colSums(data[z == c_i, ])
        } else {
          sum_data <- data[z == c_i, ]
        }
        
        mean_p <- sig_p %*% (tau_y %*% sum_data + tau0 %*% t(mu0))
        logp[c_i] <- log(n_k[c_i]) + 
          dmvnorm(data[n,], mean = mean_p, sigma = sig_p + sigma_y, log = TRUE)
      }
      logp[Nclust + 1] <- log(alpha) + 
        dmvnorm(data[n,], mean = mean_p, sigma = sigma0 + sigma_y, log = TRUE)
      
      max_logp <- max(logp)
      logp <- logp - max_logp
      loc_prob <- exp(logp)
      loc_prob <- loc_prob / sum(loc_prob)
      
      newz <- sample(1:(Nclust + 1), 1, replace = TRUE, prob = loc_prob)
      
      if (newz == Nclust + 1) {
        n_k <- c(n_k, 0)
        Nclust <- Nclust + 1
      }
      z[n] <- newz
      n_k[newz] <- n_k[newz] + 1
    }
    setTxtProgressBar(pb, iter)
    res[, iter] <- z
  }
  close(pb)
  invisible(res)
}
