library("MASS")

set.seed(11)

n <- 60

m1 <- c(1.5, 1.5)
S1 <- matrix(c(0.3, 0.05, 0.05, 0.3), ncol = 2)

clus1 <- mvrnorm(n = n, mu = m1, Sigma = S1)

m2 <- c(1.5, -1.5) 
S2 <- matrix(c(0.5, -0.08, -0.08, 0.2), ncol = 2)

clus2 <- mvrnorm(n = n, mu = m2, Sigma = S2)

m3 <- c(-1.5, 1.5)
S3 <- matrix(c(0.1, 0.03, 0.03, 0.1), ncol = 2)

clus3 <- mvrnorm(n = n, mu = m3, Sigma = S3)

m4 <- c(-1.5, -1.5)
S4 <- matrix(c(0.8, 0.5, 0.5, 0.8), ncol = 2)

clus4 <- mvrnorm(n = n, mu = m4, Sigma = S4)

datc <- rbind(clus1, clus2, clus3, clus4)



alpha <- 0.01
mu0 <- matrix(rep(0, 2), ncol = 2, byrow = TRUE)

sigma0 <- diag(2) * 3^2
sigma_y <- diag(2) * 1 

c_init <- rep(1, nrow(datc))

results <- crp_gibbs(data = datc, alpha = alpha, 
                     mu0 = mu0, sigma0 = sigma0,
                     sigma_y = sigma_y, c_init = rep(1, nrow(datc)))

tab <- apply(results, 1, FUN = function(x){
  tab <- table(x)
  ans <- names(tab[which.max(tab)])
  return(ans)
})

table(tab)
