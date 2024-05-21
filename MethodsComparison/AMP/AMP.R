generate.B.mat <- function(n, p, mu, sigma){
  
  #----------------------------------------------------------------------------#
  # INPUT: 
  # n         <- number of nodes
  # p         <- dimension of covariates
  # mu        <- signal-to-noise ratio of covariate information
  # sigma     <- community assignment vector
  #----------------------------------------------------------------------------#
  # OUTPUT: 
  # B/sqrt(p) <- the covariate information matrix scaled by sqrt(dimension)
  #----------------------------------------------------------------------------#
  
  B <- matrix(0, nrow = p, ncol = n)
  
  u <- rnorm(p, mean = 0, sd = 1)
  
  for (j in 1:n){
    
    
    
    B[ ,j] <- (sqrt(mu/n)*sigma[j]) * u + 
      rnorm(p, mean = 0, sd = 1)
    
  }
  
  return (t(B)/sqrt(p))
}

u_function <- function(x, mu, sigma){
  print(dim(x))
  tanh((mu/sigma)*x)
}

v_function <- function(v, beta, nu){
  beta*v/(beta^2+nu)
}

p_onsager <- function(x, mu, sigma){
  h <- 1/(cosh((mu/sigma)*x))
  (mu/sigma)*h^2
}

c_onsager <- function(v, beta, nu){
  h <- beta/(beta^2+nu)
}

run.amp <- function(B, c, mu, n.iter){
  
  t <- 1
  
  svd_result <- svd(B)
  
  U <- (svd_result$u[,1]) * sqrt(n)
  sd <- sqrt(mu * c)
  Sd <- diag(sd)  
  V <- (svd_result$v[,1]) * sqrt(p)
  
  
  beta <- mean(V)
  nu <- (1 + (1/c) * sd^2) / (((1/c) * sd^2) * (1 + sd^2))
  sigma <- (1 + sd^2) / ((sd^2) * (1 + (1/c) * sd^2))
  mu <- sqrt(1 - sigma^2)
  
  v <- as.matrix(V)
  f <- as.matrix(U * sqrt(nu))
  
  while(t<=n.iter){
    
    g <- as.matrix(v_function(v, beta, nu))
    c_t <- mean(c_onsager(v, beta, nu))
    u <- B %*% g - (1/c) * c_t * f
    
    sigma <- as.numeric((t(g) %*% g) / n)
    mu <- as.numeric(sigma * sd)
    
    f <- as.matrix(u_function(u, mu, sigma))
    p_t <- mean(p_onsager(u, mu, sigma))
    v <- t(B) %*% f - p_t * g
    
    nu <- as.numeric((t(f) %*% f) / n)
    beta <- as.numeric(nu * sd)
    
    estimate <- f
    t=t+1
  }
  
  return (sign(estimate))
  
}

################################################################################

overlap.amp <- function(n, gamma, sigma, mu, n.iter = 100){
  
  B.mat <- generate.B.mat(n, p, mu, sigma)
  
  sigma.hat <- run.amp(B.mat, gamma, mu, n.iter)
  
  overlap <- max(abs(sum(sigma.hat * sigma)), abs(sum(-sigma.hat * sigma)))/n
  
  return (overlap)
  
}
################################################################################

set.seed(2024)

threshold.vec <- seq(0, 4.5, 0.5)

n <- 400

p <- 500

sigma <- 2*rbinom(n, 1, 0.5)-1

gamma <- 4/5

n.resamp <- 25

m <- 2

avg.deg.vec <- sample(2:5, m, replace = T)

alpha.vec <- c(0.5, 1, 2)

overlap <- matrix(0, length(alpha.vec), length(threshold.vec))

mult <- 3

################################################################################

# Parallelising in Mac/Windows

library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

################################################################################

# Parallelising in cluster

# library(foreach)
# library(doParallel)
# library(snow)
# 
# cores <- strtoi(Sys.getenv("NSLOTS"))
# cl <- makeCluster(cores, methods = FALSE, type = "MPI")
# 
# registerDoParallel(cl)

################################################################################

start <- Sys.time()

for (k in 1:length(alpha.vec)){
  
  for (i in 1:length(threshold.vec)){
    
    alpha <- alpha.vec[k]
    
    lambda1 <- sqrt(mult * alpha * threshold.vec[i]/((1+alpha) * (m + mult-1)))
    
    lambda.vec <- c(lambda1,
                    rep(sqrt(alpha * threshold.vec[i]/((1+alpha) * (m + mult - 1)))))
    
    mu <- sqrt(gamma * threshold.vec[i]/(1+alpha))
    
    print(sum(lambda.vec^2) + ((mu^2)/gamma))
    
    writeLines(c("Starting choice ", i,"\n"), "log.txt")
    
    outputMatrix <- foreach(j=1:n.resamp, .combine=rbind) %dopar% {
      
      cat(paste("Starting iteration ", j, "\n"), file = "log.txt", append = T)
      
      output <- overlap.amp(n, gamma, sigma, mu)
      
      cat(paste("Finished iteration ", j, "\n"), file = "log.txt", append = T)
      
      output
    }
    
    overlap[k,i] <- mean(outputMatrix)
    
    print (c(sum(lambda.vec^2) + ((mu^2)/gamma),
             sum(lambda.vec^2)/((mu^2)/gamma), (mu^2)/gamma,
             overlap[k,i]))
    
  }
}

Sys.time() - start

overlap.tibble <- rbind(threshold.vec, overlap)[-1,]

overlap.tibble <- cbind(alpha.vec, overlap.tibble)

write.csv(overlap.tibble, file = "ValsCompAMP.csv")