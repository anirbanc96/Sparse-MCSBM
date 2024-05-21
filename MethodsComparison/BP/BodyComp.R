source("Functions.R")

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

power <- matrix(0, length(alpha.vec), length(threshold.vec))
overlap <- matrix(0, length(alpha.vec), length(threshold.vec))

mult <- 3

################################################################################

# Parallelising in Mac/Windows

# library(foreach)
# library(doParallel)
# 
# cores=detectCores()
# cl <- makeCluster(cores[1]-1)
# registerDoParallel(cl)

################################################################################

# Parallelising in cluster

library(foreach)
library(doParallel)
library(snow)

cores <- strtoi(Sys.getenv("NSLOTS"))
cl <- makeCluster(cores, methods = FALSE, type = "MPI")

registerDoParallel(cl)

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
      
      output <- power.comp(n, m, p, gamma, sigma, mu, lambda.vec, avg.deg.vec)
      
      cat(paste("Finished iteration ", j, "\n"), file = "log.txt", append = T)
      
      output
    }
    
    power[k,i] <- mean(outputMatrix[,1])
    overlap[k,i] <- mean(outputMatrix[,2])
    
    print (c(sum(lambda.vec^2) + ((mu^2)/gamma),
             sum(lambda.vec^2)/((mu^2)/gamma), (mu^2)/gamma,
             overlap[k,i]))
    
  }
}

Sys.time() - start

overlap.tibble <- rbind(threshold.vec, overlap)[-1,]

overlap.tibble <- cbind(alpha.vec, overlap.tibble)

write.csv(overlap.tibble, file = "ValsCompRatio.csv")