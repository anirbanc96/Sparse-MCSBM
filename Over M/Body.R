source("Functions.R")

set.seed(2024)

threshold.vec <- c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4)

n <- 400

p <- 500

sigma <- 2*rbinom(n, 1, 0.5)-1

gamma <- 4/5

n.resamp <- 25

m.vec <- c(1, 3, 5, 7, 10)

power <- matrix(0, length(m.vec), length(threshold.vec))
overlap <- matrix(0, length(m.vec), length(threshold.vec))

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

for (k in 1:length(m.vec)){
  
  m <- m.vec[k]
  
  avg.deg.vec <- sample(2:5, m, replace = T)
  
  for (i in 1:length(threshold.vec)){
    
    lambda.vec <- rep(sqrt(threshold.vec[i]/(m + 2)), m)
    
    mu <- sqrt(2) * sqrt(gamma) * sqrt(threshold.vec[i]/(m + 2))
    
    print(sum(lambda.vec^2) + (mu^2)/gamma)
    
    writeLines(c("Starting choice ", i,"\n"), "log.txt")
    
    outputMatrix <- foreach(j=1:n.resamp, .combine=rbind) %dopar% {
      
      cat(paste("Starting iteration ", j, "\n"), file = "log.txt", append = T)
      
      output <- power.comp(n, m, p, gamma, sigma, mu, lambda.vec, avg.deg.vec)
      
      cat(paste("Finished iteration ", j, "\n"), file = "log.txt", append = T)
      
      output
    }
    
    
    power[k,i] <- mean(outputMatrix[,1])
    overlap[k,i] <- mean(outputMatrix[,2])
    
    print (c(m.vec[k], threshold.vec[i], power[k,i], overlap[k,i]))
    
  }
}

Sys.time() - start

library(tidyverse)
power.overlap.tibble <- as_tibble(rbind(threshold.vec, power, overlap))

write_csv(power.overlap.tibble, file = "Vals.csv")