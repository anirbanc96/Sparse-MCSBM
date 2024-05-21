source("comdet-dcmase.R")
source("comdetmethods.R")
source("dcmase.R") 
source("getElbows.R")
source("make_dcsbm_plots.R")
source("run_graph_tool.R")
source("SpectralMethods.R")

generate.Adjacency.mat <- function(n, a, b, sigma){
  
  #----------------------------------------------------------------------------#
  # INPUT:
  # n       <- number of nodes
  # a       <- probability of intra community edge
  # b       <- probability of inter community edge
  # sigma   <- community assignment vector
  #----------------------------------------------------------------------------#
  # OUTPUT:
  # A       <- Adjacency matrix of the two community SBM
  #----------------------------------------------------------------------------#
  
  A <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      prob <- (1/n)*(a*(sigma[i] == sigma[j]) + b*(sigma[i]!=sigma[j]))
      A[i,j] <- rbinom(1, 1, prob)
    }
  }
  
  A <- A + t(A)
  return (A)
}

generate.Adjacency.list <- function(n, m, a.vec, b.vec, sigma){
  
  #----------------------------------------------------------------------------#
  # INPUT:
  # n              <- number of nodes
  # m              <- number of graphs/networks
  # a.vec          <- vector of intra community probability for m networks
  # b.vec          <- vector of inter community probability for m networks
  # sigma          <- community assignment vector
  #----------------------------------------------------------------------------#
  # OUTPUT:
  # Adjacency.list <- List containing m adjacency matrices.
  #----------------------------------------------------------------------------#
  
  Adjacency.list <- lapply(1:m, 
                           function(x){generate.Adjacency.mat(n, a.vec[x],
                                                              b.vec[x], sigma)})
  return (Adjacency.list)
  
}

overlap.comdet_dcmase <- function(n, m, gamma, sigma, lambda.vec, d.vec){

  b.vec <- d.vec - lambda.vec * sqrt(d.vec)
  a.vec <- d.vec + lambda.vec * sqrt(d.vec)
  
  Adjacency.list <- generate.Adjacency.list(n, m, a.vec, b.vec, sigma)  
  
  community.labels <- comdet_dcmase(Adjacency.list, K = 2)$community_memberships
  
  sigma.hat <- ifelse(community.labels == 2, -1, community.labels)
  
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
      
      output <- overlap.comdet_dcmase(n, m, gamma, sigma, lambda.vec, avg.deg.vec)
      
      cat(paste("Finished iteration ", j, "\n"), file = "log.txt", append = T)
      
      output
    }
    
    overlap[k,i] <- mean(outputMatrix)
    
    print (c(lambda1^2, sum(lambda.vec^2), 
             threshold.vec[i], overlap[k,i]))
    
  }
}

Sys.time() - start

overlap.tibble <- rbind(threshold.vec, overlap)

write.csv(overlap.tibble, file = "ValsCompDCMASE.csv")