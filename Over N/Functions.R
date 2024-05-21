
################################################################################

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
  
  u <- rnorm(p, mean = 0, sd = 1)
  
  B <- matrix(0, nrow = p, ncol = n)
  for (j in 1:n){
    
    B[ ,j] <- (sqrt(mu/n)*sigma[j]) * u + 
      rnorm(p, mean = 0, sd = 1)
    
  }
  
  return (B/sqrt(p))
}

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

Adjacency.Neighbourhood.list <- function(n, A){
  
  #----------------------------------------------------------------------------#
  # INPUT:
  # n                  <- number of nodes
  # A                  <- Adjacency matrix of a network
  #----------------------------------------------------------------------------#
  # OUTPUT:
  # neighbourhood.list <- list containing neighbours of all nodes in A
  #----------------------------------------------------------------------------#
  
  neighbourhood.list <- lapply(1:n, function(x){ which(A[x,]!=0)})
  
  return (neighbourhood.list)
}

Neighbourhood.list <- function(n, m, Adjacency.list){
  
  #----------------------------------------------------------------------------#
  # n                  <- number of nodes
  # m                  <- number of graphs/networks
  # Adjacency.list     <- List containing m adjacency matrices.
  #----------------------------------------------------------------------------#
  # OUTPUT:
  # list.neighbourhood <- list containing list of neighbours of nodes in each
  #                       of m networks
  #----------------------------------------------------------------------------#
  
  list.neighbourhood <- lapply(1:m, 
       function(x){Adjacency.Neighbourhood.list(n, Adjacency.list[[x]])})
  
  return(list.neighbourhood)
}

f.rho.l <- function(lambda.l, d.l){
  
  #----------------------------------------------------------------------------#
  # INPUT:
  # lambda.l  <- signal to noise ratio of lth network
  # d.l       <- average degree of lth network
  #----------------------------------------------------------------------------#
  # OUTPUT:   <- value of rho_l
  #----------------------------------------------------------------------------#
  
  return (atanh(lambda.l/sqrt(d.l)))
  
}

f.rho.n.l <- function(n, lambda.l, d.l){
  
  #----------------------------------------------------------------------------#
  # n         <- number of nodes
  # lambda.l  <- signal to noise ratio of lth network
  # d.l       <- average degree of lth network
  #----------------------------------------------------------------------------#
  # OUTPUT:   <- value of rho_{n;l}
  #----------------------------------------------------------------------------#
  
  return (atanh(lambda.l * sqrt(d.l)/(n - d.l)))
  
}

f <- function(z,rho){
  
  #----------------------------------------------------------------------------#
  # INPUT:
  # z       <- first parameter
  # rho     <- second parameter
  #----------------------------------------------------------------------------#
  # OUTPUT:
  # val     <- the value of the function f(z;rho)
  #----------------------------------------------------------------------------#
  
  if (abs(z) > 15){
    return (1/2)
  }
  else{
    val <- 0.5 * log((cosh(z+rho))/cosh(z-rho))
    return (val)
  }
}

################################################################################

eta.i.update.term3.r <- function(r, i, neigh.list, eta.mat.t,
                                 eta.t, rho, rho.n){
  
  #----------------------------------------------------------------------------#
  # INPUT:
  # r          <- index of the network
  # i          <- index of the vertex/node
  # neigh.list <- list containing list of neighbours of nodes in each
  #               of m networks
  # eta.mat.t  <- the current eta matrix
  # eta.t      <- the current eta vector
  # rho        <- vector of rho values
  # rho.n      <- vector of rho_n values
  #----------------------------------------------------------------------------#
  # OUTPUT:    <- third term in the update equations of eta
  #----------------------------------------------------------------------------#
  
  neigh.i <- neigh.list[[r]][[i]]
  
  if (length(neigh.i) > 0){
    
    term1 <- sum(sapply(neigh.i, 
                        function(k){f(eta.mat.t[[r]][k,i], rho[r])}))
    
  }
  
  else{ term1 <- 0}
  
  term2 <- sum(sapply(1:n, function(k){f(eta.t[k], rho.n[r])}))
  
  return (term1 - term2)
  
}

eta.i.update <- function(n, m, i, mu, gamma, B.mat, 
                         eta.mat.t, eta.t, eta.t_1,
                         m.t, tau.t, neigh.list,
                         lambda.vec, d.vec){
  
  #----------------------------------------------------------------------------#
  # INPUT:
  # n          <- number of nodes
  # m          <- number of graphs/networks
  # i          <- index of the vertex/node
  # mu         <- signal-to-noise ratio of covariate information
  # gamma      <- the value of n/p
  # B.mat      <- the covariate matrix
  # eta.mat.t  <- the current eta matrix
  # eta.t      <- the current eta vector
  # eta.t_1    <- the pervious eta vector
  # m.t        <- the current value of m vector
  # tau.t      <- the current value of tau vector
  # neigh.list <- list containing list of neighbours of nodes in each
  #               of m networks
  # lambda.vec <- vector of signal to noise ratio for m networks
  # d.vec      <- vector of average degrees of m networks
  #----------------------------------------------------------------------------#
  # OUTPUT:    <- updated value of eta_i
  #----------------------------------------------------------------------------#
  
  rho <- sapply(1:m, function(l){f.rho.l(lambda.vec[l], d.vec[l])})
  rho.n <- sapply(1:m, function(l){f.rho.n.l(n, lambda.vec[l], d.vec[l])})
  
  term1 <- sqrt(mu/gamma) * sum(B.mat[,i] * m.t)
  
  term2 <- (mu/gamma) * tanh(eta.t_1[i]) * sum((B.mat[,i] ^ 2)/tau.t)
  
  term3 <- sum(sapply(1:m, 
           function(r){eta.i.update.term3.r(r, i, 
                                            neigh.list, eta.mat.t, eta.t,
                                            rho, rho.n)}))
  
  return (term1 - term2 + term3)
  
}


eta.update <- function(n, m, mu, gamma, B.mat, 
                       eta.mat.t, eta.t, eta.t_1,
                       m.t, tau.t, neigh.list,
                       lambda.vec, d.vec){
  
  #----------------------------------------------------------------------------#
  # INPUT:
  # n          <- number of nodes
  # m          <- number of graphs/networks
  # mu         <- signal-to-noise ratio of covariate information
  # gamma      <- the value of n/p
  # B.mat      <- the covariate matrix
  # eta.mat.t  <- the current eta matrix
  # eta.t      <- the current eta vector
  # eta.t_1    <- the pervious eta vector
  # m.t        <- the current value of m vector
  # tau.t      <- the current value of tau vector
  # neigh.list <- list containing list of neighbours of nodes in each
  #               of m networks
  # lambda.vec <- vector of signal to noise ratio for m networks
  # d.vec      <- vector of average degrees of m networks
  #----------------------------------------------------------------------------#
  # OUTPUT:
  # update     <- the updated eta vector
  #----------------------------------------------------------------------------#
  
  update <- sapply(1:n, function(i){eta.i.update(n, m, i, mu, gamma, B.mat, 
                                                 eta.mat.t, eta.t, eta.t_1,
                                                 m.t, tau.t, neigh.list,
                                                 lambda.vec, d.vec)})
  
  return (update)
  
}

################################################################################

eta.mat.i.j.l.update <- function(n, m, i, j, l, mu, gamma, B.mat, 
                                 eta.mat.t, eta.t, eta.t_1,
                                 m.t, tau.t, neigh.list,
                                 lambda.vec, d.vec){
  
  #----------------------------------------------------------------------------#
  # INPUT:
  # n          <- number of nodes
  # m          <- number of graphs/networks
  # i          <- index of the ith vertex/node
  # j          <- index of the jth vertex/node
  # l          <- index of the lth network
  # mu         <- signal-to-noise ratio of covariate information
  # gamma      <- the value of n/p
  # B.mat      <- the covariate matrix
  # eta.mat.t  <- the current eta matrix
  # eta.t      <- the current eta vector
  # eta.t_1    <- the pervious eta vector
  # m.t        <- the current value of m vector
  # tau.t      <- the current value of tau vector
  # neigh.list <- list containing list of neighbours of nodes in each
  #               of m networks
  # lambda.vec <- vector of signal to noise ratio for m networks
  # d.vec      <- vector of average degrees of m networks
  #----------------------------------------------------------------------------#
  # OUTPUT:    <- the updated value of eta_{i->j;l}
  #----------------------------------------------------------------------------#
  
  neigh.i <- neigh.list[[l]][[i]]
  
  if (!(j %in% neigh.i)){ return (0)}
  
  rho <- sapply(1:m, function(l){f.rho.l(lambda.vec[l], d.vec[l])})
  rho.n <- sapply(1:m, function(l){f.rho.n.l(n, lambda.vec[l], d.vec[l])})
  
  term1 <- sqrt(mu/gamma) * sum(B.mat[,i] * m.t)
  
  term2 <- (mu/gamma) * tanh(eta.t_1[i]) * sum((B.mat[,i] ^ 2)/tau.t)
  
  if (m != 1){
    
    term3 <- sum(sapply(setdiff(1:m,l), 
             function(r){eta.i.update.term3.r(r, i, 
                                              neigh.list, eta.mat.t, eta.t,
                                              rho, rho.n)}))
  }
  
  else{ term3 <- 0}
  
  update.term.2 <- sum(sapply(1:n,function(x){f(eta.t[x],rho.n[l])}))
  
  effective.neigh.i <- setdiff(neigh.i, j)
  
  if (length(effective.neigh.i) > 0){
    
    update.term.1 <- sum(sapply(effective.neigh.i,
                                function(x){f(eta.mat.t[[l]][x,i],rho[l])}))
    
    
    return (term1 - term2 + term3 + update.term.1 - update.term.2)
    
  }
  
  else{
    
    return (term1 - term2 + term3 - update.term.2)
    
  }
  
}

eta.mat.l.update <- function(n, m, l, mu, gamma, B.mat, 
                             eta.mat.t, eta.t, eta.t_1,
                             m.t, tau.t, neigh.list,
                             lambda.vec, d.vec){
  
  #----------------------------------------------------------------------------#
  # INPUT:
  # n          <- number of nodes
  # m          <- number of graphs/networks
  # l          <- index of the lth network
  # mu         <- signal-to-noise ratio of covariate information
  # gamma      <- the value of n/p
  # B.mat      <- the covariate matrix
  # eta.mat.t  <- the current eta matrix
  # eta.t      <- the current eta vector
  # eta.t_1    <- the pervious eta vector
  # m.t        <- the current value of m vector
  # tau.t      <- the current value of tau vector
  # neigh.list <- list containing list of neighbours of nodes in each
  #               of m networks
  # lambda.vec <- vector of signal to noise ratio for m networks
  # d.vec      <- vector of average degrees of m networks
  #----------------------------------------------------------------------------#
  # OUTPUT:    
  # update.l <- the updated value of eta_{. -> .;l} matrix
  #----------------------------------------------------------------------------#
  
  update.l <- matrix(0, n, n)
  
  for (i in 1:n){
    
    for (j in 1:n){
      
      if (i != j){
        
        update.l[i,j] <- eta.mat.i.j.l.update(n, m, i, j, l, mu, gamma, B.mat, 
                                              eta.mat.t, eta.t, eta.t_1,
                                              m.t, tau.t, neigh.list,
                                              lambda.vec, d.vec)
        
      }
      
    }
    
  }
  
  return (update.l)
  
}

eta.mat.update <- function(n, m, mu, gamma, B.mat, 
                           eta.mat.t, eta.t, eta.t_1,
                           m.t, tau.t, neigh.list,
                           lambda.vec, d.vec){
  
  #----------------------------------------------------------------------------#
  # INPUT:
  # n          <- number of nodes
  # m          <- number of graphs/networks
  # mu         <- signal-to-noise ratio of covariate information
  # gamma      <- the value of n/p
  # B.mat      <- the covariate matrix
  # eta.mat.t  <- the current eta matrix
  # eta.t      <- the current eta vector
  # eta.t_1    <- the pervious eta vector
  # m.t        <- the current value of m vector
  # tau.t      <- the current value of tau vector
  # neigh.list <- list containing list of neighbours of nodes in each
  #               of m networks
  # lambda.vec <- vector of signal to noise ratio for m networks
  # d.vec      <- vector of average degrees of m networks
  #----------------------------------------------------------------------------#
  # OUTPUT:    
  # update     <- list containing the updated value of eta_{. -> .;l} matrix
  #               for all 1<=l<=m networks
  #----------------------------------------------------------------------------#
  
  update <- lapply(1:m, function(l){eta.mat.l.update(n, m, l, mu, gamma, B.mat, 
                                                     eta.mat.t, eta.t, eta.t_1,
                                                     m.t, tau.t, neigh.list,
                                                     lambda.vec, d.vec)})
  
  return (update)
  
}

################################################################################

tau.update <- function(mu, gamma, B.mat, eta.t){
  
  #----------------------------------------------------------------------------#
  # INPUT:
  # mu          <- signal-to-noise ratio of covariate information
  # gamma       <- the value of n/p
  # B.mat       <- the covariate matrix
  # eta.mat.t   <- the current eta matrix
  # eta.t       <- the current eta vector
  #----------------------------------------------------------------------------#
  # OUTPUT:
  # inside.term <- the updated value of tau vector
  #----------------------------------------------------------------------------#
  
  inside.term <- (1+mu) - (mu/gamma) * ((B.mat^2) %*% (1/(cosh(eta.t))^2))
  return (inside.term)
  
}

################################################################################

m.update <- function(mu, gamma, B.mat, eta.t, m.t_1, tau.t.1){
  
  #----------------------------------------------------------------------------#
  # INPUT:
  # mu         <- signal-to-noise ratio of covariate information
  # gamma      <- the value of n/p
  # B.mat      <- the covariate matrix
  # eta.t      <- the current eta vector
  # m.t_1      <- the previous value of m vector
  # tau.t.1    <- the updated value of tau vector
  #----------------------------------------------------------------------------#
  # OUTPUT:
  # m.new      <- the updated m vector
  #----------------------------------------------------------------------------#
  
  m.new <- (1/tau.t.1) * sqrt(mu/gamma) * (B.mat %*% tanh(eta.t)) -
    (mu/(gamma*tau.t.1)) * ((B.mat^2) %*% (1/(cosh(eta.t))^2)) * m.t_1
  
  return (m.new)
}

################################################################################

update.all <- function(n, m,
                       mu, gamma, lambda.vec, d.vec,
                       B.mat, neigh.list,
                       m.t, m.t_1,
                       eta.t, eta.t_1,
                       eta.mat.t,
                       tau.t){
  
  #----------------------------------------------------------------------------#
  # INPUT:
  # n          <- number of nodes
  # m          <- number of graphs/networks
  # mu         <- signal-to-noise ratio of covariate information
  # gamma      <- the value of n/p
  # lambda.vec <- vector of signal to noise ratio for m networks
  # d.vec      <- vector of average degrees of m networks
  # B.mat      <- the covariate matrix
  # neigh.list <- list containing list of neighbours of nodes in each
  #               of m networks
  # m.t        <- the current value of m vector
  # m.t_1      <- the previous value of m vector
  # eta.t      <- the current eta vector
  # eta.t_1    <- the pervious eta vector
  # eta.mat.t  <- the current eta matrix
  # tau.t      <- the current value of tau vector
  #----------------------------------------------------------------------------#
  # OUTPUT:    <- a list containing the updated and current values of eta, m
  #               vector, and the updated values of tau and eta matrix
  #----------------------------------------------------------------------------#
  
  eta.mat.new <- eta.mat.update(n, m, mu, gamma, B.mat, 
                                eta.mat.t, eta.t, eta.t_1,
                                m.t, tau.t, neigh.list,
                                lambda.vec, d.vec)
  
  eta.new <- eta.update(n, m, mu, gamma, B.mat, 
                        eta.mat.t, eta.t, eta.t_1,
                        m.t, tau.t, neigh.list,
                        lambda.vec, d.vec)
  
  eta.old <- eta.t
  
  tau.new <- tau.update(mu, gamma, B.mat, eta.t)
  m.new <- m.update(mu, gamma, B.mat, eta.t, m.t_1, tau.new)
  m.old <- m.t
  
  return (list("newEta" = eta.new, "oldEta" = eta.old,
               "newtau" = tau.new, "newm" = m.new, "oldm" = m.old,
               "newEtamat" = eta.mat.new))
}

################################################################################

power.comp <- function(n, m, p, gamma, sigma, mu, lambda.vec, d.vec,
                       n.iter = 50){
  
  #----------------------------------------------------------------------------#
  # INPUT:
  # n          <- number of nodes
  # m          <- number of graphs/networks
  # p          <- dimension of covariates
  # gamma      <- the value of n/p
  # sigma      <- community assignment vector
  # mu         <- signal-to-noise ratio of covariate information
  # lambda.vec <- vector of signal to noise ratio for m networks
  # d.vec      <- vector of average degrees of m networks
  # n.iter     <- the maximum iterations of BP
  #----------------------------------------------------------------------------#
  # OUTPUT:    <- a vector containing accept/reject decision,
  #               overlap and norm of eta^{T_max}
  #----------------------------------------------------------------------------#
  
  b.vec <- d.vec - lambda.vec * sqrt(d.vec)
  a.vec <- d.vec + lambda.vec * sqrt(d.vec)
  
  B.mat <- generate.B.mat(n, p, mu, sigma)
  
  Adjacency.list <- generate.Adjacency.list(n, m, a.vec, b.vec, sigma)
  neigh.list <- Neighbourhood.list(n, m, Adjacency.list)
  
  eta.now <- rnorm(n, mean = 0, sd = 0.1)
  eta.now_1 <- rnorm(n, mean = 0, sd = 0.1)
  
  eta.0 <- eta.now
  
  m.now <- rnorm(p, mean = 0, sd = 0.1)
  m.now_1 <- rnorm(p, mean = 0, sd = 0.1)
  
  tau.now <- tau.update(mu, gamma, B.mat, eta.now_1)
  
  print (norm(as.matrix(eta.now), type = "F"))
  
  eta.mat.now <- lapply(Adjacency.list,
                        function(x){0 * x})
  
  for (t in 1:n.iter){
    
    updated.list <- update.all(n, m,
                               mu, gamma, lambda.vec, d.vec,
                               B.mat, neigh.list,
                               m.now, m.now_1,
                               eta.now, eta.now_1,
                               eta.mat.now,
                               tau.now)
    eta.now_1 <- eta.now
    eta.now <- updated.list$newEta
    eta.mat.now <- updated.list$newEtamat
    m.now_1 <- m.now; m.now <- updated.list$newm
    tau.now <- updated.list$newtau
    
    cat(paste("End Iteration", t, "\n"), file = "log.txt", append = T)
  }
  
  reject <- (norm(as.matrix(eta.now), type = "F") > norm(as.matrix(eta.0),
                                                         type = "F"))
  
  sigma.hat <- sign(eta.now)
  overlap <- max(abs(sum(sigma.hat * sigma)), abs(sum(-sigma.hat * sigma)))/n
  
  return (c(reject, overlap, norm(as.matrix(eta.now), type = "F")))
}
