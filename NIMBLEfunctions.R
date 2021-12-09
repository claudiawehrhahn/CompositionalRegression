

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Functions:

# univariate spike an slab prior:
dspikeslab <- nimbleFunction(
  run=function(x = double(0), 
               gamma = integer(0),
               var1 = double(0),
               var2 = double(0),
               log = integer(0, default=0))
  {
    returnType(double(0))
    
    ldens <- log((1-gamma) * dnorm(x, 0, sd = sqrt(var1)) + 
                   gamma *  dnorm(x, 0, sd = sqrt(var2)))
    
    if(log) return(ldens)
    else return(exp(ldens))
  }
)
cspikeslab <- compileNimble(dspikeslab)

rspikeslab <- nimbleFunction(
  run=function(n = integer(0), 
               gamma = integer(0),
               var1 = double(0),
               var2 = double(0))
  {
    returnType(double(0))
    
    if(n != 1) {
      stop("rspikeslab: number of simulated values is not 1 (n != 1).")
    }
    
    if(gamma == 1) {
      x <- rnorm(1, mean = 0, sd = sqrt(var2))
    }
    if(gamma == 0) {
      x <- rnorm(1, mean = 0, sd = sqrt(var1))
    }
    
    return(x)
  }
)
crspikeslab <- compileNimble(rspikeslab)

# mixture of Dirichelt distributions
dDMBPP <- nimbleFunction(
  run=function(x = double(1), 
               m = integer(0),
               w = double(1),
               k = double(0),
               theta1 = double(1),
               theta2 = double(1),
               log = integer(0, default=0))
  {
    returnType(double(0))
    
    J <- length(w)
    
    dirichletDensities <- nimNumeric(length = J, value = 0)
    for( j in 1:J) {
      dirichletDensities[j] <- ddirch(x, alpha = c(ceiling(k * theta1[j]),
                                                   ceiling(k * theta2[j]),
                                                   k + m - ceiling(k * theta1[j]) - ceiling(k * theta2[j])))
    }
    
    dens <- sum( w * dirichletDensities ) 
    
    if(log) return(log(dens))
    else return(dens) 
  }
)
cdDMBPP <- compileNimble(dDMBPP)

rDMBPP <- nimbleFunction(
  run=function(n = integer(0), 
               m = integer(0),
               w = double(1),
               k = double(0),
               theta1 = double(1),
               theta2 = double(1))
  {
    returnType(double(1))
    
    index <- rcat(1, w)
    alpha0 <- c(ceiling(k * theta1[index]),
                ceiling(k * theta2[index]),
                k + m - ceiling(k * theta1[index]) - ceiling(k * theta2[index]))# nimNumeric(length = m + 1, value = 1)
    out <- rdirch(1, alpha = alpha0)
    
    return(out) 
  }
)
crDMBPP <- compileNimble(rDMBPP)

