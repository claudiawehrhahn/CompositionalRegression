
rm(list=ls())
library(nimble)

source("./NIMBLEfunctions.R")
load(file="./Data/indexes.RData")



#------------------------------------------------------------------------------
# loading data:

i <- 1 # one of the simulated data sets.
E <- indexes[[i]][1] 
prior <- indexes[[i]][2] 
nrec <- indexes[[i]][3] 
rep <- indexes[[i]][4] 

d <- 2 # dimension of simplex for simulated data set
pCov <- 2 # number of coefficients in linear predictor

Y <- scan(paste("./Data/E", E, "nrec", nrec, "rep", rep, "response.txt", sep=""))
Yaux <- matrix(Y, ncol = d, nrow = nrec, byrow = TRUE)
Ymatrix <- cbind(Yaux, 1 - rowSums(Yaux)) # in nimble argument of Dirichlet distribution need to add up to one.
X <- scan(paste("./Data/E", E, "nrec", nrec, "rep", rep, "covariate.txt", sep=""))
Xmatrix <- matrix(X, ncol = p, nrow = nrec, byrow = TRUE)
tun <- scan(paste("./Data/E", E, "nrec", nrec, "rep", rep, "tunning.txt", sep="")) # tunnig parameter for g prior

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# NIMBLE code:

# code for a m=2 and one predictor
code <- nimbleCode( {
  for(i in 1:nrec) {
    y[i, 1:(d+1)] ~ dDMBPP(d, w[i, 1:J], k, theta1[i, 1:J], theta2[i, 1:J])
    
    for(j in 1:(J-1)) {
      v[i, j] <- expit(beta1V[j] + beta2V[j] * X[i, 1])
    }
    w[i, 1:J] <- stick_breaking(v[i, 1:(J-1)])
    
    for(j in 1:J) {
      theta1Aux[i, j] <- exp(beta1T1[j] + beta2T1[j] * X[i, 1])
      theta2Aux[i, j] <- exp(beta1T2[j] + beta2T2[j] * X[i, 1])
      theta1[i, j] <- theta1Aux[i, j] / (1 + theta1Aux[i, j] + theta2Aux[i, j])
      theta2[i, j] <- theta2Aux[i, j] / (1 + theta1Aux[i, j] + theta2Aux[i, j])
    }
  }
  
  k1 ~ dpois(25)
  k <- k1 + 1
  
  for(j in 1:(J-1)) {
    beta1V[j] ~ dnorm(0, var = 100)
    beta2V[j] ~ dspikeslab(gammaV, tau1_V * xtx1[1:p, 1:p], tau2_V * xtx1[1:p, 1:p])
  }
  
  for(j in 1:J) {
    beta1T1[j] ~ dnorm(0, var = 100)
    beta1T2[j] ~ dnorm(0, var = 100)
    beta2T1[j] ~ dspikeslab(gammaT, tau1_T * xtx1[1:p, 1:p], tau2_T * xtx1[1:p, 1:p])
    beta2T2[j] ~ dspikeslab(gammaT, tau1_T * xtx1[1:p, 1:p], tau2_T * xtx1[1:p, 1:p])
  }
  
  gammaV ~ dbern(pV)
  gammaT ~ dbern(pT)
})

# initializacion, model building and mcmc samples
set.seed(1)
consts <- list(J = 25, p = pCov - 1, d = 2, nrec = nrow(Ymatrix),
               X = matrix(Xmatrix[, 2], ncol= 1),
               pV = 0.5, pT = 0.5,
               tau1_V = tun[1], tau2_V = tun[2],
               tau1_T = tun[3], tau2_T = tun[4],
               xtx1 = solve(t(Xmatrix[, 2])%*%Xmatrix[, 2]))
inits <- list(k1 = 24, gammaV = 0, gammaT = 0,
              beta1V = rep(0, consts$J - 1),
              beta2V = rep(0, consts$J - 1),
              beta1T1 = rep(0, consts$J),
              beta1T2 = rep(0, consts$J),
              beta2T1 = rep(0, consts$J),
              beta2T2 = rep(0, consts$J))
data = list(y = Ymatrix)
m <- nimbleModel(code, 
                 data = data,
                 inits = inits, 
                 constants = consts)
cm <- compileNimble(m)
mConf <- configureMCMC(m, monitors = list( 'k', 'gammaV', 'gammaT', 'beta1V', 'beta2V',
                                           'beta1T1', 'beta1T2', 'beta2T1', 'beta2T2'))
mcmc <- buildMCMC(mConf)
cmcmc <- compileNimble(mcmc, project = m)
samples <- runMCMC(cmcmc, niter = 1100, nburnin = 1000, thin=1, setSeed = TRUE) 

save(samples, file=paste("./Output/E", E, "/prior", prior, "/", nrec, "/rep", rep,"mcmcDMBPP_nimble.RData", sep=""))

