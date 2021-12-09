
rm(list = ls())

load(file="./Data/indexes.RData")

dyn.load("depSymplexCodeMCMC.so")


funMCMC <- function(i) {
    
    E <- indexes[[i]][1] 
    prior <- indexes[[i]][2] 
    nrec <- indexes[[i]][3] 
    rep <- indexes[[i]][4] 
    
    Y <- scan(paste("./Data/E", E, "nrec", nrec, "rep", rep, "response.txt", sep=""))
    X <- scan(paste("./Data/E", E, "nrec", nrec, "rep", rep, "covariate.txt", sep=""))
    tun <- scan(paste("./Data/E", E, "nrec", nrec, "rep", rep, "tunning.txt", sep=""))

    set.seed(E + prior + nrec + rep);
    result <- .C("depSymplexMCMC",      
                 as.integer(E),
                 as.integer(prior),
                 as.integer(nrec),
                 as.integer(rep),
                 as.double(Y),
                 as.double(X),
                 as.double(tun))
} 

funMCMC(1)


