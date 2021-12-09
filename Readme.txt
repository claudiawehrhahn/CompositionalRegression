

Steps 1. and 2. describe how to run the DMBPP model using the C code used in the simulation study. To see how to use the NIMBLE R package to run the model see the Notes.



1. Do R CMD SHLIB depSymplexCodeMCMC.c from console while in the directory where the depSymplexCodeMCMC.c file is. The files depSymplexCodeMCMC.o and depSymplexCodeMCMC.so should be created.

2. Open wrapperDepSymplexMCMC.R and run the code for one of the simulated data sets.


NOTES:

a. The depSymplexCodeMCMC.c file works for number of beta paramters pCov = 2 and simplex simension d = 2.

b. The files NIMBLEDepSymplex.R and NIMBLEfunctions.R allow to run the proposed compositional regression model in NIMBLE. Its usage is illustrated for one of the simulated data sets, but is simple to use for other data sets.




