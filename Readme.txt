
This Readme file describes the steps to run the Dependent multivariate Bernstein polynomial process (DMBPP) for compositional regresion. This model is proposed and described in the paper "Dependent Bayesian nonparametric modeling of compositional data using random Bernstein polynomials" by Claudia Wehrhahn, Andrés F. Barrientos and Alejandro Jara.
 
General description
- 1. and 2. describe how to run the DMBPP model using the C code for in the simulation study in the paper. In this simulation the number of beta paramters pCov = 2 and simplex simension d = 2.
- 3. described how to run the DMBPP model using the nimble R package. Its usage is illustrated for one of the simulated data sets, but is simple to adapt for other data sets.


1. Do R CMD SHLIB depSymplexCodeMCMC.c from console while in the directory where the depSymplexCodeMCMC.c file is. The files depSymplexCodeMCMC.o and depSymplexCodeMCMC.so should be created.

2. Run wrapperDepSymplexMCMC.R for one of the simulated data sets.

3. Run the NIMBLEDepSymplex.R file for one of the simulated data sets. 




