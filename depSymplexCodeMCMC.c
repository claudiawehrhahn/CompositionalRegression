//
//
//  
//
//  Created by Claudia Wehrhahn on 3/18/20.
//
//
// difference with depSymplexCode.c is where we compute post processing: here within the mcmc
// difference with depSymplexCode2.c is that now data are arguments of the function

#define d 2
#define pCov 2
#define kmax 150
#define nSB 25
#define lamK 25

#define nburn 10000  // length of burn-in in MCMC
#define njump 10      // size of thinning in MCMC
#define nsave 10000   // number of saved iterations in MCMC

#include <R.h> // to be called from R
#include <Rmath.h> // to use distributions
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "./functions.c"


void depSymplexMCMC(int *E, int *prior, int *nrec, int *rep, double *dataY, double *dataX, double *dataTunning) { // arguments in the  function need to be used with * before
    
    char buf[200];
    double tmp, tmp2, alpha[2], lpost, lpostStar;
    double yslice, L[pCov], R[pCov], wslice = 2.0;
    int i, j, l, i1, j1, jb, l1;
    int acceptK = 0;
    
    
    GetRNGstate();
    
    //----------------------------------------------------
   
    
    double tprob0, lProbMod[4];
    if(*prior == 1) {
        tprob0 = 2;
    }
    if(*prior == 2) {
        tprob0 = 10;
    }
    lProbMod[0] = log(tprob0-1) - log(tprob0);
    lProbMod[1] = log(tprob0-1) - log(2*tprob0*tprob0);
    lProbMod[2] = log(tprob0-1) - log(2*tprob0*tprob0);
    lProbMod[3] = - log(tprob0*tprob0);
    
    
    // creating file

    snprintf(buf, 100, "./Output/E%d/prior%d/%d/rep%d_k.txt" ,*E, *prior, *nrec, *rep);
    FILE *file_k = fopen(buf, "w");
    snprintf(buf, 100, "./Output/E%d/prior%d/%d/rep%d_gammaV.txt" ,*E, *prior, *nrec, *rep);
    FILE *file_gammaV = fopen(buf, "w");
    snprintf(buf, 100, "./Output/E%d/prior%d/%d/rep%d_gammaT.txt" ,*E, *prior, *nrec, *rep);
    FILE *file_gammaT = fopen(buf, "w");
    snprintf(buf, 100, "./Output/E%d/prior%d/%d/rep%d_betaV.txt" ,*E, *prior, *nrec, *rep);
    FILE *file_betaV = fopen(buf, "w");
    snprintf(buf, 100, "./Output/E%d/prior%d/%d/rep%d_betaT.txt" ,*E, *prior, *nrec, *rep);
    FILE *file_betaT = fopen(buf, "w");
    snprintf(buf, 100, "./Output/E%d/prior%d/%d/rep%d_llik.txt" ,*E, *prior, *nrec, *rep);
    FILE *file_llik = fopen(buf, "w");

    
    
    // loading data
    /*  load response, convariates and tunning params */
    double y[*nrec][d], x[*nrec][pCov];
    
    // response
    i1 = 0;
    for( i = 0; i < *nrec; i++ ) {
        for( j = 0; j < d; j++ ) {
            y[i][j] = dataY[i1];
            i1 += 1;
        }
    }
    
    
    // covariate
    i1 = 0;
    for( i = 0; i < *nrec; i++ ) {
        for( j = 0; j < pCov; j++ ) {
            x[i][j] = dataX[i1];
            i1 += 1;
        }
    }
    
    
    // scaling constants of g-prior
    // scaling constants of g-prior
    double gCntV[pCov], gCntT[pCov];
    gCntV[0] = dataTunning[0];
    gCntV[1] = dataTunning[1];
    gCntT[0] = dataTunning[2];
    gCntT[1] = dataTunning[3];
    
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    /* precomputing dirichlet densities of data points */
    // NEEDS TO BE CHANGED WHEN d > 2 (d: dimension of the symplex)
    double ****dirDenData = (double ****) malloc(*nrec * sizeof(double ***));
    for(i = 0; i < *nrec; i++) {
        dirDenData[i] = (double ***) malloc(kmax * sizeof(double **));
        for(l = 0; l < kmax; l++) { // degrees
            dirDenData[i][l] = (double **) malloc((l+1) * sizeof(double *));
            for(i1 = 0; i1 < (l+1); i1++) {
                dirDenData[i][l][i1] = (double *) malloc((l+1) *sizeof(double));
                for(j1 = 0; j1 < (l+1); j1++) {
                    if((i1 + j1) < (l+1+d-2)) {
                        alpha[0] = (double) i1 + 1;
                        alpha[1] = (double) j1 + 1;
                        dirDenData[i][l][i1][j1] = dirich_den(y[i], alpha, (int) l+1+d);
                    } else {
                        dirDenData[i][l][i1][j1] = 0.0;
                    }
                }
            }
        }
    }
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // computing (X^t)%*%X without intercept
    // xtx = sum x_i^2. THIS COMPUTATION WORKS ONLY FOR pCov=2 (pCov: number of beta parameters).
    double xtx, xtx1;
    xtx = 0.0;
    for(i = 0; i < *nrec; i++) {
        xtx += pow(x[i][pCov-1], 2);
    }
    xtx1 = 1 / xtx;
    
    
    // variance for betaV and betaT for model selection:
    double sdModSelV[2], sdModSelT[2];
    for(i = 0; i <2; i++) {
        sdModSelV[i] = sqrt(gCntV[i] * xtx1);
        sdModSelT[i] = sqrt(gCntT[i] * xtx1);
    }
    
    // var-cov matrix for betaV and betaT
    double varBetaV[2][pCov][pCov], varBetaT[2][pCov][pCov];
    for(i = 0; i < 2; i++) {
        for(j = 0; j < pCov; j++) {
            for(jb = 0; jb < pCov; jb++) {
                varBetaV[i][j][jb] = 0.0;
                varBetaT[i][j][jb] = 0.0;
            }
        }
    }
    for(i = 0; i < 2; i++) {
        varBetaV[i][0][0] = 100.0;
        varBetaT[i][0][0] = 100.0;
    }
    for(i = 0; i < 2; i++) {
        varBetaV[i][1][1] = pow(sdModSelV[i], 2);
        varBetaT[i][1][1] = pow(sdModSelT[i], 2);
    }
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    /* declaration of objects representing parameters */
    int k, k1;
    double betaV[nSB-1][pCov], v[*nrec][nSB], w[*nrec][nSB];
    double betaT[nSB][d][pCov], theta[*nrec][nSB][d];
    double betaV2[nSB-1][pCov], v2[*nrec][nSB], w2[*nrec][nSB];
    double betaT2[nSB][d][pCov], theta2[*nrec][nSB][d];
    
    /* model selection variables and probabilities */
    int gammaV = 0, gammaT = 0;
    double pGamma[4];
    
    
    /* initializingparameters */
    k = 25;
    k1 = k;
    
    for( l = 0; l < nSB - 1; l++ ) {
        for(jb = 0; jb < pCov; jb++ ) {
            betaV[l][jb] = (double) 0;
            betaV2[l][jb] = betaV[l][jb];
        }
        for( j = 0; j < d; j++ ) {
            for( jb = 0; jb < pCov; jb++ ) {
                betaT[l][j][jb] = (double) 0;
                betaT2[l][j][jb] = betaT[l][j][jb];
            }
        }
    }
    for( j = 0; j < d; j++ ) {
        for( jb = 0; jb < pCov; jb++ ) {
            betaT[nSB - 1][j][jb] = (double) 0;
            betaT2[nSB - 1][j][jb] = betaT[nSB - 1][j][jb];
        }
    }
    
    
    /* initializing weights and atoms */
    for( i = 0; i < *nrec; i++ ) {
        /* v and w */
        tmp = 0.0;
        for( jb = 0; jb < pCov; jb++ ) {
            tmp += x[i][jb] * betaV[0][jb];
        }
        v[i][0] = exp(tmp) / (1 + exp(tmp));
        v2[i][0] = v[i][0];
        w[i][0] = v[i][0];
        w2[i][0] = w[i][0];
        tmp2 = 1 - v[i][0];
        
        for( l = 1; l < nSB-1; l++ ) {
            tmp = 0.0;
            for( jb = 0; jb < pCov; jb++ ) {
                tmp += x[i][jb] * betaV[l][jb];
            }
            v[i][l] = exp(tmp) / (1 + exp(tmp));
            v2[i][l] = v[i][l];
            w[i][l] = v[i][l] * tmp2;
            w2[i][l] = w[i][l];
            tmp2 *= 1- v[i][l];
            
        }
        
        v[i][nSB-1] = 1.0;
        v2[i][nSB-1] = v[i][nSB-1];
        w[i][nSB-1] = v[i][nSB-1] * tmp2;
        w2[i][nSB-1] = w[i][nSB-1];
        
        
        /* atoms */
        for( l = 0; l < nSB; l++ ) {
            tmp2 = 1.0;
            for(j = 0; j < d; j++ ) {
                tmp = 0.0;
                for( jb = 0; jb < pCov; jb++ ) {
                    tmp += x[i][jb] * betaT[l][j][jb];
                }
                theta[i][l][j] = exp(tmp);
                tmp2 += theta[i][l][j];
            }
            for(j = 0; j < d; j++ ) {
                theta[i][l][j] /= tmp2;
                theta2[i][l][j] = theta[i][l][j];
            }
        }
        
    }
    
    // mcmc specification
    int nmcmc = nburn + njump * nsave;
    int itersave, iiter;
    int isave = 0;
    itersave = nburn + (isave) * njump;
    
    
    /* MCMC  */
    
    
    for( iiter = 0; iiter < nmcmc; iiter++ ){

        /* updating betaV */
        
        for( l = 0; l < nSB - 1; l++ ) {//
            
            /* computing likelihoods */
            lpost = 0.0;
            for( i = 0; i < *nrec; i++ ) {
                tmp = 0.0;
                for( l1 = 0; l1 < nSB; l1++ ) {
                    for(j = 0; j < d; j++) {
                        alpha[j] = (int) ceil(k * theta[i][l1][j]) - 1;
                    }
                    tmp += w[i][l1] * dirDenData[i][k-1][(int) alpha[0]][(int) alpha[1]];;
                }
                lpost += log(tmp);
            }
            tmp = ((double)(1-gammaV)) * bivariate_norm_pdf(betaV[l], varBetaV[0]);
            tmp += ((double)gammaV) * bivariate_norm_pdf(betaV[l], varBetaV[1]);
            lpost += log(tmp);

            
            /* Step a) */
            yslice = lpost - rexp((double) 1.0);
            
            
            /* Step b) */
            for(jb = 0; jb < pCov; jb++ ) {
                L[jb] = betaV[l][jb] - wslice * runif((double) 0, (double) 1);
                R[jb] = L[jb] + wslice;
            }
            
            
            /* Step c) */
            for(jb = 0; jb < pCov; jb++ ) {
                betaV2[l][jb] =  L[jb] + runif((double) 0, (double) 1) * (R[jb]-L[jb]);
            }
            
            /* computing updated likelihood */
            for( i = 0; i < *nrec; i++ ) {//
                tmp = 0.0;
                for(jb = 0; jb < pCov; jb++ ) {
                    tmp += x[i][jb] * betaV2[0][jb];
                }
                v2[i][0] = exp(tmp) / (1 + exp(tmp));
                w2[i][0] = v2[i][0];
                tmp2 = 1 - v2[i][0];
                for( l1 = 1; l1 < nSB - 1; l1++) {
                    tmp = 0.0;
                    for(jb = 0; jb < pCov; jb++ ) {
                        tmp += x[i][jb] * betaV2[l1][jb];
                    }
                    v2[i][l1] = exp(tmp) / (1 + exp(tmp));
                    w2[i][l1] = v2[i][l1] * tmp2;
                    tmp2 *= 1 - v2[i][l1];
                }
                w2[i][nSB-1] = v2[i][nSB-1] * tmp2;
            }
            
            lpostStar = 0.0;
            for( i = 0; i < *nrec; i++ ) {
                tmp = 0.0;
                for( l1 = 0; l1 < nSB; l1++ ) {
                    for(j = 0; j < d; j++) {
                        alpha[j] = (int) ceil(k * theta[i][l1][j]) - 1;
                    }
                    tmp += w2[i][l1] * dirDenData[i][k-1][(int) alpha[0]][(int) alpha[1]];
                }
                lpostStar += log(tmp);
            }
            tmp = (1-gammaV) * bivariate_norm_pdf(betaV2[l], varBetaV[0]);
            tmp += gammaV * bivariate_norm_pdf(betaV2[l], varBetaV[1]);
            lpostStar += log(tmp);
            
            
            while( yslice > lpostStar ) {
                for(jb = 0; jb < pCov; jb++ ) {
                    if( betaV2[l][jb] < betaV[l][jb] ) {
                        L[jb] = betaV2[l][jb];
                    } else  {
                        R[jb] = betaV2[l][jb];
                    }
                }
                
                for(jb = 0; jb < pCov; jb++ ) {
                    betaV2[l][jb] =  L[jb] + runif((double) 0, (double) 1) * (R[jb]-L[jb]);
                }
                
                
                /* computing updated likelihood */
                for( i = 0; i < *nrec; i++ ) {//
                    tmp = 0.0;
                    for(jb = 0; jb < pCov; jb++ ) {
                        tmp += x[i][jb] * betaV2[0][jb];
                    }
                    v2[i][0] = exp(tmp) / (1 + exp(tmp));
                    w2[i][0] = v2[i][0];
                    tmp2 = 1 - v2[i][0];
                    for( l1 = 1; l1 < nSB - 1; l1++) {
                        tmp = 0.0;
                        for(jb = 0; jb < pCov; jb++ ) {
                            tmp += x[i][jb] * betaV2[l1][jb];
                        }
                        v2[i][l1] = exp(tmp) / (1 + exp(tmp));
                        w2[i][l1] = v2[i][l1] * tmp2;
                        tmp2 *= 1 - v2[i][l1];
                    }
                    w2[i][nSB-1] = v2[i][nSB-1] * tmp2;
                }
                
                lpostStar = 0.0;
                for( i = 0; i < *nrec; i++ ) {
                    tmp = 0.0;
                    for( l1 = 0; l1 < nSB; l1++ ) {
                        for(j = 0; j < d; j++) {
                            alpha[j] = (int) ceil(k * theta[i][l1][j]) - 1;
                        }
                        tmp += w2[i][l1] * dirDenData[i][k-1][(int) alpha[0]][(int) alpha[1]];
                    }
                    lpostStar += log(tmp);
                }
                tmp = (1-gammaV) * bivariate_norm_pdf(betaV2[l], varBetaV[0]);
                tmp += gammaV * bivariate_norm_pdf(betaV2[l], varBetaV[1]);
                lpostStar += log(tmp);
            }
            
            for(jb = 0; jb < pCov; jb++ ) {//
                betaV[l][jb] = betaV2[l][jb];
            }
            for( i = 0; i < *nrec; i++ ) {
                v[i][0] = v2[i][0];
                w[i][0] = w2[i][0];
                
                for( l1 = 1; l1 < nSB -1; l1++) {
                    v[i][l1] = v2[i][l1];
                    w[i][l1] = w2[i][l1];
                }
                w[i][nSB-1] = w2[i][nSB-1];
            }
        }
        
        
        /* updating betaT */
        
        for( l = 0; l < nSB; l++ ) {
            for( j = 0; j < d; j++ ) {
                /* computing likelihoods */
                lpost = 0.0;
                for( i = 0; i < *nrec; i++ ) {
                    tmp = 0.0;
                    for( l1 = 0; l1 < nSB; l1++ ) {
                        for(j1 = 0; j1 < d; j1++) {
                            alpha[j1] = (int) ceil(k * theta[i][l1][j1]) - 1;
                        }
                        tmp += w[i][l1] * dirDenData[i][k-1][(int) alpha[0]][(int) alpha[1]];
                    }
                    lpost += log(tmp);
                }
                tmp = (1-gammaT) * bivariate_norm_pdf(betaT[l][j], varBetaT[0]);
                tmp += gammaT * bivariate_norm_pdf(betaT[l][j], varBetaT[1]);
                lpost += log(tmp);
                
                
                /* Step a) */
                yslice = lpost - rexp((double) 1.0);
                
                
                /* Step b) */
                for(jb = 0; jb < pCov; jb++ ) {
                    L[jb] = betaT[l][j][jb] - wslice * runif((double) 0, (double) 1);
                    R[jb] = L[jb] + wslice;
                }
                
                
                /* Step c) */
                for(jb = 0; jb < pCov; jb++ ) {
                    betaT2[l][j][jb] = L[jb] + runif((double) 0, (double) 1) * (R[jb]-L[jb]);
                }
                
                /* computing updated likelihood */
                /* updating the atoms */
                for( i = 0; i < *nrec; i++ ) {
                    tmp2 = 1.0;
                    for(j1 = 0; j1 < d; j1++ ) {
                        tmp = 0.0;
                        for( jb = 0; jb < pCov; jb++ ) {
                            tmp += x[i][jb] * betaT2[l][j1][jb];
                        }
                        theta2[i][l][j1] = exp(tmp);
                        tmp2 += theta2[i][l][j1];
                    }
                    for(j1 = 0; j1 < d; j1++ ) {
                        theta2[i][l][j1] /= tmp2;
                    }
                }
                
                /* computing likelihoods */
                lpostStar = 0.0;
                for( i = 0; i < *nrec; i++ ) {
                    tmp = 0.0;
                    for( l1 = 0; l1 < nSB; l1++ ) {
                        for(j1 = 0; j1 < d; j1++) {
                            alpha[j1] = (int) ceil(k * theta2[i][l1][j1]) - 1;
                        }
                        tmp += w[i][l1] * dirDenData[i][k-1][(int) alpha[0]][(int) alpha[1]];
                    }
                    lpostStar += log(tmp);
                }
                tmp = (1-gammaT) * bivariate_norm_pdf(betaT2[l][j], varBetaT[0]);
                tmp += gammaT * bivariate_norm_pdf(betaT2[l][j], varBetaT[1]);
                lpostStar += log(tmp);
                
                
                while( yslice > lpostStar ) {
                    for(jb = 0; jb < pCov; jb++ ) {
                        if( betaT2[l][j][jb] < betaT[l][j][jb] ) {
                            L[jb] = betaT2[l][j][jb];
                        } else  {
                            R[jb] = betaT2[l][j][jb];
                        }
                    }
                    
                    for(jb = 0; jb < pCov; jb++ ) {
                        betaT2[l][j][jb] = L[jb] + runif((double) 0, (double) 1) * (R[jb]-L[jb]);
                    }
                    
                    
                    /* computing updated likelihood */
                    /* updating the atoms */
                    for( i = 0; i < *nrec; i++ ) {
                        tmp2 = 1.0;
                        for(j1 = 0; j1 < d; j1++ ) {
                            tmp = 0.0;
                            for( jb = 0; jb < pCov; jb++ ) {
                                tmp += x[i][jb] * betaT2[l][j1][jb];
                            }
                            theta2[i][l][j1] = exp(tmp);
                            tmp2 += theta2[i][l][j1];
                        }
                        for(j1 = 0; j1 < d; j1++ ) {
                            theta2[i][l][j1] /= tmp2;
                        }
                    }
                    
                    /* computing likelihoods */
                    lpostStar = 0.0;
                    for( i = 0; i < *nrec; i++ ) {
                        tmp = 0.0;
                        for( l1 = 0; l1 < nSB; l1++ ) {
                            for(j1 = 0; j1 < d; j1++) {
                                alpha[j1] = (int) ceil(k * theta2[i][l1][j1]) - 1;
                            }
                            tmp += w[i][l1] * dirDenData[i][k-1][(int) alpha[0]][(int) alpha[1]];
                        }
                        lpostStar += log(tmp);
                    }
                    tmp = (1-gammaT) * bivariate_norm_pdf(betaT2[l][j], varBetaT[0]);
                    tmp += gammaT * bivariate_norm_pdf(betaT2[l][j], varBetaT[1]);
                    lpostStar += log(tmp);
                    
                }
                
                for(jb = 0; jb < pCov; jb++ ) {
                    betaT[l][j][jb] = betaT2[l][j][jb];
                }
                for( i = 0; i < *nrec; i++ ) {
                    for(j1 = 0; j1 < d; j1++ ) {
                        theta[i][l][j1] = theta2[i][l][j1];
                    }
                }
            }
        }
        
        
        /* updating the degree */
        k1 = candK(k);
        
        lpost = 0.0;
        lpostStar = 0.0;
        for( i = 0; i < *nrec; i++ ) {
            tmp = 0.0;
            tmp2 = 0.0;
            for( l = 0; l < nSB; l++ ) {
                for(j = 0; j < d; j++) {
                    alpha[j] = (int) ceil(k * theta[i][l][j]) - 1;
                }
                tmp += w[i][l] * dirDenData[i][k-1][(int) alpha[0]][(int) alpha[1]];
                for(j = 0; j < d; j++) {
                    alpha[j] = (int) ceil(k1 * theta[i][l][j]) - 1;
                }
                tmp2 += w[i][l] * dirDenData[i][k1-1][(int) alpha[0]][(int) alpha[1]];
            }
            lpost += log(tmp);
            lpostStar += log(tmp2);
        }
        lpost += pois_lden(k, lamK);
        lpostStar += pois_lden(k1, lamK);
        
        
        if( log(runif((double) 0, (double) 1)) <= lpostStar + denQk(k1) - lpost - denQk(k) ) {
            k = k1;
            acceptK += 1;
        } else {
            k1 = k;
        }
        
        
        // computing model selection probs:
       
        i = 0; // probability pGamma[0] : (gamma^V, gamma^T) = (0,0)
        pGamma[i] = 0.0;
        for( l = 0; l < nSB - 1; l++ ) {//
            pGamma[i] += dnorm((double) betaV[l][1], (double) 0.0, (double) sdModSelV[0], (int) 1);
        }
        for( l = 0; l < nSB; l++ ) {
            for( j = 0; j < d; j++ ) {
                pGamma[i] += dnorm((double) betaT[l][j][1], (double) 0.0, (double) sdModSelT[0], (int) 1);
            }
        }
        pGamma[i] += lProbMod[i];
        
        
        i = 1; // probability pGamma[1] : (gamma^V, gamma^T) = (1,0)
        pGamma[i] = 0.0;
        for( l = 0; l < nSB - 1; l++ ) {//
            pGamma[i] += dnorm((double) betaV[l][1], (double) 0.0, (double) sdModSelV[1], (int) 1);
        }
        for( l = 0; l < nSB; l++ ) {
            for( j = 0; j < d; j++ ) {
                pGamma[i] += dnorm((double) betaT[l][j][1], (double) 0.0, (double) sdModSelT[0], (int) 1);
            }
        }
        pGamma[i] += lProbMod[i];
        
        if(pGamma[0] > pGamma[1]) {
            tmp = pGamma[0];
        } else {
            tmp = pGamma[1];
        }
        
        
        i = 2; // probability pGamma[2] : (gamma^V, gamma^T) = (0,1)
        pGamma[i] = 0.0;
        for( l = 0; l < nSB - 1; l++ ) {//
            pGamma[i] += dnorm((double) betaV[l][1], (double) 0.0, (double) sdModSelV[0], (int) 1);
        }
        for( l = 0; l < nSB; l++ ) {
            for( j = 0; j < d; j++ ) {
                pGamma[i] += dnorm((double) betaT[l][j][1], (double) 0.0, (double) sdModSelT[1], (int) 1);
            }
        }
        pGamma[i] += lProbMod[i];
        
        if(pGamma[2] > tmp) {
            tmp  = pGamma[2];
        }
        
        i = 3; // probability pGamma[3] : (gamma^V, gamma^T) = (1,1)
        pGamma[i] = 0.0;
        for( l = 0; l < nSB - 1; l++ ) {//
            pGamma[i] += dnorm((double) betaV[l][1], (double) 0.0, (double) sdModSelV[1], (int) 1);
        }
        for( l = 0; l < nSB; l++ ) {
            for( j = 0; j < d; j++ ) {
                pGamma[i] += dnorm((double) betaT[l][j][1], (double) 0.0, (double) sdModSelT[1], (int) 1);
            }
        }
        pGamma[i] += lProbMod[i];
    
        if(pGamma[3] > tmp) {
            tmp  = pGamma[3];
        }
        
        tmp2 = 0.0;
        for(l = 0; l < 4; l++) {
            pGamma[l] = exp(pGamma[l] - tmp);
            tmp2 += pGamma[l];
        }
        for(l = 0; l < 4; l++) {
            pGamma[l] /= tmp2;
        }
        
        // sampling model:
        tmp2 = runif((double) 0, (double) 1);
        l = 0;
        tmp = pGamma[l];
        while(tmp2 >= tmp){
            l += 1;
            tmp += pGamma[l];
        }
        
        if(l == 0) {
            gammaV = 0; gammaT = 0;
        }
        if(l == 1) {
            gammaV = 1; gammaT = 0;
        }
        if(l == 2) {
            gammaV = 0; gammaT = 1;
        }
        if(l == 3) {
            gammaV = 1; gammaT = 1;
        }
        
        /* saving output and computing posterior quantities */
        if(iiter == itersave) {
            fprintf(file_k, "%d\n", k);
            fprintf(file_gammaV, "%d\n", gammaV);
            fprintf(file_gammaT, "%d\n", gammaT);
            for( l = 0; l < nSB-1; l++ ) {
                for( jb = 0; jb < pCov; jb++ ) {
                    fprintf(file_betaV, "%lf\n", betaV[l][jb]);
                }
            }
            for( l = 0; l < nSB; l++ ) {
                for( j = 0; j < d; j++ ) {
                    for( jb = 0; jb < pCov; jb++ ) {
                        fprintf(file_betaT, "%lf\n", betaT[l][j][jb]);
                    }
                }
            }
            lpost = 0.0;
            for( i = 0; i < *nrec; i++ ) {
                tmp = 0.0;
                for( l1 = 0; l1 < nSB; l1++ ) {
                    for(j = 0; j < d; j++) {
                        alpha[j] = (double) ceil(k * theta[i][l1][j]) - 1;
                    }
                    tmp += w[i][l1] * dirDenData[i][k-1][(int) alpha[0]][(int) alpha[1]];;
                }
                lpost += log(tmp);
            }
            fprintf(file_llik, "%lf\n", lpost);
            
            
            isave += 1;
            itersave = nburn + (isave) * njump;
            
        } // end iiter
        
    } // end mcmc
    
   
    
    PutRNGstate();
    
    
    free(dirDenData);
    
    fclose(file_k);
    fclose(file_gammaV);
    fclose(file_gammaT);
    fclose(file_betaV);
    fclose(file_betaT);
    fclose(file_llik);
    
    

}
