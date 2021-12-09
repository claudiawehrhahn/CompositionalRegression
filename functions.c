//
//  functionscopy.c
//  



// bivariate normal distribution
double bivariate_norm_pdf(double x[2], double Sigma[2][2]) {
    
    double pi = 3.141592;
    double det, iSigma[2][2], tmp[2];
    double out;
    
    det = Sigma[0][0] * Sigma[1][1] - Sigma[0][1] * Sigma[1][0];
    iSigma[0][0] = Sigma[1][1] / det;
    iSigma[1][1] = Sigma[0][0] / det;
    iSigma[0][1] = -Sigma[0][1] / det;
    iSigma[1][0] = -Sigma[1][0] / det;
    
    tmp[0] = x[0] * iSigma[0][0] + x[1] * iSigma[1][0];
    tmp[1] = x[0] * iSigma[0][1] + x[1] * iSigma[1][1];
    out = -0.5 * (tmp[0] * x[0] + tmp[1] * x[1]);
    out -= 0.5 * log(det);
    out -= log(2 * pi);
    
    return exp(out);
}


// dirichlet density
double dirich_den(double x[d], double alpha[d], int k) {
    double sumX, sumAlpha, sumAlphaX, sumGAlpha;
    int i, cond;
    
    double out; /* output */
    
    cond = 0;
    for(i = 0; i < d; i++) {
        if(x[i] == 0) {
            cond += 1;
        }
    }
    
    sumX = 0.0;
    sumAlpha = 0.0;
    if(cond == 0) {
        sumAlphaX = 0.0;
        sumGAlpha = 0.0;
        for( i = 0; i < d; i++ ) {
            sumX += x[i];
            sumAlpha += alpha[i];
            sumAlphaX += (alpha[i]-1)*log(x[i]);
            sumGAlpha += lgamma(alpha[i]);
        }
        sumAlphaX += (k-sumAlpha-1)*log(1-sumX);
        sumGAlpha += lgamma(k - sumAlpha);
        
        out = exp( lgamma(k) + sumAlphaX - sumGAlpha ); //
    } else {
        sumAlphaX = 1.0;
        sumGAlpha = 1.0;
        for( i = 0; i < d; i++ ) {
            sumX += x[i];
            sumAlpha += alpha[i];
            sumAlphaX *= pow(x[i], (alpha[i]-1));// (alpha[i]-1)*log(x[i]);
            sumGAlpha *= exp(lgamma(alpha[i]));
        }
        sumAlphaX *= pow((1-sumX), (k-sumAlpha-1));
        sumGAlpha *= exp(lgamma(k - sumAlpha));
        
        out = exp(lgamma(k)) * sumAlphaX / sumGAlpha ; //
    }
    return out;
}


/* Poisson density in log scale */
double pois_lden(int k, double lam) {
    double out; /* output */
    double param1;
    
    param1 = (double) k + 1;
    out = - lam + k * log(lam) - lgamma(param1);
    
    return out;
}

/* density of the proposal distribution for uodating k */
double denQk(int k) {
    double out; /* output */
    
    if( k == 1 ){
        out = 0;
    } else {
        out = log(0.5);
    }
    
    return out;
}


/* candidate generator for k */
int candK(int k) {
    int out; /* output */
    
    double u; /* internal variable */
    
    u = runif((double) 0, (double) 1);
    if(u <= 0.5) {
        out = k - 1;
    } else {
        out = k + 1;
    }
    
    if(out < 1 ) {
        out = 1;
    }
    
    return out;
}






