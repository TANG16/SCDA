#include "mex.h"
#include "math.h"
#include "matrix.h"

#define  PI  3.1415926535897932
#define  log2PI 1.837877066409345

/* ********************************************************************* */


void loglik_h_HMM_adapt_eff_mex(double *y, double *h, double *theta, double *mid_inv,
        mwSignedIndex T, mwSignedIndex N_q,
        double *loglik)
{
    mwSignedIndex i, j, t;
    double *bin_midpoint, *mu_bin, *exp_mu_bin, *loglik_int;
    double mu, phi, sigma2, sigma, beta, h0; 
    double Gauss_const, temp;

    mu = theta[0];
    phi = theta[1];
    sigma2 = theta[2];
    sigma = sqrt(sigma2);
    beta = theta[3];    

    Gauss_const = log2PI + log(sigma2);
    Gauss_const = - 0.5*Gauss_const;
//     mexPrintf("log2PI = %6.4f\n",log2PI);
//     mexPrintf("sigma2 = %6.4f\n",sigma2);
//     mexPrintf(" log(sigma2) = %6.4f\n", log(sigma2));
//     mexPrintf("Gauss_const = %6.4f\n",Gauss_const);
    
    h0 = mu;
    
    /* Variable size arrays */
    bin_midpoint = mxMalloc((N_q)*sizeof(double));  
    mu_bin = mxMalloc((N_q)*sizeof(double));  
    exp_mu_bin = mxMalloc((N_q)*sizeof(double));  
    loglik_int = mxMalloc((N_q)*sizeof(double));  

    
//     for (j=0; j<100; j++)
//     {
//         mexPrintf("h[%i] = %6.4f\n",j,h[j]);
//     }
    
    for (j=0; j<T/2; j++)
    {    
        t = 2*j+1;
        loglik[j] = 0;

        for (i=0; i<N_q; i++)
        {
            // determine the quantiles
            if (j == 0)
            {
                bin_midpoint[i] = phi*(h0-mu) + sigma*mid_inv[i]; 
//                 mexPrintf("bin_midpoint[%i,%i] = %6.4f\n",i,j,bin_midpoint[i]);
            }
            else
            {
                bin_midpoint[i] = phi*(h[t-2]-mu) + sigma*mid_inv[i];
//                 mexPrintf("bin_midpoint[%i,%i] = %6.4f\n",i,j,bin_midpoint[i]);
                
            }
            
//             if (j<10)
//             {
//                 mexPrintf("bin_midpoint[%i,%i] = %6.4f\n",i,j,bin_midpoint[i]);
//             }
            
            mu_bin[i] = (mu + bin_midpoint[i]);
            exp_mu_bin[i] = exp(mu + bin_midpoint[i]);        

            loglik_int[i] = Gauss_const - 0.5*(pow(h[t] - mu - phi*bin_midpoint[i],2)/sigma2);
            temp = log2PI + mu_bin[i] + pow(y[t-1] - beta*exp_mu_bin[i],2)/exp_mu_bin[i];
            loglik_int[i] = loglik_int[i] - 0.5*temp;    
            loglik_int[i] = exp(loglik_int[i]);
            
            loglik[j] = loglik[j] + loglik_int[i]; 
        }
//  loglik(1,t/2) = log(sum(exp(loglik_int)));   
        loglik[j] = log(loglik[j]);
    }
 
    /* Free allocated memory */
    mxFree(bin_midpoint); 
    mxFree(mu_bin); 
    mxFree(exp_mu_bin); 
    mxFree(loglik_int); 
}

/* ********************************************************************* */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex T, T2, N_q;                         /* size of matrix */
    double *y, *h, *theta, *mid_inv;    /* input*/
    double *loglik;                     /* output */

    /* Getting the inputs */
    y = mxGetPr(prhs[0]); 
    h = mxGetPr(prhs[1]);
    theta = mxGetPr(prhs[2]);
    mid_inv = mxGetPr(prhs[3]);
    
    T = mxGetN(prhs[0]); /* no. of observations */
    N_q = mxGetN(prhs[3]); /* no of bins  */

    T2 = T/2;
//     mexPrintf("T2 = %i\n",T2);
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,T2,mxREAL); 

    /* get a pointer to the real data in the output matrix */
    loglik = mxGetPr(plhs[0]);
  
    /* call the function */
    loglik_h_HMM_adapt_eff_mex(y, h, theta, mid_inv, T, N_q, loglik);
}
