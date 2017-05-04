# setwd("../MCMC")
######
# Main code for running MCMC algorithm for capture-recapture data
# - white stork data using MH random walk update - created by RK (Nov 08)
######

# Define the function: with input parameters:
# nt = number of iterations
# nburn  = burn-in

storkcodeMH <- function(nt,nburn){

    # Define the parameter values:
    # ni = number of release years
    # nj = number of recapture years
    # nparam = maximum number of parameters
    # ncov = maximum number of covariates
    ni = 16
    nj = 16
    nparam = 3
    ncov = 1

    # Read in the data:
    data <- matrix(c(
    19,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,
    0,33,3,0,0,0,0,0,0,0,0,0,0,0,0,0,14,
    0,0,35,4,0,0,0,0,0,0,0,0,0,0,0,0,14,
    0,0,0,42,1,0,0,0,0,0,0,0,0,0,0,0,26,
    0,0,0,0,42,1,0,0,0,0,0,0,0,0,0,0,30,
    0,0,0,0,0,32,2,1,0,0,0,0,0,0,0,0,36,
    0,0,0,0,0,0,46,2,0,0,0,0,0,0,0,0,16,
    0,0,0,0,0,0,0,33,3,0,0,0,0,0,0,0,28,
    0,0,0,0,0,0,0,0,44,2,0,0,0,0,0,0,20,
    0,0,0,0,0,0,0,0,0,43,1,0,0,1,0,0,10,
    0,0,0,0,0,0,0,0,0,0,34,1,0,0,0,0,25,
    0,0,0,0,0,0,0,0,0,0,0,36,1,0,0,0,16,
    0,0,0,0,0,0,0,0,0,0,0,0,27,2,0,0,22,
    0,0,0,0,0,0,0,0,0,0,0,0,0,22,0,0,16,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,1,17,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,8),nrow=ni,byrow=T)
    
    # Read in the covariate values:
    cov <- c(0.79,2.04,1.04,-0.15,-0.01,-0.48,1.72,-0.83,
    -0.02,0.14,-0.71,0.50,-0.62,-0.88,-1.00,-1.52)
    
    mn <- array(0,nparam)
    std <- array(0,nparam)

    # Read in the priors:
    # Beta prior on recapture probability
    alphap <- 1
    betap <- 1
    # Normal priors (mean and variance) on regression coefficients
    mu <- c(0,0)
    sig2 <- c(10,10)
    sig <- sqrt(sig2)

    # Parameters for MH updates (Uniform random walk):
    delta <- c(0.05,0.1,0.1)

    # Set initial parameter values:
    param <- c(0.9,0.7,0.1)

    # param[1] = recapture rate
    # param[2:12] = regression coefficients for survival rates

    # sample is an array in which we put the sample from the posterior distribution.
    sample <- array(0, dim=c(nt, nparam))
    sample <- data.frame(sample)

    # Label the columns of the array "sample"
    names(sample) <- c("p", "beta0", "beta1")

    # Calculate log(likelihood) for initial state using a separate function "calclikhood":
    likhood <- calclikhood(ni, nj, data, param, nparam, cov, ncov)

    # Set up a vector for parameter values and associated log(likelihood):
    output <- dim(nparam+1)

    # MCMC updates - MH algorithm: ####
    # Cycle through each iteration:

    for (t in 1:nt){
        # Update the parameters in the model using function "updateparam":
        output <- updateparam(nparam, param, ncov, cov, ni, nj, data, likhood, alphap, betap, mu, sig, delta)

        # Set parameter values and log(likelihood) value of current state to be the output from
        # the MH step:

        param <- output[1:nparam]
        likhood <- output[nparam+1]

        # Record the set of parameter values:
        for (i in 1:nparam) {
            sample[t,i] <- param[i]
        }
        if (t >= (nburn+1)) {
            for (i in 1:nparam) {
                mn[i] <- mn[i] + param[i]
                std[i] <- std[i] + param[i]**2
            }
        }
    }

    # Calculate the mean and standard deviation of the parameters
    # following burn-in:
    for (i in 1:nparam) {
        mn[i] <- mn[i]/(nt-nburn)
        std[i] <- std[i]/(nt-nburn)
        std[i] <- std[i] - (mn[i])**2
    }
    
    # Output the posterior mean and standard deviation of the parameters
    # following burn-in to the screen:
    cat("Posterior summary estimates for each parameter:  ", "\n")
    cat("\n")
    cat("mean  (SD)", "\n")
    cat("p: ", "\n")
    cat(mn[1], "   (", std[1], ")", "\n")
    for (i in 2:nparam) {
        if (mn[i] != 0) {
            cat("\n")
            cat("beta_",(i-2), "\n")
            cat(mn[i], " (", std[i], ")", "\n")
        }
    }

    # Output the sample from the posterior distribution:
    sample
}

######
# This function calculates the log likelihood of capture-recapture data
######

calclikhood <- function(ni, nj, data, param, nparam, cov, ncov){

    # First we set the survival and recapture probs:
    # Set up the size of the array containing the survival and recapture probs:
    # Set up the set of cell probabilities (q), and initially set them to be all equal to zero:

    phi <- array(0,nj)
    p <- array(0,nj)
    q <- array(0,dim=c(ni,nj+1))

    # Set the recapture and survival probs for each year of the study:

    for (i in 1:nj) {
        exprn <- param[2]+param[3]*cov[i]
        phi[i] <- 1/(1+exp(-exprn))

        p[i] <- param[1]
    }

    # Set initial value of the log(likelihood) to be zero

    likhood <- 0

    # Calculate the Multinomial cell probabilities and corresponding contribution to the
    # log(likelihood):
    # Cycle through the number of years that there are releases:
    for (i in 1:ni){
        # For diagonal elements:
        q[i,i] <- phi[i]*p[i]
        likhood <- likhood + data[i,i]*log(q[i,i])

        # Calculate the elements above the diagonal:
        if (i <= (nj-1)) {
            for (j in (i+1):nj) {
                q[i,j] <- prod(phi[i:j])*prod(1-p[i:(j-1)])*p[j]
                likhood <- likhood + data[i,j]*log(q[i,j])
            }
        }
        # Probability of an animal never being seen again
        q[i,nj+1] <- 1 - sum(q[i,i:nj])
        likhood <- likhood + data[i,nj+1]*log(q[i,nj+1])
    }

  # Output the log(likelihood) value:
  likhood
}

######
# Function for updating the parameters values:
######

updateparam <- function(nparam, param, ncov, cov, ni, nj, data, likhood, alphap, betap, mu, sig, delta){
# output <- updateparam(nparam, param, ncov, cov, ni, nj, data, likhood, alphap, betap, mu, sig, delta)
  # nparam = 3
  # ncov = 1
  # delta <- c(0.05,0.1,0.1)
  # param <- c(0.9,0.7,0.1)
  # # Beta prior on recapture probability
  # alphap <- 1
  # betap <- 1
  # # Normal priors (mean and variance) on regression coefficients
  # mu <- c(0,0)
  # sig2 <- c(10,10)
  # sig <- sqrt(sig2)
  
    # Cycle through each parameter in turn and propose to update using
    # random walk MH with Uniform proposal density:
    for (i in 1:nparam) {
        # Keep a record of the current parameter value being updated
        oldparam <- param[i]

        # Propose a new value for the parameter using a random walk with
        # Uniform proposal density
        param[i] <- runif(1, param[i]-delta[i], param[i]+delta[i])

        # Automatically reject any moves where recapture prob is outside [0,1]
        if (param[1] >= 0 & param[1] <= 1) {
            # Calculate the log(acceptance probability):
    
            # Calculate the new likelihood value for the proposed move:
            # Calculate the numerator (num) and denominator (den) in turn:
            newlikhood <- calclikhood(ni, nj, data, param, nparam, cov, ncov)

            if (i == 1) {
                # For recapture probability add in prior (Beta) terms to the acceptance probability
                num <- newlikhood + log(dbeta(param[i],alphap,betap))
                den <- likhood + log(dbeta(oldparam,alphap,betap))
            }
            else {
                # For regression coefficients add in prior (Normal) terms to the acceptance probability
                num <- newlikhood + log(dnorm(param[i],mu[i-1],sig[i-1]))
                den <- likhood + log(dnorm(oldparam,mu[i-1],sig[i-1]))
            }

            # All other prior terms (for other parameters) cancel in the acceptance probability.
            # Proposal terms cancel since proposal distribution is symmetric.
            
            # Acceptance probability of MH step:
            A <- min(1,exp(num-den))
        }
        else {
            A <- 0
        }

        # To do the accept/reject step of the algorithm:
        # Simulate a random number in [0,1]:
        u <- runif(1)
        # Accept the move with probability A:

        if (u <= A) {
            # Accept the proposed move:
            # Update the log(likelihood) value:
            likhood <- newlikhood
        }
        else {
            # Reject proposed move so parameter stays at current value:
            param[i] <- oldparam
        }
    }

    # Set the values to be outputted from the function to be the
    # set of parameter values and log(likelihood) value:
    output <- c(param, likhood)

    # Output the parameter values:
    output
}