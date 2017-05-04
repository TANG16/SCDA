#  MODEL DESCRIPTION
#  Age structured model ( 2 age classes: 1-year and 2 years or older)
#  Age at first breeding =1 year
#  Pre-breeding census
#  Little owl data (1978-2003)
#***********************************************************************

var G[(N_max+1),(N_max+1),ti], P[(N_max+1),ti], Q[(N_max+1),ti], logmu[(N_max+1),(ti-1)], loglam[(N_max+1),(ti-1)]; 
                                     
# Model
model{ 
  #*******************************************
  # Define the regression equations ####
  #*****************************************
  for (i in 1:(ti-1))
  {
    # Best model structure (phi(a2+sex+T),p(sex+t),b(t)) for little owl
    # data from Schaub et al.(2006)
    
    # Juvenile survival rate
    logit(phijM[i]) <- v[1] + v[3] + v[4]*stdT[i]  # Male
    logit(phij[i]) <- v[1] + v[4]*stdT[i]          # Female
    
    # Adult survival
    logit(phiaM[i]) <- v[1] + v[2] + v[3] +  v[4]*stdT[i]   # Male
    logit(phia[i]) <- v[1] + v[2] + v[4]*stdT[i]            # Female
    
    # Recapture rate
    logit(lambdaM[i]) <- v[5] + bp[i]     # Male
    logit(lambda[i]) <- bp[i]	      # Female
    
    # Immigration
    log(im[i]) <- v[6] + v[7]*voleH[i]     # Immigration rate as a function of 
                                           # vole abundance
                                           # voleH is a categorical variable  
                                           # with 2 levels (high,low)
  }
  
  #*******************************************
  # Define the priors for the parameters ####
  #*******************************************
  for (i in 1:7)
  {
    v[i] ~ dnorm(0,0.01)#I(-10,10)
  } #i
  
  for (i in 1:(ti-1))
  {
    bp[i] ~ dnorm(0,0.0001)I(-10,10)
    fec[i] ~ dunif(0,10)
  } 
  
  #*****************************************
  # The Integrated population model ####
  #******************************************
  
  #**************************************************
  # Likelihood for reproductive data
  #**************************************************
  for (i in 1:(ti-1))
  {
    nestlings[i] ~ dpois(rho[i])
    rho[i] <- sample.size[i]*fec[i]
  } 
  
  #***********************************************
  # Likelihood for population survey data ####
  #************************************************
  # prob1 <- rep(1/101,101)
  
  # N1prior ~ dcat(prob1)     		# 1-year
  # N1[1] <- N1prior-1
  
  # Nad_prior = rep(1/(Up+1), Up+1); entered as data
  # NadSurvprior ~ dcat(Nad_prior)  	 # Adults
  # NadSurv[1] <- NadSurvprior-1 # 20-1
  # NadImmprior ~ dcat(Nad_prior) 	# Immigrants
  # NadImm[1] <- NadImmprior-1 # 5-1
  
  # Due to the zeros trick: use a discrete uniform prior so the only influence on the posterior distr is the upper limit
  for (t in 1:ti){
    NadSurv[t] ~ dcat(Nad_prior[]) 
    NadImm[t] ~ dcat(Nad_prior[])
  }
  
  for (t in 1:ti){
    for (j in 0:N_max){
      Ntot[j+1,t] = j + NadSurv[t] + NadImm[t]
      # Q[j+1,t] <- exp(-Ntot[j+1,t] + popcount[t]*log(Ntot[j+1,t]) - logfact(popcount[t]))
      Q[j+1,t] <- exp(-Ntot[j+1,t] + popcount[t]*log(Ntot[j+1,t]) - logfact_m[(popcount[t]+1)])
     }
  }
  
  for (t in 1:(ti-1)){
    for (j in 0:N_max){ 
      logmu[j+1,t] <- log(im[t]) + log(Ntot[j+1,t])  
      loglam[j+1,t] <- log(0.5) + log(fec[t]) + log(phij[t]) + log(Ntot[j+1,t]) 
    }
  }
  
  # prior for transition probabilities
  for (j in 0:N_max){
    for (i in 0:N_max){
      G[j+1,i+1,1] <- 1/(N_max+1)
      G[j+1,i+1,2] <- 1/(N_max+1)
    }
    # prior for first augmented observation probabilities
    P[j+1,1] <- Nad_prior[1]
    # P[j+1,2] <- ifelse((Ntot[j+1,1] - NadSurv[2])>0,
    #             exp(NadSurv[2]*log(phia[1]) + (Ntot[j+1,1] - NadSurv[2])*log(1-phia[1]) # part for survivors
    #             + logfact(Ntot[j+1,1]) - logfact(abs(Ntot[j+1,1] - NadSurv[2])) - logfact(NadSurv[2])
    #             - exp(logmu[j+1,1]) + NadImm[2]*logmu[j+1,1] - logfact(NadImm[2])), # part for immigrants
    #             0)
    P[j+1,2] <- ifelse((Ntot[j+1,1] - NadSurv[2])>0,
                       exp(NadSurv[2]*log(phia[1]) + (Ntot[j+1,1] - NadSurv[2])*log(1-phia[1]) # part for survivors
                           + logfact_m[(Ntot[j+1,1]+1)] - logfact_m[(abs(Ntot[j+1,1] - NadSurv[2])+1)] - logfact_m[(NadSurv[2]+1)]
                           - exp(logmu[j+1,1]) + NadImm[2]*logmu[j+1,1] - logfact_m[(NadImm[2]+1)]), # part for immigrants
                       0)    
  }
  
  
  C <- 1000000
  loglik[1] <- sum(sum(G[,,1] %*% (P[,1] * Q[,1]))) 
  loglik[2] <- sum(sum(G[,,2] %*% (P[,2] * Q[,2]))) 
  
  for(t in 1:ti){ 
    zeros[t] ~ dpois(phi[t])
    phi[t] <- -loglik[t] + C
  }
  
  for(t in 3:ti){ 
    # Old states up to N_max
    for (j in 0:(N_max)){  # old state N1_{t-2} = j  
      # Transition  probabilities - from 0 !!!
      for (i in 0:(N_max-1)){  # new state N1_{t-1} = i <- t-1 due to the unconditional probability formula for Na_{t}
        # G[j+1,i+1,t] <- exp(-exp(loglam[j+1,t-2]) + i*loglam[j+1,t-2] - logfact(i)) # t: to match with the observation probability for Na at t
        G[j+1,i+1,t] <- exp(-exp(loglam[j+1,t-2]) + i*loglam[j+1,t-2] - logfact_m[(i+1)]) # t: to match with the observation probability for Na at t        
      }
      G[j+1,(N_max+1),t] <- max(0,1- sum(G[j+1,1:N_max,t]))
      P[j+1,t] <- ifelse((Ntot[j+1,t-1] - NadSurv[t])>0,
                      exp(NadSurv[t]*log(phia[t-1]) + (Ntot[j+1,t-1] - NadSurv[t])*log(1-phia[t-1]) # part for survivors
                      + logfact_m[(Ntot[j+1,t-1]+1)] - logfact_m[(abs(Ntot[j+1,t-1] - NadSurv[t])+1)] - logfact_m[(NadSurv[t]+1)]
                      - exp(logmu[j+1,t-1]) + NadImm[t]*logmu[j+1,t-1] - logfact_m[(NadImm[t]+1)]), # part for immigrants
                      0)  
      }
    loglik[t] <- log(sum(sum(G[,,t] %*% (P[,t] * Q[,t])))) 
  }
  

  
  # #####
  #***************************
  # System process
  #***************************
  # for (tt in 2:ti)
  # {
  #   mean1[tt] <- 0.5*fec[tt-1]*phij[tt-1]*Ntot[tt-1]
  #   N1[tt] ~ dpois(mean1[tt])
  #   
  #   mpo[tt] <- Ntot[tt-1]*im[tt-1]
  #   NadImm[tt] ~ dpois(mpo[tt])
  #   NadSurv[tt] ~ dbin(phia[tt-1],Ntot[tt-1])
  # } 
  
  #*****************************
  # Observation process
  #*****************************
  # for(tt in 1:ti)
  # {
  #   Ntot[tt] <- NadSurv[tt] + Nadimm[tt] + N1[tt]
  #   popcount[tt] ~ dpois(NadSurv[tt] + NadImm[tt] + N1[tt]) #y~Po(N1+Na) 
  # } 
  
  
  #***********************************************************
  # Likelihood for capture-recapture data : CJS models (2 age classes) ####        
  #***********************************************************

  #***********************************
  # Female capture recapture data
  #***********************************
  
  for( i in 1:(2*(ti-1))){
    m[i,1:ti] ~ dmulti(pr[i,],r[i])
  } 
  
  # m-array cell probabilities for juveniles
  for(i in 1:(ti-1))
  {
    q[i] <- 1-lambda[i]
    
    # Main diagonal
    pr[i,i]<-phij[i]*lambda[i]
    
    # above main diagonal
    for(j in (i+1):(ti-1))
    {
      pr[i,j] <- phij[i]*prod(phia[(i+1):j])*prod(q[i:(j-1)])*lambda[j]
    } 
    
    # Below main diagonal
    for( j in 1:(i-1))
    {
      pr[i,j] <- 0
    } 
    
    # Last column
    pr[i,ti] <- 1-sum(pr[i,1:(ti-1)])
  } 
  
  # m-array cell probabilities for adults
  for(i in 1:(ti-1))
  {
    # main diagonal
    pr[i+ti-1,i] <- phia[i]*lambda[i]
    
    # above main diagonal
    for(j in (i+1):(ti-1))
    {
      pr[i+ti-1,j] <- prod(phia[i:j])*prod(q[i:(j-1)])*lambda[j]
    } # j
    
    # below main diagonal
    for(j in 1:(i-1))
    {
      pr[i+ti-1,j] <- 0
    } # j
    
    # last column
    pr[i+ti-1,ti] <- 1-sum(pr[i+ti-1,1:(ti-1)])
  } 
  
  #*********************************
  # Male capture recapture data
  #*********************************
  
  for(i in 1:(2*(ti-1)))
  {
    mM[i,1:ti] ~ dmulti(prM[i,],rM[i])
  } 
  
  # m-array cell probabilities for juveniles
  for(i in 1:(ti-1))
  {
    qM[i] <- 1-lambdaM[i]
    
    # main diagonal
    prM[i,i] <- phijM[i]*lambdaM[i]
    
    # above main diagonal
    for(j in (i+1):(ti-1)) 
    {
      prM[i,j] <- phijM[i]*prod(phiaM[(i+1):j])*prod(qM[i:(j-1)])*lambdaM[j]
      
    } 
    
    # below main diagonal
    for(j in 1:(i-1)) 
    {
      prM[i,j] <- 0
    }
    
    # last column
    prM[i,ti] <- 1-sum(prM[i,1:(ti-1)])
  } 
  
  # m-array cell probabilities for adults
  for(i in 1:(ti-1)) 
  {
    # main diagonal
    prM[i+ti-1,i] <- phiaM[i]*lambdaM[i]
    
    # above main diagonal
    for(j in (i+1):(ti-1)) 
    {
      prM[i+ti-1,j] <- prod(phiaM[(i+1):j])*prod(qM[i:(j-1)])*lambdaM[j]
    } 
    
    # below main diagonal
    for(j in 1:(i-1)) 
    {
      prM[i+ti-1,j] <- 0
    } 
    
    # last column
    prM[i+ti-1,ti] <- 1-sum(prM[i+ti-1,1:(ti-1)])
  }   
}