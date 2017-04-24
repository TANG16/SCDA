#  MODEL DESCRIPTION
#  Age structured model ( 2 age classes: 1-year and 2 years or older)
#  Age at first breeding =1 year
#  Pre-breeding census
#  Little owl data (1978-2003)
#***********************************************************************

var G[(N_max+1),(N_max+1),(ti+1)], P[(N_max+1),(N_max+1),(ti+1), Q[(N_max+1),(N_max+1),(ti+1)]; 
                                     
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
    
    logit(pM[i]) <- v[5] + bp[i]     # Male
    
    logit(p[i]) <- bp[i]	      # Female
    
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
  } #i
  
  # prob1 <- rep(1/101,101)
  
  # N1prior ~ dcat(prob1)     		# 1-year
  # N1[1] <- N1prior-1
  
  # NadSurvprior ~ dcat(prob1)  	 # Adults 
  # NadSurv[1] <- NadSurvprior-1
  # 
  # Nadimmprior ~ dcat(prob1) 	# Immigrants
  # Nadimm[1] <- Nadimmprior-1
  
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

  # Due to the zeros trick: use a discrete uniform prior so the only influence on the posterior distr is the upper limit
  for (t in 1:ti){
    NadSurv[t] ~ dcat(Nad_prior[]) # Nad_prior = rep(1/(Up+1), Up+1); entered as data
    NadImm[t] ~ dcat(Nad_prior[]) # Nad_prior = rep(1/(Up+1), Up+1); entered as data
  }
  
  
  C <- 1000000
  for(t in 2:(ti+1)){ 
    zeros[t] ~ dpois(phi[t])
    phi[t] <- -loglik[t] + C
    
    Ntot[tt] <- NadSurv[tt] + NadImm[tt] + N1[tt]
    mean1[t] <- 0.5*fec[t-1]*phij[t-1]*Ntot[t-1]
    mpo[t] <- Ntot[t-1]*im[t-1]
    
    for (j in 0:(N_max-1)){  
      Ntot[j,t-1] = j + NadSurv[t-1] + NadImm[t-1]
      lam[j,t] <- 0.5*fec[t-1]*phij[t-1]*Ntot[j,t-1]
      mu[j,t] <- im[t-1]*Ntot[j,t-1]
      
      for (i in 0:(N_max-1)){  # from 0!!!
 
        # logfact is the log of the factorial: log(x!)
        G[j+1,i+1,t] <- exp(-lam[j,t] + i*log(lam[j,t]) - logfact(i))
        P[j+1,i+1,t] <- exp(Na[t-1]*log(phia[t-1]) + i*log(1-phia[t-1]) + logfact(i + Na[t-1]) - logfact(i) - logfact(Na[t-1]))  
        
        Q[j+1,i+1,t] <- exp()
      }  
      
      G[j,(N_max+1),t] <- max(0,1- sum(G[j,1:(N_max),t])0
      P[j,(N_max+1),t] <-  

    }
    
    
    loglik[t] <- sum(sum(G[,,t] %*% P[,,t] %*% Q[,,t])) 
  }
  #***************************
  # System process
  #***************************
  for (tt in 2:ti)
  {
    mean1[tt] <- 0.5*fec[tt-1]*phij[tt-1]*Ntot[tt-1]
    N1[tt] ~ dpois(mean1[tt])
    
    mpo[tt] <- Ntot[tt-1]*im[tt-1]
    NadImm[tt] ~ dpois(mpo[tt])
    NadSurv[tt] ~ dbin(phia[tt-1],Ntot[tt-1])
  } 
  
  #*****************************
  # Observation process
  #*****************************
  for(tt in 1:ti)
  {
    Ntot[tt] <- NadSurv[tt] + Nadimm[tt] + N1[tt]
    popcount[tt] ~ dpois(NadSurv[tt] + NadImm[tt] + N1[tt]) #y~Po(N1+Na) 
  } 
  
  
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
    q[i] <- 1-p[i]
    
    # Main diagonal
    pr[i,i]<-phij[i]*p[i]
    
    # above main diagonal
    for(j in (i+1):(ti-1))
    {
      pr[i,j] <- phij[i]*prod(phia[(i+1):j])*prod(q[i:(j-1)])*p[j]
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
    pr[i+ti-1,i] <- phia[i]*p[i]
    
    # above main diagonal
    for(j in (i+1):(ti-1))
    {
      pr[i+ti-1,j] <- prod(phia[i:j])*prod(q[i:(j-1)])*p[j]
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
    qM[i] <- 1-pM[i]
    
    # main diagonal
    prM[i,i] <- phijM[i]*pM[i]
    
    # above main diagonal
    for(j in (i+1):(ti-1)) 
    {
      prM[i,j] <- phijM[i]*prod(phiaM[(i+1):j])*prod(qM[i:(j-1)])*pM[j]
      
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
    prM[i+ti-1,i] <- phiaM[i]*pM[i]
    
    # above main diagonal
    for(j in (i+1):(ti-1)) 
    {
      prM[i+ti-1,j] <- prod(phiaM[(i+1):j])*prod(qM[i:(j-1)])*pM[j]
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