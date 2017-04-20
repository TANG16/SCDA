#  MODEL DESCRIPTION
#***********************************************************************
### var G[N_max,N_max,T], P[N_max,N_max,T];
var Na[T], G[N_max,(T+1)], P[N_max,(T+1)]; # time verying vectors enought in our case

# Model

model{ 
  # THE COVARIATES ####
	# Define the priors for the logistic regression parameters
	alpha1 ~ dnorm(0,0.01)
	alphaa ~ dnorm(0,0.01)
	alphar ~ dnorm(0,0.01)
	alphal ~ dnorm(0,0.01)
	beta1 ~ dnorm(0,0.01)
	betaa ~ dnorm(0,0.01)
	betar ~ dnorm(0,0.01)
	betal ~ dnorm(0,0.01)
  
	# Define the observation error prior
	sigy <- 1/tauy
	tauy ~ dgamma(0.001,0.001)

	# Define the logistic regression equations
	# for(t in 1:(T-1)){
	for(t in 1:T){
	  logit(phi1[t]) <- alpha1 + beta1*f[t] # corresponds to the year 1963
		logit(phia[t]) <- alphaa + betaa*f[t]
		
		# log(rho[t]) <- alphar + betar*t # We assume here that t=1
		log(rho[t]) <- alphar + betar*stdT[t] # We assume here that t=1

		# logit(lambda[t]) <- alphal + betal*(t+1)
		logit(lambda[t]) <- alphal + betal*stdT[t]
	}

	# THE STATE SPACE MODEL ####
	# 0-1 trick for the states process, the obervation process can stay the same
	
	# Define r[t]
	# for (t in 3:(T-1)){
	# 	r[t-2] <- (Na[t+1]+N1[t+1])/(Na[t]+N1[t])
	# }

	# Define the initial population priors
	 for(t in 1:2){
	# 	# N1[t] ~ dnorm(200,0.000001)
	# 	# Na[t] ~ dnorm(1000,0.000001) --> 1000+-1000
	#   N1[t] ~ dpois(200)
	  Na[t] ~ dbin(0.5,2000) # --> 1000+-500
	}

	#####
	# Define the system process for the census/index data using the Normal approximation
	# for(t in 3:T){
	# 	mean1[t] <- rho[t-1]*phi1[t-1]*Na[t-1]
	# 	meana[t] <- phia[t-1]*(N1[t-1]+Na[t-1])
	# 	
	# 	tau1[t] <- 1/(Na[t-1]*rho[t-1]*phi1[t-1])
	# 	taua[t] <- 1/((N1[t-1]+Na[t-1])*phia[t-1]*(1-phia[t-1]))
	# 	
	#   N1[t] ~ dnorm(mean1[t],tau1[t])
	# 	Na[t] ~ dnorm(meana[t],taua[t])
	# }
	
	# Define the system process for the census/index data using the Poisson/Binomial model 
	# NB. Need to change initial population priors as well to ensure N1 and Na take integer values
	# 
	# for(t in 3:T){
	#   bin1[t] <- N1[t-1]+Na[t-1]
	#   bin2[t] <- phia[t-1]
	#   
	#   po[t] <- Na[t-1]*rho[t-1]*phi1[t-1]
	#   
	#   N1[t] ~ dpois(po[t])
	#   Na[t] ~ dbin(bin2[t],bin1[t])
	# }
	
	# Zero trick for loglik of HMM ####
	# an observation x[i] contributes a likelihood L[i] 
	# the "zeros trick": a Poisson(phi) observation of zero has likelihood exp(-phi), 
	# if the observed data is a set of 0's, and phi[i] is set to - log(L[i]), 
	# we will obtain the correct likelihood contribution for x[i] i.e. L[i]
	# defining:
	# spy[i] ~ pdf(y[i],params)/C # scaled probability of y[i] with pdf the required formula 
	# 1 ~ dben(spy[i])
	# together yield the same thing as
	# y[i] ~ pdf(params) 
	# which essentially is the value of the pdf when y[i] has its particular value and when the params have their particular randomly generated MCMC values 
	
	for (t in 3:T){
	  # Use a discrete Uniform prior so the only influence on the posterior distr is the Upper limit
	  Na[t] ~ dcat(Na_prior[]) # Na_prior = rep(1/(Up+1), Up+1); entered as data
	}
	
	C <- 1000000
	for(t in 3:(T+1)){ 
	  zeros[t] ~ dpois(phi[t])
	  phi[t] <- -loglik[t] + C
	  
 	  lam[t] <- Na[t-1]*rho[t-1]*phi1[t-1]
 	  
	  for (i in 1:(N_max-1)){  
      # logfact is the log of the factorial: log(x!)
      G[i,t] <- exp(-lam[t] + i*log(lam[t]) - logfact(i))
#       lf[t] <- logfact(i + Na[t-1]) - logfact(i) - logfact(Na[t-1])
# 	    P[i,t] <- exp(Na[t-1]*log(phia[t-1]) + i*log(1-phia[t-1]) + lf[t])
      P[i,t] <- exp(Na[t-1]*log(phia[t-1]) + i*log(1-phia[t-1]) + logfact(i + Na[t-1]) - logfact(i) - logfact(Na[t-1]))      
	  } 
 	  
	  G[N_max,t] <- 1- sum(G[1:(N_max-1),t])
	  P[N_max,t] <- 1- sum(P[1:(N_max-1),t])
	  # loglik[t] <- sum(sum(G[,,t] %*% P[,,t] ))
	  loglik[t] <- log(sum(G[,t] * P[,t])) # piecewise multiplication enough here
	  # loglik[t] <- sum(sum(G[,,t] %*% P[,,t] %*% Q[,,t])) <--- FOR OWLS
  }
	
	# Define the observation process for the census/index data
	for(t in 3:T){
	    y[t] ~ dnorm(Na[t],tauy)
 	}

	# THE RECOVERY MODEL ####
	# Calculate the no. of birds released each year
	# for(t in 1:T1){
	#  	rel[t] <- sum(m[t,])
	# }
	 
	# Define the recovery likelihood 
	for(t in 1:T1){
		m[t, 1:(T2+1)] ~ dmulti(p[t,], rel[t])
	}

	# Calculate the cell probabilities for the recovery table 
	
	for(t1 in 1 : (T1-1)){
		# Calculate the diagonal
		p[t1, t1] <- lambda[t1]*(1-phi1[t1])
		
		# Calculate value one above the diagonal
		p[t1, t1+1] <- lambda[t1+1]* phi1[t1]*(1-phia[t1+1])

		# Calculate remaining terms above diagonal
		for(t2 in (t1+2):T2){
			for(t in (t1+1):(t2-1)){
				lphi[t1, t2, t] <- log(phia[t])
			}
			# Probabilities in table
			p[t1,t2] <- lambda[t2]*phi1[t1]*(1-phia[t2])*exp(sum(lphi[t1,t2,(t1+1):(t2-1)]))
		}
		for(t2 in 1:(t1-1)){
			# Zero probabilities in lower triangle of table
			p[t1, t2] <- 0
		}
		# Probability of an animal never being seen again
		p[t1, T2+1] <- 1 - sum(p[t1,1:T2])	
	}

	# Final row
	p[T1,T1] <- lambda[T1]*(1-phi1[T1])
	
	for(t in 1:(T1-1)){
		p[T1,t] <- 0
	}
	p[T1,T1+1] <- 1 - p[T1,T1]
}