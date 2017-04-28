# Set working directory: ####
setwd("Abadi")

save_on = TRUE
# Load required packages and fix the random seed
library(rjags)
library(coda)
library(lattice)
set.seed(134522)

# MCMC details: ####

ada=100
iter=1000
th=1
cha=2

# Load data: ###
source("Abadi_Data_HMM.R")
# Set initial parameter values: ####
source("Abadi_StartingVals_HMM.R")

if (0==1){
  # v = rnorm(7)
  # v = rep(0,7)
  v= c(-2.7, 3.4, 0.2, 0.3, -1.6, -2.0, 1.3);
  
  # bp = rnorm((ti-1))
  bp = 1.5*sin(1:(ti-1));
  # fec = runif((ti-1),0,5)  
  fec = abs(cos(1:(ti-1)))
  # Juvenile survival rate
  index = v[1] + v[3] + v[4]*stdT  # Male
  phijM <- exp(index)/(1+exp(index))
  index <- v[1] + v[4]*stdT          # Female
  phij <- exp(index)/(1+exp(index))
  # Adult survival
  index <- v[1] + v[2] + v[3] +  v[4]*stdT   # Male
  phiaM <- exp(index)/(1+exp(index))
  index <- v[1] + v[2] + v[4]*stdT            # Female
  phia <- exp(index)/(1+exp(index))
  # Recapture rate
  index <- v[5] + bp     # Male
  lambdaM <- exp(index)/(1+exp(index))
  index <- bp	      # Female
  lambda <- exp(index)/(1+exp(index))
  # Immigration
  index <- v[6] + v[7]*voleH 
  im <- exp(index)
  
  pr <- matrix(0,50,26)
  prM <- matrix(0,dim(mmal)[1],dim(mmal)[2])
  q <- rep(0,ti-1)
  qM <- rep(0,ti-1)
}


# Initialise the model: ####
tstart = proc.time()
mod <- jags.model('Abadi_Bugs_HMM.R',data,inits,n.chains=cha,n.adapt=ada)
temp = proc.time()-tstart
time_init <- temp

# Run the MCMC simulations: ####
tstart = proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
temp = proc.time()-tstart
time_sample <- temp

if (save_on) {
  save(mod, output1, time_init, time_sample, file = "Abadi_HMM_iter10000.RData")
}

summary(output1)

mat1 = as.matrix(output1[1])
mat2 = as.matrix(output1[2])

# Retrieve the names ####
mat1_names <- colnames(mat1)


frame1 <- do.call(rbind.data.frame, output1)
frame2 <- do.call(rbind.data.frame, output2)

str(output1) # List of 2
str(output2) # List of 12


plot(output1[,2])
output1[j] [,"par_name"]
# it returns the j-th chain for the parameter of interest
plot(output1[1][,2])
plot(output1[2][,2])
plot(output1[1][,2])

summary(output1)

mat1 <- as.matrix(output1)