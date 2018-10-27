ada=10000
iter=100000
th=1
cha=1
# load the full DA 
load("Results/BKM_iter10000_ada1000_linux.RData")
# load the HMM based, fixed bins
load("Results/BKM_HMM_approx_iter10000_ada1000_COMB_unifprior.RData")
# load the HMM based, adaptive bins
load("Results/BKM_HMM_approx_iter10000_ada1000_COMB_unifprior_norm.RData")


NSE_theta_DA = rep(0, times=9)
for (i in (36+seq(1:9))){
  temp = as.vector(as.matrix(output_DA[,i]))
  NSE_theta_DA[i-36] = nse.nw(temp, lag.prewhite = 0)
}

NSE_Na_DA = rep(0, times=36)
for (i in (seq(1:36))){
  temp = as.vector(as.matrix(output_DA[,i]))
  NSE_Na_DA[i] = nse.nw(temp, lag.prewhite = 0)
}



ESS_DA[1:36]/100000

std_Na_DA = apply(matDA[,1:36], 2, sd)
## ESS  ##########################################
# DA ##########################################

matDA = as.matrix(output_DA[1])
matDA <- matDA[,37:81]
mat_names_DA <- colnames(matDA) 
time_init_DA 
time_sample_DA  

ESS_DA = lapply(output_DA,effectiveSize)
ESS_DA <- matrix(unlist(ESS_DA), ncol = 1, byrow = TRUE)
ESS_DA <- ESS_DA[37:81,]
ESS_persec_DA <- ESS_DA/time_sample_DA[3]

# ADAPTIVE ##########################################

names_HMMnorm <- c("Q=5","Q=10","Q=20","Q=50") 

mat5 = as.matrix(outputHMMnorm_5[1])
mat10 = as.matrix(outputHMMnorm_10[1])
mat20 = as.matrix(outputHMMnorm_20[1])
mat50 = as.matrix(outputHMMnorm_50[1])
mat_names_HMMnorm <- colnames(mat50) 

time_HMMnorm_init <- c(time_HMMnorm_init_5[3],
                             time_HMMnorm_init_10[3],
                             time_HMMnorm_init_20[3],
                             time_HMMnorm_init_50[3])
names(time_HMMnorm_init) <- names_HMMnorm
  

time_HMMnorm_sample <- c(time_HMMnorm_sample_5[3],
                             time_HMMnorm_sample_10[3],
                             time_HMMnorm_sample_20[3],
                             time_HMMnorm_sample_50[3])
names(time_HMMnorm_sample) <- names_HMMnorm


ESS5 = lapply(outputHMMnorm_5,effectiveSize)
ESS10 = lapply(outputHMMnorm_10,effectiveSize)
ESS20 = lapply(outputHMMnorm_20,effectiveSize)
ESS50 = lapply(outputHMMnorm_50,effectiveSize)


ESS_HMMnorm <- matrix(unlist(ESS5), ncol = 1, byrow = TRUE)
ESS_HMMnorm <- cbind(ESS_HMMnorm, matrix(unlist(ESS10), ncol = 1, byrow = TRUE))
ESS_HMMnorm <- cbind(ESS_HMMnorm, matrix(unlist(ESS20), ncol = 1, byrow = TRUE))
ESS_HMMnorm <- cbind(ESS_HMMnorm, matrix(unlist(ESS50), ncol = 1, byrow = TRUE))
colnames(ESS_HMMnorm) <- names_HMMnorm
rownames(ESS_HMMnorm) <- mat_names_HMMnorm


ESS_persec_HMMnorm <- ESS_HMMnorm/time_HMMnorm_sample

# FIXED ##########################################

names_HMM <- c("B=5/169","B=15/55","B=29/29") 

mat_5_169 = as.matrix(output_HMM_5_169[1])
mat_15_55 = as.matrix(output_HMM_15_55[1])
mat_29_29 = as.matrix(output_HMM_29_29[1])
mat_names_HMM <- colnames(mat_5_169) 

time_HMM_init <- c(time_HMM_init_5_169[3],
                   time_HMM_init_15_55[3],
                   time_HMM_init_29_29[3])
names(time_HMM_init) <- names_HMM


time_HMM_sample <- c(time_HMM_sample_5_169[3],
                     time_HMM_sample_15_55[3],
                     time_HMM_sample_29_29[3])
names(time_HMM_sample) <- names_HMM


ESS_5_169 = lapply(output_HMM_5_169,effectiveSize)
ESS_15_55 = lapply(output_HMM_15_55,effectiveSize)
ESS_29_29 = lapply(output_HMM_29_29,effectiveSize)


ESS_HMM <- matrix(unlist(ESS_5_169), ncol = 1, byrow = TRUE)
ESS_HMM <- cbind(ESS_HMM, matrix(unlist(ESS_15_55), ncol = 1, byrow = TRUE))
ESS_HMM <- cbind(ESS_HMM, matrix(unlist(ESS_29_29), ncol = 1, byrow = TRUE))
colnames(ESS_HMM) <- names_HMM
rownames(ESS_HMM) <- mat_names_HMM
 
ESS_persec_HMM <- ESS_HMM/time_HMM_sample
 


# BARPLOT ESS 
barplot(t(cbind(ESS_DA[2:36],ESS_HMM[2:36,],ESS_HMMnorm[2:36,])), beside= TRUE, 
        legend =  names_all, col = myColors)


barplot(t(cbind(ESS_DA[37:45],ESS_HMM[37:45,],ESS_HMMnorm[37:45,])), beside= TRUE, 
        legend =  names_all, col = myColors)

# BARPLOT ESS PER SECOND
barplot(t(cbind(ESS_persec_DA[2:36],ESS_persec_HMM[2:36,],ESS_persec_HMMnorm[2:36,])), 
        beside= TRUE, 
        legend =  names_all, col = myColors)

barplot(t(cbind(ESS_persec_DA[37:45],ESS_persec_HMM[37:45,],ESS_persec_HMMnorm[37:45,])), 
        beside= TRUE, 
        legend =  names_all, col = myColors)
 
## PLOT MEAN NA ########################################

names_all = c('DA',names_HMM,names_HMMnorm)
# DA
plot(colMeans(matDA[,1:36]), type='l', xlab ="", ylab="", sub="Na")
# FIXED
lines(colMeans(mat_5_169[,1:36]), type='l', col ="brown")
lines(colMeans(mat_15_55[,1:36]), type='l', col ="red")
lines(colMeans(mat_29_29[,1:36]), type='l', col ="pink")
# ADAPTIVE
lines(colMeans(mat5[,1:36]), type='l', col ="blueviolet")
lines(colMeans(mat10[,1:36]), type='l', col="blue")
lines(colMeans(mat20[,1:36]), type='l', col="cyan")
lines(colMeans(mat50[,1:36]), type='l', col="green")

myColors <- c("black",
              "brown", "red", "pink", 
              "blueviolet",  "blue",  "cyan", "green")
legend(30, 1950, names_all, lty = 1, col=myColors)


## PRINT TIME ########################################
time_all = matrix(c(time_init_DA[3],time_HMM_init,time_HMMnorm_init,
                    time_sample_DA[3],time_HMM_sample,time_HMMnorm_sample),
                  byrow=FALSE,ncol= 2)
rownames(time_all) = names_all
colnames(time_all) = c("Init","Sample")


myFile <- "time_all.txt"
write.table(round(time_all, 3), file=myFile, row.names=TRUE, col.names=TRUE)

## PRINT MEANS PARAMS  ########################################

mean_DA = colMeans(matDA[,37:45])
mean_5_169 = colMeans(mat_5_169[,37:45]) 
mean_15_55 = colMeans(mat_15_55[,37:45])  
mean_29_29 = colMeans(mat_29_29[,37:45]) 
mean_5 = colMeans(mat5[,37:45]) 
mean_10 = colMeans(mat10[,37:45]) 
mean_20 = colMeans(mat20[,37:45]) 
mean_50 = colMeans(mat50[,37:45]) 


std_DA = apply(matDA[,37:45], 2, sd)
std_5_169 = apply(mat_5_169[,37:45], 2, sd) 
std_15_55 = apply(mat_15_55[,37:45], 2, sd)
std_29_29 = apply(mat_29_29[,37:45], 2, sd)
std_5 = apply(mat5[,37:45], 2, sd)
std_10 = apply(mat10[,37:45], 2, sd)
std_20 = apply(mat20[,37:45], 2, sd)
std_50 = apply(mat50[,37:45], 2, sd)

mean_all = matrix(c(mean_DA,
                    mean_5_169, mean_15_55, mean_29_29,
                    mean_5, mean_10, mean_20, mean_50),
                  byrow=TRUE,ncol= 9)
rownames(mean_all) = names_all
colnames(mean_all) = mat_names_DA[37:45]

std_all = matrix(c(std_DA,
                   std_5_169, std_15_55, std_29_29,
                   std_5, std_10, std_20, std_50),
                  byrow=TRUE,ncol= 9)
rownames(std_all) = names_all
colnames(std_all) = mat_names_DA[37:45]


myFile <- "mean_all.txt"
write.table(round(mean_all, 3), file=myFile, row.names=TRUE, col.names=TRUE)


myFile <- "std_all.txt"
write.table(round(std_all, 3), file=myFile, row.names=TRUE, col.names=TRUE)


# TRACE Na ##########################################
for (i in seq(2,36,5)){
  plotname = paste("Figures/BKM_COMB_Trace_",mat_names_DA[i],".png",seq="")
  png(filename = plotname,width = 1121, height = 797, units = "px", 
      pointsize = 16, antialias = "cleartype")
  
  par(mfrow=c(4,2))
  
  plot(matDA[,i], type='l', col =myColors[1], xlab = "", ylab="", sub = paste(names_all[1],", ", mat_names_DA[i], seq=""))
  plot(mat_5_169[,i] ,type="l", col =myColors[2], xlab = "", ylab="", sub = paste(names_all[2],", ", mat_names_DA[i], seq=""))
  plot(mat_15_55[,i], type='l', col =myColors[3], xlab = "", ylab="", sub = paste(names_all[3],", ", mat_names_DA[i], seq=""))
  plot(mat_29_29[,i], type='l', col =myColors[4], xlab = "", ylab="", sub = paste(names_all[4],", ", mat_names_DA[i], seq=""))
  
  plot(mat5[,i], type='l', col =myColors[5], xlab = "", ylab="", sub = paste(names_all[5],", ", mat_names_DA[i], seq=""))
  plot(mat10[,i] ,type="l", col =myColors[6], xlab = "", ylab="", sub = paste(names_all[6],", ", mat_names_DA[i], seq=""))
  plot(mat20[,i], type='l', col =myColors[7], xlab = "", ylab="", sub = paste(names_all[7],", ", mat_names_DA[i], seq=""))
  plot(mat50[,i], type='l', col =myColors[8], xlab = "", ylab="", sub = paste(names_all[8],", ", mat_names_DA[i], seq=""))

  dev.off()
} 


# ACF Na
for (i in seq(2,36,5)){
  plotname = paste("Figures/BKM_COMB_ACF_",mat_names_DA[i],".png",seq="")
  png(filename = plotname,width = 1121, height = 797, units = "px", 
      pointsize = 16, antialias = "cleartype")  
  
  par(mfrow=c(4,2))
  
  acf(matDA[,i], main = paste(names_all[1],", ", mat_names_DA[i], seq=""))
  acf(mat_5_169[,i], main = paste(names_all[2],", ", mat_names_DA[i], seq=""))
  acf(mat_15_55[,i], main = paste(names_all[3],", ", mat_names_DA[i], seq=""))
  acf(mat_29_29[,i], main = paste(names_all[4],", ", mat_names_DA[i], seq=""))
  
  acf(mat5[,i], main = paste(names_all[5],", ", mat_names_DA[i], seq=""))
  acf(mat10[,i], main = paste(names_all[6],", ", mat_names_DA[i], seq=""))
  acf(mat20[,i], main = paste(names_all[7],", ", mat_names_DA[i], seq=""))
  acf(mat50[,i], main = paste(names_all[8],", ", mat_names_DA[i], seq=""))
  
  dev.off()  
} 

# TRACE Param ##########################################

for (i in (36+seq(1:9))){
  plotname = paste("Figures/BKM_COMB_Trace_",mat_names_DA[i],".png",seq="")
  png(filename = plotname,width = 1121, height = 797, units = "px", 
      pointsize = 16, antialias = "cleartype")
  
  par(mfrow=c(4,2))
  
  plot(matDA[,i], type='l', col =myColors[1], xlab = "", ylab="", sub = paste(names_all[1],", ", mat_names_DA[i], seq=""))
  plot(mat_5_169[,i] ,type="l", col =myColors[2], xlab = "", ylab="", sub = paste(names_all[2],", ", mat_names_DA[i], seq=""))
  plot(mat_15_55[,i], type='l', col =myColors[3], xlab = "", ylab="", sub = paste(names_all[3],", ", mat_names_DA[i], seq=""))
  plot(mat_29_29[,i], type='l', col =myColors[4], xlab = "", ylab="", sub = paste(names_all[4],", ", mat_names_DA[i], seq=""))
  
  plot(mat5[,i], type='l', col =myColors[5], xlab = "", ylab="", sub = paste(names_all[5],", ", mat_names_DA[i], seq=""))
  plot(mat10[,i] ,type="l", col =myColors[6], xlab = "", ylab="", sub = paste(names_all[6],", ", mat_names_DA[i], seq=""))
  plot(mat20[,i], type='l', col =myColors[7], xlab = "", ylab="", sub = paste(names_all[7],", ", mat_names_DA[i], seq=""))
  plot(mat50[,i], type='l', col =myColors[8], xlab = "", ylab="", sub = paste(names_all[8],", ", mat_names_DA[i], seq=""))
  
  dev.off()
} 

# ACF Param
for (i in (36+seq(1:9))){
  plotname = paste("Figures/BKM_COMB_ACF_",mat_names_DA[i],".png",seq="")
  png(filename = plotname,width = 1121, height = 797, units = "px", 
      pointsize = 16, antialias = "cleartype")
  
  par(mfrow=c(4,2))
  
  acf(matDA[,i], main = paste(names_all[1],", ", mat_names_DA[i], seq=""))
  acf(mat_5_169[,i], main = paste(names_all[2],", ", mat_names_DA[i], seq=""))
  acf(mat_15_55[,i], main = paste(names_all[3],", ", mat_names_DA[i], seq=""))
  acf(mat_29_29[,i], main = paste(names_all[4],", ", mat_names_DA[i], seq=""))
  
  acf(mat5[,i], main = paste(names_all[5],", ", mat_names_DA[i], seq=""))
  acf(mat10[,i], main = paste(names_all[6],", ", mat_names_DA[i], seq=""))
  acf(mat20[,i], main = paste(names_all[7],", ", mat_names_DA[i], seq=""))
  acf(mat50[,i], main = paste(names_all[8],", ", mat_names_DA[i], seq=""))
  
  dev.off()
} 