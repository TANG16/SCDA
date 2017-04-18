mat1 = as.matrix(output1[1])
IF = matrix(data=NA,nrow=116,ncol=1)
for (ii in 1:116){
  acf_curr = acf(mat1[,ii],lag=iter,plot=FALSE)
  acf_curr = acf_curr$acf
  IF_curr = 1 + 2*sum(acf_curr[-1])
  IF[ii] = IF_curr
}

IF <- nrow(mat1)/ apply(mat1, 2, effectiveSize)
ESS <- apply(mat1, 2, effectiveSize)

par(mfrow=c(1,2))
barplot(ESS)
hist(ESS)


par(mfrow=c(1,2))
barplot(IF)
hist(IF)

name = "Owl_ESS.png"
savePlot("name")
, width=400,height=350)
dev.off()
dev.print(png, 'Owl_ESS.png')
acf(mat1[,115],lag=iter)


plot(mat1[,115])

plot(output1[1] [,115])

# Regression parameters for surviaval probabilities
for (ii in 2:7){
  nam = paste("Figures/Owl/Owl_acf_v",toString(ii),".png",sep="")
  png(file=nam, width=675, height=457)

  ss <- paste("v[",toString(ii),"]",sep="")
  par(mfrow=c(1,2))
  acf(mat1[,ss], main = ss)
  acf(mat1[,ss],lag=iter, main = ss)
  
  dev.off()
}

# Derived parameters
DP = c("MEPHJUF","MEPHADF","MEFE","MEIM_H","MEIM_L","MEPHJUM","MEPHADM")
for (dp in DP){
  nam = paste("Figures/Owl/Owl_acf_",dp,".png",sep="")
  png(file=nam, width=675, height=457)
  

  par(mfrow=c(1,2))
  acf(mat1[,dp], main = dp)
  acf(mat1[,dp],lag=iter, main = dp)
  
  dev.off()
}
