# params_hmm <- c('mu','sigma2','phi','h','G_odd','G_even','Q_odd','Q_even')
# params_hmm <- c('mu','sigma2','phi','G_odd','G_even')
# params_hmm <- c('mu','sigma2','phi','Q_odd','Q_even')

mat_names_HMM[1:N_bin] # G_even[1,1] - G_even[30,1]
mat_names_HMM[N_bin+(1:N_bin)] # G_even[1,2] - G_even[30,2]
mat_names_HMM[(T/2-1)*N_bin+(1:N_bin)] # G_even[1,500] - G_even[30,500]
mat_names_HMM[(T/2)*N_bin+(1:N_bin)] # G_odd[1,1] - G_odd[30,1]
mat_names_HMM[(T-1)*N_bin+(1:N_bin)] # G_odd[1,500] - G_odd[30,500]
mat_names_HMM[T*N_bin+(1:3)] # mu phi sigma2  

mat_names_HMM[1:T/2] # Q_even[1] - Q_even[500]
mat_names_HMM[T/2+(1:N_bin)] # Q_odd[1,1] - Q_odd[30,1]
mat_names_HMM[T/2+N_bin+(1:N_bin)] # Q_odd[1,2] - Q_odd[30,2]
mat_names_HMM[T/2+(T/2-1)*N_bin+(1:N_bin)] # Q_odd[1,500] - Q_odd[30,500]
mat_names_HMM[T/2+(T/2)*N_bin+(1:3)] # mu phi sigma2 
(T/2+1):(T/2+(T/2)*N_bin)
### plot G matrices
myCol = rainbow(50)
par(mfrow=c(1,1))
plot(seq(1,50),seq(1,500,by=10),col=myCol,xlab="", ylab="",main="Colors for selected times ")

plotname = "SV_HMM_G_mats.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")

par(mfrow = c(4,2))
for (j in c(100,500,750,1000)){
  G_even = matrix(mat1_HMM[j,1:((T/2-1)*N_bin+N_bin)],ncol=T/2,nrow=N_bin,byrow=FALSE)
  G_odd = matrix(mat1_HMM[j,((T/2)*N_bin+1):((T-1)*N_bin+N_bin)],ncol=T/2,nrow=N_bin,byrow=FALSE)
  plot(G_even[,1],type="l",xlab="", ylab="",main=paste("G_even, iter=", j,sep="" ))
  for (t in seq(1,500,by=10)){
    lines(G_even[,t],type="l",col=myCol[1+(t-1)/10])
  }
  plot(G_odd[,1],type="l",xlab="", ylab="",main=paste("G_odd, iter=", j,sep="" ))
  for (t in seq(1,500,by=10)){
    lines(G_odd[,t],type="l",col=myCol[1+(t-1)/10])
  }
}

dev.off()


## Q mats

plotname = "SV_HMM_Q_mats2.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")

par(mfrow = c(4,2))
for (j in c(400,600,700,900)){
  Q_even = mat1_HMM[j,1:(T/2)]
  plot(Q_even,type="l",xlab="", ylab="",main=paste("Q_even, iter=", j,sep="" ))
  
  Q_odd = matrix(mat1_HMM[j,(T/2+1):(T/2+(T/2)*N_bin)],ncol=T/2,nrow=N_bin,byrow=FALSE)
  
  plot(Q_odd[,1],type="l",xlab="", ylab="",main=paste("Q_odd, iter=", j,sep="" ))
  for (t in seq(1,500,by=10)){
    lines(Q_odd[,t],type="l",col=myCol[1+(t-1)/10])
  }
}

dev.off()


### G mats iter
plotname = "SV_HMM_G_mats_iter.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")

myCol2 = rainbow(6)
par(mfrow = c(4,2))
for (t in c(1,100,300,500)){
  plot(mat1_HMM[, (250-1)*N_bin +1],type="l",xlab="", ylab="",main="G_even,trace for time 250 bins 1:5:30")
  for (b in seq(1,30,by=5)){
    lines(mat1_HMM[, (250-1)*N_bin +b],type="l",xlab="", ylab="",col=myCol2[1+(b-1)/5])
  }
  
  plot(mat1_HMM[, (T/2)*N_bin+ (250-1)*N_bin +1],type="l",xlab="", ylab="",main="G_odd,trace for time 250 bins 1:5:30")
  for (b in seq(1,30,by=5)){
    lines(mat1_HMM[,(T/2)*N_bin+ (250-1)*N_bin +b],type="l",xlab="", ylab="",col=myCol2[1+(b-1)/5])
  }  
}

dev.off()

### G mats colsums
plotname = "SV_HMM_G_mats_colsums.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")
par(mfrow = c(4,2))
for (j in c(100,500,750,1000)){
  G_even = matrix(mat1_HMM[j,1:((T/2-1)*N_bin+N_bin)],ncol=T/2,nrow=N_bin,byrow=FALSE)
  G_odd = matrix(mat1_HMM[j,((T/2)*N_bin+1):((T-1)*N_bin+N_bin)],ncol=T/2,nrow=N_bin,byrow=FALSE)
  
  plot(colSums(G_even),type="l",xlab="",ylab="",main=paste("ColSums of G_even, iter=", j,sep="" ))
  
  plot(colSums(G_odd),type="l",xlab="",ylab="",main=paste("ColSums of G_odd, iter=", j,sep="" ))
}
dev.off()



###### ADAPTIVE

mat_names_HMM_adapt[1:(N_q*(T/2))] # "Q_odd[1,1]" - "Q_odd[10,500]"
mat_names_HMM_adapt[N_q*(T/2)+(1:(N_q*(T/2)))] # "bin_mid[1,1]" - "bin_mid[10,500]"
mat_names_HMM_adapt[T*N_q+ (1:T/2)] # "h[1]" - "h[500]"
mat_names_HMM_adapt[T*N_q+ T/2 + (1:3) ] #  "mu"     "phi"    "sigma2"


par(mfrow = c(3,1))
for (i in 1:3){
  plot(mat1_HMM_adapt[,T*N_q+T/2+i],type='l',col='blue', xlab="", ylab="", sub=mat_names_HMM_adapt[T/2+i])
  lines(mat2_HMM_adapt[,T*N_q+T/2+i],type='l',col='red')
}

par(mfrow = c(4,2))
for (j in c(100,500,750,1000)){
  Q_odd = matrix(mat1_HMM_adapt[j,1:(N_q*(T/2))],ncol=T/2,nrow=N_q,byrow=FALSE)
  bin_mid = matrix(mat1_HMM_adapt[j,N_q*(T/2)+(1:(N_q*(T/2)))],ncol=T/2,nrow=N_q,byrow=FALSE)
  
  plot(Q_odd[,1],type="l",xlab="", ylab="",main=paste("Q_odd, iter=", j,sep="" ))
  for (t in seq(1,500,by=10)){
    lines(Q_odd[,t],type="l",col=myCol[1+(t-1)/10])
  }
  plot(bin_mid[,1],type="l",xlab="", ylab="",main=paste("bin_mid, iter=", j,sep="" ))
  for (t in seq(1,500,by=10)){
    lines(bin_mid[,t],type="l",col=myCol[1+(t-1)/10])
  }
}


mat1_HMM_adapt[j,N_q*(T/2)+(1:N_q*(T/2))]
