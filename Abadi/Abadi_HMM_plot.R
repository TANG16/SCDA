par(mfrow=c(2,1),oma=c(0,0,1.5,0))
plot(colMeans(mat1[,(ti+(2:ti))]), type='l', xlab ="", ylab="", sub="NadSurv")
plot(colMeans(mat1[,(2:ti)]), type='l', xlab ="", ylab="", sub="Nadimm")
mtext("Posterior means", outer=TRUE, cex=1)


# ACF

par(mfrow=c(3,3))
acf(mat1[,ncol(mat1)-6],main = "v[1]")
acf(mat1[,ncol(mat1)-5],main = "v[2]")
acf(mat1[,ncol(mat1)-4],main = "v[3]")
acf(mat1[,ncol(mat1)-3],main = "v[4]")
acf(mat1[,ncol(mat1)-2],main = "v[5]")
acf(mat1[,ncol(mat1)-1],main = "v[6]")
acf(mat1[,ncol(mat1)],main = "v[7]")


par(mfrow=c(3,3))
acf(mat1[,1+26],main = "NadSurv[1]")
acf(mat1[,3+26],main = "NadSurv[3]")
acf(mat1[,7+26],main = "NadSurv[7]")
acf(mat1[,11+26],main = "NadSurv[11]")
acf(mat1[,14+26],main = "NadSurv[14]")
acf(mat1[,17+26],main = "NadSurv[17]")
acf(mat1[,20+26],main = "NadSurv[20]")
acf(mat1[,23+26],main = "NadSurv[23]")
acf(mat1[,26+26],main = "NadSurv[26]")



par(mfrow=c(3,3))
acf(mat1[,1],main = "NadImm[1]")
acf(mat1[,3],main = "NadImm[3]")
acf(mat1[,7],main = "NadImm[7]")
acf(mat1[,11],main = "NadImm[11]")
acf(mat1[,14],main = "NadImm[14]")
acf(mat1[,17],main = "NadImm[17]")
acf(mat1[,20],main = "NadImm[20]")
acf(mat1[,23],main = "NadImm[23]")
acf(mat1[,26],main = "NadImm[26]")


# Trace
par(mfrow=c(3,3))
plot(mat1[,ncol(mat1)-6],type='l',xlab="", ylab="",sub = "v[1]")
plot(mat1[,ncol(mat1)-5],type='l',xlab="", ylab="",sub = "v[2]")
plot(mat1[,ncol(mat1)-4],type='l',xlab="", ylab="",sub = "v[3]")
plot(mat1[,ncol(mat1)-3],type='l',xlab="", ylab="",sub = "v[4]")
plot(mat1[,ncol(mat1)-2],type='l',xlab="", ylab="",sub = "v[5]")
plot(mat1[,ncol(mat1)-1],type='l',xlab="", ylab="",sub = "v[6]")
plot(mat1[,ncol(mat1)],type='l',xlab="", ylab="",sub = "v[7]")


par(mfrow=c(3,3))
plot(mat1[,1+26],type='l',xlab="", ylab="",sub = "NadSurv[1]")
plot(mat1[,3+26],type='l',xlab="", ylab="",sub = "NadSurv[3]")
plot(mat1[,7+26],type='l',xlab="", ylab="",sub = "NadSurv[7]")
plot(mat1[,11+26],type='l',xlab="", ylab="",sub = "NadSurv[11]")
plot(mat1[,14+26],type='l',xlab="", ylab="",sub = "NadSurv[14]")
plot(mat1[,17+26],type='l',xlab="", ylab="",sub = "NadSurv[17]")
plot(mat1[,20+26],type='l',xlab="", ylab="",sub = "NadSurv[20]")
plot(mat1[,23+26],type='l',xlab="", ylab="",sub = "NadSurv[23]")
plot(mat1[,26+26],type='l',xlab="", ylab="",sub = "NadSurv[26]")



par(mfrow=c(3,3))
plot(mat1[,1],type='l',xlab="", ylab="",sub = "NadImm[1]")
plot(mat1[,3],type='l',xlab="", ylab="",sub = "NadImm[3]")
plot(mat1[,7],type='l',xlab="", ylab="",sub = "NadImm[7]")
plot(mat1[,11],type='l',xlab="", ylab="",sub = "NadImm[11]")
plot(mat1[,14],type='l',xlab="", ylab="",sub = "NadImm[14]")
plot(mat1[,17],type='l',xlab="", ylab="",sub = "NadImm[17]")
plot(mat1[,20],type='l',xlab="", ylab="",sub = "NadImm[20]")
plot(mat1[,23],type='l',xlab="", ylab="",sub = "NadImm[23]")
plot(mat1[,26],type='l',xlab="", ylab="",sub = "NadImm[26]")





# Q matrix

par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
for (it in c(1,5,10,20,50,70,80,90,100 )){
  t=1
  plot(mat1[it,(30*(t-1)+1):(t*30)],type='l',col="red", xlab ="", ylab="diag(Q)", sub=paste('MCMC iter: ',toString(it),sep=""))
  for (t in 2:(ti-1)){
    lines(mat1[it,(30*(t-1)+1):(t*30)],type='l', xlab ="", ylab="diag(Q)")
  }
  t=ti
  lines(mat1[it,(30*(t-1)+1):(t*30)],type='l',col="blue", xlab ="", ylab="diag(Q)")
}