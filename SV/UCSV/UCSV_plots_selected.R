# mat1_DA = as.matrix(output_ucsv_DA[1]) 
# mat2_DA = as.matrix(output_ucsv_DA[2]) 
# mat_names_DA <- colnames(mat1_DA)
# 
# ESS_DA = lapply(output_ucsv_DA,effectiveSize)
# ESS1_DA = as.matrix(ESS_DA[[1]])
# ESS2_DA = as.matrix(ESS_DA[[2]])
# 
# mat_names_DA[1:T] # g[1] : g[258]
# mat_names_DA[T+1] # g0
# mat_names_DA[(T+2):(2*T+1)] # h[1] : h[258]
# mat_names_DA[(2*T+2)] # h0
# mat_names_DA[(2*T+2+1)]  # omega_g
# mat_names_DA[(2*T+2+2)]  # omega_h
# mat_names_DA[(2*T+2+1):(3*T+4)] # tau[1] : tau[258]
# 
# theta1_DA = mat1_DA[,c(T+1,2*T+2,2*T+2+1,2*T+2+2)] 
# theta2_DA = mat2_DA[,c(T+1,2*T+2,2*T+2+1,2*T+2+2)]
# 
# tau1_sel_DA = mat1_DA[,seq((2*T+2+1),(3*T+4),by=50)]
# tau2_sel_DA = mat2_DA[,seq((2*T+2+1),(3*T+4),by=50)]
# 
# mean_tau1_DA = colMeans(mat1_DA[,(2*T+2+1):(3*T+4)])
# mean_tau2_DA = colMeans(mat2_DA[,(2*T+2+1):(3*T+4)])
# 
# 
# g1_sel_DA = mat1_DA[,seq(1,T,by=50)]
# g2_sel_DA = mat2_DA[,seq(1,T,by=50)]
# 
# mean_g1_DA = colMeans(mat1_DA[,1:T])
# mean_g2_DA = colMeans(mat2_DA[,1:T])
# 
# 
# h1_sel_DA = mat1_DA[,seq((T+2),(2*T+1),by=50)]
# h2_sel_DA = mat2_DA[,seq((T+2),(2*T+1),by=50)]
# 
# mean_h1_DA = colMeans(mat1_DA[,(T+2):(2*T+1)])
# mean_h2_DA = colMeans(mat2_DA[,(T+2):(2*T+1)])



T=258
# N_bin = sv_model_HMM$data()$N_bin
# N_q = sv_model_HMM_adapt$data()$N_q

plot(ESS1_DA[seq(2,T,by=2)],ylim=c(0,1000),type='l',col='red',xlab="",ylab="",
     sub=paste("ESS N_bin=",toString(N_bin),", Nq=",toString(N_q),sep=""))
lines(ESS2_DA[seq(2,T,by=2)],type='l',col='darkred')
# lines(ESS1_HMM[c(1:T/2)],type='l',col='green')
# lines(ESS2_HMM[c(1:T/2)],type='l',col='darkgreen')
# lines(ESS1_HMM_adapt[c(1:T/2)],type='l',col='blue')
# lines(ESS2_HMM_adapt[c(1:T/2)],type='l',col='darkblue')
# legend(0,1000,legend = c("DA","HMM","HMM_adapt"),lty =1,col=c("red","green","blue"))

par(mfrow=c(3,1))
plot(mean_tau1_DA,type='l',col='blue',xlab="",ylab="")
lines(mean_tau2_DA,type='l',col='darkblue')

plot(mean_h1_DA,type='l',col='red',xlab="",ylab="")
lines(mean_h2_DA,type='l',col='darkred')

plot(mean_g1_DA,type='l',col='green',xlab="",ylab="")
lines(mean_g2_DA,type='l',col='darkgreen')

write.table(mean_g1_DA,file='mean_g1_DA.csv',sep=",", col.names = FALSE,row.names = FALSE)
write.table(mean_g2_DA,file='mean_g2_DA.csv',sep=",", col.names = FALSE,row.names = FALSE)
write.table(mean_h1_DA,file='mean_h1_DA.csv',sep=",", col.names = FALSE,row.names = FALSE)
write.table(mean_h2_DA,file='mean_h2_DA.csv',sep=",", col.names = FALSE,row.names = FALSE)
write.table(mean_tau1_DA,file='mean_tau1_DA.csv',sep=",", col.names = FALSE,row.names = FALSE)
write.table(mean_tau2_DA,file='mean_tau2_DA.csv',sep=",", col.names = FALSE,row.names = FALSE)

mean_theta1_DA = colMeans(theta1_DA);
mean_theta1_DA = c(mean_theta1_DA[2],mean_theta1_DA[1],mean_theta1_DA[4],mean_theta1_DA[3])
mean_theta2_DA = colMeans(theta2_DA);
mean_theta2_DA = c(mean_theta2_DA[2],mean_theta2_DA[1],mean_theta2_DA[4],mean_theta2_DA[3])
write.table(mean_theta1_DA,file='mean_theta1_DA.csv',sep=",", col.names = FALSE,row.names = FALSE)
write.table(mean_theta2_DA,file='mean_theta2_DA.csv',sep=",", col.names = FALSE,row.names = FALSE)

par(mfrow=c(3,1))
plot(exp(mean_tau1_DA/2),type='l',col='blue',xlab="",ylab="")
lines(exp(mean_tau2_DA/2),type='l',col='darkblue')

plot(exp(mean_h1_DA/2),type='l',col='red',xlab="",ylab="")
lines(exp(mean_h2_DA/2),type='l',col='darkred')

plot(exp(mean_g1_DA/2),type='l',col='green',xlab="",ylab="")
lines(exp(mean_g2_DA/2),type='l',col='darkgreen')


colMeans(theta1_DA)
colMeans(theta2_DA)
colMeans(theta1_HMM)
colMeans(theta2_HMM)
colMeans(theta1_HMM_adapt)
colMeans(theta2_HMM_adapt)

time_sample_DA
time_sample_HMM
time_sample_HMM_adapt


time_init_DA
time_sample_DA 
mean(ESS1_DA[1:1000])
mean(ESS2_DA[1:1000])
ESS1_DA[1:3]
ESS2_DA[1:3]

time_init_HMM
time_sample_HMM 
mean(ESS1_HMM[1:1000])
mean(ESS2_HMM[1:1000])
ESS1_HMM[1:3]
ESS2_HMM[1:3]

time_init_HMM_adapt
time_sample_HMM_adapt 
mean(ESS1_HMM_adapt[1:1000])
mean(ESS2_HMM_adapt[1:1000])
ESS1_HMM_adapt[1:3]
ESS2_HMM_adapt[1:3]


par(mfrow=c(3,1))
for (i in 1:3){
  plot(theta1_HMM_adapt[,i],type= 'l',xlab='',ylab='',col='blue')
  lines(theta2_HMM_adapt[,i],type= 'l',xlab='',ylab='',col='darkblue')
  # lines(param[i] + 0*theta1_HMM_adapt[,i],type= 'l',xlab='',ylab='',col='red')
}



par(mfrow=c(2,2))
for (i in 1:4){
  plot(theta1_DA[,i],type= 'l',xlab='',ylab='',col='blue',main=colnames(theta1_DA)[i])
  lines(theta2_DA[,i],type= 'l',xlab='',ylab='',col='darkblue')
}

par(mfrow=c(2,2))
for (i in 1:4){
  acf(theta1_DA[,i],xlab='',ylab='',lag=100,main=colnames(theta1_DA)[i])
}


par(mfrow=c(1,1))
plot(mean_tau1_DA,type= 'l',xlab='',ylab='',col='blue',main='mean tau')
lines(mean_tau2_DA,type= 'l',xlab='',ylab='',col='darkblue')
lines(y,type= 'l',xlab='',ylab='',col='red')


if (D > 0){
  M = 5
} else{
  M = 6
}

par(mfrow=c(3,M))
for (i in 1:M){
  plot(tau1_sel_DA[,i],type= 'l',xlab='',ylab='',col='blue')
  lines(tau2_sel_DA[,i],type= 'l',xlab='',ylab='',col='darkblue')
}
for (i in 1:M){
  plot(h1_sel_DA[,i],type= 'l',xlab='',ylab='',col='red')
  lines(h2_sel_DA[,i],type= 'l',xlab='',ylab='',col='darkred')
}
for (i in 1:M){
  plot(g1_sel_DA[,i],type= 'l',xlab='',ylab='',col='green')
  lines(g2_sel_DA[,i],type= 'l',xlab='',ylab='',col='darkgreen')  
}
 





par(mfrow=c(3,M))
for (i in 1:M){
  acf(tau2_sel_DA[,i] ,xlab='',ylab='',lag=100, main = paste('tau1 ',toString(seq(1,T,50)[i]),sep=''))
  # lines(tau2_sel_DA[,i],type= 'l',xlab='',ylab='')
}
for (i in 1:M){
  acf(h2_sel_DA[,i] ,xlab='',ylab='',lag=100, main = paste('h1 ',toString(seq(1,T,50)[i]),sep=''))
  # lines(h2_sel_DA[,i],type= 'l',xlab='',ylab='',col='darkred')
}
for (i in 1:M){
  acf(g2_sel_DA[,i] ,xlab='',ylab='',lag=100, main = paste('g1 ',toString(seq(1,T,50)[i]),sep=''))
  # lines(g2_sel_DA[,i],type= 'l',xlab='',ylab='',col='darkgreen')  
}




par(mfrow=c(3,1))
plot(ESS1_DA[(2*T+2+1):(3*T+4)],type= 'l',xlab='',ylab='',col='blue',main='ESS tau')
lines(ESS2_DA[(2*T+2+1):(3*T+4)],type= 'l',xlab='',ylab='',col='darkblue')

plot(ESS1_DA[(T+2):(2*T+1)],type= 'l',xlab='',ylab='',col='red',main='ESS h')
lines(ESS2_DA[(T+2):(2*T+1)],type= 'l',xlab='',ylab='',col='darkred')

plot(ESS1_DA[1:T],type= 'l',xlab='',ylab='',col='green',main='ESS g')
lines(ESS2_DA[1:T],type= 'l',xlab='',ylab='',col='darkgreen')
  





par(mfrow=c(3,1))
plot(mean_tau1_DA,type= 'l',xlab='',ylab='',col='blue',main='mean tau')
lines(mean_tau2_DA,type= 'l',xlab='',ylab='',col='darkblue')

plot(mean_h1_DA,type= 'l',xlab='',ylab='',col='red',main='mean h')
lines(mean_h2_DA,type= 'l',xlab='',ylab='',col='darkred')

plot(mean_g2_DA,type= 'l',xlab='',ylab='',col='green',main='mean g')
lines(mean_g2_DA,type= 'l',xlab='',ylab='',col='darkgreen')
