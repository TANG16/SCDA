T=2000
N_bin = sv_model_HMM$data()$N_bin
N_q = sv_model_HMM_adapt$data()$N_q

plot(ESS1_DA[seq(2,T,by=2)],ylim=c(0,1000),type='l',col='red',xlab="",ylab="",
     sub=paste("ESS N_bin=",toString(N_bin),", Nq=",toString(N_q),sep=""))
lines(ESS2_DA[seq(2,T,by=2)],type='l',col='darkred')
lines(ESS1_HMM[c(1:T/2)],type='l',col='green')
lines(ESS2_HMM[c(1:T/2)],type='l',col='darkgreen')
lines(ESS1_HMM_adapt[c(1:T/2)],type='l',col='blue')
lines(ESS2_HMM_adapt[c(1:T/2)],type='l',col='darkblue')
legend(0,1000,legend = c("DA","HMM","HMM_adapt"),lty =1,col=c("red","green","blue"))

par(mfrow=c(1,1))
plot(h_true,type='l',col='black',xlab="",ylab="",
     sub=paste("posterior mean h N_bin=",toString(N_bin),", Nq=",toString(N_q),sep=""))
plot(mean_H1_DA,type='l',col='red',xlab="",ylab="")
lines(mean_H2_DA,type='l',col='darkred')
lines(seq(2,T,by=2),mean_H1_HMM,type='l',col='green')
lines(seq(2,T,by=2),mean_H2_HMM,type='l',col='darkgreen')
lines(seq(2,T,by=2),mean_H1_HMM_adapt,type='l',col='blue',xlab="",ylab="")
lines(seq(2,T,by=2),mean_H2_HMM_adapt,type='l',col='darkblue')
legend(0,-3,legend = c("true","DA","HMM","HMM_adapt"),lty =1,col=c("black","red","green","blue"))

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



par(mfrow=c(3,1))
for (i in 1:3){
  plot(theta1_DA[,i],type= 'l',xlab='',ylab='',col='blue')
  lines(theta1_DA[,i],type= 'l',xlab='',ylab='',col='darkblue')
}
 