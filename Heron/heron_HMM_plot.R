# Gamma1 with 5000 iter ####
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat1[1000,1:(80*72)], nrow = 80, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma)", sub='MCMC iter: 1000')
for (t in 2:72){
  lines(matrix(mat1[1000,1:(80*72)], nrow = 80, ncol = 72, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[2500,1:(80*72)], nrow = 80, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 2500')
for (t in 2:72){
  lines(matrix(mat1[2500,1:(80*72)], nrow = 80, ncol = 72, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[4000,1:(80*72)], nrow = 80, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 4000')
for (t in 2:72){
  lines(matrix(mat1[4000,1:(80*72)], nrow = 80, ncol = 72, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[5000,1:(80*72)], nrow = 80, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 5000')
for (t in 2:72){
  lines(matrix(mat1[5000,1:(80*72)], nrow = 80, ncol = 72, byrow = FALSE)[,t],type='l')
}
mtext("Gamma 1", outer=TRUE, cex=1)

# P2 ####
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat1[1000,(80*72+50*72)+(1:(80*72))], nrow = 80, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(P2)", sub='MCMC iter: 1000')
for (t in 2:72){
  lines(matrix(mat1[1000,(80*72+50*72)+(1:(80*72))], nrow = 80, ncol = 72, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[2500,(80*72+50*72)+(1:(80*72))], nrow = 80, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P2)", sub='MCMC iter: 2500')
for (t in 2:72){
  lines(matrix(mat1[2500,(80*72+50*72)+(1:(80*72))], nrow = 80, ncol = 72, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[4000,(80*72+50*72)+(1:(80*72))], nrow = 80, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P2)", sub='MCMC iter: 4000')
for (t in 2:72){
  lines(matrix(mat1[4000,1:(80*72)], nrow = 80, ncol = 72, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[5000,(80*72+50*72)+(1:(80*72))], nrow = 80, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P2)", sub='MCMC iter: 5000')
for (t in 2:72){
  lines(matrix(mat1[5000,(80*72+50*72)+(1:(80*72))], nrow = 80, ncol = 72, byrow = FALSE)[,t],type='l')
}
mtext("P2", outer=TRUE, cex=1)



# Gamma3 with 5000 iter ####
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat1[1000,(80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma)", sub='MCMC iter: 1000')
for (t in 2:72){
  lines(matrix(mat1[1000,(80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[2500,(80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 2500')
for (t in 2:72){
  lines(matrix(mat1[2500,(80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[4000,(80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 4000')
for (t in 2:72){
  lines(matrix(mat1[4000,(80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[5000,(80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 5000')
for (t in 2:72){
  lines(matrix(mat1[5000,(80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,t],type='l')
}
mtext("Gamma 3", outer=TRUE, cex=1)




# Q with 5000 iter ####
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat1[1000,(80*72+50*72+80*72+50*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Q)", sub='MCMC iter: 1000')
for (t in 2:72){
  lines(matrix(mat1[1000,(80*72+50*72+80*72+50*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[2500,(80*72+50*72+80*72+50*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Q)", sub='MCMC iter: 2500')
for (t in 2:72){
  lines(matrix(mat1[2500,(80*72+50*72+80*72+50*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[4000,(80*72+50*72+80*72+50*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Q)", sub='MCMC iter: 4000')
for (t in 2:72){
  lines(matrix(mat1[4000,(80*72+50*72+80*72+50*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[5000,(80*72+50*72+80*72+50*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Q)", sub='MCMC iter: 5000')
for (t in 2:72){
  lines(matrix(mat1[5000,(80*72+50*72+80*72+50*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,t],type='l')
}
mtext("Q", outer=TRUE, cex=1)




# P4 with 5000 iter ####
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat1[1000,(80*72+50*72+80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(P4)", sub='MCMC iter: 1000')
for (t in 2:72){
  lines(matrix(mat1[1000,(80*72+50*72+80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[2500,(80*72+50*72+80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P4)", sub='MCMC iter: 2500')
for (t in 2:72){
  lines(matrix(mat1[2500,(80*72+50*72+80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[4000,(80*72+50*72+80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P4)", sub='MCMC iter: 4000')
for (t in 2:72){
  lines(matrix(mat1[4000,(80*72+50*72+80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[5000,(80*72+50*72+80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P4)", sub='MCMC iter: 5000')
for (t in 2:72){
  lines(matrix(mat1[5000,(80*72+50*72+80*72)+(1:(50*72))], nrow = 50, ncol = 72, byrow = FALSE)[,t],type='l')
}
mtext("P4", outer=TRUE, cex=1)
