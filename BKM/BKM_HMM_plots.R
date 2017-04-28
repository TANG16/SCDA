load("BKM_HMM_try_iter3000_ada1000_laptop.RData")

cha = length(output1)
ss = dim(output1[[1]])
iter = ss[1]
vars = ss[2]

mat1_names <- colnames(mat1) 

mat1 = as.matrix(output1[1])
mat2 = as.matrix(output1[2])
mat3 = as.matrix(output1[3])



# Check the transition probabilities ####
Gamma_last = matrix(mat1[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)
Gamma_last_10 = diag(Gamma_last[,12])
sum(Gamma_last_10) #1
sum_Gamma_last = colSums(Gamma_last)

# 1000 iter ####
# Chain 1 with 1000 iter
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat1[100,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma)", sub='MCMC iter: 100')
for (t in 4:34){
  lines(matrix(mat1[100,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 500')
for (t in 4:34){
  lines(matrix(mat1[500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[900,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 900')
for (t in 4:34){
  lines(matrix(mat1[900,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 1000')
for (t in 4:34){
  lines(matrix(mat1[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}
mtext("Chain 1", outer=TRUE, cex=1)

# Chain 2 with 1000 iter
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat2[100,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma)", sub='MCMC iter: 100')
for (t in 4:34){
  lines(matrix(mat2[100,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat2[500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 500')
for (t in 4:34){
  lines(matrix(mat2[500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat2[900,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 900')
for (t in 4:34){
  lines(matrix(mat2[900,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat2[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 1000')
for (t in 4:34){
  lines(matrix(mat2[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}
mtext("Chain 2", outer=TRUE, cex=1)


# 2000 iter ####
# Chain 1 with 2000 iter
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat1[500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma)", sub='MCMC iter: 500')
for (t in 4:34){
  lines(matrix(mat1[500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 1000')
for (t in 4:34){
  lines(matrix(mat1[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[1500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 1500')
for (t in 4:34){
  lines(matrix(mat1[1500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[2000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 2000')
for (t in 4:34){
  lines(matrix(mat1[2000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}
mtext("Chain 1", outer=TRUE, cex=1)


# Chain 2 with 2000 iter
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat2[500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma)", sub='MCMC iter: 500')
for (t in 4:34){
  lines(matrix(mat2[500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat2[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 1000')
for (t in 4:34){
  lines(matrix(mat2[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat2[1500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 1500')
for (t in 4:34){
  lines(matrix(mat2[1500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat2[2000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 2000')
for (t in 4:34){
  lines(matrix(mat2[2000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}
mtext("Chain 2", outer=TRUE, cex=1)


# Chain 3 with 2000 iter
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat3[500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma)", sub='MCMC iter: 500')
for (t in 4:34){
  lines(matrix(mat3[500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat3[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 1000')
for (t in 4:34){
  lines(matrix(mat3[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat3[1500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 1500')
for (t in 4:34){
  lines(matrix(mat3[1500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat3[2000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 2000')
for (t in 4:34){
  lines(matrix(mat3[2000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}
mtext("Chain 3", outer=TRUE, cex=1)



# 3000 iter ####
# Chain 1 with 3000 iter
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat1[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma)", sub='MCMC iter: 1000')
for (t in 4:34){
  lines(matrix(mat1[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[2000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 2000')
for (t in 4:34){
  lines(matrix(mat1[2000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[2500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 2500')
for (t in 4:34){
  lines(matrix(mat1[2500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[3000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 3000')
for (t in 4:34){
  lines(matrix(mat1[3000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}
mtext("Chain 1", outer=TRUE, cex=1)


# Chain 2 with 3000 iter
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat2[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma)", sub='MCMC iter: 1000')
for (t in 4:34){
  lines(matrix(mat2[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat2[2000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 2000')
for (t in 4:34){
  lines(matrix(mat2[2000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat2[2500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 2500')
for (t in 4:34){
  lines(matrix(mat2[2500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat2[3000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 3000')
for (t in 4:34){
  lines(matrix(mat2[3000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}
mtext("Chain 2", outer=TRUE, cex=1)


# Chain 3 with 3000 iter
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat3[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma)", sub='MCMC iter: 1000')
for (t in 4:34){
  lines(matrix(mat3[1000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat3[2000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 2000')
for (t in 4:34){
  lines(matrix(mat3[2000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat3[2500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 2500')
for (t in 4:34){
  lines(matrix(mat3[2500,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat3[3000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 3000')
for (t in 4:34){
  lines(matrix(mat3[3000,1:(34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}
mtext("Chain 3", outer=TRUE, cex=1)


# Check the augmented observation matrix  ####
##### 1000 iter
# Chain 1 with 1000 iter
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat1[100,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(P)", sub='MCMC iter: 100')
for (t in 4:34){
  lines(matrix(mat1[100,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[500,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P)", sub='MCMC iter: 500')
for (t in 4:34){
  lines(matrix(mat1[500,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[900,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P)", sub='MCMC iter: 900')
for (t in 4:34){
  lines(matrix(mat1[900,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[1000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P)", sub='MCMC iter: 1000')
for (t in 4:34){
  lines(matrix(mat1[1000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}
mtext("Chain 1", outer=TRUE, cex=1)


# Chain 2 with 1000 iter
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat2[100,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(P)", sub='MCMC iter: 100')
for (t in 4:34){
  lines(matrix(mat2[100,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat2[500,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P)", sub='MCMC iter: 500')
for (t in 4:34){
  lines(matrix(mat2[500,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat2[900,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P)", sub='MCMC iter: 900')
for (t in 4:34){
  lines(matrix(mat2[900,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat2[1000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P)", sub='MCMC iter: 1000')
for (t in 4:34){
  lines(matrix(mat2[1000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}
mtext("Chain 2", outer=TRUE, cex=1)


##### 3000 iter
# Chain 1 with 3000 iter
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat1[1000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(P)", sub='MCMC iter: 1000')
for (t in 4:34){
  lines(matrix(mat1[1000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[2000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P)", sub='MCMC iter: 200')
for (t in 4:34){
  lines(matrix(mat1[2000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[2500,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P)", sub='MCMC iter: 2500')
for (t in 4:34){
  lines(matrix(mat1[2500,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[3000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P)", sub='MCMC iter: 3000')
for (t in 4:34){
  lines(matrix(mat1[3000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}
mtext("Chain 1", outer=TRUE, cex=1)

# Chain 2 with 3000 iter
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat2[1000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(P)", sub='MCMC iter: 1000')
for (t in 4:34){
  lines(matrix(mat2[1000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat2[2000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P)", sub='MCMC iter: 200')
for (t in 4:34){
  lines(matrix(mat2[2000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat2[2500,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P)", sub='MCMC iter: 2500')
for (t in 4:34){
  lines(matrix(mat2[2500,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat2[3000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P)", sub='MCMC iter: 3000')
for (t in 4:34){
  lines(matrix(mat2[3000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}
mtext("Chain 2", outer=TRUE, cex=1)


# Chain 3 with 3000 iter
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat3[1000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(P)", sub='MCMC iter: 1000')
for (t in 4:34){
  lines(matrix(mat3[1000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat3[2000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P)", sub='MCMC iter: 200')
for (t in 4:34){
  lines(matrix(mat3[2000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat3[2500,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P)", sub='MCMC iter: 2500')
for (t in 4:34){
  lines(matrix(mat3[2500,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat3[3000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(P)", sub='MCMC iter: 3000')
for (t in 4:34){
  lines(matrix(mat3[3000,(34*100+1+35+1):(34*100+1+35+34*100)], nrow = 100, ncol = 34, byrow = FALSE)[,t],type='l')
}
mtext("Chain 3", outer=TRUE, cex=1)




# Trace plots ####
# Chain 1
par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
plot(mat1[,"sigy"], type="l", xlab ="", ylab="", sub="sigy")
plot(mat1[,"alphar"], type="l", xlab ="", ylab="", sub="alpha r")
plot(mat1[,"betar"], type="l", xlab ="", ylab="", sub="beta r")

plot(mat1[,"Na[13]"], type="l", xlab ="", ylab="", sub="Na[13]")
plot(mat1[,"Na[23]"], type="l", xlab ="", ylab="", sub="Na[23]")
plot(mat1[,"Na[33]"], type="l", xlab ="", ylab="", sub="Na[33]")
# transition probability to N1=20
plot(mat1[,"G[20,10]"], type="l", xlab ="", ylab="", sub="Gamma[20,10], k=20,t=10")
plot(mat1[,"G[20,20]"], type="l", xlab ="", ylab="", sub="Gamma[20,20], k=20,t=20")
plot(mat1[,"G[20,30]"], type="l", xlab ="", ylab="", sub="Gamma[20,30], k=20,t=30")
mtext("Chain 1", outer=TRUE, cex=1)

# Chain 2
par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
plot(mat2[,"sigy"], type="l", xlab ="", ylab="", sub="sigy")
plot(mat2[,"alphar"], type="l", xlab ="", ylab="", sub="alpha r")
plot(mat2[,"betar"], type="l", xlab ="", ylab="", sub="beta r")

plot(mat2[,"Na[13]"], type="l", xlab ="", ylab="", sub="Na[13]")
plot(mat2[,"Na[23]"], type="l", xlab ="", ylab="", sub="Na[23]")
plot(mat2[,"Na[33]"], type="l", xlab ="", ylab="", sub="Na[33]")
# transition probability to N1=20
plot(mat2[,"G[20,10]"], type="l", xlab ="", ylab="", sub="Gamma[20,10], k=20,t=10")
plot(mat2[,"G[20,20]"], type="l", xlab ="", ylab="", sub="Gamma[20,20], k=20,t=20")
plot(mat2[,"G[20,30]"], type="l", xlab ="", ylab="", sub="Gamma[20,30], k=20,t=30")
mtext("Chain 2", outer=TRUE, cex=1)

# Chain 3
par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
plot(mat3[,"sigy"], type="l", xlab ="", ylab="", sub="sigy")
plot(mat3[,"alphar"], type="l", xlab ="", ylab="", sub="alpha r")
plot(mat3[,"betar"], type="l", xlab ="", ylab="", sub="beta r")

plot(mat3[,"Na[13]"], type="l", xlab ="", ylab="", sub="Na[13]")
plot(mat3[,"Na[23]"], type="l", xlab ="", ylab="", sub="Na[23]")
plot(mat3[,"Na[33]"], type="l", xlab ="", ylab="", sub="Na[33]")
# transition probability to N1=20
plot(mat3[,"G[20,10]"], type="l", xlab ="", ylab="", sub="Gamma[20,10], k=20,t=10")
plot(mat3[,"G[20,20]"], type="l", xlab ="", ylab="", sub="Gamma[20,20], k=20,t=20")
plot(mat3[,"G[20,30]"], type="l", xlab ="", ylab="", sub="Gamma[20,30], k=20,t=30")
mtext("Chain 3", outer=TRUE, cex=1)




# Chain 2
par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
plot(mat1[,"Na[1]"], type="l", xlab ="", ylab="", sub="Na[1] ch1")
plot(mat1[,"Na[2]"], type="l", xlab ="", ylab="", sub="Na[2] ch1")
plot(mat1[,"Na[3]"], type="l", xlab ="", ylab="", sub="Na[3] ch1")

plot(mat2[,"Na[1]"], type="l", xlab ="", ylab="", sub="Na[1] ch2")
plot(mat2[,"Na[2]"], type="l", xlab ="", ylab="", sub="Na[2] ch2")
plot(mat2[,"Na[3]"], type="l", xlab ="", ylab="", sub="Na[3] ch2")

plot(mat3[,"Na[1]"], type="l", xlab ="", ylab="", sub="Na[1] ch3")
plot(mat3[,"Na[2]"], type="l", xlab ="", ylab="", sub="Na[2] ch3")
plot(mat3[,"Na[3]"], type="l", xlab ="", ylab="", sub="Na[3] ch3")
#####


### ACF ##### burn in 500 ?
par(mfrow=c(3,1))
acf(output1[[1]][501:2000,"Na[3]"], main="Na[3], Chain 1")
acf(output1[[2]][501:2000,"Na[3]"], main="Na[3], Chain 2")
acf(output1[[3]][501:2000,"Na[3]"], main="Na[3], Chain 3")

par(mfrow=c(3,1))
acf(output1[[1]][501:2000,"Na[13]"], main="Na[13], Chain 1")
acf(output1[[2]][501:2000,"Na[13]"], main="Na[13], Chain 2")
acf(output1[[3]][501:2000,"Na[13]"], main="Na[13], Chain 3")

par(mfrow=c(3,1))
acf(output1[[1]][501:2000,"Na[23]"], main="Na[23], Chain 1")
acf(output1[[2]][501:2000,"Na[23]"], main="Na[23], Chain 2")
acf(output1[[3]][501:2000,"Na[23]"], main="Na[23], Chain 3")


par(mfrow=c(3,1))
acf(output1[[1]][501:2000,"Na[33]"], main="Na[33], Chain 1")
acf(output1[[2]][501:2000,"Na[33]"], main="Na[33], Chain 2")
acf(output1[[3]][501:2000,"Na[33]"], main="Na[33], Chain 3")


#
par(mfrow=c(4,3))
acf(output1[[1]][501:2000,"Na[3]"], main="Na[3], Chain 1")
acf(output1[[2]][501:2000,"Na[3]"], main="Na[3], Chain 2")
acf(output1[[3]][501:2000,"Na[3]"], main="Na[3], Chain 3")

acf(output1[[1]][501:2000,"Na[13]"], main="Na[13], Chain 1")
acf(output1[[2]][501:2000,"Na[13]"], main="Na[13], Chain 2")
acf(output1[[3]][501:2000,"Na[13]"], main="Na[13], Chain 3")

acf(output1[[1]][501:2000,"Na[23]"], main="Na[23], Chain 1")
acf(output1[[2]][501:2000,"Na[23]"], main="Na[23], Chain 2")
acf(output1[[3]][501:2000,"Na[23]"], main="Na[23], Chain 3")

acf(output1[[1]][501:2000,"Na[33]"], main="Na[33], Chain 1")
acf(output1[[2]][501:2000,"Na[33]"], main="Na[33], Chain 2")
acf(output1[[3]][501:2000,"Na[33]"], main="Na[33], Chain 3")


#
par(mfrow=c(4,3))
acf(output1[[1]][501:2000,"Na[19]"], main="Na[19], Chain 1")
acf(output1[[2]][501:2000,"Na[19]"], main="Na[19], Chain 2")
acf(output1[[3]][501:2000,"Na[19]"], main="Na[19], Chain 3")

acf(output1[[1]][501:2000,"Na[21]"], main="Na[21], Chain 1")
acf(output1[[2]][501:2000,"Na[21]"], main="Na[21], Chain 2")
acf(output1[[3]][501:2000,"Na[21]"], main="Na[21], Chain 3")

acf(output1[[1]][501:2000,"Na[23]"], main="Na[23], Chain 1")
acf(output1[[2]][501:2000,"Na[23]"], main="Na[23], Chain 2")
acf(output1[[3]][501:2000,"Na[23]"], main="Na[23], Chain 3")


#
par(mfrow=c(4,3))
acf(output1[[1]][501:2000,"alphar"], main="alphar, Chain 1")
acf(output1[[2]][501:2000,"alphar"], main="alphar, Chain 2")
acf(output1[[3]][501:2000,"alphar"], main="alphar, Chain 3")

acf(output1[[1]][501:2000,"betar"], main="betar, Chain 1")
acf(output1[[2]][501:2000,"betar"], main="betar, Chain 2")
acf(output1[[3]][501:2000,"betar"], main="betar, Chain 3")

acf(output1[[1]][501:2000,"alphal"], main="alphal, Chain 1")
acf(output1[[2]][501:2000,"alphal"], main="alphal, Chain 2")
acf(output1[[3]][501:2000,"alphal"], main="alphal, Chain 3")

acf(output1[[1]][501:2000,"betal"], main="betal, Chain 1")
acf(output1[[2]][501:2000,"betal"], main="betal, Chain 2")
acf(output1[[3]][501:2000,"betal"], main="betal, Chain 3")

#
par(mfrow=c(4,3))
acf(output1[[1]][501:2000,"alpha1"], main="alpha1, Chain 1")
acf(output1[[2]][501:2000,"alpha1"], main="alpha1, Chain 2")
acf(output1[[3]][501:2000,"alpha1"], main="alpha1, Chain 3")

acf(output1[[1]][501:2000,"beta1"], main="beta1, Chain 1")
acf(output1[[2]][501:2000,"beta1"], main="beta1, Chain 2")
acf(output1[[3]][501:2000,"beta1"], main="beta1, Chain 3")

acf(output1[[1]][501:2000,"alphaa"], main="alphaa, Chain 1")
acf(output1[[2]][501:2000,"alphaa"], main="alphaa, Chain 2")
acf(output1[[3]][501:2000,"alphaa"], main="alphaa, Chain 3")

acf(output1[[1]][501:2000,"betaa"], main="betaa, Chain 1")
acf(output1[[2]][501:2000,"betaa"], main="betaa, Chain 2")
acf(output1[[3]][501:2000,"betaa"], main="betaa, Chain 3")
