#graphics
dote=function(x,y) points(x,y,col="green",pch=19,cex=1)
mote=function(x,y,z,w) lines(c(x,z),c(y,w),col="gold",lwd=0.5)

cst=dbeta(0.5,0.5,0.5)*0.5 #normalising constant

#inverting f(x)=d, 2nd degree equation
hitden=function(d) 0.5+0.5*sqrt(1-4*( cst/ max(d,dbeta(0.5,0.5,0.5)))^2)*c(-1,1)
#output
curve(dbeta(x,.5,.5),0,1,ylab="density",lwd=2,col="steelblue",n=1001)
x=runif(1);
u=runif(1)*dbeta(x,0.5,0.5);
dote(x,u)
for (t in 1:100){ #100 slice steps
  bo=hitden(u)
  nx=sample(c(runif(1,0,bo[1]),runif(1,bo[2],1)),1)
  nu=runif(1)*dbeta(nx,0.5,0.5)
  mote(x,u,nx,nu)
  x=nx;
  u=nu;
  dote(x,u)
}