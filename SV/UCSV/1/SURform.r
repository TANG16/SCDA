# SURform.r
SURform<-function(X){
  r = dim(X)[1];
  c = dim(X)[2];  
  idi = kronecker(matrix(1:r),matrix(rep(1,c)));
  idj = matrix(1:(r*c));
  Xout = kronecker(diag(r),matrix(rep(1,c),1)) %*% diag(c(t(X)));
  return (Xout);
}