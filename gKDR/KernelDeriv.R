#  function KernelDeriv()
#
#  Computing effective directions for regression with RKHS.
#
#  Author: Kenji Fukumizu (Institute of Statistical Mathematics)
#  Date: July 17, 2013
#  Version: 1.00
#
#
#  Author:  Kenji Fukumizu
#  Affiliation:  The Institute of Statistical Mathematics, ROIS
#  (C) Copyright:  Kenji Fukumizu
#
#-----------------------------------------------
# KernelDeriv()
#
# Arguments
#  X:  explanatory variables (input data)
#  Y:  response variables (teaching data)
#  K:  dimension of effective subspaces
#  SGX:  bandwidth (deviation) parameter in Gaussian kernel for X
#  SGY:  bandwidth (deviation) parameter in Gaussian kernel for Y
#  EPS:  regularization coefficient
#
# Return value(s)
#  B:  orthonormal column vectors (M x K)
#  t:  value of the objective function
#
# 
#-----------------------------------------------

KernelDeriv=function(X,Y,K,SGX,SGY,EPS)
{
  X=as.matrix(X); Y=as.matrix(Y)
  
  np=dim(X);N=np[1];M=np[2];  # N: data size, M: dim of X.
  
  I=diag(N);
  
  sx2=2*SGX*SGX;
  sy2=2*SGY*SGY;
  
  # Gram matrix of X
  ab=X%*%t(X);
  aa=diag(ab);
  D=matrix(aa,N,N);
  #xx=max(D + t(D) - 2*ab, matrix(0,N,N));
  xx=D + t(D) - 2*ab;
  xx=xx*(xx>matrix(0,N,N));
  Kx=exp(-xx/sx2);  
  
  # Gram matrix of Y
  ab=Y%*%t(Y);
  aa=diag(ab);
  D=matrix(aa,N,N);
  yy=D + t(D) - 2*ab;
  yy=yy*(yy>matrix(0,N,N));
  Ky=exp(-yy/sy2);  
  
  # Derivative of k(X_i, x) w.r.t. x
  Hm=NULL
  for (k in seq(M)){
    Xij=((matrix(X[,k],N,N)-t(matrix(X[,k],N,N)))/(SGX*SGX))*Kx
    Hm=cbind(Hm,Xij)
  }
  
  # compute  sum_i H(X_i)'*Kx^-1*Ky*Kx^-1*H(X_i)
  f=(solve(Kx+N*EPS*I)%*%Ky)%*%solve(Kx+N*EPS*I);
  
  HH=t(Hm)%*%Hm
  #NNM=matrix(0,N*N,M)
  #HHm=rep(list(NNM),M)
  #for (l in seq(M)){
  #  for (k in seq(M)){
  #  HHm[[l]][,k]=as.vector(HH[(N*(k-1)+1):(N*k),(N*(l-1)+1):(N*l)])
  #  }
  #}
  R=matrix(0,M,M)
  fm=matrix(as.vector(f),N*N,M)
  for (l in seq(M)){
    HHm=matrix(0,N*N,M)
    for (k in seq(M)){
      HHm[,k]=as.vector(HH[(N*(k-1)+1):(N*k),(N*(l-1)+1):(N*l)])
    }
    R[,l]=colSums(as.matrix(HHm*fm))
  }
  
  ##???
  VL=eigen(R);
  e=VL$values;
  K = min(ncol(VL$vectors), K)
  B=VL$vectors[,1:K];
  t=sum(e[1:K]);
  
  return(list(B=B,t=t))
}
