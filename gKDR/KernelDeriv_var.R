#  function KernelDeriv_var()
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
  # KernelDeriv_var()
#
# Arguments
#  X:  explanatory variables (input data)
#  Y:  response variables (teaching data)
#  K:  dimension of effective subspaces
#  SIGX:  sigma parameter for Gaussian kernel
#  SIGY:  sigma parameter for Gaussian kernel 
#  EPS:  regularization parameter
#
# Return value(s)
#  B:  orthonormal column vectors (M x K)
#  t:  sum of eigenvalues 
#
# Description
#  This program computes the projection matrix by gradient0-based KDR method
#  The incomplete Gholesky approximation for G_Y is used for reducing the
#  required memory.
# 
#-----------------------------------------------
  
KernelDeriv_var=function(X,Y,K,SGX,SGY,EPS,NDIV=nrow(X)){
  X=as.matrix(X); Y=as.matrix(Y)
  
  np=dim(X);N=np[1];M=np[2]  # N: data size, M: dim of X.
  
  tol=0.000001;   # tolerance for incomplete cholesky approximation
  
  ridx=sample(N); 
  X=X[ridx,];
  Y=Y[ridx,];
  
  sx2=2*SGX*SGX;
  
  # Gram matrix of X
  ab=X%*%t(X);
  aa=diag(ab);
  D=matrix(aa,N,N);
  xx=D + t(D) - 2*ab;
  xx=xx*(xx>matrix(0,N,N));
  Kx=exp(-xx/sx2);  
  
  
  # incomplete cholesky approximation of Ky
  GP=chol_inc_gauss(t(Y),SGY,tol);
  Pvec=order(GP$Pvec);
  G=GP$G;
  Ry=G[Pvec,,drop=F];
  r=length(Ry[1,]);
  Ty=t(Ry)%*%solve(Kx+diag(N*EPS,N));
  rm(Ry);
  
  # random partition of data
  lx=floor(N/NDIV);
  ei=cumsum(rep(lx,NDIV));
  si=ei-rep(lx-1,NDIV);
  ei[NDIV]=N;       
  # si: staring idx, ei: ending idx
  
  
  Proj=matrix(0,M,M);
  for(i in 1:NDIV){
    # Derivative of k(X_i, x) w.r.t. x
    pi=si[i]:ei[i];
    Xp=X[pi,];
    Lp=length(pi);
    
    H=NULL
    for (k in seq(M)){
      xia=matrix(Xp[,k],N,Lp,byrow = T)
      xja=matrix(X[,k],N,Lp)
      H=cbind(H,((xia-xja)/(SGX*SGX))*Kx[,pi]); # N x Lp x M
    }
    
    # compute the matrix for gKDR
    Hy=matrix(Ty%*%H,r*Lp,M);
    R=t(Hy)%*%Hy;
    rm(Hy);
    
    # compute the first K eigenvectors
    BL=eigen(R)
    #L=BL$values[seq(K)]
    K = min(ncol(BL$vectors), K)
    B=BL$vectors[,seq(K)]
    
    Proj=Proj+B%*%t(B);
  }  
  
  BL=eigen(Proj);
  #L=BL$values[seq(K)]
  K = min(ncol(BL$vectors), K)
  B=BL$vectors[,seq(K)]
  
  return(B)
}

  