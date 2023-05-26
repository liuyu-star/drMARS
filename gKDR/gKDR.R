#  gKDR_sample
#
#  Gradient-based kernel dimension redeuction (gKDR) sample program
#
#  Author: Kenji Fukumizu
#  Affiliation: The Institute of Statistical Mathematics, ROIS
#  Date: July 17, 2013
#  Version: 1.00
#
#  (C) Copyright:  Kenji Fukumizu
#  
#  This is a R code for gKDR, gKDR-i, gKDR-v.  
#  For the details of these algorithms, see the follwoing paper
#   Gradient-based kernel dimension reduction for regression
#     by Kenji Fukumizu and Chenlei Leng
#

source('gKDR//MedianDist.R')
source('gKDR/KernelDeriv.R')
source('gKDR/KernelDeriv_var.R')
source('gKDR/chol_inc_gauss.R')

gKDR=function(X, Y, K, Method="gKDR"){

  X=as.matrix(X); Y=as.matrix(Y)
  
  #K=1;       
  # Dimension of Effective directions
  
  #Method='gKDR';
  # Data specification
  np=dim(X);
  M=np[2]; # dimensionality of X
  N=np[1]; # sample size
  AC=0.0; # parameter of the regressor 
  
  NCV=5;      # Number of cross-validation
  DEC=10;     # Number of decreasing dimensions for gKDR-i;
  
  B0=rep(0,M);
  B0[1]=1;      # B: Projection matrix onto the effective direction
  candx=c(0.25 ,0.5, 0.75 ,1 ,2);  # candidates for CV
  candy=c(0.25, 0.5 ,0.75, 1 ,2);

#  candx=c(0.75);  # candidates for CV
#  candy=c(0.75);
  
  eps = c(0.00001);
  # Gaussian kernels are used. 
  # Deviation parameter are chosen by CV with kNN classifier. 
  
  sgx0=MedianDist(X);   # Basic value for bandwidth for X
  sgy0=MedianDist(Y);   # Basic value for bandwidth for Y
  
  if (1 == 2){
    library(FNN)
    
    # For cross-validation
    ridx=sample(N);  # random order 
    Xr=X[ridx,];
    Yr=Y[ridx,];   
    lx=ceiling(N/NCV);
    ei=cumsum(rep(lx,NCV));
    si=ei-rep(lx-1,NCV);
    ei[NCV]=N;       # si: staring idx, ei: ending idx
    err_tbl=matrix(0,length(candx)*length(candy)*length(eps), NCV);
    
    for (h in 1:length(candx)) {
      sgx=sgx0*candx[h];
      for (k in 1:length(candy)) {
        sgy=sgy0*candy[k]; 
        for (ll in 1:length(eps)) {
          EPS = eps[ll];
          for (i in 1:NCV) {
            ri=si[i]:ei[i];
            Xe=Xr; Ye=Yr; 
            Xe=Xe[-ri,];
            Ye=Ye[-ri,];   # Xe, Ye: trainig sample for CV
            Xt=Xr[ri,];
            Yt=Yr[ri,];    # Xt, Yt: test sample for CV
            
            if(Method=='gKDR'){
              Bt=KernelDeriv(Xe,Ye,K,sgx,sgy,EPS);
              B=Bt$B;t=Bt$t;
            }
            
            if(Method=='gKDR-i'){
              if(M-K <=DEC){
                decd=rep(1,M-K);
              }else{
                dd=floor((M-K)/DEC);
                r=M-K-dd*DEC;
                decd=rep(dd,DEC);
                decd[1:r]=rep(dd+1,r);
              }
              B=diag(M);
              for (ii in 1:length(decd)) {
                Ze=Xe%*%B;
                Bp=B;
                dim=M-sum(decd[1:ii]);
                Bt=KernelDeriv(Ze,Ye,dim,sgx*sqrt(dim/M),sgy,EPS);
                B=Bt$B;t=Bt$t;
                B=Bp%*%B;
                #B=B/sqrtm(t(B)%*%B);
                SB= eigen(t(B)%*%B)
                B=B%*%(SB$vectors %*% diag(1/sqrt(SB$values)) %*% t(SB$vectors))
              } 
            }
            
            if(Method=='gKDR-v'){
              B=KernelDeriv_var(Xe,Ye,K,sgx,sgy,EPS,NDIV=min(floor(nrow(Xe)/2),50));
            }
            
            # kNN regression for CV
            nnidx=knnx.index(Xe%*%B,Xt%*%B,k=5,algo="kd_tree");
            Yo=matrix(0,length(ri),length(Y[1,]));
            
            for (j in 1:length(ri)) {
              ii=nnidx[j,];
              Yo[j,]=colSums(Ye[ii,,drop=F])/5;
            }
            
            dd=Yt-Yo;      
            err_tbl[(h-1)*length(candy)*length(eps)+(k-1)*length(eps)+ll,i]=sum(colSums(dd*dd))/length(ri);  
          } 
        }
      } 
    } 

midx=which.min(rowMeans(err_tbl));
opth=ceiling(midx/(length(candy)*length(eps)));
rr=midx-(opth-1)*length(candy)*length(eps);
optk=ceiling(rr/length(eps));
opte=optk*length(eps)-rr;
if(opte==0){
  opte=length(eps);
}

# Parameter
sgx=sgx0*candx[opth];
sgy=sgy0*candy[optk];
EPS=eps[opte];
}   
  
  opth = 3;
  optk = 3;
  #opte = 1;
  
  sgx=sgx0;# *candx(opth);
  sgy=sgy0; #*candy(optk);
  EPS=1.0000e-05;
  

  if(Method=='gKDR'){
    B=KernelDeriv(X,Y,K,sgx,sgy,EPS)$B;
  }
  
  if(Method=='gKDR-i'){
    if(M-K <=DEC){
      decd=rep(1,M-K);
    }else{
      dd=floor((M-K)/DEC);
      r=M-K-dd*DEC;
      decd=rep(dd,DEC);
      decd[1:r]=rep(dd+1,r);
    }
    B=diag(M);
    for (ii in 1:length(decd)) {
      Z=X%*%B;
      Bp=B;
      dim=M-sum(decd[1:ii]);
      B=KernelDeriv(Z,Y,dim,sgx*sqrt(dim/M),sgy,EPS)$B;
      B=Bp%*%B;
      #B=B/sqrtm(t(B)%*%B);
      SB= eigen(t(B)%*%B)
      B=B%*%(SB$vectors %*% diag(1/sqrt(SB$values)) %*% t(SB$vectors))
    } 
  }
  
  if(Method=='gKDR-v'){
    B=KernelDeriv_var(X,Y,K,sgx,sgy,EPS,NDIV=min(floor(N/2),50));
  }
  
  return(B)
}

