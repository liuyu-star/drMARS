library(mda)
library(earth)
source('drMARS/MARS.gradient.R')

drMARS=function(x,y,degree = NULL,Xscale=T,plus=T){
  x=as.matrix(x)
  n=nrow(x)
  p=ncol(x)
  cname=colnames(x)
  colnames(x)=c(paste0("x",seq(p)))
  
  #Xsd=diag(p)
  Xsd=x
  x <- scale(x, center = TRUE, scale = FALSE)
  if(Xscale){
    #Xmean = colMeans(x)
    #Xsd = apply(x, 2, sd)+1.0e-5
    #x = (x-matrix(Xmean, n, p, byrow=TRUE))/matrix(Xsd, n, p, byrow=TRUE)
    #Xsd=matrix(Xsd, p, p)
    
    Xsd <- crossprod(x,x)/n
    eig <- eigen(Xsd)
    Xsd <- solve(eig$vectors %*% tcrossprod(diag(sqrt(eig$values*(eig$values>0))), eig$vectors)+diag(1e-4,p))
    Xsd<-x %*% Xsd
  }
  
  qrz <- qr(x)  # QR decomp of WEIGHTED x matrix
  Q=qr.Q(qrz)[,1:qrz$rank]
  R=qr.R(qrz)[1:qrz$rank,1:qrz$rank]
  x=sqrt(n) * Q
  pp=ncol(x)
  colnames(x)=paste0("x",seq(pp))
  
  
  cv=-Inf
  sp =if(is.null(degree)) seq(min(pp,4)) else degree
  for (df in sp) {
    
    A=try(earth(x=x,y=y,degree = df),silent = T)#,minspan=50))
    if(is(A, "try-error")){
      A=mars(x=x,y=y,degree = df)#,minspan=50)
      #gcv=A$gcv
      A=mars.to.earth(A,trace=FALSE)
    }#else{
    #gcv=A$gcv
    #}
    
    
    #########################################
    B1=diag(p)
    dm=p
    if(length(A$coefficients)>1){
      dx=MARS.gradient(x,A$coefficients)$Gradient
      
      M0=crossprod(dx,dx)#t(dx)%*%dx 
      M0.eigen <- eigen(M0)##
      or <- rev(order(abs(M0.eigen$values)))
      M0.eigen$values <- M0.eigen$values[or]
      M0.eigen$vectors <- M0.eigen$vectors[,or]
      B = M0.eigen$vectors
      
      if(plus){
        M1 = M0.eigen$vectors%*%tcrossprod(diag(sqrt(M0.eigen$values*(M0.eigen$values>0))),M0.eigen$vectors)
        M0 = M1%*%crossprod(x,x)%*%M1
        M0.eigen <- eigen(tcrossprod(M0,M0))
        or <- rev(order(abs(M0.eigen$values)))
        M0.eigen$values <- M0.eigen$values[or]
        M0.eigen$vectors <- M0.eigen$vectors[,or]
        B = M1%*%M0.eigen$vectors
      }
      
      #B <- M1%*%M0.eigen$vectors
      B <- backsolve(sqrt(n)*R,B)
      B <- if (is.matrix(B)) B else matrix(B,ncol=1)
      
      #B <- B/Xsd
      #B <- Xsd%*%B
      if(nrow(B)<p){
        B1=matrix(0,p,ncol(B))
        B1[qrz$pivot[1:qrz$rank],]=B
      }else{
        B1=B
      }
      B1 = B1/matrix(sqrt(colSums(B1^2))+1.0e-5, p, ncol(B1), byrow=TRUE)
      #B <- apply(B,2,function(x) x/sqrt(sum(x^2)))
      #B = B*(abs(B)>1.0e-6)
      
      #dr.selected = which.max( eig.values[1:(max.k-1)]/(eig.values[2:max.k] + 1.0e-10))
      L = cumsum(M0.eigen$values)/sum(M0.eigen$values)
      dm = which.min(abs(L-0.95))
      dimnames(B1)<-list(cname, paste0("Dir", 1:NCOL(B1)))
    }
    
    
    #gcvi=earth(x=cbind(x%*%B1), y=y, degree=1)$gcv
    set.seed(230516)
    cvi=earth(Xsd%*%B1, y, degree=1,nfold=10)$cv.rsq.tab[11,2]
    if (cvi > cv)
    {
      cv = cvi
      B0=B1
      df0=df
      dm0=dm
    }
  }
  
  return(list(B=B0,ndir=dm0,degree=df0,cv=cv))
}
