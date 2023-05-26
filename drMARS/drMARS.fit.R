library(earth)
library(e1071)
source('drMARS/drMARS.R')
source('drMARS/drMARS.CV.R')
source('drMARS/drMARS.iter.R')
source('drMARS/VarTF.R')

drMARS.fit = function(x,y,xnew,degree = NULL,Xadd=F,Xnorm=T,Xscale=F,plus=F,iter=F,
                      ndir=c("NoPreSel","PreSel",1)[1], max.dim=min(5,ncol(x)), max.iter=50)
{
  if(Xnorm){
    VTF=VarTF(X=x, Xnew=xnew, max.kurtosis = 3)
    x=VTF$X
    xnew=VTF$Xnew
  }
  
  p = ncol(x)
  n = nrow(x)
  N = nrow(xnew)
  
  
  smp =if(is.null(degree)) seq(min(p,4)) else degree
  gcv.mars=Inf
  for (df in smp) {
    #set.seed(230516)
    #fit0=earth(x, y, degree=d,  nfold=10)
    fiti=try(earth(x, y, degree=df),silent = T)#, nfold=10
    if(is(fiti, "try-error")){
      gcvi=Inf
    }else{
      gcvi=fiti$gcv# cv.rsq.tab[11,2]
    }
    if(gcvi<gcv.mars){
      fit.mars=fiti
      gcv.mars=gcvi
      df.mars=df
    }
  }
  pred = predict(fit.mars, xnew)
  pred[which(pred<min(fit.mars$fitted))]=min(fit.mars$fitted)
  pred[which(pred>max(fit.mars$fitted))]=max(fit.mars$fitted)
  gcv.mars=log(fit.mars$gcv) + log(n)*(0+nrow(fit.mars$coefficients))/n
  #gcv.mars=log(fit.mars$rss/n) + log(n)*(0+nrow(fit.mars$coefficients))/n
  
  ################################################################
  max.dim =ifelse(is.numeric(ndir),ndir,max.dim) 
  if(iter){
    BB = drMARS.iter(x,y,degree=NULL,Xscale,plus,max.dim,max.iter)$B#[,seq(D),drop=F]
  }else{
    BB = drMARS(x,y,degree=NULL,Xscale,plus)$B[,seq(max.dim),drop=F]
  }
  #max.dim=ncol(BB)
  
  
  if(Xscale){
    #x <- scale(x, center = TRUE, scale = FALSE)
    Xmed=apply(x,2,median)
    x=x-matrix(Xmed,n,p,byrow = T)
    Sigma <- crossprod(x,x)/n
    eig <- eigen(Sigma)
    Sigma <- solve(eig$vectors %*% tcrossprod(diag(sqrt(eig$values*(eig$values>0))), eig$vectors)+diag(1e-4,p))
    x<-x %*% Sigma
    xnew=xnew-matrix(Xmed,N,p,byrow = T)
    xnew<-xnew %*% Sigma
    
    #maxy = max(y)
    #miny = min(y)
    #y = (y-miny)/(maxy-miny)
    #ynew = (ynew-miny)/(maxy-miny)
  }
  
  
  #dr.mean = mave(y~x%*%BB,max.dim = D, method="meanOPG")
  #Dr <-mave.dim(dr.mean)$dim.min
  
  cv=-Inf
  if(ndir=="PreSel"){
    for (d in seq(max.dim)) {
      B0=BB[,seq(d),drop=F]
      gcvs=sapply(seq(min(3,ncol(B0))),function(df)earth(x%*%B0, y, degree=df)$gcv)
      df=which.min(gcvs)
      
      set.seed(230516)
      cvi=earth(x%*%B0, y, degree=1,nfold=10)$cv.rsq.tab[11,2]
      if(cvi>cv){
        cv=cvi
        B=B0
        ndir=d
      }
    }
  }
  
  smd=if(ndir=="NoPreSel") seq(max.dim) else ndir 
  
  gcv0=Inf
  #gcv0=gcv.mars
  #fit0=fit.mars
  #df0=df.mars
  for (d in smd) {
    B0=BB[,seq(d),drop=F]
    if(Xadd){
      xB=cbind(x, x%*%B0)
    }else{
      xB=x%*%B0
    }
    xB=as.matrix(xB)
    
    smp =if(is.null(degree)) seq(min(ncol(xB),4)) else degree
    for (df in smp) {
      #set.seed(230516)
      fiti=try(earth(xB, y, degree=df),silent = T)#, nfold=10
      if(is(fiti, "try-error")){
        fiti=fit.mars
        gcvi=Inf
      }else{
        gcvi=log(fiti$gcv) + log(n)*(d*log(n)+nrow(fiti$coefficients))/n
        #gcvi=log(fiti$rss/n) + log(n)*(ndir*log(n)+nrow(fiti$coefficients))/n
      }
      
      if(gcvi<gcv0){
        fit0=fiti
        gcv0=gcvi
        df0=df
        B=B0
        ndir=d
      }
    }
  }
  
  if (gcv0 < gcv.mars){
    if(Xadd){
      xBnew=cbind(xnew, xnew%*%B) 
    }else{
      xBnew=xnew%*%B
    }
  }else{
    xBnew = xnew
    B=0
    #ndir=0
    gcv0=gcv.mars
    fit0=fit.mars
    df0=df.mars
  }
  pred = predict(fit0,xBnew)
  pred[which(pred<min(fit0$fitted))]=min(fit0$fitted)
  pred[which(pred>max(fit0$fitted))]=max(fit0$fitted)
  A = svm(x,y-fit0$fitted)#, gamma=L/p)
  pred = predict(A, xnew)+pred

  return(list(predicted=pred, fitted=fit0$fitted, 
              ndir=ndir, degree=df0, gcv = gcv0, B=B))
}
