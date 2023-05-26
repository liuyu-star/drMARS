library(earth)

drMARS.CV=function(BB, max.dim=5, nfold=10){
  
  cv=-Inf
  for (d in seq(max.dim)) {
    B0=BB[,seq(d),drop=F]
    gcvs=sapply(seq(min(3,ncol(B0))),function(df)earth(x%*%B0, y, degree=df)$gcv)
    df=which.min(gcvs)
    
    set.seed(230516)
    cvi=earth(x%*%B0, y, degree=1,nfold=nfold)$cv.rsq.tab[11,2]
    if(cvi>cv){
      cv=cvi
      B=B0
      ndir=d
    }
  }
  
  return(list(ndir=ndir,B=B,CV=cv))
}

