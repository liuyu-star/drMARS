drMARS.iter = function(x,y,degree = NULL,Xscale=T,plus=T,max.dim = 5,max.iter=10)
{
  x=as.matrix(x)
  p = ncol(x)
  n = nrow(x)
  
  BB0=drMARS(x,y,degree,Xscale,plus)
  df=BB0$degree
  B0=BB0$B[,1:max.dim,drop=F]
  cv0=BB0$cv
  ndir=BB0$ndir
  
  if(Xscale){
    x <- scale(x, center = TRUE, scale = FALSE)
    Xsd <- crossprod(x,x)/n
    eig <- eigen(Xsd)
    Xsd <- solve(eig$vectors %*% tcrossprod(diag(sqrt(eig$values*(eig$values>0))), eig$vectors)+diag(1e-4,p))
    x<-x %*% Xsd
  }
  
  
  cv=cv0
  dp=diag(p)
  B1 = B0
  n.iter=1;
  if(!identical(B0,dp)){
    nfail=0
    eps=1;
    Be=Bi0 = B0
    
    while ((n.iter<=max.iter)&&(eps>=1e-03)&&(nfail<=2)) 
    {
      outi = drMARS(cbind(x%*%B1, x),y,df,Xscale,plus)
      B2 = outi$B[,1:max.dim,drop=F]
      #cvi=outi$cv
      
      if(identical(B2,dp)){
        break
      }
      
      B1 = B1%*%B2[1:max.dim,,drop=F]+B2[-(1:max.dim),,drop=F]
      B1 = as.matrix(B1/matrix(sqrt(colSums(B1^2))+1.0e-5, p, max.dim, byrow=TRUE))
      eps = svd(Bi0%*%solve(t(Bi0)%*%Bi0 + diag(1.0e-5,max.dim))%*%t(Bi0) - 
                  B1%*%solve(t(B1)%*%B1+diag(1.0e-5,max.dim))%*%t(B1))$d[1]
      
      Bi0 = B1
      
      set.seed(230516)
      cvi=earth(x=x%*%B1, y=y, degree=1,nfold=10)$cv.rsq.tab[11,2]
      if(cvi>cv){
        cv=cvi
        Be = B1
        nfail=0
      }else{
        nfail=nfail+1
      }
      
      n.iter=n.iter+1
    }
    
    B1=Be
  }
  
  return(list(B0=B0, B=B1, n.iter=n.iter-1, cv=cv, degree = df))
}
