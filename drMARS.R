library(mda)
library(earth)
library(glmnet)
library(e1071)
library(moments)


MARS.gradient=function(x,Coef,Hessian=FALSE){
  x=as.matrix(x)
  np=dim(x);n=np[1];p=np[2]
  mnp=matrix(0,n,p)
  mpp=matrix(0,p,p)
  basicfun=rownames(Coef)
  
  if(length(Coef)==1){
    dX=mnp
    HX=rep(list(mpp),n)
    pks=c()
  }else{
    
    for (bi in 2:length(basicfun)) {
      bf=unlist(strsplit(basicfun[bi],split = ""))
      if("."%in%bf){
        idm=max(which(bf=="."))-1
        if("h"%in%bf){
          if(bf[idm]==")"){
            basicfun[bi]=substr(basicfun[bi],1,idm)
          }
        }else{
          basicfun[bi]=substr(basicfun[bi],1,idm)
        }
      }
    }
    
    h=function(t){ifelse(t>0,t,0)}
    
    xp=c(paste0("x",seq(p)))
    for (ip in seq(p)) {
      assign(xp[ip],x[,ip])
    }
    
    dX=mnp
    dXX=rep(list(mnp),p)
    pks=c()
    for (i in 2:length(Coef)) {
      hi=unlist(strsplit(basicfun[i],split = ""))
      if("*"%in%hi){
        hi=unlist(strsplit(basicfun[i],split = "*",fixed = T))
      }else{
        hi=basicfun[i]
      }
      
      hxj=mnp
      dhxj=mnp
      dhx=mnp
      dhxx=rep(list(mnp),p)
      pk=c()
      for (j in seq(length(hi))) {
        hj=unlist(strsplit(hi[j],split = ""))
        Xid=which("x"==hj)
        
        if("h"%in%hj){
          nk=ifelse(hj[Xid-1]=="-",nchar(hi[j]),min(which("-"==hj)))-1
          xj=eval(parse(text = substr(hi[j],3,nchar(hi[j])-1)))
          dhxjk=ifelse("-"==hj[Xid-1],-1,1)*(xj>0)
        }else{
          nk=nchar(hi[j])
          dhxjk=rep(1,n)
        }
        
        k=as.numeric(substr(hi[j],Xid+1,nk))
        pk=c(pk,k)
        hxj[,k]= eval(parse(text = hi[j]))
        dhxj[,k]= dhxjk
      }
      
      pk=sort(pk)
      for (k in pk) {
        dhxk=dhxj[,k]
        
        if(length(setdiff(pk,k))!=0){
          for (l in setdiff(pk,k)) {
            dhxk=dhxk*hxj[,l]
            ###########
            if(Hessian){
              dhxxk=dhxj[,k]*dhxj[,l]
              if(length(setdiff(pk,c(k,l)))!=0){
                for (kl in setdiff(pk,c(k,l))) {
                  dhxxk=dhxxk*hxj[,kl]
                }
              }
              dhxx[[k]][,l]=Coef[i]*dhxxk
            }
            ###########
          }
        }
        
        dhx[,k]=Coef[i]*dhxk
      }
      
      pks=c(pks,pk)
      dX=dX+dhx
      ##############
      if(Hessian){
        for (pi in pk) {
          dXX[[pi]]=dXX[[pi]]+dhxx[[pi]]
        }
      }
    }
    
    ###############
    HX=NULL
    if(Hessian){
      pks=sort(unique(pks))
      HX=rep(list(mpp),n)
      for (pi in pks) {
        for (pj in pks) {
          for (ni in seq(n)) {
            HX[[ni]][pi,pj]=dXX[[pi]][,pj][ni]
          }
        }
      }
    }
  }
  
  return(list(Gradient=dX,Hessian=HX,Var=sort(unique(pks))))
}


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


##################################################################################
varTF = function(X, Xnew=c(), max.kurtosis = 3)
{
  # variable transform to a kurtosis specified by max.kurtosis
  #  the values must be greater than 1.8
  
  # return X: transformed X
  #        Xnew: transformed Xnew
  
  max.kurtosis = max(max.kurtosis, 1.8)
  p = ncol(X)
  n = nrow(X)
  if (length(Xnew) == 0)
    n1 = 0
  else
    n1 = nrow(Xnew)
  
  N = n + n1
  
  X0 = rbind(X, Xnew)  
  for (i in 1:p)
  {
    x = X0[,i]
    I = order(x)
    x.t = x[I]
    
    xd = (x.t[2:N]-x.t[1:(N-1)])
    for (a in seq(1, 0.1, -0.1))
    {
      x.t = cumsum(c(x.t[1], xd^a))
      kur=moments::kurtosis(x.t)
      if(is.nan(kur)){
        next
      }
      if (kur <= max.kurtosis)
      {
        break
      }
    }
    x.t[I] = x.t
    X0[,i] = x.t
  }
  X = X0[1:n,]
  if (n1 > 0)
    Xnew = X0[(n+1):(n+n1),]
  
  return(list(X=X, Xnew=Xnew))
}


#main function
drMARS.fit = function(x,y,xnew,degree = NULL,Xadd=F,Xnorm=T,Xscale=F,plus=F,iter=F,
                      ndir=c("NoPreSel","PreSel",1)[1], max.dim=min(5,ncol(x)), max.iter=50)
{
  if(Xnorm){
    VTF=varTF(X=x, Xnew=xnew, max.kurtosis = 3)
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
  #pred.mars=pred
  gcv.mars=log(fit.mars$gcv) + log(n)*(0+nrow(fit.mars$coefficients))/n
  #gcv.mars=log(fit.mars$rss/n) + log(n)*(0+nrow(fit.mars$coefficients))/n
  
  ################################################################
  max.dim =ifelse(is.numeric(ndir),ndir,max.dim) 
  if(iter){
    BB = drMARS.iter(x,y,degree,Xscale,plus,max.dim,max.iter)$B#[,seq(D),drop=F]
  }else{
    BB = drMARS(x,y,degree,Xscale,plus)$B[,seq(max.dim),drop=F]
  }
  max.dim=ncol(BB)
  
  
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
  
  
  #Dr=ProjDimCV(x,y,BB)$d
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
  #err = mean(abs(pred-ynew)^2)#mean(1*(pred>0.5)!=ynew)#
  
  #return(list(pred.MARS=pred.mars, pred.drMARS=pred, fitted=fit0$fitted,cv.MARS = cv.mars,
  return(list(predicted=pred, fitted=fit0$fitted, 
              ndir=ndir, degree=df0, gcv = gcv0, B=B))
}

