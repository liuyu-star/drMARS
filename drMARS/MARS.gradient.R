
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
