
MARS.component=function(Coef,x){
  
  if(length(Coef)==1){
    return(list(fxFun=NULL,fxValue=NULL))
  }else{
    
    basicfun=rownames(Coef)
    bx=c(1:ncol(x))
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
    
    #pks=c()
    gb=as.list(rep(0,length(bx)))
    names(gb)=paste0("f(x",seq(length(bx)),")")
    #gx=Coef[1]
    for (i in 2:length(Coef)) {
      hi=unlist(strsplit(basicfun[i],split = ""))
      if("*"%in%hi){
        hi=unlist(strsplit(basicfun[i],split = "*",fixed = T))
      }else{
        hi=basicfun[i]
      }
      
      #pk=c()
      k=c()
      for (j in seq(length(hi))) {
        hj=unlist(strsplit(hi[j],split = ""))
        xid=which("x"==hj)
        
        if("h"%in%hj){
          nk=ifelse(hj[xid-1]=="-",nchar(hi[j]),min(which("-"==hj)))-1
        }else{
          nk=nchar(hi[j])
        }
        
        k=c(k,as.numeric(substr(hi[j],xid+1,nk)))
      }
      
      bf=paste(Coef[i],basicfun[i],sep = "*")
      if(length(k)>1){
        gb=c(gb,bf)
        names(gb)[length(gb)]=paste0("f(x",as.character(k)[1],as.character(k)[2],")") 
      }else{
        for (bi in bx) {
          if(k==bi){
            gb[bi]=paste(gb[bi],bf,sep = "+")
          }
        }
      }
    }
    
    
    ######################################
    
    h=function(t){ifelse(t>0,t,0)}
    
    xp=c(paste0("x",seq(ncol(x))))
    for (ip in seq(ncol(x))) {
      assign(xp[ip],x[,ip])
    }
    
    fx=matrix(0,nrow(x),length(gb)) 
    colnames(fx)=names(gb)#c(paste0("x",seq(ncol(B))))
    for (bi in seq(length(gb))) {
      fx[,bi]=eval(parse(text = gb[[bi]]))
    }
    #fx[,ncol(B)+1]=eval(parse(text = gx))
    #fxx=fx 
    
    return(list(fxFun=gb,fxValue=fx))
  }
}
