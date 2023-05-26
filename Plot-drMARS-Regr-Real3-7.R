rm(list = ls())

library(earth)
library(mda)
library(MAVE)
library(MASS)
library(glmnet)
library(e1071)
library(energy)
library(randomForest)
library(dr)
library(CVarE)
library(latex2exp)
library(plot3D)
library(gam)
source("drMARS.R")
source("D:/BaiduNetdiskWorkspace/DataSets/RegrDataSets2.R")



cores=40;Rep=120
#ratios = c(1/10,1/3,1/2,2/3)[3]
sn=c(500,1000,2000,4000)[c(3)]#[2]#
#datas=seq(20)[c(7,11,1,2,4)]#[seq(3)]#[c(8,12,1,2,4)]
datas=c(1,10,11,28)#c(1,5,10,11,28,37,38)[7]#c(1,10,11,28)
depend_packages = c("energy", "MASS", "MAVE", "earth", "glmnet","mda",'e1071','RSNNS','randomForest','dr')

max.df=6

wd="D:/BaiduNetdiskWorkspace"

fxx=function(coef,xB){
  
  if(length(coef)==1){
    return(list(fxFun=NULL,fxValue=NULL))
  }else{
    
    basicfun=rownames(coef)
    bx=c(1:ncol(xB))
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
    names(gb)=paste0("f(xB",seq(length(bx)),")")
    #gx=coef[1]
    for (i in 2:length(coef)) {
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
      
      bf=paste(coef[i],basicfun[i],sep = "*")
      if(length(k)>1){
        gb=c(gb,bf)
        names(gb)[length(gb)]=paste0("f(xB",as.character(k)[1],as.character(k)[2],")") 
      }else{
        for (bi in bx) {
          if(k==bi){
            gb[bi]=paste(gb[bi],bf,sep = "+")
          }
        }
      }
    }
    
    
    ######################################
    #xB = cbind(x,x%*%B)
    #xBnew = cbind(xnew,xnew%*%B)
    #xxBnew=rbind(xB,xBnew)
    #xB = x%*%B
    
    h=function(t){ifelse(t>0,t,0)}
    
    xp=c(paste0("x",seq(ncol(xB))))
    for (ip in seq(ncol(xB))) {
      assign(xp[ip],xB[,ip])
    }
    
    fx=matrix(0,nrow(xB),length(gb)) 
    colnames(fx)=names(gb)#c(paste0("xB",seq(ncol(B))))
    for (bi in seq(length(gb))) {
      fx[,bi]=eval(parse(text = gb[[bi]]))
    }
    #fx[,ncol(B)+1]=eval(parse(text = gx))
    #fxx=fx 
    
    return(list(fxFun=gb,fxValue=fx))
  }
}

#MARS data3
########################################################
set.seed(3)
dat=datas[3]
XY =  Data.Set(i.data=dat,wd="D:/BaiduNetdiskWorkspace/DataSets/Regression1")#,wd="/home/liuyu/files/DataSets/Regression1")#,wd="D:/BaiduNetdiskWorkspace/DataSets/Regression1")#,wd="/home/liuyu/files/DataSets/Regression1")#
X = XY$X0
Y = XY$y0
rm(XY)

if(1==2){
  #XT=varTF(X,max.kurtosis = 3)$X
  XT=X
  x_20=XT[,-seq(20)]
  x_20=varTF(x_20,max.kurtosis = 3)$X
  x_20=scale(x_20)
  x = cbind(XT[,seq(20)],x_20)# scale(X)#X#
}


x=varTF(X,max.kurtosis = 3)$X
x=scale(x)
#xvar=apply(x, 2, var)
#x=X
y=Y

#x = scale(X)
#eig = eigen(cov(x))
#V = as.matrix(eig$vectors[,cumsum(eig$values)/sum(eig$values) <= 0.99])
#x=x%*%V


n=nrow(x)
p=ncol(x)


if(1==2){
  #dr.mean = mave(y~x, method="meanOPG")
  #dr.mean.dim <-min(4,mave.dim(dr.mean)$dim.min)
  #B=dr.mean$dir[[dr.mean.dim]]
  dr.mean.dim=3
  
  
  gcv0=Inf
  #smp=seq(min(4,p))#
  
  di=1
  for (di in seq(3))
  {
    #BB = mave(y~x, method="meanOPG")$dir[[5]]
    BB =OPG.mars(x,y,degree=di,dm = 4,max.iter=10)$B
    #BB = V%*%BB
    
    n.sim=dr.mean.dim
    for (n.sim in seq(4))
    {
      Bi = BB[,1:n.sim]
      xB = cbind(x%*%Bi)
      colnames(xB)=paste0("x",seq(n.sim))
      for (d in seq(ncol(xB)))
      {
        fiti=earth(xB, y, degree=d)
        gcvi = fiti$gcv
        #gcvi = log(fiti$gcv) + log(n)*(n.sim*log(n) +nrow(fiti$coefficients))/n
        if (gcvi < gcv0)
        {
          Dr.drmars = n.sim
          gcv0 = gcvi
          B = Bi
          fit0 = fiti
          Df.drmars=d
        }
      }
    }
  }
  
  fit.drmars=fit0
  print(c(Df.drmars=Df.drmars,Dr.drmars=Dr.drmars))
}

if(1==2){
  max.dim=ncol(Bi)
  cv=-Inf
  for (d in seq(max.dim)) {
    B0=Bi[,seq(d),drop=F]
    set.seed(230516)
    cvi=earth(x%*%B0, y, degree=1,nfold=10)$cv.rsq.tab[11,2]
    if(cvi>cv){
      cv=cvi
      B=B0
      Dr=d
    }
  }
  
  gcv0=Inf
  d=ncol(B)
  xB=x%*%B
  #xB=cbind(x,x%*%B)
  sq=seq(ncol(xB))
  for (df in sq) {
    #set.seed(230516)
    fiti = earth(xB, y, degree=df)#, nfold=10)
    gcvi=fiti$gcv#log(fiti$gcv)#+log(n)*(d*log(n) + nrow(fiti$coefficients))/n
    
    if(gcvi<gcv0){
      fit0=fiti
      gcv0=gcvi
      #B=B0
      #df0=df
    }
  }
}

D=3
B = drMARS(x,y,degree = NULL,Xscale=F,plus=F)$B[,seq(D),drop=F]
xB=x%*%B
colnames(xB)=paste0("x",seq(ncol(xB)))
fit.drmars = earth(xB, y, degree=1)


#pred = fit.drmars$fitted.values
#e.drMARS = mean(abs(pred-y)^2)

#fit.mars=earth(x, y, degree=1)#
#pred = fit.mars$fitted.values
#e.MARS = mean(abs(pred-y)^2)


#######################################
coef=fit.drmars$coefficients

B0=B
B=(abs(B0)>0.01)*B0
write.csv(B, file="B-data3.csv", row.names=T)


#B[seq(2),seq(2)]=diag(2)
xb = cbind(x[,2],x[,1],x%*%B[,3])
data=data.frame(xb, y=y)
colnames(data)=c("sex","age","v3",'y')
fit.gam= gam::gam(y~ sex +age +s(v3), data=data)

win.graph()
#par(mfrow=c(np/2,2))
par(mfrow=c(1,3))
#plot(fit.gam,se=TRUE,ask=TRUE,xlab="v1",ylab="fitted function g1(v1)")#expression(paste(beta[1]^T,x))
#expression(paste("fitted function",g[1],"(",beta[1]^T,x,")")))
plot(fit.gam,ask=TRUE,residuals=FALSE, se=TRUE, cex=1,pch =20,xlab="",ylab="",xaxt="n")
1
0
axis(1, seq(min(xb[,1]),max(xb[,1]),length.out=5),round(seq(min(X[,2]),max(X[,2]),length.out=5),1),cex.axis=1)
title(xlab=TeX(r"(sex)"),ylab=TeX(r"(fitted function $g_1(\cdot)$)"),cex.lab=1)

plot(fit.gam,ask=TRUE,residuals=FALSE, se=TRUE, cex=1,pch =20,xlab="",ylab="",xaxt="n")
2
0
axis(1, seq(min(xb[,2]),max(xb[,2]),length.out=5),round(seq(min(X[,1]),max(X[,1]),length.out=5),1),cex=1)
title(xlab=TeX(r"(age)"),ylab=TeX(r"(fitted function $g_2(\cdot)$)"))

plot(fit.gam,ask=TRUE,residuals=FALSE, se=TRUE, cex=1,pch =20,xlab="",ylab="",xaxt="n")
3
0 
axis(1, seq(min(xb[,3]),max(xb[,3]),length.out=5),round(seq(min(xb[,3]),max(xb[,3]),length.out=5),1),cex=5)
title(xlab=TeX(r"($\beta_3^Tx$)"),ylab=TeX(r"(fitted function $g_3(\cdot)$)"))


if(1==2){
  win.graph()
  #par(mfrow=c(np/2,2))
  par(mfrow=c(1,3))
  #for (ip in seq(ncol(B))) {
  #fB=fit.gam$fitted.values#y
  fB=fit.drmars$fitted.values
  xB=xb[,1]
  plot(x=xB,y=fB,pch =20,xlab=TeX(r"($v_1$)"),ylab=TeX(r"(fitted values $E(Y|X)$)"))#expression(paste("fitted values G",hat(G),"(",x,")")))#,main=paste0("plot-B",ip))#ylab="fxB"
  
  xB=xb[,2]
  plot(x=xB,y=fB,pch =20,xlab=TeX(r"($v_2$)"),ylab=TeX(r"(fitted values $E(Y|X)$)"))#expression(paste("fitted values G",hat(G),"(",x,")")))#,main=paste0("plot-B",ip))#ylab="fxB"
  
  xB=xb[,3]
  plot(x=xB,y=fB,pch =20,xlab=TeX(r"($v_3$)"),ylab=TeX(r"(fitted values $E(Y|X)$)"))#expression(paste("fitted values G",hat(G),"(",x,")")))#,main=paste0("plot-B",ip))#ylab="fxB"
}

##############################################################################################

if(1==2){
  nv=300
  v1=seq(min(xb[,1]),max(xb[,1]),length.out=nv)
  v2=seq(min(xb[,2]),max(xb[,2]),length.out=nv)
  v3=seq(min(xb[,3]),max(xb[,3]),length.out=nv)
  xB=cbind(v1,v2,v3)
  
  fx=fxx(coef=fit.drmars$coefficients,xB=xB)
  fv1=fx$fxValue[,1]
  fv2=fx$fxValue[,2]
  fv3=fx$fxValue[,3]
  
  #win.graph()
  #par(mfrow=c(np/2,2))
  #par(mfrow=c(1,3))
  
  plot(x=v1,y=fv1,xlab=TeX(r"($v_1$)"),ylab=TeX(r"(fitted function $g_1(v_1)$)"),col="black",
       pch =20,type = "b")
  plot(x=v2,y=fv2,xlab=TeX(r"($v_2$)"),ylab=TeX(r"(fitted function $g_2(v_2)$)"),col="black",
       pch =20,type = "b")
  plot(x=v3,y=fv3,xlab=TeX(r"($v_3$)"),ylab=TeX(r"(fitted function $g_3(v_3)$)"),col="black",
       pch =20,type = "b")
  
  
  write.csv(B, file="B-data3.csv", row.names=FALSE)
  
}

if(1==2){
  #xbnew = xnew%*%B
  #name1 = paste("v", 1:p1, sep="")
  cnameB = paste0("B",seq(ncol(B)))
  colnames(xb)=cnameB
  
  #fm=c(paste0(paste0("s(",c("V1",cnameB)),")"))
  #fm=c(paste0(paste0("s(",cnameB),")"))
  fm=c("s(B1)", "B2" ,"s(B3)")
  fm=formula(paste("y",paste(fm,collapse = "+"),sep = "~")) 
  
  #data = data.frame(cbind(V1=z1, xb, y=y))
  data = data.frame(cbind(xb, y=y))
  fit.gam0 = mgcv::gam(fm, data=data)
  fit.gam= gam::gam(fm, data=data)
  
  
  data=data.frame(cbind(x[,c(1,2,18)], y=y))
  colnames(data)=c("age","sex","DFA",'y')
  fit.gam0 = mgcv::gam(y~age + sex +s(DFA), data=data)
  #fit.gam0 = gam::gam(y~sex + age + s(DFA), data=data)
  
  
  B=(abs(B)>0.01)*B
  B[seq(2),seq(2)]=diag(2)
  xb = x%*%B
  data=data.frame(cbind(xb, y=y))
  colnames(data)=c("age","sex","v3",'y')
  fit.gam= gam::gam(y~age + sex +s(v3), data=data)
  #fit.gam0 = mgcv::gam(y~age + sex +s(v3), data=data)
  
  ########################################
  if(ncol(B)<3){
    np=2}else if(ncol(B)<5){
      np=4}else if(ncol(B)<7){
        np=6}else{np=8}
  
  win.graph()
  #par(mfrow=c(np/2,2))
  par(mfrow=c(1,3))
  #for (ip in seq(ncol(B))) {
  #fB=c(fit.gam$smooth[,ip])
  fB=fit.gam$fitted.values#y
  xB=xb[,1]
  plot(x=xB,y=fB,pch =20,xlab="v1",ylab="fitted values E(Y|X)")#expression(paste("fitted values G",hat(G),"(",x,")")))#,main=paste0("plot-B",ip))#ylab="fxB"
  
  xB=xb[,2]
  plot(x=xB,y=fB,pch =20,xlab="v2",ylab="fitted values E(Y|X)")#expression(paste("fitted values G",hat(G),"(",x,")")))#,main=paste0("plot-B",ip))#ylab="fxB"

  xB=xb[,3]
  plot(x=xB,y=fB,pch =20,xlab="v3",ylab="fitted values E(Y|X)")#expression(paste("fitted values G",hat(G)),"(",x,")")))#,main=paste0("plot-B",ip))#ylab="fxB"
  
  #lines(xB,fB,type="o",pch=16,col="green",lty=1,lwd=0.5)
  #} 
  
  win.graph()
  #par(mfrow=c(np/2,2))
  par(mfrow=c(1,3))
  #plot(fit.gam,se=TRUE,ask=TRUE,xlab="v1",ylab="fitted function g1(v1)")#expression(paste(beta[1]^T,x))
  #expression(paste("fitted function",g[1],"(",beta[1]^T,x,")")))
  plot(fit.gam,ask=TRUE,residuals=FALSE, se=TRUE, cex=1,pch =20,xlab="",ylab="",xaxt="n")
  1
  0
  axis(1, seq(min(xb[,1]),max(xb[,1]),length.out=5),round(seq(min(X[,1]),max(X[,1]),length.out=5),1),cex.axis=1)
  title(xlab=TeX(r"(age)"),ylab=TeX(r"(fitted function $g_1(\cdot)$)"),cex.lab=1)
  
  plot(fit.gam,ask=TRUE,residuals=FALSE, se=TRUE, cex=1,pch =20,xlab="",ylab="",xaxt="n")
  2
  0
  axis(1, seq(min(xb[,2]),max(xb[,2]),length.out=5),round(seq(min(X[,2]),max(X[,2]),length.out=5),1),cex=1)
  title(xlab=TeX(r"(sex)"),ylab=TeX(r"(fitted function $g_2(\cdot)$)"))
  
  plot(fit.gam,ask=TRUE,residuals=FALSE, se=TRUE, cex=1,pch =20,xlab="",ylab="",xaxt="n")
  3
  0 
  axis(1, seq(min(xb[,3]),max(xb[,3]),length.out=5),round(seq(min(xb[,3]),max(xb[,3]),length.out=5),1),cex=5)
  title(xlab=TeX(r"($\beta_3^Tx$)"),ylab=TeX(r"(fitted function $g_3(\cdot)$)"))

  
  fB=fit.gam0$fitted.values#y
  xB=data[,1] #x[,1]
  plot(x=xB,y=fB,pch =20,xlab="age",ylab="fitted values E(Y|X)")
  
  xB=data[,2]
  plot(x=xB,y=fB,pch =20,xlab="sex",ylab="fitted values E(Y|X)")
  
  xB=data[,3]
  #plot(x=xB,y=fB,pch =20,xlab="DFA",ylab="fitted values E(Y|X)")
  plot(x=xB,y=fB,pch =20,xlab="v3",ylab="fitted values E(Y|X)")
  

  
  plot(fit.gam0,all.terms=TRUE,se=TRUE)
}




#MARS data4
################################################################
set.seed(4)
if(1==1){
  XY = as.matrix(read.csv(paste0(wd,"/DataSets/Regression1/",'Residential_Building_Data_Set.csv'),header=F))  
  X = XY[,2:103]
  if(1==2){
    #X20=matrix(0,nrow(XY),20)
    #for (ip in seq(20)) {
    #  X20[,ip]= c(XY[,1]==ip)
    #}
    X20=(matrix(XY[,1],nrow(XY),20)==matrix(seq(20),nrow(XY),20,byrow = T))+0
    X=cbind(X20,X) 
  }
  #Y = XY[, 104]  # sales price
  Y = XY[, 105]
  Y = log(Y) #N=372,p=103
}

#dat=datas[4]
#XY =  Data.Set(i.data=dat,wd="D:/BaiduNetdiskWorkspace")#,wd="/home/liuyu/files/DataSets/Regression1")#,wd="D:/BaiduNetdiskWorkspace/DataSets/Regression1")#,wd="/home/liuyu/files/DataSets/Regression1")#
#X = XY$X0
#Y = XY$y0

rm(XY)


if(1==2){
  #XT=varTF(X,max.kurtosis = 3)$X
  XT=X
  x_20=XT[,-seq(20)]
  x_20=varTF(x_20,max.kurtosis = 3)$X
  x_20=scale(x_20)
  x = cbind(XT[,seq(20)],x_20)# scale(X)#X#
}


x=varTF(X,max.kurtosis = 3)$X
x=scale(x)
#xvar=apply(x, 2, var)
#x=X
y=Y

#x = scale(X)
#eig = eigen(cov(x))
#V = as.matrix(eig$vectors[,cumsum(eig$values)/sum(eig$values) <= 0.99])
#x=x%*%V


n=nrow(x)
p=ncol(x)

if(1==2){
  #dr.mean = mave(y~x, method="meanOPG")
  #dr.mean.dim <-min(4,mave.dim(dr.mean)$dim.min)
  #B=dr.mean$dir[[dr.mean.dim]]
  dr.mean.dim=3
  
  
  gcv0=Inf
  #smp=seq(min(4,p))#
  
  di=1
  for (di in seq(3))
  {
    #BB = mave(y~x, method="meanOPG")$dir[[5]]
    BB =OPG.mars(x,y,degree=di,dm = 4,max.iter=10)$B
    #BB = V%*%BB
    
    n.sim=dr.mean.dim
    for (n.sim in seq(4))
    {
      Bi = BB[,1:n.sim]
      xB = cbind(x%*%Bi)
      colnames(xB)=paste0("x",seq(n.sim))
      for (d in seq(ncol(xB)))
      {
        fiti=earth(xB, y, degree=d)
        gcvi = fiti$gcv
        #gcvi = log(fiti$gcv) + log(n)*(n.sim*log(n) +nrow(fiti$coefficients))/n
        if (gcvi < gcv0)
        {
          Dr.drmars = n.sim
          gcv0 = gcvi
          B = Bi
          fit0 = fiti
          Df.drmars=d
        }
      }
    }
  }
  
  fit.drmars=fit0
  print(c(Df.drmars=Df.drmars,Dr.drmars=Dr.drmars))
}


D=3
B = drMARS(x,y,degree = NULL,Xscale=F,plus=F)$B[,seq(D),drop=F]
xB=x%*%B
colnames(xB)=paste0("x",seq(ncol(xB)))
fit.drmars = earth(xB, y, degree=2)

######################################################
coef=fit.drmars$coefficients


win.graph()
par(mfrow=c(1,3))

B0=B
B=(abs(B0)>0.01)*B0
xb = x%*%B
fv=fit.drmars$fitted.values
v1=xb[,1]
v2=xb[,2]
v3=xb[,3]

scatter3D(v1,v2, fv, phi = 0, theta = 40,xlab="v1",ylab="v2",zlab="fitted values E(Y|X)",#main=expression(paste("b1=",beta[1],"  ","b2=",beta[2])),
          colkey = FALSE,col="black",#xaxt="n",yaxt="n", #bty = "g",
          pch =20, ticktype = "detailed",cex.axis=2,cex.lab=2.5)
#axis(1, seq(min(v1),max(v1),length.out=5),round(seq(min(X[,1]),max(X[,1]),length.out=5),1),cex=1)
#axis(2, seq(min(v2),max(v2),length.out=5),round(seq(min(X[,2]),max(X[,2]),length.out=5),1),cex=1)

scatter3D(v1,v3, fv, phi = 0, theta = 40,xlab="v1" ,ylab="v3",zlab="fitted values E(Y|X)",#main=expression(paste("b1=",beta[1],"  ","b3=",beta[3])),
          colkey = FALSE,col="black",# bty = "g",
          pch =20, ticktype = "detailed",cex.axis=2,cex.lab=2.5)

scatter3D(v2,v3, fv, phi = 0, theta = 40,xlab="v2" ,ylab="v3",zlab="fitted values E(Y|X)",#main=expression(paste("b2=",beta[2],"  ","b3=",beta[3])),
          colkey = FALSE,col="black", #bty = "g",
          pch =20, ticktype = "detailed",cex.axis=2,cex.lab=2.5)

#######################


nv=30
v1=seq(min(xb[,1]),max(xb[,1]),length.out=nv)
v2=seq(min(xb[,2]),max(xb[,2]),length.out=nv)
v3=seq(min(xb[,3]),max(xb[,3]),length.out=nv)
v12 <- mesh(v1,v2)
v13 <- mesh(v1,v3)
v23 <- mesh(v2,v3)


xB=cbind(c(v12$x),c(v12$y),c(v13$y))
fx=fxx(coef=fit.drmars$coefficients,xB=as.matrix(xB))
fv12=matrix(rowSums(fx$fxValue[,c(1,2,4)]),ncol = nv)
fv13=matrix(rowSums(fx$fxValue[,c(1,3,5)]),ncol = nv)


xB=cbind(c(v12$x),c(v23$x),c(v23$y))
fx=fxx(coef=fit.drmars$coefficients,xB=as.matrix(xB))
fv23=matrix(rowSums(fx$fxValue[,c(2,3)]),ncol = nv)


win.graph()
par(mfrow=c(1,3))

persp(v1, v2, fv12, xlab= "v1",ylab="v2",zlab="fitted function g1(v1)+g2(v2)+g12(v1,v2)",
      phi = 0, theta = 40,ticktype = "detailed", col = "lightblue",cex.axis=2,cex.lab=2.5)


persp(v1, v3, fv13, xlab="v1" ,ylab="v3",zlab="fitted function g1(v1)+g3(v3)+g13(v1,v3)",
      phi = 0, theta = 40,ticktype = "detailed", col = "lightblue",cex.axis=2,cex.lab=2.5)

persp(v2, v3, fv23, xlab="v2" ,ylab="v3",zlab="fitted function g2(v2)+g3(v3)",#+g23(v2,v3)
      phi = 0, theta = 40,ticktype = "detailed", col = "lightblue",cex.axis=2,cex.lab=2.5)



B0=B
B=(abs(B0)>0.1)*B0
write.csv(B, file="B-data4.csv", row.names=FALSE)


###################################
if(1==2){
  x <- seq(0, 4, length.out=100)
  alpha <- 1:5
  
  plot(x, xlim=c(0, 4), ylim=c(0, 10), 
       xlab='x', ylab=TeX(r'($\alpha  x^\alpha$, where $\alpha \in \{1 \ldots 5\}$)'), 
       type='n', main=TeX(r'(Using $\LaTeX$ for plotting in base graphics!)', bold=TRUE))
  
  for (a in alpha) {
    lines(x, a*x^a, col=a)
  }
  
  legend('topleft', 
         legend=TeX(sprintf(r'($\alpha = %d$)', alpha)), 
         lwd=1, 
         col=alpha)
  
  
  x <- seq(0, 4, length.out=100)
  alpha <- 1:5
  data <- map_df(alpha, ~ tibble(v=.*x^., x=x, alpha=.))
  
  p <- ggplot(data, aes(x=x, y=v, color=as.factor(alpha))) +
    geom_line() + 
    ylab(TeX(r'($\alpha  x^\alpha$, where $\alpha \in 1\ldots 5$)')) +
    ggtitle(TeX(r'(Using $\LaTeX$ for plotting in ggplot2. I $\heartsuit$ ggplot!)')) +
    coord_cartesian(ylim=c(-1, 10)) +
    guides(color=guide_legend(title=NULL)) +
    scale_color_discrete(labels=lapply(sprintf(r'($\alpha = %d$)', alpha), TeX)) 
  # Note that ggplot2 legend labels must be lists of expressions, not vectors of expressions
  
  print(p)
  
  latex2exp_examples(cex=0.9)
}


if(1==2){
  for (ip in 2:ncol(B)) {
    #xb[,ip] = lm(xb[,ip]~xb[,ip-1])$residuals
    xb[,ip] = lm(xb[,ip]~xb[,1:(ip-1)])$residuals
  } 
}


if(1==2){
  xxb=x#cbind(x,x%*%B)#
  for (ib in seq(ncol(B))) {
    Bi = OPG.mars(xxb,y,degree=NULL,dm = ib,max.iter=10, method=c("S1","S2","S12")[1])$B
    xbi=xxb%*%Bi
    xxb = cbind(x,xbi)
  }
  xb=xxb[,-seq(ncol(x))]
}


if(1==2){
  #xbnew = xnew%*%B
  #name1 = paste("v", 1:p1, sep="")
  cnameB = paste0("B",seq(ncol(B)))
  colnames(xb)=cnameB
  
  #fm=c(paste0(paste0("s(",c("V1",cnameB)),")"))
  fm=c(paste0(paste0("s(",cnameB),")"))
  fm=formula(paste("y",paste(fm,collapse = "+"),sep = "~")) 
  
  #data = data.frame(cbind(V1=z1, xb, y=y))
  data = data.frame(cbind(xb, y=y))
  fit.gam0 = mgcv::gam(fm, data=data)
  fit.gam= gam::gam(fm, data=data) 
}



########################################
if(1==2){
  if(ncol(B)<3){
    np=2}else if(ncol(B)<5){
      np=4}else if(ncol(B)<7){
        np=6}else{np=8}
  
  plot(fit.gam,se=TRUE,ask=TRUE,xlab="v1",ylab="fitted function g1(v1)")#expression(paste(beta[1]^T,x))
  #expression(paste("fitted function",g[1],"(",beta[1]^T,x,")")))
  1
  0
  plot(fit.gam,se=TRUE,ask=TRUE,xlab="v2",ylab="fitted function g2(v2)")
  2
  0
  plot(fit.gam,se=TRUE,ask=TRUE,xlab="v3",ylab="fitted function g3(v3)")
  3
  0
  
  #for (i in 1:3) {
  #  plot(b,select=i,rug=FALSE,col="green",
  ##       col.axis="white",col.lab="white",all.terms=TRUE)
  #  for (j in 1:2) axis(j,col="white",labels=FALSE)
  #  box(col="white")
  #  eval(parse(text=paste("fx <- f",i,"(x)",sep="")))
  #  fx <- fx-mean(fx)
  #  lines(x,fx,col=2) ## overlay `truth' in red
  #}
}

