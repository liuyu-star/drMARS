Data.Set = function(i.data,dir=paste0(getwd(),"/DataSets"))
{
  # Regression datasets.
  if(i.data==1){
    data = read.table(paste0(dir,'/Concrete.data'))
    data = data.matrix(data)
    N = nrow(data)
    p = ncol(data)-1
    X0 = data[,1:p]
    y0 = data[,p+1] #N=1030,p=8
  }
  if(i.data==2){
    data = read.csv(paste0(dir,'/kc_house_data.csv'))
    N = nrow(data)
    p = ncol(data)
    X0 = as.matrix(data[,4:p])
    y0 = log(as.matrix(data[,3])) #N=21613,p=18
  }
  if(i.data==3){
    data = read.csv(paste0(dir,'/parkinsons.csv'))
    p = ncol(data)
    data = as.matrix(data[,2:p])
    y0 = data[,5]
    X0 = data[,-c(4, 5)]
    N = nrow(data)
    p = ncol(data)
    #  X0 = as.matrix(data[,1:p-1]);
    #  y0 = as.matrix(data[,p]);
    I = (apply(X0, 2, sd)>0)
    X0 = X0[,I]
    p = ncol(X0) #N=5875,p=19
  }
  if(i.data==4){
    XY = as.matrix(read.csv(paste0(dir,'/Residential_Building_Data_Set.csv',header=F)))  
    X0 = XY[,2:103]
    y = XY[, 105] # construction cost
    N = length(y)
    p = ncol(X0)
    y0 = log(y) #N=372,p=102
  }
  
  # Classification datasets.
  if(i.data==5){
    yx=read.csv(paste0(dir,"/Pistachio_28_Features_Dataset.csv"))
    p=ncol(yx)-1
    n=nrow(yx) #n=2148,p=28
    yx[,p+1]=as.integer(as.factor(yx[,p+1]))
    y0=rep(0,n)
    y0[yx[,p+1]==2]=1
    X0=as.matrix(yx[,1:p])
  }
  if(i.data==6){
    yx1=read.table(paste0(dir,"/data6/Hill_Valley_without_noise_Testing.txt"),header=TRUE,sep=",")
    yx2=read.table(paste0(dir,"/data6/Hill_Valley_without_noise_Training.txt"),header=TRUE,sep=",")
    yx=rbind(yx1,yx2)
    n=dim(yx)[1]
    p=dim(yx)[2]-1 #n=1212,p=100
    yx[,p+1]=as.integer(as.factor(yx[,p+1]))
    X0=as.matrix(yx[,1:p])
    y0=rep(0,n)
    y0[yx[,p+1]==2]=1
  }
  if(i.data==7){
    om<-function(a) return(sum(is.na(a)))
    yx=read.csv(paste0(dir,"/2018_Financial_Data.csv"))
    yx=yx[apply(yx,1,om)<=3,]
    yx=yx[,apply(yx,2,om)==0]
    yx=yx[,-c(1,218)]
    p=ncol(yx)-1
    n=nrow(yx) #n=986,p=217
    y0=yx[,p+1]
    X0=as.matrix(yx[,1:p])
  }
  
  ####################################################
  #setwd(wd0)
  n = dim(X0)[1]
  p = dim(X0)[2]
  
  #I = union(c(which(rowSums(is.na(X0))>0),which(rowSums(X0==0)==p)),which(is.na(y0)))
  I = which((rowSums(is.na(X0))>0)|is.na(y0))#,which(rowSums(X==0)==p))
  if(1==2){
    
    if(n>1000){
      #X = X[seq(2000),]X
      #Y = Y[seq(2000)]
      
      for (ip in seq(p)) {
        id=which(X0[,ip]!=0)
        Xi=X0[id,ip]
        boxs=boxplot.stats(Xi)
        I=union(I,id[which(Xi%in%boxs$out)])
        #I=union(I,id[which((Xi<boxs$stats[2])|(Xi>boxs$stats[4]))])
      }
    }
  }
  if(length(I)>0){
    X0=X0[-I,]
    y0=y0[-I]
  }
  
  if(1==2){
    for (i in 1:p)
    {
      xi = X0[,i]
      xi = xi - min(xi) + (max(xi)-min(xi))/100
      result = boxcox(xi~1, lambda = seq(-5,5,0.5))
      mylambda = result$x[which.max(result$y)]
      xi = (xi^mylambda-1)/mylambda
      I = (xi-mean(xi) > 3*sd(xi))
      if (length(I) > 0)
        xi[I] = mean(xi) + 3*sd(xi)
      
      I = (xi-mean(xi) < -3*sd(xi))
      if (length(I) > 0)
        xi[I] = mean(xi) - 3*sd(xi)
      
      X0[,i] = xi
    }
  }
  
  #I = which(apply(X0, 2, sd)>0)
  #X0 = X0[,I]
  return(list(X = X0, y=y0))
}

