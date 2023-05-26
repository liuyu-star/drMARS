chol_inc_gauss=function(x,sigma,tol){

# CHOL_INC_FUN - incomplete Cholesky decomposition of the Gram matrix defined
#                by data x, with the Gaussiab kernel with width sigma
#                Symmetric pivoting is used and the algorithms stops 
#                when the sum of the remaining pivots is less than TOL.
# 

# CHOL_INC returns returns an uPvecer triangular matrix G and a permutation 
# matrix P such that P'*A*P=G*G'.

# P is ONLY stored as a reordering vector PVEC such that 
#                    A(Pvec,Pvec)= G*G' 
# consequently, to find a matrix R such that A=R*R', you should do
# [a,Pvec]=sort(Pvec); R=G(Pvec,:);

# Copyright (c) Francis R. Bach, 2002.

########################################################################################

sqdist=function(a,b){
  # SQDIST - computes squared Euclidean distance matrix
  #          computes a rectangular matrix of pairwise distances
  # between points in A (given in columns) and points in B
  
  # NB: very fast implementation taken from Roland Bunschoten
  a=as.matrix(a);b=as.matrix(b)
  aa = colSums(as.matrix(a*a));
  bb = colSums(as.matrix(b*b));
  #bb =apply(as.matrix(b*b), 2, sum); 
  ab = t(a)%*%b; 
  d = abs(matrix(aa,length(aa),length(bb)) + matrix(bb,length(aa),length(bb),by=T) - 2*ab);
   
  return(d)
 }


x=as.matrix(x)
p=ncol(x);
Pvec= 1:p;
I = NULL;
#calculates diagonal elements (all equal to 1 for gaussian kernels)
diagG=rep(1,p);
i=1;
G=NULL;

while ((sum(diagG[i:p])>tol)&&(i<=p)){
  G=cbind(G,rep(0,p));
  # find best new element
  if(i>1){
    jast=which.max(diagG[i:p]);
    jast=jast+i-1;
    #updates permutation
    Pvec[c(i,jast)] = Pvec[c(jast,i)];
    # updates all elements of G due to new permutation
    G[c(i,jast),1:i]=G[c(jast,i),1:i];
    # do the cholesky update
  }else{
    jast=1;
  }

  
  #A(Pvec(i),Pvec(i));
  G[i,i]=sqrt(diagG[jast]);
  if (i<p){
    #calculates newAcol=A(Pvec((i+1):n),Pvec(i))
    newAcol = exp(-0.5/sigma^2*sqdist(x[,Pvec[(i+1):p],drop=F],x[,Pvec[i],drop=F]));
    
    if (i>1){
      G[(i+1):p,i]=1/G[i,i]*(newAcol - G[(i+1):p,1:(i-1),drop=F]%*%t(G[i,1:(i-1),drop=F]));
    }else{
      G[(i+1):p,i]=1/G[i,i]*newAcol;
    }
    
    # updates diagonal elements
    diagG[(i+1):p]=rep(1,p-i) - rowSums(G[(i+1):p,1:i,drop=F]^2);
  }

   i=i+1;
}

return(list(G=G,Pvec=Pvec))
}
  

