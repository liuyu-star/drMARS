MedianDist=function(X){
  n=nrow(X)
  ab=X%*%t(X)
  aa=diag(ab)
  Dx=matrix(aa,n,n)+matrix(aa,n,n,byrow = T)-2*ab
  dx=as.vector(Dx-diag(diag(Dx)))
  dx=dx[-which(dx==0)]
  dx=sqrt(median(dx))
  return(dx)
}