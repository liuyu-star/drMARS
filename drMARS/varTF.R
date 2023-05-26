library(moments)

VarTF = function(X, Xnew=c(), max.kurtosis = 3)
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
