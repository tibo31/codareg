est_var <-
function(X,Y,B){ # page 146 kevin
  if (missing(B)){B=mlm(Y,X)}
  n=dim(Y)[1]
  sig = (1/n)*t(Y-X%*%B)%*%(Y-X%*%B)
  return(sig)
}
