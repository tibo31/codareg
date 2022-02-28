b_cov <-
function(X,Y,B){
  sig = est_var(X,Y,B)
  cov = kronecker(sig,solve(t(X)%*%X))
  return(cov)
}
