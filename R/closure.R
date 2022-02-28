closure <-
function(x, k){ # x the NxD matrix, k the constant (k=1 if missing)
  if (missing(k)){
    newx <- (x)/rowSums(x)
    return(newx)
  }
  else{
    newx <- (k*x)/rowSums(x)
    return(newx)
  }
}
