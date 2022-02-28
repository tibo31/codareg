share_ratio_summary <-
function(results){
  Bs = results$B_simplex
  xnames = rownames(Bs)
  constant = results$constant
  if (constant){
    if (dim(Bs)[1]==2){
      Bs = t(as.matrix(Bs[-1,]))
      rownames(Bs) = c(xnames[[2]])}
    else{Bs = as.matrix(Bs[-1,])}
  }
  res=list()
  P = (dim(Bs)[1])
  for (k in 1:P){
    dimy = dim(Bs)[2]
    matres = matrix(0,dimy,dimy)
    for (i in 1:dimy){
      for (j in 1:dimy){
        matres[i,j] = log(Bs[k,j]/Bs[k,i])
      }
    }
    colnames(matres) = colnames(Bs)
    rownames(matres) = colnames(Bs)
    res[[xnames[[k+1]]]] = matres
  }
  return(res)
}
