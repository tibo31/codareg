inv_permutation_D <-
function(mat,P){ # undoes the permutation
  D = dim(mat)[2]
  if (D==P){
    return(mat)
  }
  if (1==P){
    return(cbind(mat[,D],mat[,1:(D-1)]))
  }
  else{
    return(cbind(mat[,1:(P-1)],mat[,D],mat[,(P):(D-1)]))
  }
}
