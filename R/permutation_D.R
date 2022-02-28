permutation_D <-
function(mat,P){ # this function takes as an input a matrix (mat) and position (P) and returns the same matrix but with the column P as the last column, i.e: [x1,x2...,xP,...,xN] -> [x1,x2...,xN,xP]
  D = dim(mat)[2]
  if (D==P){
    return(mat)
  }
  if (1==P){
    return(cbind(mat[,2:D],mat[,1]))
  }
  else{
    return(cbind(mat[,1:(P-1)],mat[,(P+1):D],mat[,P]))
  }
}
