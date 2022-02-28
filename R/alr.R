alr <-
function(x,R=dim(x)[2]){ # x: N x D matrix, R : Reference (last term by default)
  if (is.null(attr(x, "space"))){attr(x, "space") <- "simplex"}
  if (attr(x, "space")=="alr_coord"){stop("The matrix is already in the ALR coordinates space.")}
  if (attr(x, "space")=="ilr_coord"){stop("The matrix is in the ILR coordinates space. Matrix should be in the simplex space.")}
  copy = x
  for (i in 1:(dim(x)[1])){
    trans = numeric(dim(x)[2])
    for (j in 1:(dim(x)[2])){
      trans[j]= log(x[i,j]/x[i,R])
    }
    copy[i,] = trans
  }
  copy = as.matrix(copy[,-R])
  attr(copy, "space") <- "alr_coord"
  attr(copy, "reference") <- R
  return(copy)
}
