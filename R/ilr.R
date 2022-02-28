ilr <-
function(x, V){ # x matrix N x D
  if (is.null(attr(x, "space"))){attr(x, "space") <- "simplex"}
  if (attr(x, "space")=="ilr_coord"){stop("The matrix is already in the ILR coordinates space.")}
  if (attr(x, "space")=="alr_coord"){stop("The matrix is in the ALR coordinates space. Matrix should be in the simplex space.")}
  if (missing(V)){V=V_D(dim(x)[2])}
  x = closure(x)
  ilr_mat = t(t(V)%*%t(log(x)))
  attr(ilr_mat, "space") <- "ilr_coord"
  attr(ilr_mat, "contrast") <- V
  return(ilr_mat)
}
