inversealr <-
function(x,R=dim(x)[2]+1){ # x = NxD-1 matrix, R original reference variable
  if (is.null(attr(x, "space"))){attr(x, "space") <- "alr_coord"}
  if (attr(x, "space")=="simplex"){stop("The matrix is already in the simplex space.")}
  D_1 = dim(x)[2]
  if (R==(D_1+1)){aumented = cbind(x,numeric(dim(x)[1]))}
  if (R==1){aumented = cbind(numeric(dim(x)[1]),x)}
  if ((R!=(D_1+1)&(R!=1))){aumented = cbind(x[,1:R-1],numeric(dim(x)[1]),x[,R:D_1])}
  res = closure(exp(aumented))
  attr(res,"space") = "simplex"
  attr(res,"reference") = NULL
  return(res)
}
