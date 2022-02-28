inverseilr <-
function(x,V){
  if (is.null(attr(x, "space"))){attr(x, "space") <- "ilr_coord"}
  if (attr(x, "space")=="simplex"){stop("The matrix is already in the simplex space.")}
  x = t(x)
  if (missing(V)){
    if (is.null(attr(x,"contrast"))){
      V=V_D(dim(x)[1]+1)
    }
    else{V = attr(x,"contrast")}
  }
  inv_mat = closure(t(exp(V%*%x)))
  attr(inv_mat,"space") = "simplex"
  attr(inv_mat,"contrast") = NULL
  return(inv_mat)
}
