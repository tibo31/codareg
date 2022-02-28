var2ilr <-
function(Y,X,Z=list(matrix()),V=list(matrix()),constant=TRUE){ 
  
  # determine if variables are in the simplex or in ILR coordinates
  n_simplex = 0
  n_ilr     = 0
  
  if (is.null(attr(Y,"space"))){
    n_simplex = 1
  }
  
  if (!is.null(attr(Y,"space"))){
    if (attr(Y,"space")=="simplex"){
      n_simplex=1
    }
    else{
      n_ilr=1
    }
  }
  
  for (k in 1:length(X)){
    if (is.null(attr(X[[k]],"space"))){
    n_simplex = n_simplex+1}
    if (!is.null(attr(X[[k]],"space"))){
      if (attr(X[[k]],"space")=="simplex"){
        n_simplex=n_simplex+1
        }
      else{
        n_ilr=n_ilr+1
      }
    }
  }
  
  # case when X, Y is in the simplex
  
  # if V specified
  
  if (!is.na(V[[1]][1,1])){
    if (length(V)!=(1+length(X))){stop("When using the optional argument 'V' it is needed to specify one base matrix for the independent variable and one for each of the q dependent variable in a list that takes the form V = list(V_Y,V_X1...,V_Xq).")}
    }
    
  # if V not specified
  
  if (n_ilr==0){
    if (is.na(V[[1]][1,1])){
      V=list(V_D(dim(Y)[2])) # this the 'standard' contrast matrix for Y
      for (i in 1:length(X)){ # here 'standard' contrast matrices for X are added
        V=c(V,list(V_D(dim(X[[i]])[2])))
      }
    }
    
  # now we have V for the simplex variables case
  
  # transform Y into ilr coordinates

  ty = ilr(Y, V[[1]])
    
  # transform X to the ilr coordinates
  tx = matrix(numeric(dim(X[[1]])[1]))+1 # here it is created a matrix which is a vector of ones (should be removed afterwords if c=False)
  for (i in 1:length(X)){
    tx = cbind(tx,ilr(X[[i]],V[[i+1]])) # the transformed X is created as a unique matrix
    }
  if (constant==FALSE){ # should remove the vector of zeros if no constant is required
    tx=tx[,-1]
    }
  if (!is.na(Z[[1]][1,1])){
    for (z in Z){
      tx=cbind(tx,z)
    }
  }
  }
  
  
  # case where some variables in the simplex, some in ilr
  
  if (min(c(n_simplex,n_ilr))!=0){
    stop("Every variable should be in the same space, whether ilr coordinates or the simplex.")
  }
    
  # every variable in ilr coordinates: V should be specified in the arguments or in the attributes of the data Y,X
  
  if (!is.na(V[[1]][1,1])){
    if (!is.null(attr(Y,"contrast"))){
      equal = TRUE
      for (k in 1:(length(X)+1)){
        if (k==1){
          vk = V[[1]]
          va = attr(Y,"contrast")
          tof = all(dim(vk)==dim(va)) && all(vk==va)
          if (!tof){
            equal = FALSE
          }
        }
        else{
          vk = V[[k]]
          va = attr(X[[k-1]],"contrast")
          tof = all(dim(vk)==dim(va)) && all(vk==va)
          if (!tof){
            equal = FALSE
          }
        }
      }
      if (!equal){warning("At least one contrast matrix given in V does not coincide with its analogous in the attributes of the data.\n")}
    }
  }
  
  if (n_simplex==0){
    # construct ty (and recover V if not specified)
    if (is.na(V[[1]][1,1])){
      if (is.null(attr(Y,"contrast"))){
        stop("When using data in the coordinates space, their associated contrast matrices are needed. They can be attached as an attribute of each matrix of data or directly given as an input in V.")}
      else{
        V[[1]]=attr(Y,"contrast")
        for (k in 1:length(X)){V[[k+1]]=attr(X[[k]],"contrast")
        }
      }
    }
    ty = Y
    # construct tx (and recover V)
    tx = matrix(numeric(dim(X[[1]])[1]))+1
    for (k in 1:length(X)){
      tx = cbind(tx,X[[k]])
      }
    if(constant==FALSE){ # should remove the vector of zeros if no constant is required
      tx=tx[,-1]
      }
    if (!is.na(Z[[1]][1,1])){
      for (z in Z){
        tx=cbind(tx,z)
      }
    }
  }
  
  return (list(Y=ty,X=tx,V=V,constant=constant))
}
