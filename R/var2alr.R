var2alr <-
function(Y,X,Z=list(matrix()),R=list(matrix()),constant=TRUE){ 
  
  # determine if variables are in the simplex or in ALR coordinates
  n_simplex = 0
  n_alr     = 0
  
  if (is.null(attr(Y,"space"))){
    n_simplex = 1
  }
  
  if (!is.null(attr(Y,"space"))){
    if (attr(Y,"space")=="simplex"){
      n_simplex=1
    }
    else{
      n_alr=1
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
        n_alr=n_alr+1
      }
    }
  }
  
  # case when X, Y is in the simplex
  
  # if R specified
  
  if (!is.na(R[[1]])){
    if (length(R)!=(1+length(X))){stop("When using the optional argument 'R' it is needed to specify one reference for the independent variable and one for each of the q dependent variable in a list that takes the form R = list(R_Y,R_X1...,R_Xq).")}
    }
    
  # if R not specified
  
  if (n_alr==0){
    if (is.na(R[[1]])){
      R=list(dim(Y)[2])
      for (i in 1:length(X)){ # here 'standard' contrast matrices for X are added
        R[[i+1]]=dim(X[[i]])[2]
      }
    }
    
  # now we have R for the simplex variables case
  
  # transform Y into alr coordinates

  ty = alr(Y, R[[1]])
    
  # transfrom X to the alr coordinates
  tx = matrix(numeric(dim(X[[1]])[1]))+1 # here it is created a matrix which is a vector of ones (should be removed afterwords if c=False)
  for (i in 1:length(X)){
    tx = cbind(tx,alr(X[[i]],R[[i+1]])) # the transformed X is created as a unique matrix
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
  
  
  # case where some variables in the simplex, some in alr
  
  if (min(c(n_simplex,n_alr))!=0){
    stop("Every variable should be in the same space, whether alr coordinates or the simplex.")
  }
    
  # every variable in alr coordinates: R should be specified in the arguments or in the attributes of the data Y,X
  
  if (!is.na(R[[1]])){
    if (!is.null(attr(Y,"reference"))){
      equal = TRUE
      for (k in 1:(length(X)+1)){
        if (k==1){
          vk = R[[1]]
          va = attr(Y,"reference")
          tof = all(dim(vk)==dim(va)) && all(vk==va)
          if (!tof){
            equal = FALSE
          }
        }
        else{
          vk = R[[k]]
          va = attr(R[[k-1]],"reference")
          tof = all(dim(vk)==dim(va)) && all(vk==va)
          if (!tof){
            equal = FALSE
          }
        }
      }
      if (!equal){warning("At least one alr reference given does not coincide with its analogous in the attributes of the data.\n")}
    }
  }
  
  if (n_simplex==0){
    # construct ty (and recover R if not specified)
    if (is.na(R[[1]])){
      if (is.null(attr(Y,"reference"))){
        stop("When using data in the coordinates space, their associated references list is needed. It can be attached as an attribute of each matrix of data or directly given as an input in R.")}
      else{
        R[[1]]=attr(Y,"reference")
        for (k in 1:length(X)){R[[k+1]]=attr(X[[k]],"reference")
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
  
  return (list(Y=ty,X=tx,R=R,constant=constant))
}
