alr_x_reg <-
function(Y,X,Z=list(matrix()),R=list(NULL),constant=TRUE, verbose=FALSE, pres=FALSE){ # Y dependent, X indep compositional, Z indep. not compositional, c if constant is included, V list of contrast matrices. order: (Y, X1, X2, ...) 
# pres = TRUE to print with a better presentation (name of variables and sd if it is the case) of B and B*.
  
  if (constant){x_col_names = c("intercept")}
  else{x_col_names = c()}
  
  for (i in 1:length(X)){
    if (is.null(colnames(X[[i]]))){
      for (j in 1:(dim(X[[i]])[2])){
        x_col_names = c(x_col_names, paste("X_",i,"_",j,sep=""))
      }
    }
    else{x_col_names = c(x_col_names,colnames(X[[i]]))}
  }
  if (!is.na(Z[[1]][1,1])){
    for (i in 1:length(Z)){
      if (is.null(colnames(Z[[i]]))){
      for (j in 1:(dim(Z[[i]])[2])){
        x_col_names = c(x_col_names, paste("Z_",i,"_",j,sep=""))
      }
    }
      else{x_col_names = c(x_col_names,colnames(Z[[i]]))
      }
    }
  }
  
  if (is.null(colnames(Y))){
    y_col_names = c()
    for (i in 1:(dim(Y)[2])){
      y_col_names = c(y_col_names, paste("Y_",i,sep=""))
    }
  }
  else{y_col_names = colnames(Y)}
  
  ty = Y #already in coordinates
  if (is.null(R[[1]])){
    for (i in 1:length(X)){
      R[[i]] = dim(X[[i]])[2]
      X[[i]]=permutation_D(X[[i]],R[[i]])
    }
  }
  
  else{
    for (i in 1:length(X)){
      X[[i]]=permutation_D(X[[i]],R[[i]])
      }
    }
  
  tx = matrix(numeric(dim(X[[1]])[1]))+1 # here it is created a matrix which is a vector of ones (should be removed afterwords if c=False)
  for (i in 1:length(X)){
    tx = cbind(tx,alr(X[[i]],R[[i]])) # the transformed X is created as a unique matrix
    }
  if (constant==FALSE){ # should remove the vector of zeros if no constant is required
    tx=tx[,-1]
    }
  if (!is.na(Z[[1]][1,1])){
    for (z in Z){
      tx=cbind(tx,z)
    }
  }
  
  if (sum(is.na(ty))>0){stop("There are NAs in the Y variable.")}
  if (sum(is.na(tx))>0){stop("There are NAs in the X variable.")}
  
  
  if (verbose){
  print("yx REGRESSION")
  print("===========================")
  
  print("ALR Y-Matrix")
  print(ty)
  print("ALR X-Matrix")
  print(tx)}

  Bstar = mlm(ty,tx)
  bcov = b_cov(tx,ty,Bstar)
  Bs = t(matrix(numeric(dim(Y)[2])))
  
  if (verbose){
  print("B* in the coordinates space:")
  print(Bstar)
  print("---------------------------")}
  if (constant==TRUE){
    b0=Bstar[1,]
    Bs=rbind(Bs, Bstar[1,])
    if (verbose){
    print("b0 in the simplex:")
    print(b0)}
  }
  h=2
  xdimlist=list()
  if (verbose){print("---------------------------")}
  for (i in 1:length(X)){
    if (verbose){print(paste("B (in the simplex) for the comp. variable, X", i, sep=""))}
    # number of parts?
    if (is.null(attr(X[[i]],"space"))){dimx = dim(X[[i]])[2]} #number of parts of this variable i of X D parts --> D-1 param in the ilr 
    if (!is.null(attr(X[[i]],"space"))){
      if (attr(X[[i]],"space")=="simplex"){
        dimx = dim(X[[i]])[2]
      }
      else{dimx = dim(X[[i]])[2]+1}
    }
    xdimlist[i] = dimx
    if (h==(h+dimx-2)){matrix_x = t(matrix(Bstar[h,]))}
    if (h!=(h+dimx-2)){matrix_x = Bstar[h:(h+dimx-2),] }
    B_q = t(matrix_x)%*%F_D(dimx)
    B_q = inv_permutation_D(B_q,R[[i]])
    Bs= rbind(Bs, as.matrix(t(B_q)))
    h=h+dimx-1
    if (verbose){print(t(B_q))}
  }
  if (verbose){print("---------------------------")}
  
  if (!is.na(Z[[1]][1,1])){
    for (i in 1:length(Z)){
      if (verbose){print(paste("B (in the simplex) for the non-comp. variable, Z", i, sep=""))}
      dimz = dim(Z[[i]])[2]
      if (h==(h+dimz-1)){matrix_z = t(matrix(Bstar[h,]))}
      if (h!=(h+dimz-1)){matrix_z = Bstar[h:(h+dimz-1),]} # it is non comp, it does not lose a dimension in the transformation
      h=h+dimz
      c_k = as.matrix(matrix_z)
      Bs = rbind(Bs, c_k)
      if (verbose){print(c_k)}
    }
  }
  Bs = as.matrix(Bs[-1,])
  attr(Bstar,"space")="alr_coord"
  attr(Bstar,"reference") = R
  
  attr(tx,"space") = "alr_coord"
  attr(tx,"reference") = R
  
  Bs = as.matrix(Bs)
  
  colnames(Bs) = y_col_names
  rownames(Bs) = x_col_names
  attr(Bs,"space")="simplex"
  
  fittedvalues = tx%*%Bstar
  resid = ty- fittedvalues
  attr(resid,"space") = NULL
  attr(resid,"contrast") = NULL
  
  reslist = list(Y_coord = ty, constant = constant, X_coord = tx, B_coord = Bstar, B_simplex = Bs , Bcoord_cov=bcov, residuals_coord = resid, fitted_v_coord = fittedvalues, xdim=xdimlist)
  attr(reslist,"reg_type") = "alr_x_reg"
  if (!pres){return(reslist)} # return results if not pres

  
  # part of PRES = TRUE
  
  if (pres){
    dfsimplex = data.frame(Bs)
    name=c()
    for (i in 1:(dim(Y)[2])){
      name = c(name, paste("Y_", i, sep=""))
    }
    colnames(dfsimplex)=name
    
    rname = c()
    if (constant){rname = c(rname, "intercept")}
    for (i in 1:length(X)){
      for (j in 1:(dim(X[[i]])[2])){
        rname = c(rname, paste("X",i,j,sep="_"))
      }
    }
    
    if (!is.na(Z[[1]][1,1])){
      for (i in 1:length(Z)){
        for (j in 1:(dim(Z[[i]])[2])){
           rname = c(rname, paste("Z",i,j,sep="_"))
        }
      }
    }
    
    rownames(dfsimplex)=rname
    
    cnames = c()
    for (i in (1:(dim(Y)[2]))){
      cnames = c(cnames, paste("Yalr_",i,sep=""))
      cnames = c(cnames, paste("(sd_Yalr_",i,")",sep=""))
    }
    
    rname = c()
    if (constant){rname = c(rname, "intercept")}
    for (i in 1:length(X)){
      for (j in 1:(dim(X[[i]])[2]-1)){
        rname = c(rname, paste("X",i,j,sep="_"))
          }
        }
    
    if (!is.na(Z[[1]][1,1])){
      for (i in 1:length(Z)){
        for (j in 1:(dim(Z[[i]])[2])){
           rname = c(rname, paste("Z",i,j,sep="_"))
        }
      }
    }

    
    # B star and sds
    dfalr = data.frame(Bstar)
    
    
    
    # sds
    dimalr = dim(Y)[2]
    nc = dim(dfalr)[1]
    for (i in dimalr:1){
      if (i==dimalr){dfalr = cbind(dfalr,matrix(nrow=nc))}
      else{
        dfalr = cbind(dfalr[,1:i],matrix(nrow=nc),dfalr[,(i+1):(dim(dfalr)[2])])
      }
    }
    
    dfalr = t(dfalr)
    
    for (i in 1:(dim(bcov)[1])){
      dfalr[2*i] = sqrt(bcov[i,i]) # standard deviation as sqrt of var
    }
    dfalr = t(dfalr)
    
    rownames(dfalr) = rname
    colnames(dfalr) = cnames
    
    
    print("      B (COEFFICIENTS IN THE SIMPLEX):      ")
    print(as.matrix(dfsimplex))
    print("--------------------------------------------")
    print("     B* (COEF. and SD IN THE ALR SPACE):    ")
    print(as.matrix(dfalr))
    print("--------------------------------------------")
    print(paste("n obs  : ", dim(tx)[1], sep="")) 
    print(paste("SSR    : ", sum(resid^2), sep="")) # sum of squared residuals
  }
}
