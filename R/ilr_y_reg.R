ilr_y_reg <-
function(Y,X,V=NaN,constant=TRUE, verbose=FALSE, pres=FALSE){ # Y dependent, X indep. not compositional, c if constant is included, V list of contrast matrices. order: (Y, X1, X2, ...) 
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
  
  if (is.null(colnames(Y))){
    y_col_names = c()
    for (i in 1:(dim(Y)[2])){
      y_col_names = c(y_col_names, paste("Y_",i,sep=""))
    }
  }
  else{y_col_names = colnames(Y)}
  
  tx = matrix(numeric(dim(X[[1]])[1]))+1 # X = tx as it is already in coord.
  for (i in 1:length(X)){
    tx = cbind(tx,X[[i]])
    }
  if (constant==FALSE){ # should remove the vector of ones if no constant is required
    tx=tx[,-1]
  }
  
  # transform Y to coordinates
  if (!is.na(V)){ # if V specified as argument
    ty = ilr(Y,V)
  }
  
  else{
    V=V_D(dim(Y)[2]) # if V not specified as argument
    ty=ilr(Y,V)}
  
  if (sum(is.na(ty))>0){stop("There are NAs in the Y variable.")}
  if (sum(is.na(tx))>0){stop("There are NAs in the X variable.")}
  
  if (verbose){
  print("yx REGRESSION")
  print("===========================")
  
  print("ILR Y-Matrix")
  print(ty)
  print("ILR X-Matrix")
  print(tx)}
  
  Bstar = mlm(ty,tx)
  bcov = b_cov(tx,ty,Bstar) 
  Bs = t(matrix(numeric(dim(Y)[2])))
  
  if (verbose){
  print("B* in the coordinates space:")
  print(Bstar)
  print("---------------------------")}
  if (constant==TRUE){
    b0=inverseilr(t(Bstar[1,]),V)
    Bs = rbind(Bs, b0)
    if (verbose){
    print("b0 in the simplex:")
    print(b0)}
  }
  h=2
  if (verbose){print("---------------------------")}
  if (!is.na(X[[1]][1,1])){
    for (i in 1:length(X)){
      if (verbose){print(paste("B (in the simplex) for the non-comp. variable, X", i, sep=""))}
      dimx = dim(X[[i]])[2]
      matrix_x = Bstar[h:(h+dimx-1),]
      if (h==(h+dimx-1)){matrix_x = t(matrix(Bstar[h,]))}
      if (h!=(h+dimx-1)){matrix_x = Bstar[h:(h+dimx-1),]} # it is non comp, it does not lose a dimension in the transformation
      h=h+dimx
      c_k = inverseilr(matrix_x, V)
      Bs = rbind(Bs, c_k)
      if (verbose){print(c_k)}
    }
  }
  Bs = Bs[-1,]
  attr(Bstar,"space")="ilr_coord"
  attr(Bstar,"contrast") = V
  
  attr(tx,"space") = "ilr_coord"
  attr(tx,"contrast") = NaN
  

  
  fittedvalues = tx%*%Bstar
  resid = ty- fittedvalues
  attr(resid,"space") = NULL
  attr(resid,"contrast") = NULL
  colnames(Bs) = y_col_names
  rownames(Bs) = x_col_names
  attr(Bs,"space")="simplex"
  
  list_return = list(Y_coord = ty, X_coord = tx, constant = constant, B_coord = Bstar, B_simplex = Bs, Bcoord_cov=bcov, residuals_coord = resid, fitted_v_coord = fittedvalues)
  attr(list_return, "reg_type") = "ilr_y_reg"
  
  if (!pres){return(list_return)} 
# return results if not pres + Bcoord_cov=bcov
  
  # part of PRES = TRUE
  
  if (pres=="TRUE"){
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
    
    rownames(dfsimplex)=rname
    
    cnames = c()
    for (i in (1:(dim(Y)[2]-1))){
      cnames = c(cnames, paste("Yilr_",i,sep=""))
      cnames = c(cnames, paste("(sd_Yilr_",i,")",sep=""))
    }
    
    rname = c()
    if (constant){rname = c(rname, "intercept")}
    for (i in 1:length(X)){
      for (j in 1:(dim(X[[i]])[2])){
        rname = c(rname, paste("X",i,j,sep="_"))
          }
        }

    
    # B star and sds
    dfilr = data.frame(Bstar)
    
    
    # sds
    dimilr = dim(Y)[2]-1
    nc = dim(dfilr)[1]
    for (i in dimilr:1){
      if (i==dimilr){dfilr = cbind(dfilr,matrix(nrow=nc))}
      else{
        dfilr = cbind(dfilr[,1:i],matrix(nrow=nc),dfilr[,(i+1):(dim(dfilr)[2])])
      }
    }
    
    dfilr = t(dfilr)
    
    for (i in 1:(dim(bcov)[1])){
      dfilr[2*i] = sqrt(bcov[i,i]) # standard deviation as sqrt of var
    }
    dfilr = t(dfilr)
    
    rownames(dfilr) = rname
    colnames(dfilr) = cnames
    
    
    print("      B (COEFFICIENTS IN THE SIMPLEX):      ")
    print(as.matrix(dfsimplex))
    print("--------------------------------------------")
    print("     B* (COEF. and SD IN THE ILR SPACE):    ")
    print(as.matrix(dfilr))
    print("--------------------------------------------")
    print(paste("n obs  : ", dim(tx)[1], sep="")) 
    print(paste("SSR    : ", sum(resid^2), sep="")) # sum of squared residuals
  }
}
