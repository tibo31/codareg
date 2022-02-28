alr_yx_reg <-
function(Y,X,Z=list(matrix()),R=list(NULL),constant=TRUE, verbose=FALSE, pres=FALSE){ # Y dependent, X indep compositional, Z indep. not compositional, c if constant is included, V list of reference variables order: (Y, X1, X2, ...) 
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
  
  
  if (!is.null(R[[1]])){
    Y=permutation_D(Y,R[[1]])
    for (j in 1:length(X)){
      X[[j]]=permutation_D(X[[j]],R[[j+1]])
    input = var2alr(Y,X,Z,constant=constant)
    }
  }
  else{
    input = var2alr(Y,X,Z,constant=constant)
    R=input$R
  }
  
  ty = input$Y
  tx = input$X
  if (sum(is.na(ty))>0){stop("There are NAs in the Y variable.")}
  if (sum(is.na(tx))>0){stop("There are NAs in the X variable.")}
  dimy = dim(ty)[2]+1
  constant = input$constant
  
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
    b0=inversealr(t(Bstar[1,]),R[[1]])
    Bs = rbind(Bs, b0)
    if (verbose){
    print("b0 in the simplex:")
    print(b0)}
  }
  h=2
  if (verbose){print("---------------------------")}
  for (i in 1:length(X)){
    if (verbose){print(paste("B (in the simplex) for the comp. variable, X", i, sep=""))}
    # number of parts?
    if (is.null(attr(X[[i]],"space"))){dimx = dim(X[[i]])[2]} #number of parts of this variable i of X D parts --> D-1 param in the alr 
    if (!is.null(attr(X[[i]],"space"))){
      if (attr(X[[i]],"space")=="simplex"){
        dimx = dim(X[[i]])[2]
      }
      else{dimx = dim(X[[i]])[2]+1}
    }
    if (h==(h+dimx-2)){matrix_x = t(matrix(Bstar[h,]))}
    if (h!=(h+dimx-2)){matrix_x = Bstar[h:(h+dimx-2),] }
    B_q = K_D(dimy)%*%t(matrix_x)%*%F_D(dimx)
    B_q = inv_permutation_D(B_q,R[[i+1]])
    Bs= rbind(Bs, t(B_q))
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
      c_k = inversealr(matrix_z, R[[1]]) # change R[[1]] with the other...
      Bs = rbind(Bs, c_k)
      if (verbose){print(c_k)}
    }
  }
  Bs = Bs[-1,]
  attr(Bstar,"space")="alr_coord"
  attr(Bstar,"reference") = R
  
  attr(tx,"space") = "alr_coord"
  attr(tx,"reference") = R[-1]
  
  fittedvalues = tx%*%Bstar
  resid = ty- fittedvalues
  attr(resid,"space") = NULL
  attr(resid,"contrast") = NULL
  
  rownames(Bs) = x_col_names
  colnames(Bs) = y_col_names
  attr(Bs,"space")="simplex"
  
  reslist = list(Y_coord = ty, X_coord = tx, constant = constant, B_coord = Bstar, B_simplex = Bs , Bcoord_cov=bcov, residuals_coord = resid, fitted_v_coord = fittedvalues)
  attr(reslist,"reg_type") = "alr_yx_reg"
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
    for (i in (1:(dim(Y)[2]-1))){
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
    dimalr = dim(Y)[2]-1
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
    print("     B* (COEF. and SD IN THE COORDINATES):    ")
    print(as.matrix(dfalr))
    print("--------------------------------------------")
    print(paste("n obs  : ", dim(tx)[1], sep="")) 
    print(paste("SSR    : ", sum(resid^2), sep="")) # sum of squared residuals
  }
}
