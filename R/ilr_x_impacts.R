ilr_x_impacts <-
function(results){ # the unique argument is the named list that comes from the ilr_x_reg function
    B = results$B_coord
    xnames = rownames(results$B_simplex)
    ynames = colnames(results$B_simplex)
    constant = results$constant
    if (constant){xnames = xnames[-1]}
    Vx = attr(B,"contrast")
    if (constant){
      if (dim(B)[2]==1){B=matrix(B[-1,])}
      else{
        B = B[-1,]
      }
    }
    h=1
    impacts=list()
    for (i in 1:length(Vx)){
      dimx = dim(Vx[[i]])[1]
      if (is.null(dim(B))){B=t(matrix(B))}
      bstar = B[h:(h+dimx-2),]
      h=h+dimx-1
      impacts[[i]]=Vx[[i]]%*%bstar
    }
    cc=1
    for (k in (1:length(impacts))){
      val = impacts[[k]]
      ldimx = dim(val)[1]
      colnames(val) = ynames
      rownames(val) = xnames[cc:(cc+ldimx-1)]
      cc = cc + ldimx
      impacts[[k]] = val
    }
    ## part of pooling impacts into a single matrix 
    impacts_mat = impacts[[1]]
    if (length(impacts)>1){
      for (j in (2:length(impacts))){
        impacts_mat = rbind(impacts_mat, impacts[[j]])
      }
    }
    return(impacts_mat)
    return(impacts)
}
