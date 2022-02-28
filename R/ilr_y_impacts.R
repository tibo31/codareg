ilr_y_impacts <-
function(results){ #results argument regers to the output of ilr_y_reg
  # first we need EY*
  ccnames = colnames(results$B_simplex)
  rcnames = rownames(results$B_simplex)
  constant = results$constant
  Bcoord = results$B_coord
  Vy = attr(results$Y_coord, "contrast")
  xcoord = results$X_coord
  EYcoord = xcoord%*%Bcoord
  EYsimpl = inverseilr(EYcoord,Vy)
  D = dim(Vy)[1]
  impacts=list()
  for (i in 1:(dim(EYsimpl)[1])){
    Wi = diag(D) - matrix(1,D,1)%*%t(EYsimpl[i,])
    imp = t(Wi%*%Vy%*%t(Bcoord))
    colnames(imp) = ccnames
    rownames(imp) = rcnames
    if (constant){
        if (dim(imp)[1]==2){ 
          imp = t(as.matrix(imp)[-1,])
          rownames(imp) = c(rcnames[[2]])
        }
        else{
          imp = imp[-1,]
        }
      }
    impacts[[i]] = imp
  }
  return(impacts)
}
