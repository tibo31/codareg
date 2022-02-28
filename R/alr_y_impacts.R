alr_y_impacts <-
function(results){ #results argument regers to the output of ilr_y_reg
  # first we need EY*
  constant = results$constant
  Bcoord = results$B_coord
  coln = colnames(results$B_simplex)
  rown = rownames(results$B_simplex)
  R = attr(results$Y_coord, "reference")
  xcoord = results$X_coord
  EYcoord = xcoord%*%Bcoord
  EYsimpl = inversealr(EYcoord,R)
  D = dim(EYcoord)[2]+1
  impacts=list()
  for (i in 1:(dim(EYsimpl)[1])){
    Wi = diag(D) - matrix(1,D,1)%*%t(EYsimpl[i,]) # Wi are the same in both cases
    Wistar = Wi%*%P_D(D) #this is the only important thing that changes
    tempres = t(Wistar%*%t(Bcoord))
    colnames(tempres) = coln
    rownames(tempres) = rown
    if (constant){
      if (dim(tempres)[1]==2){ 
          tempres = t(as.matrix(tempres)[-1,])
          rownames(tempres) = c(rown[[2]])
        }
        else{
          tempres = tempres[-1,]
        }
      }
    
    impacts[[i]] = tempres
  }
  return(impacts)
}
