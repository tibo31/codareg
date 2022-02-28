P_D <-
function(D){
  return(rbind(diag(D-1),numeric(D-1)))
}
