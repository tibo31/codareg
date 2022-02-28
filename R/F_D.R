F_D <-
function(D){
  fd = diag(D-1)
  fd = cbind(fd, numeric(D-1)-1)
  return(fd)
}
