K_D <-
function(D){
  upm = diag(D-1) - (numeric(D-1)+1)%*%t(numeric(D-1)+1)/D
  down = -t(numeric(D-1)+1)/D
  upm = rbind(upm,down)
  return (upm)
}
