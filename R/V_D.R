V_D <-
function(D){ # D is the dimension (parts)
  V_D = matrix(0L, nrow =D , ncol = D-1)
  for (i in 1:(D-1)){  # pseudodiagonal elements
    V_D[i,i] = (D-i)/(sqrt((D-i+1)*(D-i)))
  }
  for (i in 2:D){ # below pseudodiagonal elements
    for (j in 1:i-1){
      V_D[i,j]=-1/(sqrt((D-j+1)*(D-j)))
    }
  } 
  return(V_D)
}
