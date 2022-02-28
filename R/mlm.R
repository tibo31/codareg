mlm <-
function(Y,X){ # only the dependent and independent matrices are needed
  ## use solve(t(X),Y) or QR algorithm qr.solve(), cholesky (I HAVE USED QR)
  B = solve(crossprod(X),crossprod(X,Y),tol=1e-21)
  return(B)
}
