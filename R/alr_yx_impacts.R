alr_yx_impacts <-
function(results){ # the unique argument is the named list that comes from the alr_yx_reg function
    x_names = rownames(results$B_simplex)
    y_names = colnames(results$B_simplex)
    constant = results$constant
    tx = results$X_coord
    B = results$B_coord
    exp_Yilr = tx%*%B
    Ry = attr(B,"reference")[[1]]
    exp_Y_simplex = inversealr(exp_Yilr,Ry)
    I_D = diag(dim(exp_Y_simplex)[2])
    N=dim(exp_Y_simplex)[1]
    Bsimplex = results$B_simplex
    impacts = list()
    for (i in (1:N)){
      W_i = I_D - matrix(1,dim(exp_Y_simplex)[2],1)%*%exp_Y_simplex[i,]
      imp_i = t(W_i%*%t(Bsimplex))
      colnames(imp_i) = y_names
      rownames(imp_i) = x_names
      imp_i = imp_i[-1,]
      impacts[[i]] = imp_i
      }
    return(impacts)
}
