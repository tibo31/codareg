reg <-
function(dataset, formula, transformation = "ILR", pres=FALSE, V=list(matrix()), R){
  input = check_formula(dataset, formula)
  Y = input$Y
  X = input$X
  Z = input$Z
  constant = input$constant
  model_type = input$model_type
  if (transformation == "ILR"){
    if (model_type == "yx_reg"){
      if (pres){
      print("YX-compositional model. ILR transformation. ")
      print("------------------------------------------- ")
      }
      return(ilr_yx_reg(Y=Y,X=X,Z=Z,V=V,constant = constant, pres = pres))
    }
    if (model_type == "x_reg"){
      if (pres){
      print("X-compositional model. ILR transformation.  ")
      print("------------------------------------------- ")
      }
      return(ilr_x_reg(Y=Y,X=X,Z=Z,V=V,constant=constant,pres=pres))
    }
    if (model_type == "y_reg"){
      if (pres){
      print("Y-compositional model. ILR transformation.  ")
      print("------------------------------------------- ")
      }
      return(ilr_y_reg(Y=Y,X=X,V=V,constant=constant,pres=pres))
    }
  }
  if (transformation == "ALR"){
    if (model_type == "yx_reg"){
      if (pres){
      print("YX-compositional model. ALR transformation. ")
      print("------------------------------------------- ")
      }
      if (missing(R)){R=list(NULL)}
      return(alr_yx_reg(Y=Y,X=X,Z=Z,R=R,constant=constant,pres=pres))
    }
    if (model_type == "x_reg"){
      if (pres){
      print("X-compositional model. ALR transformation.  ")
      print("------------------------------------------- ")
      }
      if (missing(R)){R=list(NULL)}
      return(alr_x_reg(Y=Y,X=X,Z=Z,R=R,constant=constant,pres=pres))
    }
    if (model_type == "y_reg"){
      if (pres){
      print("Y-compositional model. ALR transformation.  ")
      print("------------------------------------------- ")
      }
      if (missing(R)){R=NaN}
      return(alr_y_reg(Y=Y,X=X,R=R,constant=constant,pres=pres))
    }
  }
}
