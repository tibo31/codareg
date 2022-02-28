check_formula <-
function(dataset, formula){ # returns a list of X, a list of Y, a list of Z and the type of model that is specified. It also checks that every value is positive, etc (ask Thibault)
  constant = as.logical(attr(terms(formula),"intercept")) 
  compositional_y = FALSE
  ytext = as.character(formula)[2]
  if (substring(ytext,1,6) == "cbind("){
    compositional_y = TRUE
  }
  vars_y = all.vars(formula[[2]])
  Y = as.matrix(dataset[vars_y[[1]]])
  if (length(vars_y)>1){
    for (i in 2:length(vars_y)){
      Y = cbind(Y, as.matrix(dataset[vars_y[[i]]]))  
    }
  } 
  
  compositional_x = FALSE
  xtext = as.character(formula)[3]
  splitted = strsplit(xtext, "\\+|\\-")[[1]]
  numx=1
  numz=1
  X = list() #compo vars
  Z = list() #standard vars
  neg = sum(Y<=0)
  nonumeric = sum(!(is.numeric(Y)))
  nonumeric = nonumeric + sum((is.nan(Y)))
  nonumeric = nonumeric + sum((is.na(Y)))
  
  for (i in splitted){
    temp_text = trimws(i)
    if (!(temp_text=="1" | temp_text=="0")){
        if (!(substring(temp_text,1,6)=="cbind(")){
          Z[[numz]]=as.matrix(dataset[temp_text])
          nonumeric = nonumeric + sum(!(is.numeric(Z[[numz]])))
          nonumeric = nonumeric + sum((is.na(Z[[numz]])))
          nonumeric = nonumeric + sum((is.nan(Z[[numz]])))
          numz = numz + 1
        }
        else{
          compo_text = strsplit(substring(temp_text,7,nchar(temp_text)-1),",")[[1]]
          compo_mat = as.matrix(dataset[trimws(compo_text[1])])
          for (k in 2:length(compo_text)){
            compo_mat = cbind(compo_mat, as.matrix(dataset[trimws(compo_text[k])]))
          }
          X[[numx]] = compo_mat
          neg = neg + sum(X[[numx]]<=0)
          nonumeric = nonumeric + sum(!(is.numeric(X[[numx]])))
          nonumeric = nonumeric + sum((is.nan(X[[numx]])))
          nonumeric = nonumeric + sum((is.na(X[[numx]])))
          numx = numx + 1 
        }
    }
  }
  
  # no negative or 0 values
  if (nonumeric>0){stop("Variables should be numeric and missing values are not admitted.")}
  if (neg>0){stop("Values of endogenous and exogenous variables which are compositional should be strictly positive.")}
  
  # which compositional model?
  
  if (compositional_y == TRUE){
    if (length(X)>0){
      cmodel = "yx_reg"
      if (length(Z)==0){Z=list(matrix())}
    }
    else{
      cmodel = "y_reg"
      X=Z
      Z=list(matrix())
    }
  }
  else{
    if (length(X)>0){
      cmodel = "x_reg"
      if (length(Z)==0){Z=list(matrix())}
    }
    else{
      stop("At least one compositional variable should be used in the model.")
    }
  }
  
  return(list(Y=Y,X=X,Z=Z, model_type = cmodel, constant = constant))
}
