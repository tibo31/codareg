impacts <-
function(results){
  if (attr(results,"reg_type")=="alr_x_reg"){return(alr_x_impacts(results))}
  if (attr(results,"reg_type")=="alr_y_reg"){return(alr_y_impacts(results))}
  if (attr(results,"reg_type")=="alr_yx_reg"){return(alr_yx_impacts(results))}
  if (attr(results,"reg_type")=="ilr_x_reg"){return(ilr_x_impacts(results))}
  if (attr(results,"reg_type")=="ilr_y_reg"){return(ilr_y_impacts(results))}
  if (attr(results,"reg_type")=="ilr_yx_reg"){return(ilr_yx_impacts(results))}
  else{stop("Input in 'impacts()' function is not a valid output of one of the regression functions.")}
}
