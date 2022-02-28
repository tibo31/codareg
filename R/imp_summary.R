imp_summary <-
function(impact_results){
  ynames = colnames(impact_results[[1]])
  xnames = rownames(impact_results[[1]])
  ncols = length(ynames)*length(xnames)
  xaxis = numeric(ncols)
  cc=1
  for (i in xnames){
    for (j in ynames){
      xaxis[cc] = paste(i,j,sep=" -> ") #to?
      cc=cc+1
    }
  }
  resmat = matrix(0,ncols,6)
  rownames(resmat) = xaxis
  colnames(resmat) = c("Mean", "Min.", "q1 (25%)", "q2 (50%)", "q3 (75%)", "Max.")
  qq=1
  for (i in 1:length(xnames)){
    for (j in 1:length(ynames)){
      values = numeric(length(impact_results))
      for (n in 1:length(impact_results)){
        values[n] = impact_results[[n]][i,j]
      }
      resmat[qq,1] = mean(values)
      resmat[qq,2] = min(values)
      resmat[qq,3] = quantile(values,0.25)
      resmat[qq,4] = quantile(values,0.5)
      resmat[qq,5] = quantile(values,0.75)
      resmat[qq,6] = max(values)
      qq=qq+1
    }
  }
  return(resmat)
}
