\name{share_ratio_summary}
\alias{share_ratio_summary}
\title{
Share Ratios. 
}
\description{
This function computes the share ratios for each explanatory variable in the Y-compositional as shown in (...).}
\usage{
share_ratio_summary(results)
}
\arguments{
  \item{results}{Results of the concerned regression.
}
}

\value{
It is returned a list of matrices. Each list corresponds to one explanatory, non-compositional variable. Each term s_ij of the k-th matrix is computed as log(b_jk/b_ik). Where b_ik and b_jk are the coefficients associated to explanatory variable k and the j-th and i-th part of the exogenous variable, respectively.
}
\references{
London Slides (??)}
\author{
Ivan Rendo Barreiro
}

\examples{
# compositional endogenous variable
y_var = rbind(c(0.1,0.4,0.25,0.25),
    c(0.3,0.3,0.3,0.1),
    c(0.9,0.05,0.02,0.03),
    c(0.1,0.7,0.15,0.05),
    c(0.8,0.1,0.05,0.05),
    c(0.35,0.25,0.3,0.1),
    c(0.4,0.1,0.1,0.4)) 
    
# standard exogenous variable
x_var1 = matrix(c(23,1,32,12,9,11,3)) 
# standard exogenous variable, but multivariate
x_var2 = rbind(c(3,-4.3),c(0,-2),c(4,5),c(-2,3),
    c(4,2), c(1,11), c(23,2)) 

results = ilr_y_reg(y_var,list(x_var1,x_var2))
share_ratio_summary(results)
}

