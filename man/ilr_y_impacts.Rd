\name{ilr_y_impacts}
\alias{ilr_y_impacts}
\title{
Compute impacts of covariates in the y-compositional model through an ILR transformation.}
\description{
This function takes the results of the already computed regression as an input and computes the impacts of each covariate. This function is for the ILR transformation and y-compositional model case. 
}
\usage{
ilr_y_impacts(results)
}

\arguments{
  \item{results}{'results' object is a named list that is the output of ilr_y_reg(...).
}
}
\details{
It is important to use the 'results' object obtained from the 'ilr_y_reg(...)' function and not from another. }
\value{
It is returned a list of matrices. Each matrix corresponds to one statistical unit. The element from row 'i' and column 'j' in the matrix which is in position 'k' in the list corresponds to the computed impact that the i-th exogenous variable has on the j-th endogenous variable for the k-th statistical unit. The first row of each matrix corresponds to the impacts that derives from the intercept, if it exists. (Y should be compositional and X should be non-compositional, but may be multivariable). 
}
\references{
Morais, J., & Thomas-Agnan, C. (2021). Impact of covariates in compositional models and simplicial derivatives. Austrian Journal of Statistics, 50(2), 1-15.}
\author{
Ivan Rendo Barreiro}


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
x_var2 = rbind(c(3,-4.3),c(0,-2),c(4,5),
    c(-2,3), c(4,2), c(1,11), c(23,2)) 

results = ilr_y_reg(y_var,list(x_var1, x_var2))

ilr_y_impacts(results)}


