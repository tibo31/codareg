\name{ilr_x_impacts}
\alias{ilr_x_impacts}

\title{
Compute impacts of covariates in the x-compositional model through an ILR transformation.
}
\description{
This function takes the results of the already computed regression as an input and computes the impacts of each covariate. This function is for the ILR transformation and x-compositional model case.}
\usage{
ilr_x_impacts(results)
}

\arguments{
  \item{results}{'results' object is a named list that is the output of ilr_x_reg
}
}
\details{
It is important to use the 'results' object obtained from the 'ilr_x_reg(...)' function and not from another. }
\value{
It is returned a list of matrices. Each matrix corresponds to each compositional variable in X. The element from row 'i' and column 'j' in the matrix which is in position 'k' in the list corresponds to the computed impact that the j-th part of the k-th compositional variable has on the i-th variable of the endogenous. (Y should be non-compositional, but may be multivariable). 
}
\references{
Morais, J., & Thomas-Agnan, C. (2021). Impact of covariates in compositional models and simplicial derivatives. Austrian Journal of Statistics, 50(2), 1-15.}
\author{
Ivan Rendo Barreiro
}

\examples{
# standard endogenous variable
y_var = matrix(c(23,1,32,12,9,11,3)) 
# compositional explicative variable
x_var = rbind(c(0.1,0.4,0.25,0.25),
    c(0.3,0.3,0.3,0.1),
    c(0.9,0.05,0.02,0.03),
    c(0.1,0.7,0.15,0.05),
    c(0.8,0.1,0.05,0.05),
    c(0.35,0.25,0.3,0.1),
    c(0.4,0.1,0.1,0.4)) 

# results from the x-compositional regression
results = ilr_x_reg(y_var,list(x_var)) 

ilr_x_impacts(results)}

