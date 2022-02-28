\name{var2ilr}
\alias{var2ilr}
\title{
Transform every variable to ILR coordinates.
}
\description{
This function is just an intermediate function that takes every variable to be used in a regression and transform it to the coordinates, if the variable is not already in coordinates. Moreover, bundles into a single matrix all the exogenous variables to be able to use it as an input for the mlm(...) function inside ilr_yx_reg(...).
}
\usage{
var2ilr(Y, X, Z = list(matrix()), V = list(matrix()), constant = TRUE)
}
\arguments{
  \item{Y}{
  It is a matrix with data from the endogenous variable, which should be compositional. Each column represents one part of the variable and rows its different observations. For example, the element in the row i and column j of Y represents the value observed for the statistical unit i with respect to part j. This argument cannot be missing.}
  \item{X}{
  It is a list of matrices that contains the explanatory variables that are compositional. Columns represent parts of the variables and rows observations. This argument cannot be missing.}
  \item{Z}{
 It is a list of matrices that contains the explanatory variables that are non-compositional (may or not be multivariable). Columns represent parts of the variables and rows observations. This argument may be missing.
}
  \item{V}{
 It is a list. Each element of V designs which contrast matrices are desired to be employed in the ILR transformation. They should be declared in order, starting with the contrast matrix for Y. For example, if X=list(X1,X2,X3) and V=list(V1,V2,V3,V4), the contrast matrix for Y is V1, for X1 is V2, for X2 is V3, and  for X3 is V4. By default, if no list V is determined, contrast matrices constructed using the Helmert transformation are used. 
}
  \item{constant}{
  It is a logical variable. 'constant'=TRUE if an intercept is desired to be added to the regression, 'constant'=FALSE if not.
}
}

\value{
It returns a named list with the following elements:

\item{Y}{A matrix with the endogenous variable transformed into the ALR coordinates.}
\item{X}{One matrix with the exogenous variables bundled (both compositional and non-compositional) in the ILR coordinates.}
\item{V}{It is a list with the contrast matrices that are eventually used inside the function.}
\item{constant}{A logical variable that returns TRUE if an intercept is considered and FALSE otherwise.}
}

\author{
Ivan Rendo Barreiro
}

\examples{
var_y = rbind(c(0.25,0.3,0.2),
    c(0.5,0.2,0.3),
    c(0.9,0.05,0.05),
    c(0.33,0.4,0.24),
    c(0.2,0.3,0.5),
    c(0.91,0.05,0.04),
    c(0.1,0.7,0.2)) 
    
var_x = rbind(c(0.1,0.4,0.25,0.25),
    c(0.3,0.3,0.3,0.1),
    c(0.9,0.05,0.02,0.03),
    c(0.1,0.7,0.15,0.05),
    c(0.8,0.1,0.05,0.05),
    c(0.35,0.25,0.3,0.1),
    c(0.4,0.1,0.1,0.4)) 
    
var_z = matrix(c(23,1,32,12,9,11,3)) 

var2ilr(var_y,list(var_x),list(var_z))
}

