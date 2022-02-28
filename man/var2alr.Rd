\name{var2alr}
\alias{var2alr}
\title{
Transform every variable to ALR coordinates.
}
\description{
This function is just an intermediate function that takes every variable to be used in a regression and transform it to the coordinates, if the variable is not already in coordinates. Moreover, bundles into a single matrix all the exogenous variables to be able to use it as an input for the mlm(...) function inside alr_yx_reg(...).}
\usage{
var2alr(Y, X, Z = list(matrix()), R = list(matrix()), constant = TRUE)
}
\arguments{
  \item{Y}{
  It is a matrix with data from the endogenous variable, which should be compositional. Each column represents one part of the variable and rows its different observations. For example, the element in the row i and column j of Y represents the value observed for the statistical unit i with respect to part j. This argument cannot be missing.}
  \item{X}{
  It is a list of matrices that contains the explanatory variables that are compositional. Columns represent parts of the variables and rows observations. This argument cannot be missing.}
  \item{Z}{
 It is a list of matrices that contains the explanatory variables that are non-compositional (may or not be multivariable). Columns represent parts of the variables and rows observations. This argument may be missing.
}
  \item{R}{
  It is a list. Each element of R designs which part of the compositional variable is desired to be taken in the ALR transformation as referece for each variable in Y and the list X, in order. For example, if X=list(X1,X2,X3) and R=list(4,2,3,1), the reference variable for Y is the 4th, whereas for X1 is the second one, for X2 is the third one and, for X3, the first one. By default, the last part of each variable is taken as reference. 
}
  \item{constant}{
  It is a logical variable. 'constant'=TRUE if an intercept is desired to be added to the regression, 'constant'=FALSE if not.
}
}
\value{
It returns a named list with the following elements:

\item{Y}{A matrix with the endogenous variable transformed into the ALR coordinates.}
\item{X}{One matrix with the exogenous variables bundled (both compositional and non-compositional) in the ALR coordinates.}
\item{R}{It is a list with integers that denote the references that are eventually used inside the function.}
\item{constant}{A logical variable that returns TRUE if an intercept is considered and FALSE otherwise.}
}

\author{
Ivan Rendo Barreiro}

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

var2alr(var_y,list(var_x),list(var_z))

}

