\name{alr_yx_reg}
\alias{alr_yx_reg}
\title{
Compute the yx-compositional regression using an ALR tranformation.
}
\description{
Estimate the coefficients of the parameters associated to the yx-compositional regression both in the simplex space (invariant to the transformation) and in the ALR coordinates space. It also returns the variance-covariance matrix of the parameters in the ALR coordinates space.}
\usage{
alr_yx_reg(Y, X, Z = list(matrix()), R = list(NULL),
    constant = TRUE, verbose = FALSE, pres = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{Y}{
  It is a matrix with data from the endogenous variable, which is compositional. Each column represents one part of the variable and rows its different observations. For example, the element in the row i and column j of Y represents the value observed for the statistical unit i with respect to part j. This argument cannot be missing.
}
  \item{X}{
  It is a list of matrices that contains the explanatory variables that are compositional. Columns represent parts of the variables and rows observations. This argument cannot be missing.
}
  \item{Z}{
  It is a list of matrices that contains the explanatory variables that are non-compositional (may or not be multivariable). Columns represent parts of the variables and rows observations. This argument may be missing.
}
  \item{R}{
  It is a list. Each element of R designs which part of the compositional variable is desired to be taken in the ALR transformation as referece for each variable in Y and the list X, in order. For example, if X=list(X1,X2,X3) and R=list(4,2,3,1), the reference variable for Y is the 4th, whereas for X1 is the second one, for X2 is the third one and, for X3, the first one. By default, the last part of each variable is taken as reference. 
}
  \item{constant}{
  It is a logical variable. 'constant'=TRUE if an intercept is desired to be added to the regression, 'constant'=FALSE if not.
}
  \item{verbose}{
  It is a logical variable. 'verbose' = TRUE and some intermediate prints are displayed. 
}
  \item{pres}{
  It is a logical variable. 'pres' = TRUE and the results of the regression are displayed in a nice way, instead of returning the named list of results that are used to compute impacts.
}
}

\value{
If pres=TRUE a named list is returned with the following components:

\item{Y_coord}{Matrix of the endogeneous variable data in the ALR coordinate space.}
\item{X_coord}{Matrix of the exogenous variables (both X and Z) data (only one matrix) in the ALR coordinate space.}
\item{constant}{TRUE if an intercept is added, FALSE otherwise.}
\item{B_coord}{Matrix of estimated coefficients of the regression in the ALR coordinate space.}
\item{B_simplex}{Matrix of estimated coefficients of the regression in the simplex space.}
\item{Bcoord_cov}{Variance-covariance matrix of the estimation of the parameters.}
\item{residuals_coord}{It is a matrix that contains the residuals of the regression in the ALR coordinates.}
\item{fitted_v_coord}{It is a matrix that contains the fitted values of the regression in the ALR coordinates.}

This returned list has also an attribute "type_reg" which is "alr_yx_reg" in order to keep the traceability for the main 'impacts' function.

If pres=FALSE nothing is returned, only printed.
}
\references{
Morais, J., & Thomas-Agnan, C. (2021). Impact of covariates in compositional models and simplicial derivatives. Austrian Journal of Statistics, 50(2), 1-15.}
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

# compositional exogenous variable
x_var = rbind(c(0.25,0.3,0.2),
    c(0.5,0.2,0.3),
    c(0.9,0.05,0.05),
    c(0.33,0.4,0.24),
    c(0.2,0.3,0.5),
    c(0.91,0.05,0.04),
    c(0.1,0.7,0.2))  

# standard exogenous variable
z_var = matrix(c(23,1,32,12,9,11,3)) 

alr_yx_reg(y_var,list(x_var),list(z_var))
}


