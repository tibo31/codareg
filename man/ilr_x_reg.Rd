\name{ilr_x_reg}
\alias{ilr_x_reg}
\title{
Compute the x-compositional regression using an ILR tranformation.
}
\description{
Estimate the coefficients of the parameters associated to the x-compositional regression both in the simplex space (invariant to the transformation) and in the ILR coordinates space. It also returns the variance-covariance matrix of the parameters in the ILR coordinates space. }
\usage{
ilr_x_reg(Y, X, Z = list(matrix()), V = list(matrix()),
    constant = TRUE, verbose = FALSE, pres = FALSE)
}
\arguments{
  \item{Y}{
  It is a matrix with data from the standard variable (i.e. non-compositional) where columns represent different variables and rows different observations. For example, the element in the row i and column j of Y represents the value observed for the statistical unit i with respect to variable j. This argument cannot be missing.}
  \item{X}{
  It is a list of matrices that contains the explanatory variables that are compositional. Columns represent parts of the variables and rows observations. This argument cannot be missing.}
  \item{Z}{
  It is a list of matrices that contains the explanatory variables that are non-compositional (may or not be multivariable). Columns represent parts of the variables and rows observations. This argument can be missing.}
  \item{V}{
 It is a list. Each element of V designs which contrast matrix is desired to be employed in the ILR transformation. For example, if X=list(X1,X2,X3) and V=list(V1,V2,V3), the contrast matrix for X1 is V1, for X2 is V2, and  for X3 is V3. By default, if no list V is determined, contrast matrices constructed using the Helmert transformation are used.}
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
\details{
Note that it is not important which transformation is used in order to compute impacts or with respects to the results in the simplex space. However, as ALR and ILR coordiante spaces differ, results in this aspect are different. Results in the simplex or their impacts may differ if a severe problem of imperfect multicollinearity exists.}
\value{
If pres=TRUE a named list is returned with the following components:

\item{Y_coord}{Matrix of the endogeneous variable data in the ILR coordinate space.}
\item{X_coord}{Matrix of the exogenous variables (both X and Z) data (only one matrix) in the ILR coordinate space.}
\item{constant}{TRUE if an intercept is added, FALSE otherwise.}
\item{B_coord}{Matrix of estimated coefficients of the regression in the ILR coordinate space.}
\item{B_simplex}{Matrix of estimated coefficients of the regression in the simplex space.}
\item{Bcoord_cov}{Variance-covariance matrix of the estimation of the parameters.}
\item{residuals_coord}{It is a matrix that contains the residuals of the regression in the ILR coordinates.}
\item{fitted_v_coord}{It is a matrix that contains the fitted values of the regression in the ILR coordinates.}

This returned list has also an attribute "type_reg" which is "ilr_x_reg" in order to keep the traceability for the main 'impacts' function.

If pres=FALSE nothing is returned, only printed.
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

ilr_x_reg(y_var,list(x_var))
}
