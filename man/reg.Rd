\name{reg}
\alias{reg}
\title{
General regression function.
}
\description{
This function performs the estimation of coefficients in models where compositional variables are involved (X-compositional, Y-compositional and YX-compositional).}
\usage{
reg(dataset, formula, transformation = "ILR", pres = FALSE, V = list(matrix()), R)
}
\arguments{
  \item{dataset}{Dataset in dataframe format.
}
  \item{formula}{Formula of the desired model. Compositional variables should be inserted as a concatenation of their parts thorugh 'cbind()' and standard ones, although multivariate, with a plus (+) sign between all of them.
}
  \item{transformation}{
This argument should be "ILR" (by default) or "ALR" and indicates which type of transformation should be employed.
}
  \item{pres}{
  It is a logical variable. 'pres' = TRUE and the results of the regression are displayed in a nice way, instead of returning the named list of results that are used to compute impacts.
}
  \item{V}{
 (Only ised if transformation="ILR"). It is a list. Each element of V designs which contrast matrices are desired to be employed in the ILR transformation. They should be declared in order, starting with the contrast matrix for Y. For example, if X=list(X1,X2,X3) and V=list(V1,V2,V3,V4), the contrast matrix for Y is V1, for X1 is V2, for X2 is V3, and  for X3 is V4. By default, if no list V is determined, contrast matrices constructed using the Helmert transformation are used.
}
  \item{R}{
(Only used if transformation="ALR"). It is a list. Each element of R designs which part of the compositional variable is desired to be taken in the ALR transformation as referece for each variable in Y and the list X, in order. For example, if X is composed by X1, X2 and X3 and R=list(4,2,3,1), the reference variable for Y is the 4th, whereas for X1 is the second one, for X2 is the third one and, for X3, the first one. By default, the last part of each variable is taken as reference.
}
}

\value{
If pres=TRUE a named list is returned with the following components:

\item{Y_coord}{Matrix of the endogeneous variable data in the coordinate space.}
\item{X_coord}{Matrix of the exogenous variables (both X and Z) data (only one matrix) in the ILR coordinate space.}
\item{constant}{TRUE if an intercept is added, FALSE otherwise.}
\item{B_coord}{Matrix of estimated coefficients of the regression in the ALR coordinate space.}
\item{B_simplex}{Matrix of estimated coefficients of the regression in the simplex space.}
\item{Bcoord_cov}{Variance-covariance matrix of the estimation of the parameters.}
\item{residuals_coord}{It is a matrix that contains the residuals of the regression in the coordinates.}
\item{fitted_v_coord}{It is a matrix that contains the fitted values of the regression in the coordinates.}

This returned list has also an attribute "type_reg" which indicates which type of model has been used in order to keep the traceability for the 'impacts' function.

If pres=FALSE nothing is returned, only printed.
}

\author{
Ivan Rendo Barreiro
}

\seealso{
See the vignette.
}
