\name{b_cov}
\alias{b_cov}
\title{
Compute the variance-covariance matrix of the estimation of the parameters.
}
\description{
Compute the variance-covariance matrix of the estimation of the parameters for any compositional model.}
\usage{
b_cov(X, Y, B)
}
\arguments{
  \item{X}{Matrix with the exogenous data where columns are variables and rows are observations.}
  \item{Y}{
Matrix with the endogenous data where columns are variables and rows are observations.}
  \item{B}{Matrix of estimated coefficients of the regression.}
}
\details{
This is an internal function and it is better to use the correspondant regression function.}
\value{
Returns the usual symmetric variance-covariance matrix where the element (i,j) corresponds to the covariance between parameter "i" estimation and parameter "j" estimation.
}
\references{
Kim, K., & Timm, N. (2006). Univariate and multivariate general linear models: theory and applications with SAS. CRC Press.}
\author{
Ivan Rendo Barreiro
}

\examples{
dataX = rbind(c(0.1,0.4,0.25,0.25),
    c(0.3,0.3,0.3,0.1),
    c(0.9,0.05,0.02,0.03),
    c(0.1,0.7,0.15,0.05),
    c(0.8,0.1,0.05,0.05),
    c(0.35,0.25,0.3,0.1),
    c(0.4,0.1,0.1,0.4)) 

dataY = rbind(c(0.2,0.8),c(0.3,0.7),
    c(0.9,0.1),c(0.75,0.25),
    c(0.3,0.7),c(0.15,0.85),c(0.9,0.1))

# this is to obtain the matrix of est. coefficients
B = mlm(dataY, dataX) 

b_cov(dataX,dataY,B)}
