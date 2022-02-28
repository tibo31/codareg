\name{mlm}
\alias{mlm}
\title{
Computes the estimated coefficients of the multivariate GLM.
}
\description{
This function computes the BLUE estimator of the coefficients of the parameters for the multivariate general model as B_hat = (X'X)^(-1)X'Y as shown in the reference in equation (5.13).
}
\usage{
mlm(Y, X)
}
\arguments{
  \item{Y}{
It is a matrix that contains the endogenous variable data, which could be multivariate. Each row is an observation and each column a variable.}
  \item{X}{
  It is a matrix that contains the exogenous variables data, which could be multivariate. Each row is an observation and each column a variable.
}
}
\details{ If Y is univariate, this model corresponds to the usual multiple linear regression one.
}
\value{
The function returns the matrix B of estimated coefficients of the parameters.
}
\references{ Kim, K., & Timm, N. (2006). Univariate and multivariate general linear models: theory and applications with SAS. CRC Press.}

\author{
Ivan Rendo Barreiro}

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

mlm(dataY,dataX) 
}
