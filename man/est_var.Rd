\name{est_var}
\alias{est_var}
\title{
Estimation of the common covariance matrix.}
\description{
It is estimated the common covariance matrix variance following equation (5.24) in Kim & Timm (2006). This function is used as an intermediate function in order to obtain the variance-covariance matrix of the parameter estimation.}
\usage{
est_var(X, Y, B)
}
\arguments{
  \item{X}{Matrix with the exogenous variables of the regression. 
}
  \item{Y}{Matrix with the endogenous variable of the regression. 
}
  \item{B}{Matrix of coefficients of the regression. If missing, it is used B=mlm(Y,X), the one obtained in through the usual regression. 
}
}

\value{
It is returned a DxD matrix (D is the number of parts of the endogenous variable), where the element (i,j) desings the covariance between parts i and j.
}

\references{
Kim, K., & Timm, N. (2006). Univariate and multivariate general linear models: theory and applications with SAS. CRC Press.}

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

dataY = rbind(c(0.2,0.8),
    c(0.3,0.7),c(0.9,0.1),
    c(0.75,0.25),c(0.3,0.7),
    c(0.15,0.85),c(0.9,0.1))

est_var(dataX,dataY)}

