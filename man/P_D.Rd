\name{P_D}
\alias{P_D}
\title{
Computation of matrix P_D.}
\description{
This Dx(D-1) matrix 'P_D' is needed to carry out some computations with respect the inverse ALR transformation and the impacts. This function just creates this matrix, that only depends on its dimension.
}
\usage{
P_D(D)
}
\arguments{
  \item{D}{
Number of rows of the matrix, which is Dx(D-1).}
}

\value{
It is returned the matrix P_D which is the identity of dimension D-1 plus a last row of zeros which is a row vector 1x(D-1).
}

\author{
Ivan Rendo Barreiro}

\examples{
P_D(5)}
