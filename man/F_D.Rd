\name{F_D}
\alias{F_D}
\title{
Computation of matrix F_D.}
\description{
This (D-1)xD matrix 'F_D' is needed to carry out some computations with respect the ALR transformation. This function just creates this matrix, that only depends on its dimension.}
\usage{
F_D(D)
}
\arguments{
  \item{D}{Number of columns of the returned matrix, which is a (D-1)xD matrix.}
}

\value{
It is returned the matrix F_D = (I_D|-1) where I_D is the identity of dimension D and -1 denotes a column vector (D-1)x1.
}
\references{
WHICH?
}
\author{
Ivan Rendo Barreiro}

\examples{
F_D(4)
}

