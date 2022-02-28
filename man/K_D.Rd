\name{K_D}
\alias{K_D}
\title{
Computation of matrix K_D.}
\description{
This Dx(D-1) matrix 'K_D' is needed to carry out some computations with respect the inverse ALR transformation and the impacts. This function just creates this matrix, that only depends on its dimension.
}
\usage{
K_D(D)
}
\arguments{
  \item{D}{
Number of rows of the matrix, which is Dx(D-1).}
}

\value{
It is returned the matrix K_D which is defined in the reference.
}
\references{
Chakir, R., Laurent, T., Ruiz-Gazen, A., Thomas-Agnan, C. & Vignes, C. (2017). Land Use Predictions on a Regular Grid at Different Scales and with Easily Accessible Covariates. Revue economique, 68, 435-469. https://doi.org/10.3917/reco.683.0435}
\author{Ivan Rendo Barreiro
}

\examples{
K_D(5)
}

