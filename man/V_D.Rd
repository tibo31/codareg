\name{V_D}
\alias{V_D}
\title{
Construction of the helmertian contrast matrices}
\description{
This function just construct the Dx(D-1) contrast matrix as shown in the equation (2.11) of the reference, following the helmertian approach.}
\usage{
V_D(D)
}
\arguments{
  \item{D}{ Number of rows of the desired contrast matrix, whcih is Dx(D-1).
}
}

\value{
Returns the Dx(D-1) contrast matrix based on the helmertian approach.
}
\references{
van den Boogaart, K. G., & Tolosana-Delgado, R. (2013). Fundamental concepts of compositional data analysis. In Analyzing Compositional Data with R (pp. 13-50). Springer, Berlin, Heidelberg.}
\author{
Ivan Rendo Barreiro
}

\examples{
V_D(4)
}

