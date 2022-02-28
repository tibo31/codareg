\name{inverseilr}
\alias{inverseilr}

\title{
Inverse of the ILR transformation.
}
\description{
This function computes the ILR inverse transformation that brings to the simplex any data in the ILR coordinates.}
\usage{
inverseilr(x, V)
}
\arguments{
  \item{x}{
Matrix with data in the ILR coordinates. Each row is a statistical unit whereas each column is a different variable in the coordinates. This matrix should feature the attributte "space" with value "ilr_coord" or not to feature the attributte at all.}
  \item{V}{
It is the contrast matrix that has been taken for the first ILR transformation. If V is not given and the matrix 'x' features one matrix as the attribute 'contrast', that matrix is used. Otherwise, a contrast matrix build through the Helmert matrices is used.}
}

\details{
Of course, inverseilr(ilr(X)) == closure(X) for any matrix X. 
}
\value{
It is returned a matrix that lies in the simplex that is the ILR inverse transformation of the input data. This matrix features one attribute "space", whose value is "simplex".}
\references{
van den Boogaart, K. G., & Tolosana-Delgado, R. (2013). Fundamental concepts of compositional data analysis. In Analyzing Compositional Data with R (pp. 13-50). Springer, Berlin, Heidelberg.}
\author{
Ivan Rendo Barreiro
}

\examples{
simplex_data = rbind(c(0.1,0.4,0.25,0.25),
    c(0.3,0.3,0.3,0.1),
    c(0.9,0.05,0.02,0.03),
    c(0.1,0.7,0.15,0.05),
    c(0.8,0.1,0.05,0.05),
    c(0.35,0.25,0.3,0.1),
    c(0.4,0.1,0.1,0.4)) 
    
inverseilr(simplex_data)
}


