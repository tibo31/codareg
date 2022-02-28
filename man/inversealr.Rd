\name{inversealr}
\alias{inversealr}
\title{
Inverse of the ALR transformation.}
\description{
This function computes the ALR inverse transformation that brings to the simplex any data in the ALR coordinates.
}
\usage{
inversealr(x, R = dim(x)[2] + 1)
}
\arguments{
  \item{x}{
Matrix with data in the ALR coordinates. Each row is a statistical unit whereas each column is a different variable in the coordinates. This matrix should feature the attributte "space" with value "alr_coord" or not to feature the attributte at all.}
  \item{R}{
It is an intenger that designs which variable has been taken as reference for the ALR transformation. By default it is understood that the last variable (of the simplex) has been used.
}
}
\details{Of course, inversealr(alr(X)) == closure(X) for any matrix X. 
}
\value{
It is returned a matrix that lies in the simplex that is the ALR inverse transformation of the input data. This matrix features one attribute "space", whose value is "simplex".
}
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
    
inversealr(simplex_data)
}
