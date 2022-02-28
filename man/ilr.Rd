\name{ilr}
\alias{ilr}
\title{
ILR transformation
}
\description{
This function computes the ILR transformation to coordinates for data in the simplex space.}
\usage{
ilr(x, V)
}
\arguments{
  \item{x}{
Matrix with data in the simplex. Each row is a statistical unit whereas each column is a different part of the compositional variable. This matrix should feature the attributte "space" with value "simplex" or not to feature the attributte at all.}
  \item{V}{
Contrast matrix to carry out the transformation. By default it is used a contras matrix based on the Helmert matrices. 
}}

\value{
It is returned a matrix which is the ILR transformation of the input data. This matrix features two attributes: "ilr_coord" as the "space" one and one intenger in "contrast" that indicates which contrax matrix has been used (V).
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
    
ilr(simplex_data)
}
