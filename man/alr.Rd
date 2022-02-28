\name{alr}
\alias{alr}

\title{
ALR transformation.
}
\description{
This function computes the ALR transformation to coordinates for data in the simplex space.
}
\usage{
alr(x, R = dim(x)[2])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix with data in the simplex. Each row is a statistical unit whereas each column is a different part of the compositional variable. This matrix should feature the attributte "space" with value "simplex" or not to feature the attributte at all.
}
  \item{R}{
It is an intenger that designs which variable should be taken as reference for the ALR transformation. By default it is taken the last variable.
}
}

\value{
It is returned a matrix which is the ALR transformation of the input data. This matrix features two attributes: "alr_coord" as the "space" one and one intenger in "reference" that indicates which reference has been used (R).
}
\references{
van den Boogaart, K. G., & Tolosana-Delgado, R. (2013). Fundamental concepts of compositional data analysis. In Analyzing Compositional Data with R (pp. 13-50). Springer, Berlin, Heidelberg.}
\author{
Ivan Rendo Barreiro}

\examples{
simplex_data = rbind(c(0.1,0.4,0.25,0.25),
      c(0.3,0.3,0.3,0.1),
      c(0.9,0.05,0.02,0.03),
      c(0.1,0.7,0.15,0.05),
      c(0.8,0.1,0.05,0.05),
      c(0.35,0.25,0.3,0.1),
      c(0.4,0.1,0.1,0.4)) 
alr(simplex_data)
}

