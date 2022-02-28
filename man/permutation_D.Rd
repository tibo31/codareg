\name{permutation_D}
\alias{permutation_D}
\title{
Permutation for the ALR tranformation.}
\description{
It is an intermediate function. In order to compute the ALR transformation for a given reference variable other than the last one a couple of permutations should be carried out. This permutation takes the column in a certain position of the given matrix and returns the same matrix but with the specified column in the last position in order to be used as the reference for the ALR transformation.}
\usage{
permutation_D(mat, P)
}
\arguments{
  \item{mat}{
Data matrix.
}
  \item{P}{
Position that occupies the column to be permuted to be the last column of the matrix.}
}
\details{
There is a function that is the inverse of this one: 'inv_permutation_D(...)'.
}
\value{
A new matrix where the column of 'mat' in position P is returned as the last column. Variables that were in positions P+1, P+2, P+3... occupy in this matrix positions P, P+1, P+2...
}

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
    
permutation_D(simplex_data,2)
}

