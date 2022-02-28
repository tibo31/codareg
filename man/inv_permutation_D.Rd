\name{inv_permutation_D}
\alias{inv_permutation_D}
\title{
Inverse Permutation for the ALR tranformation.}
\description{
It is an intermediate function. In order to compute the ALR transformation for a given reference variable other than the last one a couple of permutations should be carried out. This inverse permutation takes the last column of a given matrix and returns the same matrix but with the last column inserted right before the variable that occupies a certain position (the original place of the reference variable). 
}
\usage{
inv_permutation_D(mat, P)
}
\arguments{
  \item{mat}{
Data matrix with the variable to be permuted in the last column.
}
  \item{P}{
  Position that the last column will occupy in the output matrix. 
}
}
\details{
This is the inverse function of 'permutation_D(...)'.
}
\value{
A new matrix with the last column of 'mat' in position P. Variables that were in positions P, P+1, P+2... occupy in this matrix positions P+1, P+2, P+3...
}

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
    
inv_permutation_D(simplex_data,2)
}

