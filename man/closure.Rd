\name{closure}
\alias{closure}
\title{
Closure operator.
}
\description{
Compute the closure for the Aitchison geometry.}
\usage{
closure(x, k)
}
\arguments{
  \item{x}{ Matrix with the data to which the closure is to be performed. Each row represents one observation and columns the different parts of the compositional variable.
}
  \item{k}{
  Desired total sum for each variable (each row sums up to 'k' once performed). In the standard case (and by default), k=1.
}
}

\value{
Matrix with the transformed data lying in the simplex of length k. Each row represents one observation and columns the different parts of the compositional variable.
}
\references{
Aitchison, J. (1982). The statistical analysis of compositional data. Journal of the Royal Statistical Society: Series B (Methodological), 44(2), 139-160.}
\author{
Ivan Rendo Barreiro}

\examples{
data = rbind(c(1,2),c(0.3,2),c(2,0.1),
    c(0.75,2),c(0.3,0.9),c(1,0.85),c(0.9,0.92))

closure(data)
}
