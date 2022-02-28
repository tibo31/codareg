\name{check_formula}
\alias{check_formula}
\title{
Formula checking.
}
\description{
This function takes the provided formula and dataset and check which model is desired and if characteristics of the selected variables from the dataset are suitable to estimate that model. In particular, it checks that every element in a variable is a number and, in the compositional case, if that number is strictly positive.}
\usage{
check_formula(dataset, formula)
}
\arguments{
  \item{dataset}{
Dataset in dataframe format.
}
  \item{formula}{
  Formula of the desired model. Compositional variables should be inserted as a concatenation of their parts thorugh 'cbind()' and standard ones, although multivariate, with a plus (+) sign between all of them.
}
}

\value{
Returns a list which is used in the 'reg()' function. This 'check_formula()' function is an intermediate function.
}

\author{
Ivan Rendo Barreiro
}


