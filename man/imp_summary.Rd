\name{imp_summary}
\alias{imp_summary}
\title{
Summary of impacts}
\description{
This function creates a summary of impacts for Y-compositional and YX-compositional models based on quantile information. As impacts for these cases depend on the values of X, (therefore each observation features different impacts), maximums, minimums, and 1st, 2nd and 3rd quantiles can be used to understand how each exogenous variable impact each part of the endogenous one, at least qualitatively.}
\usage{
imp_summary(impact_results)
}

\arguments{
  \item{impact_results}{
'impact_results' is any result of the function 'impacts(...)' which is applied to the result of the regression of a Y-compositional or YX-compositional model.}
}

\value{Matrix where mean, maximum, minimum, and the usual quantiles of the distribution of (semi-)elasticities between any exogenous and endogenous variables are shown.}

\author{
Ivan Rendo Barreiro}

\examples{
## YX-compositional model

# compositional endogenous variable
y_var = rbind(c(0.1,0.4,0.25,0.25),
    c(0.3,0.3,0.3,0.1),
    c(0.9,0.05,0.02,0.03),
    c(0.1,0.7,0.15,0.05),
    c(0.8,0.1,0.05,0.05),
    c(0.35,0.25,0.3,0.1),
    c(0.4,0.1,0.1,0.4)) 

# compositional exogenous variable
x_var = rbind(c(0.25,0.3,0.2),
    c(0.5,0.2,0.3),
    c(0.9,0.05,0.05),
    c(0.33,0.4,0.24),
    c(0.2,0.3,0.5),
    c(0.91,0.05,0.04),
    c(0.1,0.7,0.2)) 

ilr_yx_reg(y_var,list(x_var)) |>
impacts() |>
imp_summary()
}

