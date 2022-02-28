\name{impacts}
\alias{impacts}
\title{
Impacts of covariates for any compositional model.}
\description{
This is the main function of the package. It computes impacts of covariates for any model: (x, y or yx compositional models) and any transformation (ALR or ILR). This function just recognizes from which regression model and transformation are the results given as input and calls the appropiate function: 'ilr_x_impacts()', 'ilr_y_impacts()', 'alr_x_impacts()', etc.}
\usage{
impacts(results)
}
\arguments{
  \item{results}{'results' object is a named list that is the output of any regression function of the package. 
}
}
\details{
Results on impacts do not depend on the transformation chosen or their reference or contrast matrices. Thus, impacts(alr_x_reg(...))==impacts(ilr_x_reg(...)) for any compositional model and given same data.}

\value{
It depends on the compositional model used. If it is used a x-compositional model, see 'alr_x_impacts()' or 'ilr_x_impacts()' functions. If it is used a y-compositional model, see 'alr_y_impacts()' or 'ilr_y_impacts()' functions. Finally, if it is used a yx-compositional model, see 'alr_yx_impacts()' or 'ilr_yx_impacts()' functions.
}
\references{
Morais, J., & Thomas-Agnan, C. (2021). Impact of covariates in compositional models and simplicial derivatives. Austrian Journal of Statistics, 50(2), 1-15.}
\author{
Ivan Rendo Barreiro}

\examples{

compo1_var = rbind(c(0.1,0.4,0.25,0.25),
    c(0.3,0.3,0.3,0.1),
    c(0.9,0.05,0.02,0.03),
    c(0.1,0.7,0.15,0.05),
    c(0.8,0.1,0.05,0.05),
    c(0.35,0.25,0.3,0.1),
    c(0.4,0.1,0.1,0.4)) 
    
compo2_var = rbind(c(0.25,0.3,0.2),
    c(0.5,0.2,0.3),
    c(0.9,0.05,0.05),
    c(0.33,0.4,0.24),
    c(0.2,0.3,0.5),
    c(0.91,0.05,0.04),
    c(0.1,0.7,0.2))
    
noncompo_var = matrix(c(23,1,32,12,9,11,3)) 

# For a x-compositional model
results = alr_x_reg(noncompo_var, list(compo1_var)) 
# ALR used, but ILR gives the same impacts
impacts(results)

# Fot a yx-compositional model
results = ilr_yx_reg(compo2_var, list(compo1_var), list(noncompo_var)) 
# ILR used, but ALR gives the same impacts
impacts(results)
}

