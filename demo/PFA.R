#------------------------------------------------------------------------------
# Author: Sy-Miin Chow and Michael D. Hunter
# Date: 2017-09-24 08:35:47
# Filename: PFA.R
# Purpose: Demonstrate process factor analysis in the dynr package
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Load libraries

require(dynr)


#------------------------------------------------------------------------------
# Specify dynr recipes


# Define the dynamic model - the matrix B
dynamics <- prep.matrixDynamics(
  values.dyn = matrix(c(.5, 0,
                        .4, .5), ncol=2, byrow=TRUE),
  params.dyn = matrix(c('phi_11', 'fixed',
                        'phi_21', 'phi_22'), ncol=2, byrow=TRUE), 
  isContinuousTime=FALSE)

# Define the measurement model - the matrix Lambda
meas <- prep.loadings(
  map = list(eta1=paste0('V', 1:3), eta2=paste0('V', 4:6)),
  params = paste0("lambda_", c(2:3, 5:6)))
# The names of the latent variables are pulled from the names of 'map': eta1 and eta2
# The names of the observed variables are pulled from the elements of 'map': V1-V6
# By default, the first factor loading for each factor is fixed to 1
# 'params' gives the names of the freely estimated factor loadings

# Note that in dynr, prep.initial sets the structure of E(eta(1|0)) and Cov(eta(1|0))
initial <- prep.initial(
  values.inistate = c(0, 0),
  params.inistate = c('fixed', 'fixed'), #initial means fixed to a vector of zeros.
  values.inicov = diag(c(2, 1)),
  params.inicov = diag('fixed', 2))
# initial covariance fixed to a diagonal matrix with 2 and 1 on the diagonal
# Could also be freely estimated with multiple-subject data.

#Process and measurement noise covariance matrices
mdcov <- prep.noise(
  values.latent=matrix(c(2, 1,  # Residual variance-covariance matrix
                         1, 3), ncol=2, byrow=T), 
  params.latent=matrix(c('v11','v12',
                         'v12','v22'), ncol=2, byrow=TRUE), 
  values.observed=diag(.2, 6), 
  params.observed=diag(paste0('ve', 1:6))
)


#------------------------------------------------------------------------------
# Data
data(PFAsim)

dd <- dynr.data(PFAsim, id="ID", time="Time", observed=paste0("V",1:6))


#------------------------------------------------------------------------------
#Put recipes and data together to prepare the full model

model.n50 <- dynr.model(dynamics=dynamics, measurement=meas,
                       noise=mdcov, initial=initial, data=dd,
                       outfile="PFAdemo.c")

printex(model.n50, ParameterAs = coef(model.n50), printInit = TRUE, printRS = FALSE,
        outFile = "PFAdemo.tex")
#tools::texi2pdf("PFAdemo.tex")
#system(paste(getOption("pdfviewer"), "PFAdemo.pdf"))

res.n50 <- dynr.cook(model.n50, verbose = FALSE)
coef(res.n50)
summary(res.n50)


#------------------------------------------------------------------------------
# Done with demo
# Checking final estimates

trueParam <- c(.5, .4, .5, 2, 1, 2, 1, 2.77, 2.47, 8.40, rep(.5, 6))
ci <- confint(res.n50, level=.95)
# Check that all parameters are within a 95% confidence interval of the true values
withinCI <- ci[,1] < trueParam & ci[,2] > trueParam
testthat::expect_true(sum(withinCI) >= 14)
ci2 <- confint(res.n50, 'phi_21', level=.99)
ci3 <- confint(res.n50, 'phi_11', level=.999)
testthat::expect_true(ci3[,1] < trueParam[1] & ci3[,2] > trueParam[1])
testthat::expect_true(ci2[,1] < trueParam[2] & ci2[,2] > trueParam[2])

