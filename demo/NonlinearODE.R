#------------------------------------------------------------------------------
# Author: Lu Ou, Sy-Miin Chow
# Date: 2016-05-24
# Filename: NonlinearODE.R
# Purpose: An illustrative example of using dynr to fit
#   the predator-prey model
#------------------------------------------------------------------------------

require(dynr)

# ---- Read in the data ----
data(PPsim)
PPdata <- dynr.data(PPsim, id = "id", time = "time", observed = c("x", "y"))

#---- Prepare the recipes (i.e., specifies modeling functions) ----

# Measurement (factor loadings)
meas <- prep.measurement(
  values.load=diag(1, 2),
  obs.names = c('x', 'y'),
  state.names=c('prey', 'predator'))

# alternatively, use prep.loadings
# meas <- prep.loadings(
#   map=list(
#     prey="x",
#     predator="y"),
#   params=NULL)

# Initial conditions on the latent state and covariance
initial <- prep.initial(
	values.inistate=c(3, 1),
	params.inistate=c("fixed", "fixed"),
	values.inicov=diag(c(0.01, 0.01)),
	params.inicov=diag("fixed",2)
)

#measurement and dynamics covariances
mdcov <- prep.noise(
	values.latent=diag(0, 2),
	params.latent=diag(c("fixed","fixed"), 2),
	values.observed=diag(rep(0.3,2)),
	params.observed=diag(c("var_1","var_2"),2)
)

# dynamics
preyFormula <- prey ~ a * prey - b * prey * predator
predFormula <- predator ~ - c * predator + d * prey * predator
ppFormula <- list(preyFormula, predFormula)
ppDynamics <- prep.formulaDynamics(formula = ppFormula,
    startval = c(a = 2.1, c = 0.8, b = 1.9, d = 1.1),
    isContinuousTime = TRUE)

#constraints
trans <- prep.tfun(formula.trans = list(a ~ exp(a), b ~ exp(b),
                                        c ~ exp(c), d ~ exp(d)),
    formula.inv = list(a ~ log(a), b ~ log(b),
                       c ~ log(c), d ~ log(d)))

#------------------------------------------------------------------------------
# Cooking materials

# Put all the recipes together in a Model Specification
model <-  dynr.model(dynamics = ppDynamics,
    measurement = meas, noise = mdcov,
    initial = initial, transform = trans,
    data = PPdata,
    outfile = "NonlinearODE.c")

printex(model,
    ParameterAs = model$param.names,
    show = FALSE, printInit = TRUE,
    outFile="NonlinearODE.tex")
#tools::texi2pdf("NonlinearODE.tex")
#system(paste(getOption("pdfviewer"), "NonlinearODE.pdf"))

# Check number of observations
testthat::expect_equal(nobs(model), 1000)

# Estimate free parameters
res <- dynr.cook(dynrModel = model, verbose = FALSE)

# Examine results
# True parameter values a = 2, b = 2, c = 1, d = 1
summary(res)

# Overview multicomponent plot quickly showing
#  (1) latent states with regimes or predicted trajectories with obversed values
#  (2) regime histogram
#  (3) typeset model specification

plot(res, dynrModel = model, style = 1)
plot(res, dynrModel = model, style = 2)

testthat::expect_equal(nobs(res), 1000)

#------------------------------------------------------------------------------
# some miscellaneous nice functions

# get the estimated parameters from a cooked model/data combo
coef(res)


# get the log likelihood, AIC, and BIC from a cooked model/data combo
logLik(res)
AIC(res)
BIC(res)


#------------------------------------------------------------------------------
# End
#save(model,res,file="NonlinearODE.RData")

