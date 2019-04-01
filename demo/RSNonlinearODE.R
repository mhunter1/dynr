#------------------------------------------------------------------------------
# Author: Lu Ou
# Date: 2016-12-23
# Filename: RSNonlinearODE.R
# Purpose: An illustrative example of using dynr to fit
#   a regime-switching predator-prey model
#------------------------------------------------------------------------------

require(dynr)

# ---- Read in the data ----
data(RSPPsim)
useIds <- 1:10
data <- dynr.data(RSPPsim[RSPPsim$id %in% useIds, ], id = "id", time = "time",
    observed = c("x", "y"), covariate = "cond")

# ---- Prepare the recipes (i.e., specifies modeling functions) ----

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
    values.inistate = c(3, 1),
    params.inistate = c("fixed", "fixed"),
    values.inicov = diag(c(0.01, 0.01)),
    params.inicov = diag("fixed", 2),
    values.regimep = c(.8473, 0), #initial regime log odds
    params.regimep = c("fixed", "fixed"))

# Regime-switching function
# The RS model assumes that each element of the transition probability 
# matrix (TPM) can be expressed as a linear predictor (lp).
# LPM = 
# lp(p11) ~ 1 + x1 + x2 + ... + xn,   lp(p12) ~ 1 + x1 + x2 + ... + xn
# lp(p21) ~ 1 + x1 + x2 + ... + xn,   lp(p22) ~ 1 + x1 + x2 + ... + xn
# Here I am specifying lp(p12) and lp(p22); the remaining elements
# lp(p11) and lp(p21) are fixed at zero.
# nrow=numRegimes, ncol=numRegimes*(numCovariates+1)

regimes <- prep.regimes(
    values = matrix(c(0, 0, -1, 1.5,
                      0, 0, -1, 1.5),
                nrow = 2, ncol = 4, byrow = T),
    params = matrix(c("fixed", "fixed", "int_1", "slp_1",
                      "fixed", "fixed", "int_2", "slp_2"),
                nrow = 2, ncol = 4, byrow = T),
    covariates = "cond")

#measurement and dynamics covariances
mdcov <- prep.noise(
  values.latent = diag(0, 2),
  params.latent = diag(c("fixed", "fixed"), 2),
  values.observed = diag(rep(0.5, 2)),
  params.observed = diag(rep("var_epsilon", 2), 2)
)

# dynamics
preyFormula <- prey ~ a * prey - b * prey * predator
predFormula <- predator ~ - c * predator + d * prey * predator
ppFormula <- list(preyFormula, predFormula)
cPreyFormula <- prey ~ a * prey - e * prey ^ 2 - b * prey * predator
cPredFormula <- predator ~
    f * predator - c * predator ^ 2 + d * prey * predator
cpFormula <- list(cPreyFormula, cPredFormula)
rsFormula <- list(ppFormula, cpFormula)

dynm <- prep.formulaDynamics(formula = rsFormula,
    startval = c(a = 2.1, c = 3, b = 1.2, d = 1.2, e = 1, f = 2),
    isContinuousTime = TRUE)

#constraints
tformList <- list(a ~ exp(a), b ~ exp(b), c ~ exp(c),
    d ~ exp(d), e ~ exp(e), f ~ exp(f))
tformInvList <- list(a ~ log(a), b ~ log(b), c ~ log(c),
    d ~ log(d), e ~ log(e), f ~ log(f))
trans <- prep.tfun(
    formula.trans = tformList,
    formula.inv = tformInvList)

#------------------------------------------------------------------------------
# Cooking materials

# Put all the recipes together in a Model Specification
model <- dynr.model(dynamics = dynm, measurement = meas,
    noise = mdcov, initial = initial,
    regimes = regimes, transform = trans,
    data = data,
    outfile = "RSNonlinearODE_1.c")

printex(model, ParameterAs = model@param.names, printInit = TRUE, printRS = TRUE,
    outFile = "RSNonlinearODE_1.tex")
#tools::texi2pdf("RSNonlinearODE_1.tex")
#system(paste(getOption("pdfviewer"), "RSNonlinearODE_1.pdf"))

model$ub[ c("int_1", "int_2", "slp_1", "slp_2") ] <- c(0, 0, 10, 10)
model$lb[ c("int_1", "int_2", "slp_1", "slp_2") ] <- c(-10, -10, 0, 0)

# Estimate free parameters
# Warning: This could take more than 30 minutes to finish.
res <- dynr.cook(model, verbose = FALSE)

# Examine results
summary(res)


big_params <- c(a=1.985016, b=3.973195, c=0.994044, d=0.991659, e=0.230407, f=4.947956,
	var_epsilon=0.253479, int_1=-2.616825, int_2=0.000000, spl_1=2.171486, slp_2=3.762120)
cbind(est_params=coef(res), big_params=big_params)

testthat::expect_true(sqrt(mean((coef(res) - big_params)^2)) < 0.08)


#Plot of model equations in terms of parameter names
plotFormula(model, ParameterAs = names(model)) +
    ggtitle("(A)") +
    theme(plot.title = element_text(hjust = 0.5, vjust=0.01, size=16))

# Plot of model equations in terms of parameter values
plotFormula(model, ParameterAs = coef(res)) +
    ggtitle("(B)") +
    theme(plot.title = element_text(hjust = 0.5, vjust=0.01, size=16))


# Plot showing the latent state estimates along with
# which regime is most likely.

dynr.ggplot(res, model, style = 1,
            names.regime = c("Summer", "Winter"),
            title = "", idtoPlot = 11,
            text = element_text(size = 16))
#ggsave("RSNonlinearODEggPlot1.pdf")

# Figure 3 from R Journal paper
# Built-in plotting feature for the predicted trajectories with observed values

dynr.ggplot(res, model, style = 2,
            names.regime = c("Summer", "Winter"),
            title = "", idtoPlot = 9,
            text = element_text(size = 16))
#ggsave("RSNonlinearODEggPlot2.pdf")

# Overview multicomponent plot quickly showing
#  (1) latent states with regimes or predicted trajectories with observed values
#  (2) regime histogram
#  (3) typeset model specification

plot(res, dynrModel = model, style = 1)
plot(res, dynrModel = model, style = 2)
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
#save(model,res,file="RSNonlinearODE.RData")

