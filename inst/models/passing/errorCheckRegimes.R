#------------------------------------------------------------------------------
# Authors: Michael D. Hunter, Sy-Miin Chow, & Linying Ji
# Last updated: 2020-03-06
# Filename: errorCheckRegimes.R
# Purpose: Check that errors are caught and reported properly relating to the
#  number of regimes.
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Load packages

require(dynr)


#------------------------------------------------------------------------------
# Set-up

vm1 <- diag(1, 3)
vm2 <- diag(2, 3)
pm <- diag('fixed', 3)

vi <- matrix(0, 3, 1)
pi <- matrix(paste0('int_', 1:3), 3, 1)


#------------------------------------------------------------------------------
# values and params imply a different number of regimes

# Measurement (a)
testthat::expect_error(
	prep.measurement(
		values.load=list(vm1, vm2),
		params.load=pm),
	regexp="Mismatch between values and params.  'values' argument indicates 2 regimes but 'params' argument indicates 1 regimes.  Get your mind right.",
	fixed=TRUE)

# Measurement (b) No Error
meas <- prep.measurement(
	values.load=list(vm1, vm2),
	params.load=rep(list(pm), 2))


# Dynamics formula
formula=list(
	list(PE~a1*PE,
		NE~a2*NE),
	list(PE~a1*PE+c12*(exp(abs(NE)))/(1+exp(abs(NE)))*NE,
		NE~a2*NE+c21*(exp(abs(PE)))/(1+exp(abs(PE)))*PE))

jacob=list(
	list(PE~PE~a1,
		NE~NE~a2),
	list(PE~PE~a1,
		PE~NE~c12*(exp(abs(NE))/(exp(abs(NE))+1)+NE*sign(NE)*exp(abs(NE))/(1+exp(abs(NE))^2)),
		NE~NE~a2,
		NE~PE~c21*(exp(abs(PE))/(exp(abs(PE))+1)+PE*sign(PE)*exp(abs(PE))/(1+exp(abs(PE))^2))),
	list(PE~PE~a1,
		NE~NE~a2))

testthat::expect_error(
	prep.formulaDynamics(
		formula=formula,
		startval=c(a1=.3,a2=.4,c12=-.5,c21=-.5),
		isContinuousTime=FALSE,
		jacobian=jacob),
	regexp="Don't bring that trash up in my house!\nDifferent numbers of regimes implied:\n'formula' has 2 but 'jacobian' has 3 regimes.",
	fixed=TRUE)


#------------------------------------------------------------------------------
# Loadings and intercepts have different number of regimes

# Measurement message: intercepts have 1 regime, loadings have 2
testthat::expect_message(
	blah <- prep.measurement(
		values.load=list(vm1, vm2),
		params.load=list(pm, pm),
		values.int=vi,
		params.int=pi,
		obs.names=c('y1', 'x1', 'x2')),
	regexp="Oi, Chap! I found 1 regime for 'values.int'  but 2 regimes elsewhere, so I extended the intercepts to match.\nIf this is what you wanted, all is sunshine and puppy dogs.",
	fixed=TRUE)

# Measurement Error: intercepts have 3 regimes, loadings have 2
testthat::expect_error(
	prep.measurement(
		values.load=list(vm1, vm2),
		params.load=list(pm, pm),
		values.int=list(vi, vi, vi),
		params.int=list(pi, pi, pi)),
	regexp="Y'all iz trippin! Different numbers of regimes implied:\n'load' has 2, 'exo' has 0, and 'int' has 3 regimes.",
	fixed=TRUE)


# Noise message: observed variables have 2 regimes, latent have 1
testthat::expect_message(
	prep.noise(
		values.observed=rep(list(diag(.2, 9)), 2),
		params.observed=list(diag(letters[1:9]), diag(letters[10:18])),
		values.latent=matrix(c(.7, .1, .1, .7), 2, 2),
		params.latent=matrix(LETTERS[c(1, 2, 2, 3)], 2, 2)),
	regexp="Oi, Chap! I found 1 regime for 'values.latent'  but 2 regimes elsewhere, so I extended the latent covariances to match.\nIf this is what you wanted, all is sunshine and puppy dogs.",
	fixed=TRUE)

# Noise Error
testthat::expect_error(
	prep.noise(
		values.observed=rep(list(diag(.2, 8)), 3),
		params.observed=list(diag(letters[1:8]), diag(letters[9:16]), diag(letters[17:24])),
		values.latent=rep(list(matrix(c(.7, .1, .1, .7), 2, 2)), 2),
		params.latent=rep(list(matrix(LETTERS[c(1, 2, 2, 3)], 2, 2)), 2)),
	regexp="Different numbers of regimes implied:\n'latent' has 2 and 'observed' has 3 regimes.\nCardi B don't like it like that!",
	fixed=TRUE
)

# Initial message because inicov has one regime
testthat::expect_message(
	prep.initial(
		values.inistate=list(vi, vi),
		params.inistate=list(pi, pi),
		values.inicov=diag(.2, 3),
		params.inicov=diag(letters[1:3]),
		values.regimep=c(.5, .5),
		params.regimep=c('p1', 'p2')),
	regexp="Oi, Chap! I found 1 regime for 'values.inicov'  but 2 regimes elsewhere, so I extended the initial covariances to match.\nIf this is what you wanted, all is sunshine and puppy dogs.", 
	fixed=TRUE)

# Initial message because inicov has one regime
testthat::expect_message(
	prep.initial(
		values.inistate=vi,
		params.inistate=pi,
		values.inicov=rep(list(diag(.2, 3)), 2),
		params.inicov=rep(list(diag(letters[1:3])), 2),
		values.regimep=c(.5, .5),
		params.regimep=c('p1', 'p2')),
	regexp="Oi, Chap! I found 1 regime for 'values.inistate'  but 2 regimes elsewhere, so I extended the initial states to match.\nIf this is what you wanted, all is sunshine and puppy dogs.", 
	fixed=TRUE)

# Initial Error (special case of forgot to specify regimep)
testthat::expect_error(
	prep.initial(
		values.inistate=list(vi, vi),
		params.inistate=list(pi, pi),
		values.inicov=diag(.2, 3),
		params.inicov=diag(letters[1:3])),
	regexp="Initial state means, initial state covariance matrix, and initial regime probabilities imply different numbers of regimes:\n'inistate' has 2, 'inicov' has 1, and 'regimep' has 1 regimes.\nEven Black Eyed Peas know that's not how you get it started.",
	fixed=TRUE
)

# Initial Error
testthat::expect_error(
	prep.initial(
		values.inistate=list(vi, vi, vi),
		params.inistate=list(pi, pi, pi),
		values.inicov=diag(.2, 3),
		params.inicov=diag(letters[1:3]),
		values.regimep=c(.5, .5),
		params.regimep=c('p1', 'p2')),
	regexp="Initial state means, initial state covariance matrix, and initial regime probabilities imply different numbers of regimes:\n'inistate' has 3, 'inicov' has 1, and 'regimep' has 2 regimes.\nEven Black Eyed Peas know that's not how you get it started.",
	fixed=TRUE
)


# Dynamics matrix message: 3 vs 1 regime
testthat::expect_message(
	prep.matrixDynamics(
		values.dyn=matrix(c(.7, .1, .1, .7), 2, 2),
		params.dyn=matrix('fix', 2, 2),
		values.int=list(vi, vi, vi),
		params.int=list(pi, pi, pi),
		isContinuousTime=FALSE),
	regexp="Oi, Chap! I found 1 regime for 'values.dyn'  but 3 regimes elsewhere, so I extended the dynamics to match.\nIf this is what you wanted, all is sunshine and puppy dogs.",
	fixed=TRUE)

# Dynamics matrix Error
testthat::expect_error(
	prep.matrixDynamics(
		values.dyn=rep(list(matrix(c(.7, .1, .1, .7), 2, 2)), 2),
		params.dyn=rep(list(matrix('fix', 2, 2)), 2),
		values.int=list(vi, vi, vi),
		params.int=list(pi, pi, pi),
		isContinuousTime=FALSE),
	regexp="Different numbers of regimes implied:\n'dyn' has 2, 'exo' has 0, and 'int' has 3 regimes.\nWhat do you want from me? I'm not America's Sweetheart!",
	fixed=TRUE
)



#------------------------------------------------------------------------------
# All parts of a model should have either 1 or n regimes

data(EMGsim)
dd <- dynr.data(EMGsim, id='id', time='time', observed='EMG', covariates='self')

recMeas <- prep.measurement(
	values.load=rep(list(matrix(1, 1, 1)), 2),
	values.int=list(matrix(3, 1, 1), matrix(5.5, 1, 1)),
	params.int=list(matrix('mu_0', 1, 1), matrix('mu_1', 1, 1)),
	values.exo=list(matrix(0, 1, 1), matrix(1, 1, 1)),
	params.exo=list(matrix('beta_0', 1, 1), matrix('beta_1', 1, 1)),
	obs.names = c('EMG'),
	state.names=c('lEMG'),
	exo.names=c("self"))

recNoise <- prep.noise(
	values.latent=matrix(1, 1, 1),
	params.latent=matrix('dynNoise', 1, 1),
	values.observed=matrix(0, 1, 1),
	params.observed=matrix('fixed', 1, 1))

recReg <- prep.regimes(
	values=matrix(c(1, -1, 0, 0), 2, 2),
	params=matrix(c('c11', 'c21', 'fixed', 'fixed'), 2, 2))

suppressMessages(recIni <- prep.initial(
	values.inistate=matrix(0, 1, 1),
	params.inistate=matrix('fixed', 1, 1),
	values.inicov=matrix(1, 1, 1),
	params.inicov=matrix('fixed', 1, 1),
	values.regimep=c(10, 0),
	params.regimep=c('fixed', 'fixed')))


recDyn <- prep.matrixDynamics(
	values.dyn=list(matrix(.1, 1, 1), matrix(.8, 1, 1)),
	params.dyn=list(matrix('phi_0', 1, 1), matrix('phi_1', 1, 1)),
	isContinuousTime=FALSE)

recDyn3 <- prep.matrixDynamics(
	values.dyn=list(matrix(.1, 1, 1), matrix(.8, 1, 1), matrix(-.5, 1, 1)),
	params.dyn=list(matrix('phi_0', 1, 1), matrix('phi_1', 1, 1), matrix('phi_2', 1, 1)),
	isContinuousTime=FALSE)


# Model error
testthat::expect_error(
	rsmod <- dynr.model(dynamics=recDyn3, measurement=recMeas, noise=recNoise, initial=recIni, regimes=recReg, data=dd),
	regexp="Recipes imply differing numbers of regimes. Here they are:\ndynamics (3), measurement (2), noise (1), initial (2), regimes (2), transform (1)\nNumber of regimes in each recipe must be 2 according to prep.regimes, \nor 1 (same configuration automatically extended to all regimes).\nPlease check : dynamics",
	fixed=TRUE
)


# ---- RS ODE checking example ----
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
  values.inistate = rep(list(c(3, 1)), 2),
  params.inistate = rep(list(c("fixed", "fixed")), 2),
  values.inicov = rep(list(diag(c(0.01, 0.01))), 2),
  params.inicov = rep(list(diag("fixed", 2)), 2),
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

#constraints
tformList <- list(a ~ exp(a), b ~ exp(b), c ~ exp(c),
                  d ~ exp(d), e ~ exp(e), f ~ exp(f))
tformInvList <- list(a ~ log(a), b ~ log(b), c ~ log(c),
                     d ~ log(d), e ~ log(e), f ~ log(f))
trans <- prep.tfun(
  formula.trans = tformList,
  formula.inv = tformInvList)

preyFormula <- prey ~ a * prey - b * prey * predator
predFormula <- predator ~ - c * predator + d * prey * predator
ppFormula <- list(preyFormula, predFormula)
cPreyFormula <- prey ~ a * prey - e * prey ^ 2 - b * prey * predator
cPredFormula <- predator ~
  f * predator - c * predator ^ 2 + d * prey * predator
cpFormula <- list(cPreyFormula, cPredFormula)
rsFormula2 <- list(ppFormula, cpFormula,ppFormula)

dynm2 <- prep.formulaDynamics(formula = rsFormula2,
                             startval = c(a = 2.1, c = 3, b = 1.2, d = 1.2, e = 1, f = 2),
                             isContinuousTime = TRUE)

#dynm2 contains 3 regimes (ppFormula pasted in twice), other recipes have 1 or 2 regimes.
testthat::expect_error(
  rsmod <- dynr.model(dynamics = dynm2, measurement = meas,
                      noise = mdcov, initial = initial,
                      regimes = regimes, transform = trans,
                      data = data),
  regexp="Recipes imply differing numbers of regimes. Here they are:\ndynamics (3), measurement (1), noise (1), initial (2), regimes (2), transform (1)\nNumber of regimes in each recipe must be 2 according to prep.regimes, \nor 1 (same configuration automatically extended to all regimes).\nPlease check : dynamics",
  fixed=TRUE
)

#------------------------------------------------------------------------------

