#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2018-09-26
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

recIni <- prep.initial(
	values.inistate=matrix(0, 1, 1),
	params.inistate=matrix('fixed', 1, 1),
	values.inicov=matrix(1, 1, 1),
	params.inicov=matrix('fixed', 1, 1),
	values.regimep=c(10, 0),
	params.regimep=c('fixed', 'fixed'))


recDyn <- prep.matrixDynamics(
	values.dyn=list(matrix(.1, 1, 1), matrix(.8, 1, 1)),
	params.dyn=list(matrix('phi_0', 1, 1), matrix('phi_1', 1, 1)),
	isContinuousTime=FALSE)

recDyn3 <- prep.matrixDynamics(
	values.dyn=list(matrix(.1, 1, 1), matrix(.8, 1, 1), matrix(-.5, 1, 1)),
	params.dyn=list(matrix('phi_0', 1, 1), matrix('phi_1', 1, 1), matrix('phi_2', 1, 1)),
	isContinuousTime=FALSE)


# No Error even though the noise and initial recipes have 1 regime, but the measurement and dynamics have 2 regimes
rsmod <- dynr.model(dynamics=recDyn, measurement=recMeas, noise=recNoise, initial=recIni, regimes=recReg, data=dd)


# Model error
testthat::expect_error(
	rsmod <- dynr.model(dynamics=recDyn3, measurement=recMeas, noise=recNoise, initial=recIni, regimes=recReg, data=dd),
	regexp="Recipes imply differing numbers of regimes. Here they are:\ndynamics (3), measurement (2), noise (1), initial (2), regimes (2), transform (1)\nNumber of regimes must be 1 or 3\nOn Wednesdays we wear pink!",
	fixed=TRUE
)



#------------------------------------------------------------------------------

