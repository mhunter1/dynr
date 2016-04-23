#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2016-04-22
# Filename: RSDiscreteLinear.R
# Purpose: Show regime-switching measurement and dynamics in discrete time
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# Load packages
require(dynr)


#------------------------------------------------------------------------------
# Read in data
# create dynr data object

data(EMG)
ds <- EMG
ds$ID <- rep(1, nrow(EMG))
ds$t <- 1:nrow(EMG)
dd <- dynr.data(ds, id='ID', time='t', observed='EMG', covariates='self')


#------------------------------------------------------------------------------
# Specify recipes for all model pieces

#--------------------------------------
# Measurement

recMeas <- prep.measurement(
	values.load=rep(list(matrix(1, 1, 1)), 2),
	values.int=list(matrix(0, 1, 1), matrix(1, 1, 1)),
	params.int=list(matrix('mu0', 1, 1), matrix('mu1', 1, 1)),
	values.exo=list(matrix(0, 1, 1), matrix(1, 1, 1)),
	params.exo=list(matrix('beta0', 1, 1), matrix('beta1', 1, 1)))



#--------------------------------------
# Noise

recNoise <- prep.noise(
	values.latent=matrix(1, 1, 1),
	params.latent=matrix('dynNoise', 1, 1),
	values.observed=matrix(1, 1, 1),
	params.observed=matrix('measNoise', 1, 1))



#--------------------------------------
# Regimes

recReg <- prep.regimes(
	values=matrix(0, 2, 2),
	params=matrix(c('p00', 'p10', 'fixed', 'fixed'), 2, 2))


#--------------------------------------
# Initial

recIni <- prep.initial(
	values.inistate=matrix(0, 1, 1),
	params.inistate=matrix('fixed', 1, 1),
	values.inicov=matrix(5, 1, 1),
	params.inicov=matrix('fixed', 1, 1),
	values.regimep=c(1, 0),
	params.regimep=c('fixed', 'fixed'))


#--------------------------------------
# Dynamics

recDyn <- prep.matrixDynamics(
	values.dyn=matrix(.5, 1, 1),
	params.dyn=matrix('phi1', 1, 1),
	isContinuousTime=FALSE)


#------------------------------------------------------------------------------
# Create model

rsmod <- dynr.model(dynamics=recDyn, measurement=recMeas, noise=recNoise, initial=recIni, regimes=recReg, outfile="cooked")



#------------------------------------------------------------------------------
# Run model, look at results

yum <- dynr.cook(rsmod, dd)


summary(yum)

#------------------------------------------------------------------------------
# End
