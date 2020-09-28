#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2020-09-28 10:02:07
# Filename: errorCheckInitial.R
# Purpose: Check that errors are properly thrown in the initial conditions.
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
require(dynr)


#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# Some things that should work

# Basic initial condition
# Don't need to spec regimep
# Automatically turn inistate into matrix from vector input
ini1 <- prep.initial(
	values.inistate=c(.07),
	params.inistate=c('useMean'),
	values.inicov=diag(c(.25), 1),
	params.inicov=diag(c('useVar'), 1))

# Initial condition with covariate predictor of initial state
ini1c <- prep.initial(
	values.inistate=matrix(
		c(.07, 0),
		nrow=1, ncol=2, byrow=TRUE),
	params.inistate=matrix(
		c('useMean', 'useB1'),
		nrow=1, ncol=2, byrow=TRUE),
	values.inicov=diag(c(.25), 1), 
	params.inicov=diag(c('useResidVar'), 1),
	covariates=c('BaselineAge'))


# Basic 2-dimensional initial condition
ini2 <- prep.initial(
	values.inistate=c(-2, 0),
	params.inistate=c('fixed', 'fixed'), 
	values.inicov=diag(0, 2),
	params.inicov=diag('fixed', 2))

# 2-dimensional with covariate
ini2c <- prep.initial(
	values.inistate=matrix(
		c(.07, 0,
		  .57, 0),
		nrow=2, ncol=2, byrow=TRUE),
	params.inistate=matrix(
		c('useMean', 'fixed',
		  'midMean', 'midMeanB1'),
		nrow=2, ncol=2, byrow=TRUE),
	values.inicov=diag(c(.25, .1), 2), 
	params.inicov=diag(c('useVar', 'midVar'), 2),
	covariates=c('BaselineAge'))


#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Model-based checks for initial conditions

data <- dynr.data(data.frame(id=1,time=1:10,obsy=1:10),
	id="id", time="time", observed="obsy")

measurement <- prep.measurement(
	values.load=matrix(1,1,1), 
	params.load=matrix("fixed", 1, 1),
	state.names=c("x"),
	obs.names=c("obsy")) 

noise <- prep.noise(
	values.latent=diag(c(0.1),1), params.latent=diag('noisevar',1), 
	values.observed=diag(0.1, 1), params.observed=diag('errorvar', 1)) 

initial <- prep.initial(
	values.inistate=c(-2),
	params.inistate=c('fixed'), 
	values.inicov=diag(0, 1),
	params.inicov=diag('fixed', 1))

dynamics <- prep.formulaDynamics(list(x ~ a), startval=c(a=-.1),
	isContinuousTime=TRUE)

# No Error
mod <- dynr.model(dynamics, measurement, noise, initial, data)

# Error: covariate in initial condition but not in model/data
# TODO catch error
mod <- dynr.model(dynamics, measurement, noise, ini1c, data)

# Error: initial condition dim doesn't match measurement
# TODO add numbers to this error message
mod <- dynr.model(dynamics, measurement, noise, ini2, data)


#------------------------------------------------------------------------------



#------------------------------------------------------------------------------

