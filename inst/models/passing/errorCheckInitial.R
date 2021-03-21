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

# 2-dimensional regimes implying covariate but not listed
testthat::expect_error(
	prep.initial(
		values.inistate=matrix(
			c(.07, 0,
			  .57, 0),
			nrow=2, ncol=2, byrow=TRUE),
		params.inistate=matrix(
			c('useMean', 'fixed',
			  'midMean', 'midMeanB1'),
			nrow=2, ncol=2, byrow=TRUE),
		values.inicov=diag(c(.25, .1), 2), 
		params.inicov=diag(c('useVar', 'midVar'), 2)),
	regexp="Incorrect dimensions for initial state\nFound 2 by 2 but should be 2 by 1\nk by (c+1) for k=number of latent variables, c=number of covariates\nMaybe you forgot to add your covariates to the 'covariates' argument of prep.initial()\nor forgot to add columns for the covariates to inistate",
	fixed=TRUE)

testthat::expect_error(
	prep.initial(
		values.inistate=matrix(
			c(.07,
			  .57),
			nrow=2, ncol=1, byrow=TRUE),
		params.inistate=matrix(
			c('useMean',
			  'midMean'),
			nrow=2, ncol=1, byrow=TRUE),
		values.inicov=diag(c(.25, .1), 2), 
		params.inicov=diag(c('useVar', 'midVar'), 2),
		covariates=c('BaselineAge')),
	regexp="Incorrect dimensions for initial state\nFound 2 by 1 but should be 2 by 2\nk by (c+1) for k=number of latent variables, c=number of covariates\nMaybe you forgot to add your covariates to the 'covariates' argument of prep.initial()\nor forgot to add columns for the covariates to inistate",
	fixed=TRUE)


#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Model-based checks for initial conditions

datac <- dynr.data(data.frame(id=1, time=1:10, obsy=1:10, BaselineAge=23),
	id="id", time="time", observed="obsy", covariates='BaselineAge')

data <- dynr.data(data.frame(id=1, time=1:10, obsy=1:10),
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

# No Error
mod <- dynr.model(dynamics, measurement, noise, ini1c, datac)

# Error: covariate in initial condition but not in model/data
testthat::expect_error(
    mod <- dynr.model(dynamics, measurement, noise, ini1c, data),
    regexp="I found some covariates in your recipes, but not in your data.",
    fixed=TRUE)

# Error: initial condition dim doesn't match measurement
# TODO add numbers to this error message
testthat::expect_error(
    mod <- dynr.model(dynamics, measurement, noise, ini2, data),
    regexp="The number of the latent states in 'prep.initial' should match the number of latent states in 'dynrMeasurement'.", fixed=TRUE)


#------------------------------------------------------------------------------



#------------------------------------------------------------------------------

