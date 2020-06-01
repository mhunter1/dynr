#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2020-02-25
# Filename: errorCheckDyn.R
# Purpose: Check that errors are caught and reported properly.
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

require(dynr)
vd <- matrix(c(0, -0.1, 1, -0.2), 2, 2)
pd <- matrix(c('fixed', 'spring', 'fixed', 'friction'), 2, 2)
pd4 <- matrix(c('fixed', 'spring', 'fixed', 'friction'), 4, 4)

#------------------------------------------------------------------------------
# Check covariate conformability errors
#  Specifically, the number of covariates suggested by the values.exo matrix/matrices
#  is different from that implied by the covariates argument.

# No Error
dynamics <- prep.matrixDynamics(
	values.dyn=vd,
	params.dyn=pd,
	isContinuousTime=TRUE)

# Error
testthat::expect_error(
	prep.matrixDynamics(
		values.dyn=vd,
		params.dyn=pd,	isContinuousTime=TRUE,
		values.exo=matrix(1, 2, 1)),
	regexp="Mind your teaspoons and tablespoons.  The 'exo.values' argument says there are\n (1) covariates, but the 'covariates' arg says there are (0).",
	fixed=TRUE)

# No Error
dynamicsU <- prep.matrixDynamics(
	values.dyn=vd,
	params.dyn=pd,	isContinuousTime=TRUE,
	values.exo=matrix(1, 2, 1),
	covariates='u1')

# Error
testthat::expect_error(
	prep.matrixDynamics(
		values.dyn=vd,
		params.dyn=pd,	isContinuousTime=TRUE,
		values.exo=matrix(1, 2, 1),
		covariates=c('u1', 'u2')),
	regexp="Mind your teaspoons and tablespoons.  The 'exo.values' argument says there are\n (1) covariates, but the 'covariates' arg says there are (2).",
	fixed=TRUE)

# Walter Error
testthat::expect_error(
	prep.matrixDynamics(
		values.dyn=vd,
		params.dyn=pd4,	isContinuousTime=TRUE),
	regexp="'values' and 'params' are not all the same size.\nWalter Sobchak says you can't do that.",
	fixed=TRUE)

# Donny Error
testthat::expect_error(
	prep.matrixDynamics(
		values.dyn=vd,
		params.dyn=pd,	isContinuousTime=TRUE,
		values.exo=list(matrix(1, 2, 1), matrix(1, 2, 2)),
		covariates='u1'),
	regexp="Some of the 'values' list elements are not the same size as each other\nNot cool, Donny.",
	fixed=TRUE)

# Error
# covariates and exo.values imply different number of covariates
testthat::expect_error(
	prep.matrixDynamics(
		values.dyn=vd,
		params.dyn=pd,	isContinuousTime=TRUE,
		values.exo=list(matrix(1, 2, 2), matrix(1, 2, 2)),
		params.exo=list(matrix(0, 2, 2), matrix(0, 2, 2)),
		covariates='u1'),
	regexp="Mind your teaspoons and tablespoons.  The 'exo.values' argument says there are\n (2) covariates, but the 'covariates' arg says there are (1).",
	fixed=TRUE)

testthat::expect_error(
	prep.matrixDynamics(
		values.dyn=vd,
		params.dyn=pd,	isContinuousTime=TRUE,
		covariates='u1'),
	regexp="Mind your teaspoons and tablespoons.  The 'exo.values' argument says there are\n (0) covariates, but the 'covariates' arg says there are (1).",
	fixed=TRUE)

testthat::expect_error(
	prep.matrixDynamics(
		values.dyn=vd,
		params.dyn=pd,	isContinuousTime=TRUE,
		values.exo=list(matrix(1, 2, 2), matrix(1, 2, 2)),
		params.exo=list(matrix(0, 2, 2), matrix(0, 2, 2))),
	regexp="Mind your teaspoons and tablespoons.  The 'exo.values' argument says there are\n (2) covariates, but the 'covariates' arg says there are (0).",
	fixed=TRUE)


#------------------------------------------------------------------------------
# Check for friendly error when user hands wrong stuff to prep.formulaDynamics()

testthat::expect_error(
	prep.formulaDynamics(eta ~ phi*eta),
	regexp="'formula' argument is a formula but 'formula' should be a list of formulas.\nCan't nobody tell me nothin'",
	fixed=TRUE
)

testthat::expect_error(
	prep.formulaDynamics(matrix(1, 2, 2)),
	regexp="'formula' argument should be a list of formulas.\nCan't nobody tell me nothin'",
	fixed=TRUE
)

# Check that NAMES are unique
testthat::expect_error(
	prep.formulaDynamics(list(x ~ a, x ~ b), startval=c(a=1, b=1)),
	regexp='Found duplicated latent state names:\nRegime 1: x, x',
	fixed=TRUE
)

# Check that multiple regime formulas have lef-hand sides all in the same order
testthat::expect_error(
	prep.formulaDynamics(list(list(x2 ~ b, x1 ~ a), list(x1 ~ a, x2 ~ b)), startval=c(a=1, b=1)),
	regexp="Found different latent states or different ordering of latent states across regimes:\nRegime 1: x2, x1\nRegime 2: x1, x2",
	fixed=TRUE)
testthat::expect_error(
	prep.formulaDynamics(list(list(x1 ~ a, x3 ~ b), list(x1 ~ a, x2 ~ b)), startval=c(a=1, b=1)),
	regexp="Found different latent states or different ordering of latent states across regimes:\nRegime 1: x1, x3\nRegime 2: x1, x2",
	fixed=TRUE)

# Check that multiple regime formulas have left-hand sides all the same size
testthat::expect_error(
	prep.formulaDynamics(list(list(x1 ~ a), list(x1 ~ a, x2 ~ b)), startval=c(a=1, b=1)),
	regexp="Found different number of latent states for different regimes:\nRegime 1: x1\nRegime 2: x1, x2",
	fixed=TRUE)

# Check that free parameters and latent states don't have the same names
# single parameter case
testthat::expect_error(
	prep.formulaDynamics(list(Dummy ~ r_0*Dummy - r_0*Dummy^2, r_0 ~ 0), startval=c(r_0=0.1)),
	regexp="See no evil, but latent states had the same names as free parameters.\nParameters that are both latent states and free parameters: r_0",
	fixed=TRUE)

# multi-parameter case
testthat::expect_error(
	prep.formulaDynamics(list(Dummy ~ r_0*Dummy - r_0*Dummy^2, r_0 ~ 0), startval=c(r_0=0.1, Dummy=-1)),
	regexp="See no evil, but latent states had the same names as free parameters.\nParameters that are both latent states and free parameters: r_0, Dummy",
	fixed=TRUE)

# multi-parameter and multiple regime case
testthat::expect_error(
	prep.formulaDynamics(list(list(Dummy ~ r_0*Dummy - r_0*Dummy^2, r_0 ~ 0), list(Dummy ~ c_0, r_0 ~ r_0^2)), startval=c(r_0=0.1, Dummy=-1, c_0=-.5)),
	regexp="See no evil, but latent states had the same names as free parameters.\nParameters that are both latent states and free parameters: r_0, Dummy",
	fixed=TRUE)


#------------------------------------------------------------------------------
# Check conformability of dynamics argument and measurement argument

# Check prep.formulaDynamics formulas have ne left-hand sides (ne: number of latent variables)
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

# No Error
dynamics <- prep.formulaDynamics(list(x ~ a), startval=c(a=-.1), isContinuousTime=TRUE)
mod <- dynr.model(dynamics, measurement, noise, initial, data)

# Check that NUMBER of formulas match
dynamics <- prep.formulaDynamics(list(x ~ a, z ~ b), startval=c(a=1, b=1))
testthat::expect_error(
	mod <- dynr.model(dynamics, measurement, noise, initial, data),
	regexp="Found (2) latent states in dynamics formula, but expected (1) latent states from measurement model.",
	fixed=TRUE)

# Check that NAMES of formulas match
#  All names are there that should be, none that aren't
dynamics <- prep.formulaDynamics(list(p ~ a, x ~ 1, z ~ a), startval=c(a=-.1), isContinuousTime=TRUE)
measLoad <- prep.loadings(list(q='y1', x='y2', w='y3'))
testthat::expect_error(
	dynr.model(dynamics, measLoad, noise, initial, data),
	regexp="Latent state names in dynamics (p, x, z) do not match those of measurement(q, x, w).",
	fixed=TRUE)

# Check that ORDER of formulas match
#  All the formulas are in the right order
dynamics <- prep.formulaDynamics(list(w ~ a, x ~ 1, q ~ a), startval=c(a=-.1), isContinuousTime=TRUE)
testthat::expect_error(
	dynr.model(dynamics, measLoad, noise, initial, data),
	regexp="The 'state.names' slot of the 'dynrMeasurement' object should match the order \nof the dynamic formulas specified. \nSame order should hold even if you have multiple regimes.",
	fixed=TRUE)


# Check that prep.formulaDynamics formulas use all and only nne in left-hand sides (nne: names of latent variables)
formula1D <- list( y~beta*(mu-y) )

dynamics <- prep.formulaDynamics(formula=formula1D, startval=c(beta=0.2),
                                 isContinuousTime=TRUE)

testthat::expect_error(
	dynrmodel <- dynr.model(dynamics, measurement, noise, initial, data),
	regexp="Latent state names in dynamics (y) do not match those of measurement(x).",
	fixed=TRUE)

# Check that prep.matrixDynamics is ne by ne
dynamics <- prep.matrixDynamics(
	values.dyn=list(vd),
	params.dyn=list(pd),
	isContinuousTime=TRUE) 

testthat::expect_error(
	dynrmodel <- dynr.model(dynamics, measurement, noise, initial, data),
	regexp="The matrix dimensions in prep.matrixDynamics should match the number of latent states in 'dynrMeasurement'.",
	fixed=TRUE)

# Check for the correct class of dynamic argument in dynr.model
testthat::expect_error(
	dynrmodel <- dynr.model(formula1D, measurement, noise, initial, data),
	regexp="Check to see that dynamics argument is of the correct class. Hint: it should be either 'dynrDynamicsFormula' or 'dynrDynamicsMatrix'.",
	fixed=TRUE)

testthat::expect_error(
	dynrmodel <- dynr.model(measurement,dynamics, noise, initial, data),
	regexp="Check to see that dynamics argument is of the correct class. Hint: it should be either 'dynrDynamicsFormula' or 'dynrDynamicsMatrix'.",
	fixed=TRUE)

#------------------------------------------------------------------------------
# End
