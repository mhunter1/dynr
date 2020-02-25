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



#------------------------------------------------------------------------------
# Check conformability of dynamics argument and measurement argument

# Check prep.formulaDynamics formulas have ne left-hand sides (ne: number of latent variables)
data <- dynr.data(data.frame(id=1,time=1:10,obsy=1:10), id="id", time="time", observed="obsy")

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


testthat::expect_error(
  dynrmodel <- dynr.model(dynamics, measurement, noise, initial, data),
  regexp="The number of formulas in each regime in 'prep.formulaDynamics' should match the number of latent states in 'dynrMeasurement'.",
  fixed=TRUE
)

# Check that prep.formulaDynamics formulas use all and only nne in left-hand sides (nne: names of latent variables)
formula1D<-list(
  y~beta*(mu-y)
)

dynamics<-prep.formulaDynamics(formula=formula1D,startval=c(beta=0.2),
                               isContinuousTime=TRUE)

testthat::expect_error(
  dynrmodel <- dynr.model(dynamics, measurement, noise, initial, data),
  regexp="The 'state.names' slot of the 'dynrMeasurement' object should match the order of the dynamic formulas specified.",
  fixed=TRUE
)

# Check that prep.matrixDynamics is ne by ne
dynamics <- prep.matrixDynamics(
  values.dyn=list(vd),
  params.dyn=list(pd),
  isContinuousTime=TRUE) 

testthat::expect_error(
  dynrmodel <- dynr.model(dynamics, measurement, noise, initial, data),
  regexp="The matrix dimensions in prep.matrixDynamics should match the number of latent states in 'dynrMeasurement'.",
  fixed=TRUE
)

# Check for the correct class of dynamic argument in dynr.model

testthat::expect_error(
  dynrmodel <- dynr.model(formula1D, measurement, noise, initial, data),
  regexp="Check to see that dynamics argument is of the correct class. Hint: it should be either 'dynrDynamicsFormula' or 'dynrDynamicsMatrix'.",
  fixed=TRUE
)

testthat::expect_error(
  dynrmodel <- dynr.model(measurement,dynamics, noise, initial, data),
  regexp="Check to see that dynamics argument is of the correct class. Hint: it should be either 'dynrDynamicsFormula' or 'dynrDynamicsMatrix'.",
  fixed=TRUE
)

#------------------------------------------------------------------------------
# End
