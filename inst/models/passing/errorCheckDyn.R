#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2017-10-30
# Filename: errorCheckDyn.R
# Purpose: Check that errors are caught and reported properly.
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

require(dynr)
vd <- matrix(c(0, -0.1, 1, -0.2), 2, 2)
pd <- matrix(c('fixed', 'spring', 'fixed', 'friction'), 2, 2)


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
		params.dyn=pd,	isContinuousTime=TRUE,
		values.exo=list(matrix(1, 2, 1), matrix(1, 2, 2)),
		covariates='u1'),
	regexp="'values' and 'params' are not all the same size.\nWalter Sobchak says you can't do that.",
	fixed=TRUE)

# Error
testthat::expect_error(
	prep.matrixDynamics(
		values.dyn=vd,
		params.dyn=pd,	isContinuousTime=TRUE,
		values.exo=list(matrix(1, 2, 1), matrix(1, 2, 2)),
		params.exo=list(matrix(0, 2, 1), matrix(0, 2, 2)),
		covariates='u1'),
	regexp="Mind your teaspoons and tablespoons.  The 'exo.values' argument says there are\n (1, 2) covariates, but the 'covariates' arg says there are (1).",
	fixed=TRUE)

#------------------------------------------------------------------------------
# End
