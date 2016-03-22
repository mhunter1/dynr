#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2016-02-25
# Filename: LinearSDE.R
# Purpose: Translate LinearSDE.R into the user-spec for the same model
#  in dynr.
#------------------------------------------------------------------------------



require(dynr)
require(KFAS)

meas <- dynr.matrixLoadings(
	values=matrix(c(1,0), 1, 2),
	params=matrix(0, 1, 2))

ecov <- dynr.matrixErrorCov(
	values.latent=diag(c(0.00001,1)), params.latent=diag(c(0, 3)),
	values.observed=diag(1.5,1), params.observed=diag(4, 1))
ecov$c.string
ecov$startval

initial<-dynr.initial(c(0,1), c(5,0), diag(1,2), diag(0,2))
writeLines(initial$c.string)
initial$startval

dynamics <- dynr.linearDynamics(
	params.dyn=matrix(c(0, 1, 0, 2), 2, 2),
	values.dyn=matrix(c(0, 1, 1, 1), 2, 2),
	time="contin")

regimes <- dynr.regimes()

# Proto-example of cooking
dynr.cook(file=stdout(), meas, ecov$c.string, initial$c.string, dynamics, regimes)


