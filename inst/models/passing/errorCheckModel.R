#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2021-03-19 17:35:47
# Filename: errorCheckModel.R
# Purpose: Check that we properly catch and throw errors from dynr.model()
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# Load packages
require(dynr)


#------------------------------------------------------------------------------
# Minimal Model

meas <- prep.measurement(
	values.load=matrix(c(1, 0), 1, 2),
	params.load=matrix(c('fixed', 'fixed'), 1, 2),
	state.names=c("Position","Velocity"),
	obs.names=c("y1"))

ecov <- prep.noise(
	values.latent=diag(c(0, 1), 2),
	params.latent=diag(c('fixed', 'dnoise'), 2),
	values.observed=diag(1.5, 1),
	params.observed=diag('mnoise', 1))

initial <- prep.initial(
	values.inistate=c(0, 1),
	params.inistate=c('inipos', 'fixed'),
	values.inicov=diag(1, 2),
	params.inicov=diag('fixed', 2))

dynamics <- prep.matrixDynamics(
	values.dyn=matrix(c(0, -0.1, 1, -0.2), 2, 2),
	params.dyn=matrix(c('fixed', 'spring', 'fixed', 'friction'), 2, 2),
	isContinuousTime=TRUE)

data(Oscillator)
data <- dynr.data(Oscillator, id="id", time="times", observed="y1")

model <- dynr.model(dynamics=dynamics, measurement=meas,
	noise=ecov, initial=initial, data=data)


#------------------------------------------------------------------------------
# Check the type of all arguments

testthat::expect_error(
	dynr.model(dynamics=meas, measurement=meas,
		noise=ecov, initial=initial, data=data),
	regexp="Check that the 'dynamics' argument is of the correct class.\nHint: it should be either a 'dynrDynamicsFormula' or 'dynrDynamicsMatrix'.\nCreate it with prep.formulaDynamics() or prep.matrixDynamics().",
	fixed=TRUE)

testthat::expect_error(
	dynr.model(dynamics=dynamics, measurement=dynamics,
		noise=ecov, initial=initial, data=data),
	regexp="Check that the 'measurement' argument is of the correct class.\nHint: it should be a 'dynrMeasurement'.\nCreate it with prep.measurement() or prep.loadings().",
	fixed=TRUE)

testthat::expect_error(
	dynr.model(dynamics=dynamics, measurement=meas,
		noise=dynamics, initial=initial, data=data),
	regexp="Check that the 'noise' argument is of the correct class.\nHint: it should be a 'dynrNoise'.\nCreate it with prep.noise().",
	fixed=TRUE)

testthat::expect_error(
	dynr.model(dynamics=dynamics, measurement=meas,
		noise=ecov, initial=dynamics, data=data),
	regexp="Check that the 'initial' argument is of the correct class.\nHint: it should be a 'dynrInitial'.\nCreate it with prep.initial().",
	fixed=TRUE)

testthat::expect_error(
	dynr.model(dynamics=dynamics, measurement=meas,
		noise=ecov, initial=initial, data=dynamics),
	regexp="Check that the 'data' argument is of the correct class.\nHint: it should be a 'list'.\nCreate it with dynr.data().",
	fixed=TRUE)

testthat::expect_error(
	dynr.model(dynamics=dynamics, measurement=meas,
		noise=ecov, initial=initial, data=Oscillator),
	regexp="Check that the 'data' argument is of the correct class.\nHint: it should be a 'list'.\nCreate it with dynr.data().",
	fixed=TRUE)


#------------------------------------------------------------------------------
# Check whether NAs are inserted correctly with discrete-time models
data(Oscillator)
Oscillator$id[401:800] <- 2
Oscillator$id[801:1000] <- 3
Oscillator$y1[c(10,799)] <- NA
Oscillator <- Oscillator[!is.na(Oscillator$y1),]
data <- dynr.data(Oscillator, id="id", time="times", observed="y1")
Oscillator_id1 <- Oscillator[Oscillator$id==1,]
data_id1 <- dynr.data(Oscillator_id1, id="id", time="times", observed="y1")

dynamics <- prep.matrixDynamics(
  values.dyn=matrix(c(0, -0.1, 1, -0.2), 2, 2),
  params.dyn=matrix(c('fixed', 'spring', 'fixed', 'friction'), 2, 2),
  isContinuousTime=FALSE)

testthat::expect_invisible({
  model <- dynr.model(dynamics=dynamics, measurement=meas,
             noise=ecov, initial=initial, data=data)
}
)

testthat::expect_true(
  {model <- dynr.model(dynamics=dynamics, measurement=meas,
                      noise=ecov, initial=initial, data=data_id1);
  all(diff(model$data$time)==1)}
)

#------------------------------------------------------------------------------



#------------------------------------------------------------------------------



#------------------------------------------------------------------------------



#------------------------------------------------------------------------------

