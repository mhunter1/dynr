#------------------------------------------------------------------------------
# Author: Yanling Li
# Date: 2020-03-21
# Filename: errorCheckNoise.R
# Purpose: Check that errors are caught and reported properly.
#------------------------------------------------------------------------------
# Load packages
require(dynr)
vd <- matrix(c(0.4, 0.3, -0.3, -0.2), 2, 2)
pd <- matrix(c('a', 'a1', 'b', 'b1'), 2, 2)

#------------------------------------------------------------------------------
# Check dimension conformability errors
# Specifically, the dimension of the noise covariance matrix in prep.noise() is different from the 
# number of latent/observed variables in prep.measurement()

dynamics <- prep.matrixDynamics(
  values.dyn=vd,
  params.dyn=pd,
  isContinuousTime=TRUE)

data <- dynr.data(data.frame(id=1,time=1:10,wp=1:10,hp=1:10), id="id", time="time", observed=c("wp","hp"))

measurement <- prep.measurement(
  values.load=matrix(c(1,0,
                       0,1),ncol=2,byrow=T), 
  params.load=matrix(rep("fixed",4),ncol=2),
  state.names=c("wp","hp"),
  obs.names=c("wp","hp")
)


initial <- prep.initial(
  values.inistate=c(0,0),
  params.inistate=c('fixed', 'fixed'),
  values.inicov=matrix(diag(1,2),byrow=T,ncol=2),
  params.inicov=matrix(c("fixed","fixed",
                         "fixed","fixed"),byrow=T,ncol=2))


noise1 <- prep.noise(
  values.latent=diag(0.1,2), 
  params.latent=matrix(c("v_wp","c_hw",
                         "c_hw","v_hp"),byrow=T,ncol=2), 
  values.observed=diag(0,1), 
  params.observed=diag('errorvar', 1)) 

noise2 <- prep.noise(
  values.latent=diag(0.1,1), 
  params.latent=diag('noisevar', 1), 
  values.observed=diag(0,2), 
  params.observed=diag(c('fixed','fixed'),2))

testthat::expect_error(
  dynrmodel <- dynr.model(dynamics, measurement, noise1, initial, data),
  regexp="The dimension of the measurement noise convariance matrix in prep.noise should match the number of observed variables in 'dynrMeasurement'.",
  fixed=TRUE
)

testthat::expect_error(
  dynrmodel <- dynr.model(dynamics, measurement, noise2, initial, data),
  regexp="The dimension of the dynamic noise covariance matrix in prep.noise should match the number of latent states in 'dynrMeasurement'.",
  fixed=TRUE
)


#------------------------------------------------------------------------------
# End
