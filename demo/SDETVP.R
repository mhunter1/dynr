#--------------------------------------------------------
# Filename: LogisticSetpoint.R
# Purpose: An illustrative example of using dynr to fit
#   a stochastic damped oscillator model with time-
#   varying (logistic growth) setpoint parameter
#   (Chen, Chow & Hunter, expected 2018)
#--------------------------------------------------------

#---- (1) Load packages ----
library(dynr)
#---- (2) Read in data ----
data(LogisticSetPointSDE)
dataLog <- dynr.data(LogisticSetPointSDE, id="id", time="times", observed="obsy")
#---- (3) Specify recipes for all model pieces ----
#---- (3a) Measurement ----
measLog <- prep.measurement(
  values.load=matrix(c(1, 0, 0), 1, 3), # starting values and fixed values
  params.load=matrix("fixed", 1, 3),
  state.names=c("x","y","z"),
  obs.names=c("obsy")) # parameter numbers or indication that parameter is fixed

#---- (3b) Dynamic and measurement noise cov structures----
ecovLog <- prep.noise(
  values.latent=diag(c(0,0.1,0), 3), params.latent=diag(c('fixed', 'pnoise','fixed'), 3), 
  values.observed=diag(0.1, 1), params.observed=diag('mnoise', 1)) 

#---- (3c) Initial condition specification ----
initialLog <- prep.initial(
  values.inistate=c(-5, 0 ,.1),
  params.inistate=c('fixed', 'fixed','fixed'), 
  values.inicov=diag(0, 3),
  params.inicov=diag('fixed', 3)) 


#---- (3d) Dynamic model ----
formulaLogistic<-list(
  x~y,
  y~eta*(x-z)+damp*y,
  z~r*z*(1-b*z)
)

dynamicsLog<-prep.formulaDynamics(formula=formulaLogistic,startval=c(eta=-.1,damp=-.05,b=0.5,r=0.3),
                                  isContinuousTime=TRUE)

#---- (4) Create model and cook it all up ----

modelLog <- dynr.model(dynamics=dynamicsLog, measurement=measLog, 
                       noise=ecovLog, initial=initialLog, #transform=trans,
                       data=dataLog, outfile="setpointLog.c")

modelLog$ub['pnoise'] <- 10


resLog <- dynr.cook(modelLog)
#---- (5) Serve it! ----

summary(resLog)

# True parameters: freq=-1; damp=-.1; b=0.1; r=.5; sigma1=0.1, mnoise=1
