#------------------------------------------------------------------------------
# Author: Lu Ou
# Date: 2016-04-13
# Filename: RSNonlinearDFADiscreteTime.R
# Purpose: An illustrative example of using dynr to fit
#   a regime-switching nonlinear dynamic factor analysis (discrete-time) model 
#   with nonlinear vector autoregressive (VAR) relation at the factor level
#------------------------------------------------------------------------------

rm(list=ls(all=TRUE))
require(dynr)

#------------------------------------------------------------------------------
# Preparation

# Define all the model components via the RECIPE functions

# Measurement (factor loadings)
meas <- prep.loadings(
  map=list(
    eta1=paste0('y', 1:3),
    eta2=paste0('y', 4:6)),
  params=1:4)


# Initial conditions on the latent state and covariance
initial <- prep.initial(
	values.inistate=c(0, 0),
	params.inistate=c(0, 0),
	values.inicov=diag(1, 2), 
	params.inicov=diag(1, 2),
	values.regimep=c(.8824, 1-.8824),
	params.regimep=c(0, 0)
)

# Regime-switching function
regimes <- prep.regimes(
	values=matrix(c(0, 0, 0, 0), 2, 2), #nrow=numRegimes, ncol=numRegimes*(numCovariates+1)
	params=matrix(c(5, 0, 0, 6), 2, 2))
# Self-transition (diagonals) are estimated
# Off-diagonal elements are fixed by the softmax scaling


#measurement and dynamics covariances
mdcov <- prep.noise(
	values.latent=diag(0, 2),
	params.latent=diag(17:18, 2),
	values.observed=diag(0, 6),
	params.observed=diag(11:16, 6))

# dynamics
formula=list(
  list(x1~param[6]*x1,
       x2~param[7]*x2),
  list(x1~param[6]*x1+param[8]*(exp(abs(x2))/(1+exp(abs(x2))))*x2,
       x2~param[7]*x2+param[9]*(exp(abs(x1))/(1+exp(abs(x1))))*x1) 
)

jacob=list(
  list(x1~x1~param[6],
       x2~x2~param[7]),
  list(x1~x1~param[6],
       x1~x2~param[8]*(exp(abs(x2))/(exp(abs(x2))+1)+x2*sign(x2)*exp(abs(x2))/(1+exp(abs(x2))^2)),
       x2~x2~param[7],
       x2~x1~param[9]*(exp(abs(x1))/(exp(abs(x1))+1)+x1*sign(x1)*exp(abs(x1))/(1+exp(abs(x1))^2))))

dynm<-prep.formulaDynamics(formula=formula,startval=c(.3,.4,-.5,-.5),isContinuousTime=FALSE,jacobian=jacob)
#cat(writeCcode(dynm)$c.string)

#------------------------------------------------------------------------------
# Cooking materials

# Model
# Put all the recipes together in a Model Specification
model <- dynr.model(dynamics=dynm, measurement=meas, noise=mdcov, initial=initial,regimes=regimes, outfile="RSNonlinearDFA")
# View specified model in latex
printex(model)

# Data
thedata <- read.table(paste0("./data/NonlinearVARsimT300n10NoMissing.txt"),na.strings = "NaN",sep=",")
discreteNA <- t(thedata)
n = 10; nT = 300
discreteNA <- cbind(rep(1:n,each=nT),rep(1:nT,n),discreteNA)
colnames(discreteNA)<-c("id", "Time", "y1", "y2", "y3", "y4", "y5", "y6")
missingprop=.1
partialmiss=sample(1:nT,nT*missingprop)
fullmiss=sample(setdiff(1:nT,partialmiss),nT*missingprop)
discreteNA[partialmiss, c("y1","y5","y6")]<-NA
discreteNA[fullmiss, c("y1","y2","y3","y4","y5","y6")]<-NA
head(discreteNA)
data <- dynr.data(discreteNA, id="id", time="Time",observed=colnames(discreteNA)[c(3:8)])



# Estimate free parameters
tfun <- function(x){c(x[1:4],
                      exp(x[5])/(1+exp(x[5])), exp(x[6])/(1+exp(x[6])),
                      x[7:10],
                      exp(x[11:18]))}
res <- dynr.cook(model=model, data=data,transformation=tfun,debug_flag=FALSE)

# Examine results
summary(res)


#plot(res, data=data, graphingPar=list(cex.main=1, cex.axis=1, cex.lab=1.2), numSubjDemo=2)

#dynr.ggplot(res, data.dynr=data, states=c(1,2), names.regime=1:2,title="Smoothed State Values", numSubjDemo=2)

# Examine results
summary(res)

#------------------------------------------------------------------------------
# some miscellaneous nice functions

# get the estimated parameters from a cooked model/data combo
coef(res)


# get the log likelihood, AIC, and BIC from a cooked model/data combo
logLik(res)
AIC(res)
BIC(res)


# compare true parameters to estimated ones
truepar <- c(1.2, 1.2, 1.1, .95, 
             log(.98/(1-.98)), log(.85/(1-.85)), 
             .2, .25, -.6, -.8,
             log(c(.28, .10, .12, .13, .12, .11)),
             log(c(.35, .3)))
#1.200000  1.200000  1.100000  0.950000  3.891820  1.734601  
#0.200000  0.250000 -0.600000 -0.800000 
#-1.272966 -2.302585 -2.120264-2.040221 -2.120264 -2.207275 -1.049822 -1.203973
data.frame(name=c('Spring', 'Damping', 'DynVar', 'MeasVar', 'IniPos'), true=truepar, estim=coef(res))


# compare estimated smoothed latent states to true
# simulated ones
sm <- data.frame(t(res@eta_smooth_final))
cor(sm, t(tx)[-1,])


#------------------------------------------------------------------------------
# End


