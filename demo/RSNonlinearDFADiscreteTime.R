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
  params=c("lambda_{21}","lambda_{31}","lambda_{52}","lambda_{62}"))


# Initial conditions on the latent state and covariance
initial <- prep.initial(
	values.inistate=c(0, 0),
	params.inistate=c("fixed", "fixed"),
	values.inicov=diag(1, 2), 
	params.inicov=diag("fixed", 2),
	values.regimep=c(.8824, 1-.8824),
	params.regimep=c("fixed", "fixed")
)

# Regime-switching function
regimes <- prep.regimes(
	values=matrix(c(.9, 0, 0, .9), 2, 2), #nrow=numRegimes, ncol=numRegimes*(numCovariates+1)
	params=matrix(c("p11", 0, 0, "p22"), 2, 2))
# Self-transition (diagonals) are estimated
# Off-diagonal elements are fixed by the softmax scaling


#measurement and dynamics covariances
mdcov <- prep.noise(
	values.latent=diag(0.3, 2),
	params.latent=diag(paste0("zeta_",1:2), 2),
	values.observed=diag(0.1, 6),
	params.observed=diag(paste0("epsilon_",1:6), 6))

# dynamics
formula=list(
  list(x1~a1*x1,
       x2~a2*x2),
  list(x1~a1*x1+c12*(exp(abs(x2))/(1+exp(abs(x2))))*x2,
       x2~a2*x2+c21*(exp(abs(x1))/(1+exp(abs(x1))))*x1) 
)

jacob=list(
  list(x1~x1~a1,
       x2~x2~a2),
  list(x1~x1~a1,
       x1~x2~c12*(exp(abs(x2))/(exp(abs(x2))+1)+x2*sign(x2)*exp(abs(x2))/(1+exp(abs(x2))^2)),
       x2~x2~a2,
       x2~x1~c21*(exp(abs(x1))/(exp(abs(x1))+1)+x1*sign(x1)*exp(abs(x1))/(1+exp(abs(x1))^2))))

dynm<-prep.formulaDynamics(formula=formula,startval=c(a1=.3,a2=.4,c12=-.5,c21=-.5),isContinuousTime=FALSE,jacobian=jacob)
#cat(writeCcode(dynm)$c.string)

trans<-prep.tfun(formula.trans=list(p11~exp(p11)/(1+exp(p11)), p22~exp(p22)/(1+exp(p22))), formula.inv=list(p11~log(p11/(1-p11)),p22~log(p22/(1-p22))), transCcode=FALSE)
#cat(writeCcode(trans)$c.string)
#------------------------------------------------------------------------------
# Cooking materials

# Model
# Put all the recipes together in a Model Specification
model <- dynr.model(dynamics=dynm, measurement=meas, noise=mdcov, 
                    initial=initial, regimes=regimes, transform=trans, 
                    outfile="RSNonlinearDFA")
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
res <- dynr.cook(model, data=data, debug_flag=FALSE)

# Examine results
summary(res)

# get the estimated parameters from a cooked model/data combo
coef(res)

# get the log likelihood, AIC, and BIC from a cooked model/data combo
logLik(res)
AIC(res)
BIC(res)


# compare true parameters to estimated ones
truepar <- c(
  .2, .25, -.6, -.8,
  1.2, 1.2, 1.1, .95, 
  c(.35, .3),
  c(.28, .10, .12, .13, .12, .11),
  0.98,0.85)
data.frame(name=res@param.names , true=truepar, estim=coef(res))

dynr.ggplot(res, data.dynr=data, states=c(1,2), names.regime=1:2,title="Smoothed State Values", numSubjDemo=2)

#------------------------------------------------------------------------------
# End


