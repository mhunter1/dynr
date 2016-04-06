rm(list=ls(all=TRUE))
# Example for fitting a nonlinear dynamic factor analysis model with
# nonlinear vector autoregressive (VAR) relation at the factor level

require(dynr)

#------------------------------------------------------------------------------
# Define all the model components via the RECIPE functions
# formula=list(x1~param[6]*x1+param[8]*(exp(abs(x2))/(1+exp(abs(x2))))*x2,
#              x2~param[7]*x2+param[9]*(exp(abs(x1))/(1+exp(abs(x1))))*x1)
# jacob=list(x1~x1~param[6],
#            x1~x2~param[8]*(exp(abs(x2))/(exp(abs(x2))+1)+x2*sign(x2)*exp(abs(x2))/pow(1+exp(abs(x2)),2)),
#            x2~x2~param[7],
#            x2~x1~param[9]*(exp(abs(x1))/(exp(abs(x1))+1)+x1*sign(x1)*exp(abs(x1))/pow(1+exp(abs(x1)),2)))
# dynm<-dynr.nonlindynamics(formula,jacob,isContinuosTime=FALSE)
# writeLines(dynm)

#------------------------------------------------------------------------------
# Create recipes for all the model pieces

# # Measurement (factor loadings)
# loads <- dynr.loadings(
# 	map=list(
# 		eta1=paste0('y', 1:3),
# 		eta2=paste0('y', 4:6)),
# 	params=1:4)
# loads
# meas <- dynr.matrixLoadings(loads$values, loads$params)
# 
# 
# # Initial conditions on the latent state and covariance
# initial <- dynr.initial(
# 	values.inistate=c(0, 0),
# 	params.inistate=c(0, 0),
# 	values.inicov=diag(exp(0), 2), #Should these be left without the exp?
# 	params.inicov=diag(exp(0), 2),
# 	values.regimep=c(.8824, 1-.8824),
# 	params.regimep=c(0, 0)
# )
# 
# cat(initial$c.string)
# # Eek.  Lu, check out the empty for loop!
# 
# 
# # Regime-switching function
# regimes <- dynr.regimes(
# 	values=matrix(c(0, 0, 0, 0), 2, 2), #nrow=numRegimes, ncol=numRegimes*(numCovariates+1)
# 	params=matrix(c(5, 0, 0, 6), 2, 2))
# # Self-transition (diagonals) are estimated
# # Off-diagonal elements are fixed by the softmax scaling
# 
# 
# #measurement and dynamics covariances
# mdcov <- dynr.matrixErrorCov(
# 	values.latent=diag(0, 2),
# 	params.latent=diag(17:18, 2),
# 	values.observed=diag(0, 6),
# 	params.observed=diag(11:16, 6))
# 
# 
# # dynamics
# dynamics <- ""
# 
# 
# #--------------------------------------
# # Cook the recipes together
# fname <- "./demo/CookedRSNonlinearDiscrete.c"  #NOTE: USE MUST BE IN THE dynr DIRECTORY FOR THIS LINE
# dynr.prep(file=fname, meas, mdcov$c.string, initial$c.string, dynamics, regimes)



#------------------------------------------------------------------------------
# Put the cooked recipes together in a Model Specification


#------------------------------------------------------------------------------
# Define all the model components via the RECIPE functions

# # measurement
# meas <- dynr.matrixLoadings(
#   values=matrix(c(1,0), 1, 2),
#   params=matrix(0, 1, 2))
# 
# # observation and dynamic noise components
# ecov <- dynr.matrixErrorCov(
#   values.latent=diag(c(0.00001,1)), params.latent=diag(c(0, 3)),
#   values.observed=diag(1.5,1), params.observed=diag(4, 1))
# ecov$c.string
# ecov$startval
# 
# # initial covariances and latent state values
# initial <- dynr.initial(
#   values.inistate=c(0,1),
#   params.inistate=c(5,0),
#   values.inicov=diag(1,2),
#   params.inicov=diag(0,2))
# writeLines(initial$c.string)
# initial$startval
# 
# # define the differential equation
# dynamics <- dynr.linearDynamics(
#   params.dyn=matrix(c(0, 1, 0, 2), 2, 2),
#   values.dyn=matrix(c(0, 1, 1, 1), 2, 2),
#   time="contin")
# 
# # null regimes, i.e. non-regime switching model
# regimes <- dynr.regimes()
# 
# # Proto-example of cooking
# # put all the strings together
# fname <- "./demo/CookedLinearSDE.c"  #NOTE: USE MUST BE IN THE dynr DIRECTORY FOR THIS LINE
# dynr.prep(file=fname, meas, ecov$c.string, initial$c.string, dynamics, regimes)
# 
# formula=list(x1~param[6]*x1+param[8]*(exp(abs(x2))/(1+exp(abs(x2))))*x2,
#              x2~param[7]*x2+param[9]*(exp(abs(x1))/(1+exp(abs(x1))))*x1)
# jacob=list(x1~x1~param[6],
#            x1~x2~param[8]*(exp(abs(x2))/(exp(abs(x2))+1)+x2*sign(x2)*exp(abs(x2))/pow(1+exp(abs(x2)),2)),
#            x2~x2~param[7],
#            x2~x1~param[9]*(exp(abs(x1))/(exp(abs(x1))+1)+x1*sign(x1)*exp(abs(x1))/pow(1+exp(abs(x1)),2)))
# dynm<-dynr.nonlindynamics(formula,jacob,isContinuosTime=FALSE)
# writeLines(dynm)

#--------------------------------------
# Put the cooked recipes together in a Model Specification

# Data

#Reading in simulated data
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
data <- dynr.data(discreteNA, id="id", time="Time",observed=colnames(discreteNA)[c(3:8)])

pstart <- log(.9/(1-.9))

# Data

#Reading in simulated data
# thedata <- read.table(paste0("./data/NonlinearVARsimT300n10.txt"),na.strings = "NaN",sep=",")
# discreteNA <- t(thedata)
# n = 10; nT = 300
# discreteNA <- cbind(rep(1:n,each=nT),rep(1:nT,n),discreteNA)
# 
# colnames(discreteNA)<-c("id", "Time", "y1", "y2", "y3", "y4", "y5", "y6")
# data <- dynr.data(discreteNA, id="id", time="Time",observed=colnames(discreteNA)[c(3:8)])

#Im = ones(1,InfDS2.nl);
#Eyem = eye(InfDS2.nl);
#Am = [Eyem - (p_Mat);
#      Im        ];
#End = [zeros(InfDS2.nl,1);
#       1      ];
#St0 = (Am' * Am)\Am' * End;       %initial values of pr(St|y_t), Eq. 4.49, p. 71
       
model <- dynr.model(
  num_regime=2,
  dim_latent_var=2,
  xstart=c(1,1,1,1,pstart,pstart,.3,.4,-.5,-.5,
           rep(log(.1),6),
           rep(log(.3),2)),
  ub=rep(9999,18),
  lb=rep(9999,18),
  #ub=c(rep(3,4),5,5,1.5,1.5,1.5,1.5,rep(log(5),6),rep(log(5),2)),
  #lb=c(rep(0,4),-5,-5,-1.5,-1.5,-1.5,-1.5,rep(log(.001),6),rep(log(.001),2)),
  options=list(maxtime=60*60, maxeval=1000,ftol_rel=as.numeric(1e-10),
               xtol_rel=as.numeric(1e-7)),
  isContinuousTime=FALSE,
  infile="./demo/RSNonlinearDiscrete.c", 
  outfile="./demo/RSNonlinearDiscrete2", 
  verbose=TRUE,
  compileLib=TRUE
)

# Estimate free parameters
tfun <- function(x){c(x[1:4],
                      exp(x[5])/(1+exp(x[5])), exp(x[6])/(1+exp(x[6])),
                      x[7:10],
                      exp(x[11:18]))}
res <- dynr.cook(model=model, data=data,transformation=tfun,debug_flag=FALSE)

# Examine results
summary(res)


#True values = 
#c(1.2, 1.2, 1.1, .95, .98, .85, .2, .25, -.6, -.8,
# .28, .10, .12, .13, .12, .11,
# .35, .3)

truepar <- c(1.2, 1.2, 1.1, .95, 
             log(.98/(1-.98)), log(.85/(1-.85)), 
             .2, .25, -.6, -.8,
             log(c(.28, .10, .12, .13, .12, .11)),
             log(c(.35, .3)))
#1.200000  1.200000  1.100000  0.950000  3.891820  1.734601  
#0.200000  0.250000 -0.600000 -0.800000 
#-1.272966 -2.302585 -2.120264-2.040221 -2.120264 -2.207275 -1.049822 -1.203973
#------------------------------------------------------------------------------
# End


#plot(res, data=data, graphingPar=list(cex.main=1, cex.axis=1, cex.lab=1.2), numSubjDemo=2)

#dynr.ggplot(res, data.dynr=data, states=c(1,2), names.regime=1:2,title="Smoothed State Values", numSubjDemo=2)

