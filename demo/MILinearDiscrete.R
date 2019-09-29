#---------------------------------------------------------------------
# Author: Yanling Li and Linying Ji
# Date: 2019-09-29
# Filename: MILinearDiscrete.R
# Purpose: An illustrative example of using dynr.mi to implement
# multiple imputation with a vector autoregressive model
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# Load packages

require(dynr)


#---------------------------------------------------------------------
# The data were generated from a vector autogressive (VAR) model 
# with two observed dependent variables (i.e., wp and hp) and 
# two covariates (i.e., ca and cn).
# x1 and x2 are two auxiliary varialbes used in the generation of missing data.
# The data contain 100 subjects (N = 100) and 100 time points (T = 100) for each subject.
# Missing data were under missing at random (MAR) condtion with around 30% missing rate.
# The true parameters used to generate the data:
# auto-regression parameters
a=0.4; a1=0.3; 
# cross-regression parameters
b=-0.3; b1=-0.2; 

# coefficients of covariates
c=0.3; c1=0.3;
d=-0.5; d1=-0.4;

# noise variance and covariance
v_wp = 1
c_hw = 0.3
v_hp = 1

# Load data 
data(VARsim)
VARsim$ca <- as.factor(VARsim$ca)

# Declare the data with the dynr.data() function
rawdata <- dynr.data(VARsim, id="ID", time="Time", 
                     observed=c("wp","hp"),covariates=c("ca","cn"))

# Define elements of the measurement model
meas <- prep.measurement(
  values.load=matrix(c(1,0,
                       0,1),ncol=2,byrow=T), 
  params.load=matrix(rep("fixed",4),ncol=2),
  state.names=c("wp","hp"),
  obs.names=c("wp","hp")
)

# Define elements of the dynamic model
formula =list(wp ~ a*wp + b*hp + c*ca + d*cn,
       hp ~ a1*hp + b1*wp +c1*ca + d1*cn)

dynm  <- prep.formulaDynamics(formula=formula,
                              startval=c(a = .4, b = -.3, b1=-.2, a1=.3, 
                                         c = .3, c1=.3, d=-.5, d1=-.4
                              ), isContinuousTime=FALSE)

# Define the initial conditions of the model
initial <- prep.initial(
  values.inistate=c(.15,.15),
  params.inistate=c('mu_wp', 'mu_hp'),
  values.inicov=matrix(c(1,.1,
                         .1,1),byrow=T,ncol=2),
  params.inicov=matrix(c("v_11","c_21",
                         "c_21","v_22"),byrow=T,ncol=2))

# Define the covariance structures of the measurement noise
# covariance matrix and the dynamic noise covariance matrix
mdcov <- prep.noise(
  values.latent=matrix(c(1,.3,
                         .3,1),byrow=T,ncol=2), 
  params.latent=matrix(c("v_wp","c_hw",
                         "c_hw","v_hp"),byrow=T,ncol=2), 
  values.observed=diag(rep(0,2)), 
  params.observed=diag(c('fixed','fixed'),2)) 


# Pass data and submodels to dynrModel object
model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=rawdata,
                    outfile=paste("trial.c",sep=""))

# Plot the Formula 
printex(model, ParameterAs = model$param.names, printInit = TRUE, printRS = FALSE,
        outFile = "MILinearDiscrete.tex")

plotFormula(model, ParameterAs = model$param.names, printDyn = TRUE, printMeas = TRUE)

# An example of using dynr.mi() function to implement multiple imputation and parameter estimation procedures
result <- dynr.mi(model, which.aux=c("x1","x2"), 
  which.lag=c("wp","hp"), lag=1,
  which.lead=NULL, lead=0,
  m=5, iter=10, 
  imp.obs=FALSE, imp.exo=TRUE,
  diag = TRUE, Rhat=1.1,
  conf.level=0.95,
  verbose=FALSE, seed=12345)

# Compare true parameters to estimated ones
truep <- c(a=0.4, b=-0.3, b1=-0.2, a1=0.3,
           c=0.3, c1=0.3, d=-0.5, d1=-0.4,
           v_wp = 1, c_hw = 0.3,v_hp = 1)
estp <- result$estimation.result[1:11, ]
data.frame(truep, estp)

# Convergence diagnostic check 
# Rhat plot
result$Rhat.plot
