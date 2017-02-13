#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2016-05-24
# Filename: LinearSDE.R
# Purpose: An illustrative example of using dynr to fit
#   a linear stochastic differential equation model
#------------------------------------------------------------------------------

#rm(list=ls(all=TRUE))

#------------------------------------------------------------------------------
# Load packages
require(dynr)
require(mvtnorm)
require(Matrix)

#------------------------------------------------------------------------------
# Example 2
# Damped linear oscillator example
#	There is measurement noise and there are unmeasured dynamic disturbances.
#	These disturbances are called dynamic noise.
#	There is a single indicator.
# TODO Figure out what combination of variables can really be estimated for 
#	this kind of model.  It is clear that measurement noise and SOME dynamic
#	noise can be estimated, but not all.


#--------------------------------------
# Data Generation

xdim <- 2
udim <- 1
ydim <- 1
tdim <- 1000
set.seed(315)
tA <- matrix(c(0, -.3, 1, -.7), xdim, xdim)
tB <- matrix(c(0), xdim, udim)
tC <- matrix(c(1, 0), ydim, xdim)
tD <- matrix(c(0), ydim, udim)
tQ <- matrix(c(0), xdim, xdim); diag(tQ) <- c(0, 2.2)
tR <- matrix(c(0), ydim, ydim); diag(tR) <- c(1.5)

x0 <- matrix(c(0, 1), xdim, 1)
P0 <- diag(c(1), xdim)
tdx <- matrix(0, xdim, tdim+1)
tx <- matrix(0, xdim, tdim+1)
tu <- matrix(0, udim, tdim)
ty <- matrix(0, ydim, tdim)

tT <- matrix(0:tdim, nrow=1, ncol=tdim+1)

tI <- diag(1, nrow=xdim)

tx[,1] <- x0
for(i in 2:(tdim+1)){
	q <- t(rmvnorm(1, rep(0, xdim), tQ))
	tdx[,i] <- tA %*% tx[,i-1] + tB %*% tu[,i-1] + q
	expA <- as.matrix(expm(tA * (tT[,i]-tT[,i-1])))
	intA <- solve(tA) %*% (expA - tI)
	tx[,i] <- expA %*% tx[, i-1] + intA %*% tB %*% tu[,i-1] + intA %*% q
	ty[,i-1] <- tC %*% tx[,i] + tD %*% tu[,i-1] + t(rmvnorm(1, rep(0, ydim), tR))
}



rownames(ty) <- paste('y', 1:ydim, sep='')
rownames(tx) <- paste('x', 1:xdim, sep='')


#plot(tx[1,], type='l')
#plot(tT[,-1], ty[1,], type='l')


#------------------------------------------------------------------------------
# Define all the model components via the RECIPE functions

# measurement
# this is the factor loadings matrix, Lambda in SEM notation or C in OpenMx notation
meas <- prep.measurement(
	values.load=matrix(c(1, 0), 1, 2), # starting values and fixed values
	params.load=matrix(c('fixed', 'fixed'), 1, 2),
	state.names=c("Position","Velocity"),
	obs.names=c("y1")) # parameter numbers or indication that parameter is fixed
# Look
meas
# no free parameters in the factor loadings


# observation and dynamic noise components
# the latent noise is the dynamic noise, Psi in SEM notation or Q in OpenMx notation
# the observed noise is the measurement noise, Theta in SEM notation or R in OpenMx notation
ecov <- prep.noise(
	values.latent=diag(c(0, 1), 2), params.latent=diag(c('fixed', 'dnoise'), 2), # uses free parameter 3
	values.observed=diag(1.5, 1), params.observed=diag('mnoise', 1)) # uses free parameter 4
# Look
ecov
# dynr takes steps to make sure covariance matrices are positive definite

#ecov <- prep.noise(
#	values.latent=rep(list(diag(c(0, 1), 2)), 2), params.latent=rep(list(diag(c('fixed', 'dnoise'), 2)), 2), # uses free parameter 3
#	values.observed=list(diag(0.5, 1), diag(1.5, 1)), params.observed=list(diag('mnoise1', 1), diag('mnoise2', 1))) # uses free parameter 4


# initial covariances and latent state values
# These initialize the recursive algorithm (extended Kalman filter) that dynr uses
# These are x0 and P0 in OpenMx notation
initial <- prep.initial(
	values.inistate=c(0, 1),
	params.inistate=c('inipos', 'fixed'), #initial position is free parameter 5, initial slope is fixed at 1
	values.inicov=diag(1, 2),
	params.inicov=diag('fixed', 2)) #initial covariance is fixed to a diagonal matrix of 1s.

#initial <- prep.initial(
#	values.inistate=list(c(0, 1), c(0, 1)),
#	params.inistate=list(c('inipos', 'fixed'), c('fixed', 'inivel')), #initial position is free parameter 5, initial slope is fixed at 1
#	values.inicov=list(diag(1, 2), diag(1, 2)),
#	params.inicov=list(diag('fixed', 2), diag('fixed', 2))) #initial covariance is fixed to a diagonal matrix of 1s.


# define the differential equation
dynamics <- prep.matrixDynamics(
	values.dyn=matrix(c(0, -0.1, 1, -0.2), 2, 2),
	params.dyn=matrix(c('fixed', 'spring', 'fixed', 'friction'), 2, 2), #uses parameters 1 and 2
	isContinuousTime=TRUE)

# Data
simdata <- cbind(id=rep(1,100), t(ty), times=tT[,-1])
data <- dynr.data(simdata, id="id", time="times", observed="y1")

# Prepare for cooking
# put all the recipes together
model <- dynr.model(dynamics=dynamics, measurement=meas, noise=ecov, initial=initial, data=data, outfile="LinearSDE.c")

# set upper bounds, if you want
model$ub <- c(100, 100, 100, 100, 100)
model$ub <- c(friction=101, spring=100, inipos=103, 100, 100)
model$ub['dnoise'] <- 99

# Developer
# Check that an error message is thrown
testthat::expect_error(model$ub <- c(friction=101, spring=100, inipos=103, 100, 100, 4), regexp="I'm going over my borders. You gave me 6 things, but I need 5 (the number of free parameters).", fixed=TRUE)


printex(model,ParameterAs=model$param.names,show=FALSE,printInit=TRUE,
        outFile="LinearSDE.tex")
#tools::texi2pdf("LinearSDE.tex")
#system(paste(getOption("pdfviewer"), "LinearSDE.pdf"))


# Estimate free parameters
res <- dynr.cook(model)
# Examine results
summary(res)

plotFormula(model, res@transformed.parameters)
#ggsave("LinearSDEPlotFml.pdf")
#------------------------------------------------------------------------------
# some miscellaneous nice functions

plot(res, dynrModel=model, textsize=6, style = 1)
plot(res, dynrModel=model, textsize=6, style = 2)


autoplot(res, model, numSubjDemo=1)


# get the estimated parameters from a cooked model/data combo
coef(res)

# get the log likelihood, AIC, and BIC from a cooked model/data combo
logLik(res)
AIC(res)
BIC(res)


# compare true parameters to estimated ones
trueParams <- c(-.3, -.7, 2.2, 1.5, 0)
data.frame(name=c('Spring', 'Damping', 'DynVar', 'MeasVar', 'IniPos'), true=trueParams, estim=coef(res))

(CI <- confint(res))

# Check that all true parameters are within the confidence intervals of the estimated params
withinIntervals <- CI[,1] < trueParams & trueParams < CI[,2]
testthat::expect_true(all(withinIntervals))

# compare estimated smoothed latent states to true
# simulated ones
sm <- data.frame(t(res@eta_smooth_final))
cor(sm, t(tx)[-1,])

#------------------------------------------------------------------------------
# End
#save(model,res,file="LinearSDE.RData")
