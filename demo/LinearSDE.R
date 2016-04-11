#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2016-02-25
# Filename: LinearSDE.R
# Purpose: Translate LinearSDE.R into the user-spec for the same model
#  in dynr.
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Load packages

require(mvtnorm)
require(Matrix)
require(dynr)

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
	values=matrix(c(1, 0), 1, 2), # starting values and fixed values
	params=matrix(c('fixed', 'fixed'), 1, 2)) # parameter numbers or indication that parameter is fixed
# Look
meas
# no free parameters in the factor loadings


# observation and dynamic noise components
# the latent noise is the dynamic noise, Psi in SEM notation or Q in OpenMx notation
# the observed noise is the measurement noise, Theta in SEM notation or R in OpenMx notation
ecov <- prep.noise(
	values.latent=diag(c(0, 1), 2), params.latent=diag(c('fixed', 3), 2), # uses free parameter 3
	values.observed=diag(1.5, 1), params.observed=diag(4, 1)) # uses free parameter 4
# Look
ecov
# ... wait what happened to the values?
log(c(1e-6, 1))
log(1.5)
# dynr takes steps to make sure covariance matrices are positive definite


# initial covariances and latent state values
# These initialize the recursive algorithm (extended Kalman filter) that dynr uses
# These are x0 and P0 in OpenMx notation
initial <- prep.initial(
	values.inistate=c(0, 1),
	params.inistate=c(5, 'fixed'), #initial position is free parameter 5, initial slope is fixed at 1
	values.inicov=diag(1, 2),
	params.inicov=diag('fixed', 2)) #initial covariance is fixed to a diagonal matrix of 1s.


# define the differential equation
dynamics <- prep.linearDynamics(
	values.dyn=matrix(c(0, -0.1, 1, -0.2), 2, 2),
	params.dyn=matrix(c('fixed', 1, 'fixed', 2), 2, 2), #uses parameters 1 and 2
	time="continuous")


# Prepare for cooking
# put all the recipes together
model <- dynr.prep(dynamics=dynamics, measurement=meas, noise=ecov, initial=initial, outfile="cooked")


# Data
simdata <- cbind(id=rep(1,100), t(ty), times=tT[,-1])
data <- dynr.data(simdata, id="id", time="times", observed="y1")

#define a transformation function for those pesky variances
tfun <- function(x){return(c(x[1:2], exp(x[3:4]), x[5]))}

# Estimate free parameters
res <- dynr.cook(model, data, tfun)

# Examine results
summary(res)

#------------------------------------------------------------------------------
# some miscellaneous nice functions

# print recipes in forms that will look nice in LaTeX
printex(meas)

printex(ecov)


# get the estimated parameters from a cooked model/data combo
coef(res)


# get the log likelihood, AIC, and BIC from a cooked model/data combo
logLik(res)
AIC(res)
BIC(res)


# compare true parameters to estimated ones
trueParams <- c(-.3, -.7, 2.2, 1.5, 0)
data.frame(name=c('Spring', 'Damping', 'DynVar', 'MeasVar', 'IniPos'), true=trueParams, estim=coef(res))


# compare estimated smoothed latent states to true
# simulated ones
sm <- data.frame(t(res@eta_smooth_final))
cor(sm, t(tx)[-1,])


#------------------------------------------------------------------------------
# End

