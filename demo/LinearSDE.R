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
	#ty[,i-1] <- tC %*% tx[,i-1] + tD %*% tu[,i-1] + t(rmvnorm(1, rep(0, ydim), tR))
	ty[,i-1] <- tC %*% tx[,i] + tD %*% tu[,i-1] + t(rmvnorm(1, rep(0, ydim), tR))
}

#plot(tx[1,], type='l')

rownames(ty) <- paste('y', 1:ydim, sep='')
rownames(tx) <- paste('x', 1:xdim, sep='')

#plot(tT[,-1], ty[1,], type='l')


#------------------------------------------------------------------------------
# Define all the model components via the RECIPE functions

# measurement
meas <- prep.matrixLoadings(
	values=matrix(c(1, 0), 1, 2),
	params=matrix(0, 1, 2))

# observation and dynamic noise components
ecov <- prep.matrixErrorCov(
	values.latent=diag(c('Free', 1)), params.latent=diag(c('fixed', 3)),
	values.observed=diag(1.5,1), params.observed=diag(4, 1))
ecov$c.string
ecov$startval

# initial covariances and latent state values
initial <- prep.initial(
	values.inistate=c('freed',1),
	params.inistate=c(5,'fix'),
	values.inicov=diag(1,2),
	params.inicov=diag('fix',2))
writeLines(initial$c.string)
initial$startval

# define the differential equation
dynamics <- prep.linearDynamics(
	params.dyn=matrix(c(0, 1, 0, 2), 2, 2),
	values.dyn=matrix(c(0, 1, 1, 1), 2, 2),
	time="contin")


# Proto-example of cooking
# put all the strings together
fname <- "./demo/CookedLinearSDE.c"  #NOTE: USE MUST BE IN THE dynr DIRECTORY FOR THIS LINE
dynr.prep(file=fname, meas, ecov$c.string, initial$c.string, dynamics)


#--------------------------------------
# Put the cooked recipes together in a Model Specification

# Data
data <- dynr.data(cbind(id=rep(1,100),t(ty), times=tT[,-1]), id="id", time="times", observed="y1")

# Model
model <- dynr.model(
              num_regime=1,
              dim_latent_var=2,
              xstart=c(-0.1,-0.2,log(1),log(1.5),0),
              ub=c(rep(9999,4),10),lb=c(-10,-10,log(10^(-6)),9999,-10),
              options=list(maxtime=30*60, 
                           maxeval=5000,
                           ftol_rel=as.numeric(1e-8),
                           xtol_rel=as.numeric(1e-8)),
              isContinuousTime=TRUE,
              infile=fname, #Cooked recipes go here
              outfile="./demo/LinearSDE2", 
              verbose=TRUE,
              compileLib=TRUE
)

# Estimate free parameters
res <- dynr.cook(model, data)

# Examine results
summary(res)

#------------------------------------------------------------------------------
# End

