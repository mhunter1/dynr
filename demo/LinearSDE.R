#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2016-02-25
# Filename: LinearSDE.R
# Purpose: Translate LinearSDE.R into the user-spec for the same model
#  in dynr.
#------------------------------------------------------------------------------



require(dynr)
require(KFAS)

meas <- dynr.matrixLoadings(
	values=matrix(c(1,0), 1, 2),
	params=matrix(0, 1, 2))

ecov <- dynr.matrixErrorCov(
	values.latent=diag(c(0.00001,1)), params.latent=diag(c(0, 3)),
	values.observed=diag(1.5,1), params.observed=diag(4, 1))
ecov$c.string
ecov$startval

initial<-dynr.initial(c(0,1), c(5,0), diag(1,2), diag(0,2))
writeLines(initial$c.string)
initial$startval

dynamics <- dynr.linearDynamics(
	params.dyn=matrix(c(0, 1, 0, 2), 2, 2),
	values.dyn=matrix(c(0, 1, 1, 1), 2, 2),
	time="contin")

regimes <- dynr.regimes()

# Proto-example of cooking
fname <- "./demo/CookedLinearSDE.c"
dynr.cook(file=fname, meas, ecov$c.string, initial$c.string, dynamics, regimes)

#--------------------------------------
# Model Specification
data <- dynr.data(cbind(id=rep(1,100),t(ty), times=tT[,-1]), id="id", time="times", observed="y1")
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
              infile=fname, 
              outfile="./demo/LinearSDE2", 
              verbose=TRUE,
              compileLib=TRUE
)

res <- dynr.run(model, data)

summary(res)




