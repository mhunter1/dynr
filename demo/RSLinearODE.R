#------------------------------------------------------------------------------
# Author: Sy-Miin Chow
# Date: 2016-04-14
# Filename: RSLinearODE.R
# Purpose: An illustrative example of using dynr to fit
#   a regime-switching linear ODE
#------------------------------------------------------------------------------
#TO DO: 
#1. Suppress warnings about "In matrix(as.numeric(x), numRow, numCol) : NAs introduced by coercion"
#2. Add intercepts and allow for regime-dependent intercepts and Lambda
#3. prep.matrixDynamics should be allowed to be regime-dependent
#4. prep.tfun then cat(writeCcode(trans)$c.string) generated error
#   "Error in .local(object) : object 'formula.trans' not found"
#5. Get rid of tfun - still need to do trans to tfun conversion
#6. Remove "transformation" from dynr.cook

rm(list=ls(all=TRUE))
require(dynr)
options(scipen=999)

# ---- Read in the data ----
nT = 500; n = 10; batch = 1
thedata = read.table(paste0("./data/New2CovT",nT,"n",n,"batch",batch,"ODEsimData.txt"))
thedata$V6 <- as.numeric(thedata$V6)
colnames(thedata) = c("ID","Time","y1","y2","x1","x2")
data <- dynr.data(thedata, id="ID", time="Time",observed=paste0('y', 1:2), 
                  covariates=paste0('x', 1:2))

#---- Prepare the recipes (i.e., specifies modeling functions) ----

# Measurement (factor loadings)
meas <- prep.loadings(
  map=list(
    eta1=paste0('y', 1),
    eta2=paste0('y', 2)),
    params = NULL)

# Initial conditions on the latent state and covariance
initial <- prep.initial(
	values.inistate=c(70, 40),
	params.inistate=c("fixed", "fixed"),
	values.inicov=diag(1, 2), 
	params.inicov=diag(c("P_11", "P22"),2),
	values.regimep=c(225, 100),
	params.regimep=c("fixed", "fixed")
)

# Regime-switching function
# The RS model assumes that each element of the transition probability 
# matrix (TPM) can be expressed as a linear predictor (lp).
# LPM = 
# lp(p11) ~ 1 + x1 + x2 + ... + xn,   lp(p12) ~ 1 + x1 + x2 + ... + xn
# lp(p21) ~ 1 + x1 + x2 + ... + xn,   lp(p22) ~ 1 + x1 + x2 + ... + xn
# Here I am specifying lp(p11) and lp(p21); the remaining elements
# lp(p12) and lp(p22) are fixed at zero.

regimes <- prep.regimes(
  values=matrix(c(6,.5,-.3,rep("fixed",3),
                  rep("fixed",3),-3,-1.5,-1), 
                nrow=2, ncol=6,byrow=T), # nrow=numRegimes, ncol=numRegimes*(numCovariates+1)
  params=matrix(c("a_{11}","d_{11,1}","d_{11,2}",rep("fixed",3),
                  rep("fixed",3),"a_{21}","d_{21,1}","d_{21,2}"), 
                nrow=2, ncol=6,byrow=T), covariates=c('x1', 'x2'))

#measurement and dynamics covariances
mdcov <- prep.noise(
	values.latent=diag(1e-6, 2),
	params.latent=diag(c("fixed","fixed"), 2),
	values.observed=diag(c(10,10)),
	params.observed=diag(c("sigma^2_{epsilon,1}","sigma^2_{epsilon,2}"),2))

# dynamics
formula=list(
  list(x1~ -r1 * x1,
       x2~ -r2 * (x2 - base)),
  list(x1~ a12 * (x2 - x1),
       x2~ - a21 * (x2 - x1)))

jacob=list(
  list(x1~x1~-r1,x1~x2~0,
       x2~x1~0, x2~x2~-r2),
  list(x1~x1~-a12,x1~x2~a12,
       x2~x1~a21, x2~x2~-a21)
  )


dynm<-prep.formulaDynamics(formula=formula,jacobian=jacob,
                           startval=c(r1=.1,r2=.1,a12=.1,a21=.1),
                           isContinuousTime=TRUE)
#cat(writeCcode(dynm)$c.string)

trans<-prep.tfun(formula.trans=list(r1~exp(r1), 
                                    r2~exp(r2),
                                    a12~exp(a12),
                                    a21~exp(a21)),
                 formula.inv=list(r1~log(r1),
                                  r2~log(r2),
                                  a12~log(a12),
                                  a21~log(a21)))
#cat(writeCcode(trans)$c.string)
#------------------------------------------------------------------------------
# Cooking materials

# Model
# Put all the recipes together in a Model Specification
model <- dynr.model(dynamics=dynm, measurement=meas, 
                    noise=mdcov, initial=initial, 
                    regimes=regimes, transform=trans, 
                    outfile="RSODEmodelRecipe.c")
# View specified model in latex
printex(model)

# Estimate free parameters
#Need to remove transformation in dynr.cook
res <- dynr.cook(model, data=data,transformation=tfun,debug_flag=FALSE)

# Examine results
summary(res)

#True values should be
#c(log(.2), log(.1), log(.3), log(.2),  100, log(9.0), log(9.0), 
# -4, 8.5, -1, 1,-2, -1)


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


#------------------------------------------------------------------------------
# End


