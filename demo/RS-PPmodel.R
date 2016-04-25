#------------------------------------------------------------------------------
# Author: Sy-Miin Chow
# Date: 2016-04-14
# Filename: RSLinearODE.R
# Purpose: An illustrative example of using dynr to fit
#   a regime-switching linear ODE
#------------------------------------------------------------------------------


rm(list=ls(all=TRUE))
require(dynr)
options(scipen=999)

# ---- Read in the data ----
thedata = read.table(paste0("./data/PPsimData.txt"))
thedata$V6 <- as.numeric(thedata$V10)
colnames(thedata) = c("ID","Time",paste0("y",1:6),"x1","x2")
data <- dynr.data(thedata, id="ID", time="Time",observed=paste0('y', 1:6), 
                  covariates=paste0('x', 1:2))

#---- Prepare the recipes (i.e., specifies modeling functions) ----

# Measurement (factor loadings)
meas <- prep.loadings(
  map=list(
    eta1=paste0('y', 1:3),
    eta2=paste0('y', 4:6)),
  params=c("lambda_{21}","lambda_{31}","lambda_{52}","lambda_{62}"))

#cat(writeCcode(meas)$c.string) #Can't write C code yet

# Initial conditions on the latent state and covariance
initial <- prep.initial(
	values.inistate=c(10, 25),
	params.inistate=c("fixed", "fixed"),
	values.inicov=diag(c(9,4)), 
	params.inicov=diag("fixed",2),
	values.regimep=c(.5, 0),
	params.regimep=c("p0", "fixed")
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
  values=matrix(c(6,.5,-.3,rep(0,3),
                  -3,-1.5,-1,rep(0,3)), 
                nrow=2, ncol=6,byrow=T), # nrow=numRegimes, ncol=numRegimes*(numCovariates+1)
  params=matrix(c("a_{11}","d_{11,1}","d_{11,2}",rep("fixed",3),
                  "a_{21}","d_{21,1}","d_{21,2}",rep("fixed",3)), 
                nrow=2, ncol=6,byrow=T), covariates=c('x1', 'x2'))

#measurement and dynamics covariances
mdcov <- prep.noise(
	values.latent=diag(1e-6, 2),
	params.latent=diag(c("fixed","fixed"), 2),
	values.observed=diag(rep(5,6)),
	params.observed=diag(c(paste0("sigma_e",1:6)),6)
)

# dynamics
formula=list(
  list(prey~ r1*prey - a12*prey*predator - a11*prey^2,
       predator~ -r2*predator + a21*prey*predator - a22*predator^2),
  list(prey~ r1*prey - a12*prey*predator,
       predator~ -r2*predator + a21*prey*predator))

jacob=list(
  list(prey~prey~r1 - a12*predator - 2*a11*prey,
       prey~predator~ -a12*prey,
       predator~prey~ a21*predator,
       predator~predator~-r2 + a21*prey - 2*a22*predator),
  list(prey~prey~r1 - a12*predator,
       prey~predator~ -a12*prey,
       predator~prey~ a21*predator,
       predator~predator~-r2 + a21*prey))


#Starting values are on constrained scale
dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(r1=3, r2=2, a12 = 1, a21 = 1, 
                                      a11 =.5, a22 = .5),
                           isContinuousTime=TRUE,jacobian=jacob)

#cat(writeCcode(dynm)$c.string)

trans<-prep.tfun(formula.trans=list(r1~exp(r1), 
                                    r2~exp(r2),
                                    a12~exp(a12),
                                    a21~exp(a21),
                                    a11~exp(a11),
                                    a22~exp(a22)),
                 formula.inv=list(r1~log(r1),
                                  r2~log(r2),
                                  a12~log(a12),
                                  a21~log(a21),
                                  a11~log(a11),
                                  a22~log(a22)))

#cat(writeCcode(trans)$c.string)
#------------------------------------------------------------------------------
# Cooking materials

# Model
# Put all the recipes together in a Model Specification
model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, 
                    regimes=regimes, transform=trans,
                    outfile="RSPPmodelRecipe.c")

#cat(writeCcode(model$dynamics)$c.string)

# View specified model in latex
#printex(model)

# Estimate free parameters
#Need to remove transformation in dynr.cook
res <- dynr.cook(model, data=data,debug_flag=FALSE)

# Examine results
summary(res)

#True values should be
#c(log(.2), log(.1), log(.3), log(.2),  100, log(9.0), log(9.0), 
# 4.5, -4, 1,-1,-1, -2)




p1 = dynr.ggplot(res, data.dynr=data, states=c(1:2), 
            names.regime=c("Exploration","Proximity-seeking"),
            names.state=c("Mom","Infant"),
            title="Results from RS-linear ODE model", numSubjDemo=2,idtoPlot=c(1,2),
            shape.values = c(1,2),
            text=element_text(size=16))

print(p1)
plot(res,data.dynr = data,model)

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


