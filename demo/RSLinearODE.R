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

#cat(writeCcode(meas)$c.string) #Can't write C code yet

# Initial conditions on the latent state and covariance
initial <- prep.initial(
	values.inistate=c(70, 40),
	params.inistate=c("fixed", "fixed"),
	values.inicov=diag(c(225,100)), 
	params.inicov=diag("fixed",2),
	values.regimep=c(1, 0),
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
	values.observed=diag(c(10,10)),
	params.observed=diag(c("sigmasq_e1","sigmasq_e2"),2))

# dynamics
formula=list(
  list(eta1~ -r1 * eta1,
       eta2~ -r2 * (eta2 - base)),
  list(eta1~ a12 * (eta2 - eta1),
       eta2~ - a21 * (eta2 - eta1)))

#jacob=list(
#  list(x1~x1~-r1,x1~x2~0,
#       x2~x1~0, x2~x2~-r2),
#  list(x1~x1~-a12,x1~x2~a12,
#       x2~x1~a21, x2~x2~-a21)
#  )


dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(r1=.1,r2=.1,a12=.1,a21=.1,base=95),
                           isContinuousTime=TRUE) #,jacobian=jacob

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
#Can extract model$xstart to see the vector of starting values
model$xstart; model$param.names

#Use the `@' sign to set upper and lower boundaries for the parameters
model@ub=c(rep(1.5, 4), 200, 5, 5, rep(30, 6))
model@lb=c(rep(-20, 4), 50, -10, -10, rep(-30, 6))
model@dynamics

#cat(writeCcode(model$dynamics)$c.string)

# View specified model in latex
printex(model$dynamics)
printex(dynm)

# Estimate free parameters
#Need to remove transformation in dynr.cook
res <- dynr.cook(model, data=data)

# Examine results
summary(res)

#True values should be
#c(log(.2), log(.1), log(.3), log(.2),  100, log(9.0), log(9.0), 
# 4.5, -4, 1,-1,-1, -2)

printex(model@dynm)


p1 = dynr.ggplot(res, data.dynr=data, states=c(1:2), 
            names.regime=c("Exploration","Proximity-seeking"),
            names.state=c("Mom","Infant"),
            title="Results from RS-linear ODE model", numSubjDemo=2,idtoPlot=c(1,2),
            shape.values = c(1,2),
            text=element_text(size=16))

print(p1)
plot(res,data.dynr = data,model)


#par(cex.axis=1.3,cex.lab=1.5,mgp=c(2,.8,0))
#plot(thedata$y1, res@residuals,type="n",
#     xlab="y1(t)/dt", ylab=expression(paste("Residuals y1(t)")))
#grid(5, 7, lwd = 2,lty="solid") # grid lines
#points(thedata$y1, res@residuals)
#abline(0,coef(lm(coef(res@residuals~y1))[2],col="red",lty=2,lwd=3)
##lines(thedata$y1, coef(gall)["dx"]*dxall$dx,lwd=2,col="red",lty=2)
#lines(loess.smooth(thedata$y1, res@residuals, 
#                   span = 2/3, degree = 2),
#      col="green",lwd=3)
#
##loessLine(dxall$dx, residuals(gall)+coef(gall)["dx"]*dxall$dx, 
##          col="green",smoother.args=list(),log.x=FALSE,log.y=FALSE)
#legend("topright",c("Linear least sqs","Loess"),lty=c(2,1),lwd=c(3,1),col=c("red","green"),
#       bty="n",cex=1.3)
#title(main="Component-Plus-Residual Plot: dx(t)/dt",cex.main=1.5) 

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


