#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Author: Sy-Miin Chow
# Date: 2016-04-14
# Filename: RSLinearODE.R
# Purpose: An illustrative example of using dynr to fit
#   a regime-switching linear ODE
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#rm(list=ls(all=TRUE))
require(dynr)
options(scipen=999,digits=2)

# ---- Read in the data ----
nT = 500; n = 10; batch = 1
thedata = read.table(paste0("./data/New2CovT",nT,"n",n,"batch",batch,"ODEsimData.txt"))
thedata$V6 <- as.numeric(thedata$V6)
colnames(thedata) = c("ID","Time","y1","y2","x1","x2")
data <- dynr.data(thedata, id="ID", time="Time",observed=paste0('y', 1:2), 
                  covariates=paste0('x', 1:2))

#---- Prepare the recipes (i.e., specifies modeling functions) ----

# Measurement (factor loadings)
meas <- prep.measurement(
  values.load=diag(c(1,1)),
  obs.names = c('y1','y2'),
  state.names=c('Mom','Baby'))

#cat(writeCcode(meas)$c.string) #Can't write C code yet

#measequ<-paste0(paste0("$\\frac{",paste(dyn[[2]]$left[1],"}{dt}",collapse="\\\\"), "$"),
#               " = ",paste0("$",paste(dyn[[2]]$right[1],collapse="\\\\"), "$"))
#plot(TeX(dynequ))

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
  params=matrix(c("a_11","d_111","d_112",rep("fixed",3),
                  "a_21","d_211","d_212",rep("fixed",3)), 
                nrow=2, ncol=6,byrow=T), covariates=c('x1', 'x2'))

#measurement and dynamic noise covariance structurres
mdcov <- prep.noise(
	values.latent=diag(0, 2),
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
#  list(eta1~eta1~-r1,x1~x2~0,
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
#---- Cooking it up! ----

# Model
# Put all the recipes together in a Model Specification
model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, 
                    regimes=regimes, transform=trans,
                    data=data,
                    outfile="RSODEmodelRecipe.c")
#Extract parameter names to set ub and lb (optional)
model$param.names

#Use the `$' sign to set upper and lower boundaries for the parameters
model$ub=c(rep(4.5, 4), 200, 150, 150, rep(30, 6))
model$lb=c(rep(2.1e-9, 4), 50, 4.5e-5, 4.5e-5, rep(-30, 6))

#Check model by printing out LaTeX code
printex(model,show=FALSE,printInit = TRUE, printProb=TRUE,outFile="./demo/RSLinearODE.tex")

TeXed <- printFormula(model,namestoPop = model$param.names)
#Also can pop= signif(res@transformed.parameters,digits=2))

#See dyn example for everything gotten substituted away
p3 <-plotFormula(TeXed,model,toPlot="dyn")
p3 <-plotFormula(TeXed,model,toPlot="meas")
plotFormula(TeXed,model,toPlot="both")

# Estimate free parameters
res <- dynr.cook(model)

#---- Serve it! ----
p1 = dynr.ggplot(res, data.dynr=data, states=c(1:2), 
                 names.regime=c("Exploration",
                                "Proximity-seeking"),
                 names.state=c("Mom","Infant"),
                 title="Results from RS-linear ODE model", 
                 numSubjDemo=2,idtoPlot=c(1,2),
                 shape.values = c(1,2),
                 text=element_text(size=16))

print(p1)
plot(res,data.dynr = data,model=model)
summary(res)

# get the estimated parameters from a cooked model/data combo
coef(res)
#True values should be
#c(log(.2), log(.1), log(.3), log(.2),  100, log(9.0), log(9.0), 
# 4.5, -4, 1,-1,-1, -2)

# get the log likelihood, AIC, and BIC from a cooked model/data combo
logLik(res)
AIC(res)
BIC(res)

#---- End of demo ----


