#------------------------------------------------------------------------------
# Author: Lu Ou, Sy-Miin Chow
# Date: 2016-05-05
# Filename: RS-PPmodel.R
# Purpose: An illustrative example of using dynr to fit
#   a regime-switching predator-prey model
#------------------------------------------------------------------------------

require(dynr)

# ---- Read in the data ----
data(RSPPsim)
data <- dynr.data(RSPPsim, id="id", time="time",observed=c("x","y"),covariate="cond")

#---- Prepare the recipes (i.e., specifies modeling functions) ----

# Measurement (factor loadings)
meas <- prep.loadings(
  map=list(
    prey=c("x"),
    predator=c("y")),
  params=NULL)

# Initial conditions on the latent state and covariance
initial <- prep.initial(
  values.inistate=c(3, 1),
  params.inistate=c("fixed", "fixed"),
  values.inicov=diag(c(0.01,0.01)), 
  params.inicov=diag("fixed",2),
	values.regimep=c(.5, .5),
	params.regimep=c("fixed", "fixed")
)


# Regime-switching function
# The RS model assumes that each element of the transition probability 
# matrix (TPM) can be expressed as a linear predictor (lp).
# LPM = 
# lp(p11) ~ 1 + x1 + x2 + ... + xn,   lp(p12) ~ 1 + x1 + x2 + ... + xn
# lp(p21) ~ 1 + x1 + x2 + ... + xn,   lp(p22) ~ 1 + x1 + x2 + ... + xn
# Here I am specifying lp(p11) and lp(p22); the remaining elements
# lp(p11) and lp(p21) are fixed at zero.

regimes <- prep.regimes(
  values=matrix(c(0,0,-1,1,
                  0,0,-1,1),
                nrow=2, ncol=4,byrow=T), # nrow=numRegimes, ncol=numRegimes*(numCovariates+1)
  params=matrix(c("fixed","fixed","int","slp",
                  "fixed","fixed","int","slp"), 
                nrow=2, ncol=4,byrow=T), 
  covariates="cond")

#measurement and dynamics covariances
mdcov <- prep.noise(
  values.latent=diag(0, 2),
  params.latent=diag(c("fixed","fixed"), 2),
  values.observed=diag(rep(0.5,2)),
  params.observed=diag(rep("var",2),2)
)

# dynamics
formula=list(
  list(prey~ r1*prey - a12*prey*predator,
       predator~ -r2*predator + a21*prey*predator),
  list(prey~ r1*prey - a1*prey^2 - a12*prey*predator,
       predator~ a2*predator - r2*predator^2 + a21*prey*predator ))

dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(r1=2.1, r2=0.8, a12 = 1.9, a21 = 1.1,
                                      a1 =1, a2 = 1),
                           isContinuousTime=TRUE)

#constraints
trans<-prep.tfun(formula.trans=list(r1~exp(r1), 
                                    r2~exp(r2),
                                    a12~exp(a12),
                                    a21~exp(a21),
                                    a1~exp(a1),
                                    a2~exp(a2)),
                 formula.inv=list(r1~log(r1),
                                  r2~log(r2),
                                  a12~log(a12),
                                  a21~log(a21),
                                  a1~log(a1),
                                  a2~log(a2))
                 )

#------------------------------------------------------------------------------
# Cooking materials

# Put all the recipes together in a Model Specification
model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial,
                    regimes=regimes, transform=trans,
                    data=data,
                    outfile="RSPPmodelRecipe.c")

printex(model, ParameterAs = model@param.names, printInit=TRUE, printRS=TRUE,
        outFile="demo/RSNonlinearODE.tex")
tools::texi2pdf("demo/RSNonlinearODE.tex")
system(paste(getOption("pdfviewer"), "RSNonlinearODE.pdf"))

model@ub[model@param.names%in%c("int","slp")]<-c(0,10)
model@lb[model@param.names%in%c("int","slp")]<-c(-10,0)
# Estimate free parameters
res <- dynr.cook(model)

# Examine results
summary(res)

#plotFormula(model, ParameterAs=signif(res@transformed.parameters,2)) 

p1 = dynr.ggplot(res, data.dynr=data, states=c(1:2), 
            names.regime=c("Free","Constrained"),
            names.state=c("Prey","Predator"),
            title="Results from RS-nonlinear ODE model", numSubjDemo=2,
            shape.values = c(1,2),
            text=element_text(size=16))

#plot(res,data.dynr = data,model)
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


