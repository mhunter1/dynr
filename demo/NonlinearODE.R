#------------------------------------------------------------------------------
# Author: Lu Ou, Sy-Miin Chow
# Date: 2016-05-02
# Filename: PPmodel.R
# Purpose: An illustrative example of using dynr to fit
#   the predator-prey model
#------------------------------------------------------------------------------


#rm(list=ls(all=TRUE))
require(dynr)
options(scipen=999)

# ---- Read in the data ----
data(PPsim)
data <- dynr.data(PPsim, id="id", time="time",observed=c("x","y"))

#---- Prepare the recipes (i.e., specifies modeling functions) ----

# Measurement (factor loadings)
meas <- prep.loadings(
  map=list(
    eta1="x",
    eta2="y"),
  params=NULL)

# Initial conditions on the latent state and covariance
initial <- prep.initial(
	values.inistate=c(3, 2),
	params.inistate=c("fixed", "fixed"),
	values.inicov=diag(c(0.05,0.05)), 
	params.inicov=diag("fixed",2)
)

#measurement and dynamics covariances
mdcov <- prep.noise(
	values.latent=diag(0, 2),
	params.latent=diag(c("fixed","fixed"), 2),
	values.observed=diag(rep(0.02,2)),
	params.observed=diag(c("var_1","var_2"),2)
)

# dynamics
formula=list(list(prey~ r1*prey - a12*prey*predator,
             predator~ -r2*predator + a21*prey*predator))
dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(r1=2.1, r2=0.8, a12 = 1.9, a21 = 1.1),
                           isContinuousTime=TRUE)
#constraints
trans<-prep.tfun(formula.trans=list(r1~exp(r1), 
                                    r2~exp(r2),
                                    a12~exp(a12),
                                    a21~exp(a21)),
                 formula.inv=list(r1~log(r1),
                                  r2~log(r2),
                                  a12~log(a12),
                                  a21~log(a21)))

#------------------------------------------------------------------------------
# Cooking materials

# Put all the recipes together in a Model Specification
model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial,
                    transform=trans, data=data,
                    outfile="PPmodelRecipe.c")

# Estimate free parameters
res <- dynr.cook(dynrModel=model)

# Examine results
summary(res)

p1 = dynr.ggplot(res, data.dynr=data, states=c(1:2), 
            names.regime=c("Exploration","Proximity-seeking"),
            names.state=c("Mom","Infant"),
            title="Results from RS-linear ODE model", numSubjDemo=2,idtoPlot=c(1,2),
            shape.values = c(1,2),
            text=element_text(size=16))

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


