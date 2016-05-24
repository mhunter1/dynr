#------------------------------------------------------------------------------
# Author: Lu Ou, Sy-Miin Chow
# Date: 2016-05-02
# Filename: NonlinearODE.R
# Purpose: An illustrative example of using dynr to fit
#   the predator-prey model
#------------------------------------------------------------------------------

require(dynr)

# ---- Read in the data ----
data(PPsim)
data <- dynr.data(PPsim, id="id", time="time",observed=c("x","y"))

#---- Prepare the recipes (i.e., specifies modeling functions) ----

# Measurement (factor loadings)
meas <- prep.loadings(
  map=list(
    prey="x",
    predator="y"),
  params=NULL)

# Initial conditions on the latent state and covariance
initial <- prep.initial(
	values.inistate=c(3, 1),
	params.inistate=c("fixed", "fixed"),
	values.inicov=diag(c(0.01,0.01)), 
	params.inicov=diag("fixed",2)
)

#measurement and dynamics covariances
mdcov <- prep.noise(
	values.latent=diag(0, 2),
	params.latent=diag(c("fixed","fixed"), 2),
	values.observed=diag(rep(0.3,2)),
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

printex(model, ParameterAs=model$param.names,show=FALSE,printInit=TRUE,
        outFile="./demo/NonlinearODE.tex")
tools::texi2pdf("demo/NonlinearODE.tex")
system(paste(getOption("pdfviewer"), "NonlinearODE.pdf"))

# Estimate free parameters
res <- dynr.cook(dynrModel=model)

# Examine results
summary(res)
plot(res,dynrModel=model)

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
#save(model,res,file="NonlinearODE.RData")


