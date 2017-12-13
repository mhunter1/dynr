#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2017-10-02 10:08:35
# Filename: VanDerPol.R
# Purpose: SAEM demo for trying the Van der Pol oscillator
#------------------------------------------------------------------------------

require(dynr)

# ---- Read in the data ----
# TODO fill in arguments to read in data file
#vdpData <- read.table() 
vdpData <- data.frame(id=1, time=1:1000, y1=rnorm(1000), y2=rnorm(1000), y3=rnorm(1000))
data <- dynr.data(vdpData, id="id", time="time", observed=c('y1', 'y2', 'y3'))

#---- Prepare the recipes (i.e., specifies modeling functions) ----

# Measurement (factor loadings)
meas <- prep.measurement(
  values.load=matrix(c(1, 1, 1, 0, 0, 0), 3, 2),
  params.load=matrix(c('fixed', 'lambda_21', 'lambda_31', rep('fixed', 3)), 3, 2),
  values.int=matrix(0, 3, 1),
  params.int=matrix(paste0('mu_', 1:3), 3, 1),
  obs.names = c('y1', 'y2', 'y3'),
  state.names=c('x1', 'x2'))


#TODO adjust initial condition
# Initial conditions on the latent state and covariance
initial <- prep.initial(
	values.inistate=c(3, 1),
	params.inistate=c("fixed", "fixed"),
	values.inicov=diag(c(0.01,0.01)), 
	params.inicov=diag("fixed",2)
)

# TODO adjust noise
#measurement and dynamics covariances
mdcov <- prep.noise(
	values.latent=diag(0, 2),
	params.latent=diag(c("fixed","fixed"), 2),
	values.observed=diag(rep(0.3,2)),
	params.observed=diag(c("var_1","var_2"),2)
)

# dynamics
#formula=list(list(prey~ a*prey - b*prey*predator,
#             predator~ -c*predator + d*prey*predator))
formula <- list(
	x1 ~ x2,
	x2 ~ -omega*x1 + zeta*(1-x1^2)*x2 ) # omega = (2*pi/gamma)^2
# TODO put in a numeric constant for omega
dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(zeta=2.1),
                           isContinuousTime=TRUE)

#------------------------------------------------------------------------------
# Cooking materials

# Put all the recipes together in a Model Specification
model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data,
                    outfile="VanDerPol.c", armadillo=TRUE)
# Note: with armadillo=TRUE this will write armadillo code using
#  the writeArmadilloCode methods.  As of 2017-10-30 10:29:47 this
#  will do almost nothing because the armadillo methods all write "".

printex(model, ParameterAs=model$param.names,show=FALSE,printInit=TRUE,
        outFile="VanDerPol.tex")
#tools::texi2pdf("VanDerPol.tex")
#system(paste(getOption("pdfviewer"), "VanDerPol.pdf"))

# Look at VanDerPol.c to see if the dynamics and jacobian functions are written in Armadillo.
# In the branch armadillo, modify R/dynrRecipe
#setMethod("writeCcode", "dynrDynamicsFormula",
# to write armadillo code instead of C code.

# Later Michael Hunter will change this to a writeArmadilloCode method
#  and connect it to dynr.model with appropriate arguments to run either the 
#  Ccode method or the ArmadilloCode method.

#------------------------------------------------------------------------------
# End
