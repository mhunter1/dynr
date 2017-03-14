#------------------------------------------------------------------------------
# Author: Lu Ou
# Date: 2016-05-24
# Filename: RSNonlinearDiscrete.R
# Purpose: An illustrative example of using dynr to fit
#   a regime-switching nonlinear dynamic factor analysis (discrete-time) model 
#   with nonlinear vector autoregressive (VAR) relation at the factor level
#------------------------------------------------------------------------------

#rm(list=ls(all=TRUE))
require(dynr)

#---- (0) Read in data ----
# Data
data(NonlinearDFAsim)
data <- dynr.data(NonlinearDFAsim, id="id", time="time",observed=colnames(NonlinearDFAsim)[c(3:8)])

#--- (1) Prepare the recipes ----

#---- (1a) Measurement (factor loadings) -----
meas <- prep.measurement(
  values.load=matrix(c(1, .8, .8, rep(0, 3),
                       rep(0, 3), 1, .8, .8), ncol=2),
  params.load=matrix(c("fixed", "lambda_21", "lambda_31", rep("fixed",3),
                       rep("fixed",3), "fixed", "lambda_52","lambda_62"), ncol=2),
  state.names=c('PE', 'NE'))

# alternatively, use prep.loadings
# meas <- prep.loadings(
#   map=list(
#     PE=paste0('y', 1:3),
#     NE=paste0('y', 4:6)),
#   params=c("lambda_21","lambda_31","lambda_52","lambda_62"))


#---- (1b) Initial conditions ----
initial <- prep.initial(
	values.inistate=c(0, 0),
	params.inistate=c("fixed", "fixed"),
	values.inicov=diag(1, 2), 
	params.inicov=diag("fixed", 2),
	values.regimep=c(1.3865, 0),
	params.regimep=c("fixed", "fixed")
)

#---- (1c) Regime-switching function ----
regimes <- prep.regimes(
	values=matrix(c(.9, 0, 0, .9), 2, 2), #nrow=numRegimes, ncol=numRegimes*(numCovariates+1)
	params=matrix(c("p11", 0, 0, "p22"), 2, 2))
# Self-transition (diagonals) are estimated
# Off-diagonal elements are fixed for identification purposes


#---- (1d) measurement and dynamics noise covariance structure ----
mdcov <- prep.noise(
	values.latent=diag(0.3, 2),
	params.latent=diag(paste0("zeta_",1:2), 2),
	values.observed=diag(0.1, 6),
	params.observed=diag(paste0("epsilon_",1:6), 6))

#---- (1e) dynamics ----
formula=list(
  list(PE~a1*PE,
       NE~a2*NE),
  list(PE~a1*PE+c12*(exp(abs(NE)))/(1+exp(abs(NE)))*NE,
       NE~a2*NE+c21*(exp(abs(PE)))/(1+exp(abs(PE)))*PE) 
)

#For nonlinear functions you need to specify the jacobian
#matrix containing the differentiation of each dynamic function
#above with respect to each latent variable (e.g., PE and NE
#in this case)
jacob=list(
  list(PE~PE~a1,
       NE~NE~a2),
  list(PE~PE~a1,
       PE~NE~c12*(exp(abs(NE))/(exp(abs(NE))+1)+NE*sign(NE)*exp(abs(NE))/(1+exp(abs(NE))^2)),
       NE~NE~a2,
       NE~PE~c21*(exp(abs(PE))/(exp(abs(PE))+1)+PE*sign(PE)*exp(abs(PE))/(1+exp(abs(PE))^2))))

dynm<-prep.formulaDynamics(formula=formula,startval=c(a1=.3,a2=.4,c12=-.5,c21=-.5),isContinuousTime=FALSE,jacobian=jacob)

trans<-prep.tfun(formula.trans=list(p11~exp(p11)/(1+exp(p11)), p22~exp(p22)/(1+exp(p22))), formula.inv=list(p11~log(p11/(1-p11)),p22~log(p22/(1-p22))), transCcode=FALSE)

#----(2) Put together the model and cook it! ----

# Put all the recipes together in a Model Specification
model <- dynr.model(dynamics=dynm, measurement=meas, noise=mdcov, 
                    initial=initial, regimes=regimes, transform=trans, 
                    data=data, 
                    outfile="RSNonlinearDiscrete")

printex(model,ParameterAs=model$param.names,printInit=TRUE, printRS=TRUE,
        outFile="RSNonlinearDiscrete.tex")
#tools::texi2pdf("RSNonlinearDiscrete.tex")
#system(paste(getOption("pdfviewer"), "RSNonlinearDiscrete.pdf"))

res <- dynr.cook(model)

#---- Examine and "serve" the results
summary(res)

# get the estimated parameters from a cooked model/data combo
coef(res)

# get the log likelihood, AIC, and BIC from a cooked model/data combo
logLik(res)
AIC(res)
BIC(res)


# compare true parameters to estimated ones
truepar <- c(
  .2, .25, -.6, -.8,
  1.2, 1.2, 1.1, .95, 
  c(.35, .3),
  c(.28, .10, .12, .13, .12, .11),
  0.98,0.85)
data.frame(name=res@param.names , true=truepar, estim=coef(res))

dynr.ggplot(res, dynrModel=model, style=1,
            names.regime=c("Decoupled (linear)","Coupled (nonlinear)"),
            title="Results from RS Nonlinear DFA model", numSubjDemo=1,idtoPlot=c(1),
            shape.values = c(1,2),
            text=element_text(size=16))
#ggsave("RSNonlinearDiscreteggPlot1.pdf")
dynr.ggplot(res, dynrModel=model, style=2,
            names.observed=c("y1","y4"),
            names.regime=c("Decoupled (linear)","Coupled (nonlinear)"),
            title="Results from RS Nonlinear DFA model", numSubjDemo=1,idtoPlot=c(1),
            text=element_text(size=16))
#ggsave("RSNonlinearDiscreteggPlot2.pdf")

plotFormula(dynrModel=model, ParameterAs=model@param.names, printDyn=TRUE, printMeas=TRUE)
plot(res, dynrModel=model, style = 1)
plot(res, dynrModel=model, style = 2, names.observed=c("y1","y4"))
#---- Done ----
#save(model,res,file="RSNonlinearDiscrete.RData")

