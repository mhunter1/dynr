# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Author: Hui-Ju Hung
# Date: 2018-05-15
# Filename: vanderpol_dynr.R
# Purpose: Fit van der Pol oscillator model (cook workable)
#          For Sy-Miin to get good initial estimates for SAEM. 
# Note: Workable for developer dynr on master branch 0.1.14-17
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library('dynr')
library('dplyr')
library('magrittr')

nPeople = 200
nTimes = 300
vdpData <- read.csv("../data/TrueInit_Y14.txt", header=FALSE)
colnames(vdpData) <- c('batch', 'kk', 'trueInit', 'time', 'y1','y2','y3', 'co1', 'co2')
vdpData$id <- rep(1:nPeople, each=nTimes)
data <- dynr.data(vdpData, id="id", time="time",
                  observed=c('y1','y2','y3'),
                  covariates=c("co1","co2"))

# ---- Get starting values for parameters in Beta ----

mdcov <- prep.noise(
    values.latent=diag(0, 2),
    params.latent=diag(c("fixed","fixed"), 2),
    values.observed=diag(rep(0.3,3)),
    params.observed=diag(c("var_1","var_2","var_3"),3)
)

meas <- prep.measurement(
    values.load = matrix(c(1, 0), 3, 2,byrow=TRUE),
    params.load = matrix(c(
                  "fixed","fixed",
                  "lambda1","fixed",
                  "lambda2","fixed"),3,2,byrow=TRUE),
    obs.names = c('y1', 'y2', 'y3'),
    state.names = c('x1', 'x2'))


initial <- prep.initial(
    values.inistate=c(1, 1),
    params.inistate=c("mu_x1", "mu_x2"),
    values.inicov=matrix(c(1,.3,
                           .3,1),ncol=2,byrow=TRUE), 
    params.inicov=matrix(c("sigma2_bx1","sigma_bx1x2",
                           "sigma_bx1x2","sigma2_bx2"),2,2,byrow=TRUE)
)
if (length(unlist(initial$params.inistate[!initial$params.inistate== "fixed"]))>0)
    initial$params.inistate

formula=
    list(x1 ~ x2,
         x2 ~ -61.68503 * x1 + (zeta0  + co1 * zeta1 + co2 * zeta2) * (1 - x1^2) * x2
    )


dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(zeta0=2,
                                      zeta1=.3,
                                      zeta2=.3),
                           isContinuousTime=TRUE)

model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data, armadillo=FALSE,
                    outfile="VanDerPol_.c")
#model$xstart<-c(3,0.5,0.3,0,0,0.5,0.5,0.5)
fitted_model <- dynr.cook(model, optimization_flag = TRUE, hessian_flag = FALSE)


# ---- Get starting values for random effects ---- 
# Take output from this initial model, particularly estimates for initial condition parameters
# (if applicable), put random effects as additional latent variables in dyn model. Also using
# parameter estimates from previously cooked model, assuming no random effects, as the starting
# values for the SAEM estimation

coef(fitted_model)
#Hui-Ju: you can use this to put the final estimated parameter values back into the starting
#values in model.
#Bring up in next week's dynr meeting about how to call this function. 
#Right now this function is hidden. And it is not working correctly with the LDL
#source('../R/dynrCook.R')
#model.step1 <<-PopBackModel(model, fitted_model@transformed.parameters)

mdcov2 <- prep.noise(
  values.latent=diag(0, 3),
  params.latent=diag(c("fixed","fixed","fixed"), 3),
  values.observed=diag(coef(fitted_model)[c("var_1","var_2","var_3")]),
  params.observed=diag(c("var_1","var_2","var_3"),3)
)

meas2 <- prep.measurement(
  values.load = matrix(c(1, 0, 0,
                         1, 0, 0,
                         1, 0, 0), 3, 3,byrow=TRUE),
  params.load = matrix(c(
    "fixed","fixed","fixed",
    "lambda1","fixed","fixed",
    "lambda2","fixed","fixed"),3,3,byrow=TRUE),
  obs.names = c('y1', 'y2', 'y3'),
  state.names = c('x1', 'x2',"bzeta"))


# Try to populate estimates from initial in the previous chunk of code into initial2
if (length(unlist(initial$params.inistate[!initial$params.inistate== "fixed"]))>0)
  initial$params.inistate

v1 = matrix(coef(fitted_model)[c("sigma2_bx1","sigma_bx1x2",
                                 "sigma_bx1x2","sigma2_bx1")],2,2)

initial2 <- prep.initial(
  values.inistate=c(coef(fitted_model)[c("mu_x1","mu_x2")],0),
  params.inistate=c("mu_x1", "mu_x2",0),
  values.inicov=matrix(c(1,.3,0,
                         .3,1,0,
                         0,0,1),ncol=3,byrow=TRUE), 
  params.inicov=matrix(c("sigma2_bx1","sigma_bx1x2","fixed",
                         "sigma_bx1x2","sigma2_bx2","fixed",
                         "fixed","fixed","sigma2_bzeta"),3,3,byrow=TRUE)
)

if (length(unlist(initial2$params.inistate[!initial2$params.inistate== "fixed"]))>0)
  initial2$params.inistate

formula2=
  list(x1 ~ x2,
       x2 ~ -61.68503 * x1 + (zeta0  + co1 * zeta1 + co2 * zeta2 + bzeta) * (1 - x1^2) * x2,
       bzeta ~ 0
  )


dynm2<-prep.formulaDynamics(formula=formula2,
                           startval=c(zeta0=2,
                                      zeta1=.3,
                                      zeta2=.3),
                           isContinuousTime=TRUE)

model2 <- dynr.model(dynamics=dynm2, measurement=meas2,
                    noise=mdcov2, initial=initial2, data=data, armadillo=FALSE,
                    outfile="VanDerPol2_.c")


# ---- This portion sets the parameter starting values to those estimated in part 1
test = function(x,y){
  pos = grep(x,
             y)
  return(pos)
}

pos=unlist(lapply(names(coef(fitted_model)[1:5]),test,names(model2$xstart)))
model2@xstart[pos] = coef(fitted_model)[1:5] 
#Need to resolve the LDL issues. Cannot set parameter values in cov matrices correctly
#Also still cannot use model2$xstart

# ---- End of parameter value setting

fitted_model2 <- dynr.cook(model2, optimization_flag = TRUE, 
                           hessian_flag = FALSE, verbose=FALSE, debug_flag=FALSE)


locc=plyr::ddply(data.frame(id=data$id,time=data$time,index=1:length(data$time)), 
           .(id), function(x){x$index[which(x$time==max(x$time))]})[,2]

#Get estimates for bzeta from fitted_model2 and use them as estimates for b
bEst = fitted_model2@eta_smooth_final[3,locc] #Use these as the starting values for InfDS.b
coef.Est = coef(fitted_model2) #Use these as the starting values for Beta


