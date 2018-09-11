#---------------------------------------------------------------------
# Author: Dongjun You
# Date: 2018-05-15
# Filename: OutlierDetection.R
# Purpose: An illustrative example of using dynr.taste to detect
#   outliers from a state space model
#---------------------------------------------------------------------
# Load packages
library(dynr)

# The data were generated from a discrete time state-space model with six observed variables and two latent variables
# 100 participants, 60 time-points for each participant
# For each participant, 3 innovative outliers and 6 additive outliers 
data("Outliers")

# The true parameters used to generate the data
true_params <- c(.6, -.2, -.2, .5,# beta
                 .9, .8, .9, .8,# lambda
                 .3, -.1, .3,# psi
                 .2, .2, .2, .2, .2, .2)# theta

data_shk <- dynr.data(Outliers, id='id', time='time',
                      observed=c('V1','V2','V3','V4','V5','V6'))

# measurement
# this is the factor loadings matrix, Lambda in SEM notation
meas_shk <- prep.measurement(
  values.load=matrix(c(1.0, 0.0,
                       0.9, 0.0,
                       0.8, 0.0,
                       0.0, 1.0,
                       0.0, 0.9,
                       0.0, 0.8), ncol=2, byrow=TRUE),
  params.load=matrix(c('fixed','fixed',
                       'l_21','fixed',
                       'l_31','fixed',
                       'fixed','fixed',
                       'fixed','l_52',
                       'fixed','l_62'), ncol=2, byrow=TRUE),
  state.names=c('eta_1','eta_2'),
  obs.names=c('V1','V2','V3','V4','V5','V6') )

# observation and dynamic noise components
# the latent noise is the dynamic noise, Psi in SEM notation
# the observed noise is the measurement noise, Theta in SEM
nois_shk <- prep.noise(
  values.latent=matrix(c(0.3, -0.1,
                         -0.1, 0.3), ncol=2, byrow=TRUE),
  params.latent=matrix(c('psi_11','psi_12',
                         'psi_12','psi_22'), ncol=2, byrow=TRUE),
  values.observed=diag(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2), ncol=6, nrow=6),
  params.observed=diag( paste0('e_', 1:6), 6) )

# initial covariances and latent state values
init_shk <- prep.initial(
  values.inistate=c(0,0),
  params.inistate=c('mu_1','mu_2'),
  values.inicov=matrix(c(0.3, -0.1,
                        -0.1,  0.3), ncol=2, byrow=TRUE),
  params.inicov=matrix(c('c_11','c_12',
                         'c_12','c_22'), ncol=2, byrow=TRUE) )

# define the differential equation
dynm_shk <- prep.matrixDynamics(
  values.dyn=matrix(c(0.8, -0.2,
                      -0.2,  0.7), ncol=2, byrow=TRUE),
  params.dyn=matrix(c('b_11','b_12',
                      'b_21','b_22'), ncol=2, byrow=TRUE),
  isContinuousTime=FALSE)

formula_shk <- list(list(eta_1 ~ b_11*eta_1 - b_12*eta_2,
                         eta_2 ~ b_21*eta_1 + b_22*eta_2))
dynm_shk2 <- prep.formulaDynamics(formula=formula_shk,
                                  startval=c(b_11=.6, b_12=-.2, b_21=-.2, b_22=.6),
                                 isContinuousTime=FALSE)
# Prepare for cooking
# put all the recipes together
model_shk <- dynr.model(dynamics=dynm_shk, measurement=meas_shk,
                        noise=nois_shk, initial=init_shk,
                        data=data_shk, outfile=paste0("model_shk.c"))
model_shk2 <- dynr.model(dynamics=dynm_shk2, measurement=meas_shk,
                         noise=nois_shk, initial=init_shk,
                         data=data_shk, outfile=paste0("model_shk2.c"))

# Estimate free parameters
# SHOULD cook with 'debug_flag=TRUE'
cook_shk <- dynr.cook(model_shk, debug_flag=TRUE)
cook_shk2 <- dynr.cook(model_shk2, debug_flag=TRUE)
cocook <- coef(cook_shk2)
cocook[model_shk2$dynamics@paramnames]
plotFormula(model_shk2, cook_shk2@transformed.parameters)

#------------------------------------------------------------
# detect outliers
# output of 'dynr.taste' is a list containing lists of results from the outlier detection process,
#   which will be input of the autoplot.dynrTaste and dynr.taste2
taste_shk <- dynr.taste(model_shk, cook_shk, conf.level=.99)
taste_shk2 <- dynr.taste(model_shk2, cook_shk2, dynm_shk2, conf.level=.99)

## create ggplot objects
# for 2 participants (default) who have largest chi-square values,
# ggplot objects for chi-square, and t statistics for all variables
taste_plot1 <- autoplot(taste_shk)
taste_plot1

# ggplot objects for chi-square, and t statistics for state variable `eta_1' and observed variables c('V1', 'V3')
taste_plot2 <- autoplot(taste_shk, 
                        names.state='eta_1', names.observed=c('V1', 'V3'))
taste_plot2

# for id number c(3, 5), 
# ggplot objects for chi-square, and t statistics for state variable `eta_1' and observed variables c('V1', 'V3')
taste_plot3 <- autoplot(taste_shk, 
                        names.state='eta_1', names.observed=c('V3', 'V1'), 
                        idtoPlot=c(3,5))
taste_plot3

# apply the detected outliers, and re-fit the model
taste2_shk <- dynr.taste2(model_shk, cook_shk, taste_shk,
                          newOutfile="taste2_shk.c")
taste2_shk2 <- dynr.taste2(model_shk2, cook_shk2, taste_shk2,
                          newOutfile="taste2_shk2.c")

# compare true parameter and estimated ones
#  cook: parameter estimates when the data include true outliers
#  taste2: parameter estimates after applying the estimated outliers by the function 'dynr.taste2'
# check the how closely recover the true parameter from 'dynr.taste2', i.e., after applying outliers
data.frame(true=true_params, cook=coef(cook_shk)[1:17],
           taste2=coef(taste2_shk$dynrCook_new)[1:17])

