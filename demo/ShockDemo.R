#---------------------------------------------------------------------
# Author: Dongjun You
# Date: 2018-05-15
# Filename: ShockDemo.R
# Purpose: An illustrative example of using dynr.taste to detect
#   outliers from a state space model
#---------------------------------------------------------------------
library(dynr)

beta2 <- matrix(c(0.6, -0.2,
                 -0.2,  0.5), ncol=2, byrow=TRUE)
lambda2 <- matrix(c(1.0, 0.0,
                    0.9, 0.0,
                    0.8, 0.0,
                    0.0, 1.0,
                    0.0, 0.9,
                    0.0, 0.8), ncol=2, byrow=TRUE)
psi2 <- matrix(c(0.3, -0.1,
                -0.1,  0.3), ncol=2, byrow=TRUE)
theta2 <- diag(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2), ncol=6, nrow=6)

true_params <- c(.6, -.2, -.2, .5,# beta
                 .9, .8, .9, .8,# lambda
                 .3, -.1, .3,# psi
                 .2, .2, .2, .2, .2, .2)# theta

data(ShockData)

data_shk <- dynr.data(ShockData, id='id', time='time',
                      observed=c('V1','V2','V3','V4','V5','V6'))

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
  state.names=c("eta_1",'eta_2'),
  obs.names=c('V1','V2','V3','V4','V5','V6') )

nois_shk <- prep.noise(
  values.latent=matrix(c(0.3, -0.1,
                         -0.1, 0.3), ncol=2, byrow=TRUE),
  params.latent=matrix(c('psi_11','psi_12',
                         'psi_12','psi_22'), ncol=2, byrow=TRUE),
  values.observed=diag(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2), ncol=6, nrow=6),
  params.observed=diag( paste0('e_', 1:6), 6) )

init_shk <- prep.initial(
  values.inistate=c(0,0),
  params.inistate=c('mu_1','mu_2'),
  values.inicov=matrix(c(0.3, -0.1,
                        -0.1,  0.3), ncol=2, byrow=TRUE),
  params.inicov=matrix(c('c_11','c_12',
                         'c_12','c_22'), ncol=2, byrow=TRUE) )

dynm_shk <- prep.matrixDynamics(
  values.dyn=matrix(c(0.8, -0.2,
                      -0.2,  0.7), ncol=2, byrow=TRUE),
  params.dyn=matrix(c('b_11','b_12',
                      'b_21','b_22'), ncol=2, byrow=TRUE),
  isContinuousTime=FALSE)

model_shk <- dynr.model(dynamics=dynm_shk, measurement=meas_shk,
                        noise=nois_shk, initial=init_shk,
                        data=data_shk, outfile=paste0("model_shk.c"))

cook_shk <- dynr.cook(model_shk, debug_flag=TRUE)
saveRDS(cook_shk, "cook_shk.rds")
cook_shk <- readRDS("cook_shk.rds")

# detect outliers
taste_shk <- dynr.taste(model_shk, cook_shk, conf.level=.99)

## create ggplot objects
# for 2 participants (default) who have largest chi-square values,
# ggplot objects for chi-square, and t statistics for all variables
taste_plot1 <- autoplot(taste_shk)
# ggplot objects for chi-square, and t statistics for state variable `eta_1' and observed variables c('V1', 'V3')
taste_plot2 <- autoplot(taste_shk, 
                        names.state='eta_1', names.observed=c('V1', 'V3'))
# for id number c(3, 5), 
# ggplot objects for chi-square, and t statistics for state variable `eta_1' and observed variables c('V1', 'V3')
taste_plot3 <- autoplot(taste_shk, 
                        names.state='eta_1', names.observed=c('V3', 'V1'), 
                        idtoPlot=c(3,5))

taste_shk11 <- dynr.taste(model_shk, cook_shk, 
                          conf.level=.99, debug_flag=TRUE)

# apply the detected outliers, and re-fit the model
taste2_shk <- dynr.taste2(model_shk, cook_shk, taste_shk,
                          newOutfile="taste2_shk.c")
saveRDS(taste2_shk, "taste2_shk.rds")
taste2_shk <- readRDS("taste2_shk.rds")
# (optional) detect outliers from the re-fitted model
taste_taste2_shk <- dynr.taste(taste2_shk$dynrModel_new,
                               taste2_shk$dynrCook_new,
                               conf.level=.99)

# compare true parameter and estimated ones
#  cook: parameter estimates when the data include true outliers
#  taste2: parameter estimates after applying the estimated outliers by the function 'dynr.taste2'
# check the how closely recover the true parameter from 'dynr.taste2', i.e., after applying outliers
data.frame(true=true_params, cook=coef(cook_shk)[1:17],
           taste2=coef(taste2_shk$dynrCook_new)[1:17])

