#---------------------------------------------------------------------
# Author: Dongjun You
# Date: 2018-05-15
# Filename: OutlierDetectionLinear.R
# Purpose: An illustrative example of using dynr.taste to detect
#   outliers from a state space model.
# Note: This is the small version of OutlierDetection.R. For the full
#  version see inst/models/passing/OutlierDetection.R
#---------------------------------------------------------------------
# Load packages
library(dynr)

# The data were generated from a discrete time state-space model with six observed variables and two latent variables
# 100 participants, 100 time-points for each participant
# For each participant, 3 innovative outliers and 6 additive outliers 
data("Outliers")
outdata <- Outliers$generated$y
shockO <- Outliers$generated$shocks$shockO
shockL <- Outliers$generated$shocks$shockL

useIds <- 1:20
outdata <- outdata[outdata$id %in% useIds, ]

# The true parameters used to generate the data
true_params <- c(.6, -.2, -.2, .5,# beta
                 .9, .8, .9, .8,# lambda
                 .3, -.1, .3,# psi
                 .2, .2, .2, .2, .2, .2)# theta

data_shk <- dynr.data(outdata, id='id', time='time',
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

# define the transition matrix
dynm_shk <- prep.matrixDynamics(
  values.dyn=matrix(c(0.6, -0.2,
                      -0.2,  0.5), ncol=2, byrow=TRUE),
  params.dyn=matrix(c('b_11','b_12',
                      'b_21','b_22'), ncol=2, byrow=TRUE),
  isContinuousTime=FALSE)

# Prepare for cooking
# put all the recipes together
model_shk <- dynr.model(dynamics=dynm_shk, measurement=meas_shk,
                        noise=nois_shk, initial=init_shk,
                        data=data_shk, outfile=paste0("model_shk.c"))

# Estimate free parameters
# SHOULD cook with 'debug_flag=TRUE' to extract information 
#   used for outlier detection
cook_shk <- dynr.cook(model_shk, debug_flag=TRUE, verbose=FALSE)

#------------------------------------------------------------
# detect outliers
# output of 'dynr.taste' is a list containing lists of results from the outlier detection process,
#   which will be input of the autoplot.dynrTaste and dynr.taste2
taste_shk <- dynr.taste(model_shk, cook_shk, conf.level=.99)

# apply the detected outliers, and re-fit the model
taste2_shk <- dynr.taste2(model_shk, cook_shk, taste_shk,
                          newOutfile="taste2_shk.c")

# compare true parameter and estimated ones
#  cook: parameter estimates when the data include true outliers
#  taste2: parameter estimates after applying the estimated outliers by the function 'dynr.taste2'
# check the how closely recover the true parameter
# from 'dynr.taste2', i.e., after applying outliers
comp_param <- data.frame(true=true_params, 
           cook=coef(cook_shk)[1:17],
           taste2=coef(taste2_shk$dynrCook_new)[1:17])
# The bigg improvement here is really the noise parameters
# psi_* and e_*.
# Compare parameter estimates before removing outliers (cook)
#  to those after removing outliers (taste2).
rms <- function(x, y){sqrt(mean((x-y)^2))}
testthat::expect_true(
	rms(comp_param$cook, comp_param$true) > rms(comp_param$taste2, comp_param$true)
)
# The parameter estimates after removing outliers are closer to the true
#  values.


##---------------------
## Checking detected outliers, comparing it with true outliers
# A helper function to extract detected outliers
shock_check <- function(dynrModel, dynrTaste, shocks) {
  tstart <- dynrModel$data$tstart
  shock_inn <- dynrTaste$t.inn.shock
  shock_add <- dynrTaste$t.add.shock
  nlat <- nrow(shock_inn)
  nobs <- nrow(shock_add)
  id <- dynrModel$data$id
  IDs <- unique(id)
  nID <- length(IDs)
  
  shocked_L <- vector('list', nID*nlat)
  shocked_O <- vector('list', nID*nobs)
  for(j in 1:nID) {
    beginTime <- tstart[j] + 1
    endTime <- tstart[j+1]
    
    for(l in 1:nlat) {
      time_L <- which(shock_inn[l, beginTime:endTime])
      if (length(time_L) > 0) {
        shocked_L[[(j-1)+l]] <- data.frame(id=j, time_L=time_L, lat=l)
      } else {
        shocked_L[[(j-1)+l]] <- NULL
      }
    }
    for(o in 1:nobs) {
      time_O <- which(shock_add[o, beginTime:endTime])
      if (length(time_O) > 0) {
        shocked_O[[(j-1)+o]] <- data.frame(id=j, time_O=time_O, obs=o)
      } else {
        shocked_O[[(j-1)+o]] <- NULL
      }
    }
  }
  list(sh_L = do.call("rbind", shocked_L), 
       sh_O = do.call("rbind", shocked_O))
}

# pre-saved detected outliers
detect_L <- Outliers$detect_L
detect_O <- Outliers$detect_O
detect_L <- detect_L[detect_L$id %in% useIds, ]
detect_O <- detect_O[detect_O$id %in% useIds, ]

# detected outliers
detected <- shock_check(model_shk, taste_shk)
detected_L <- detected$sh_L
detected_O <- detected$sh_O

#concatenate three variables into one.  check for hit rates acorss data sets
trueLatentOutliers <- apply(detect_L, 1, function(x){paste(unlist(x), sep='', collapse='_')})
foundLatentOutliers <- apply(detected_L, 1, function(x){paste(unlist(x), sep='', collapse='_')})
trueObsOutliers <- apply(detect_O, 1, function(x){paste(unlist(x), sep='', collapse='_')})
foundObsOutliers <- apply(detected_O, 1, function(x){paste(unlist(x), sep='', collapse='_')})

# percent of found outliers that were in the true set
tpr_l <- sum(foundLatentOutliers %in% trueLatentOutliers)/length(foundLatentOutliers)
testthat::expect_true(tpr_l >= 0.93)
tpr_o <- sum(foundObsOutliers %in% trueObsOutliers)/length(foundObsOutliers)
testthat::expect_true(tpr_o >= 0.71)

# percent of true outliers that were found
pf_l <- sum(trueLatentOutliers %in% foundLatentOutliers)/length(trueLatentOutliers)
testthat::expect_identical(pf_l, 1)
pf_o <- sum(trueObsOutliers %in% foundObsOutliers)/length(trueObsOutliers)
testthat::expect_true(pf_o >= .90)

