#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2017-09-14
# Filename: dynrTaste.R
# Purpose: Compute shocks and chi-squared diagnostics following
#  Chow, Hamaker, and Allaire (2009).  Using Innovative Outliers to
#    Detect Discrete Shifts in Dynamics in Group-Based State-Space Models.
#    Multivariate Behavioral research.  44:465–496.
#    DOI: 10.1080/00273170903103324
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
##' Detect discrete shifts in in time series
##'
##' @param dynrModel an object of `dynrModel' class.
##' @param dynrCook the `dynrCook' object fitted with `debug_flag=TRUE' for the `dynrModel' object. The default is NULL. 
##' If the dynrCook object were not provided, or the object were cooked
##' with `debug_flag=FALSE', 
##' \code{"dynr.taste"} will fit the dynrModel object with `debug_flag=TRUE' internally.
##'
##' @return an object of `dynr.taste' class.
##' The function summary is used to obtain a shock detection summary and optional plots.
##' 
##' @references 
##' Chow, S.-M., Hamaker, E. L., & Allaire, J. C. (2009). 
##' Using innovative outliers to detectdiscrete shifts in dynamics in group-based state-space models. _Multivariate BehavioralResearch_, 44, 465–496. 
dynr.taste <- function(dynrModel, dynrCook=NULL) {
  # check for non-regime switching
  if (dynrModel$num_regime > 2) {
    stop("This test is for non-regime switching models.")
  }
  if ( is.null(dynrCook) ||
    !("residual_cov" %in% slotNames(dynrCook)) ) {# re-cooking
    # N.B. 'optimization_flag=FALSE' produces different logLik with 'TRUE'
    dynrCook <- dynr.cook(dynrModel, verbose=FALSE, debug_flag=TRUE)
  }
  stateName <- dynrModel$measurement$state.names
  latentDim <- length(stateName)
  observedDim <- length(dynrModel$measurement$obs.names)
  timeDim <- length(dynrModel$data$time)
  personID <- unique(dynrModel$data$id)
  nID <- length(personID)
  
  # fitted parameters
  coefx <- coef(dynrCook)
  coef(dynrModel) <- coefx # modifies all values throughout dynrModel
  Lambda <- dynrModel$measurement$values.load[[1]]
  B <- dynrModel$dynamics$values.dyn[[1]]
  
  # Array of inverse covariance matrices (i.e. information matrices) for the observed variables
  # F^-1  in Chow, Hamaker, and Allaire
  F_inv <- array(apply(dynrCook$residual_cov, 3, solve),
                 c(observedDim, observedDim, timeDim))
  
  # Compute Kalman gain for every person/time combination
  P_pred <- dynrCook$error_cov_predicted
  t_Lambda <- t(Lambda)
  v <- dynrCook$innov_vec
  K <- array(NA, c(latentDim, observedDim, timeDim))
  obsChi <- numeric(timeDim)
  
  for(i in 1:timeDim){
    F_inv_i <- F_inv[,,i]
    K[,,i] <- P_pred[,,i] %*% t_Lambda %*% F_inv_i # [lat, obs]
    v_i <- v[,i] # [obs, 1]
    obsChi[i] <- t(v_i) %*% F_inv_i %*% v_i
  }
  # N.B. the first element in P_pred is the initial, predicted latent covariance matrix, i.e. from dynrInitial
  #  The second element is the predicted cov for time=2 given the updated cov from time=1.
  
  # For each person, loop backward from their final time to their first time
  # computing r and N
  r <- matrix(NA, latentDim, timeDim)
  N <- array(NA, c(latentDim, latentDim, timeDim))
  u <- matrix(0, observedDim, timeDim)
  tstart <- dynrModel$data$tstart
  latChi <- numeric(timeDim)
  
  for(j in 1:nID){
    beginTime <- tstart[j] + 1 # Add 1 because model_tstart indexes from 0 rather than 1 for interface to C
    endTime <- tstart[j+1] # Use the 0-indexed beginning position of the next person as the 1-indexed end of the current person
    
    # set endTime r and N to 0
    r[,endTime] <- 0
    N[,,endTime] <- 0
    
    for(i in endTime:(beginTime+1)){
      # Probably could/should compute shocks in this same loop!!!
      ri <- matrix(r[,i], latentDim, 1)
      Ni <- matrix(N[,,i], latentDim, latentDim)
      Ki <- matrix(K[,,i], latentDim, observedDim)
      Li <- B - Ki %*% Lambda
      obsInfI <- matrix(F_inv[,,i], observedDim, observedDim) #TODO check for off by one index
      vi <- matrix(v[,i], nrow=observedDim, ncol=1) #TODO check for off by one index
      ui <- obsInfI %*% vi - t(Ki) %*% ri 
      u[, i] <- ui
      # Could also just compute shocks in this loop
      rnew <- t_Lambda %*% ui + t(B) %*% ri
      r[,i-1] <- rnew
      Nnew <- t_Lambda %*% obsInfI %*% Lambda + t(Li) %*% Ni %*% Li
      N[,,i-1] <- Nnew
      Ninv <- try(solve(Nnew), silent=TRUE)
      if(class(Ninv) == "try-error"){Ninv <- MASS::ginv(Nnew)}
      latChi[i-1] <- t(rnew) %*% Ninv %*% rnew
    }
  }

  # t-values
  t_value <- matrix(NA, latentDim, timeDim)
  for (t in 1:timeDim) {# CHOW, HAMAKER, ALLAIRE (2009) eq. (27)
    t_value[,t] <- r[,t] / sqrt(diag(N[,,t]))  
  }
  
  ################ delta estimate #############################
  delta <- matrix(NA, latentDim, timeDim)
  rownames(delta) <- paste0("delta_", stateName)
  # assume independent shock indicator, TODO: function arguments??
  # assume delta [latentDim,1]
  # p 469. delta = [q, 1]. how to choose q??
  W_t <- diag(1, latentDim)
  # no measurement shock assumed, TODO: function arguments??
  X_t <- matrix(0, observedDim, latentDim)
  
  for(i in 1:nID){
    beginTime <- tstart[i] + 1
    endTime <- tstart[i+1]
    time_i <- endTime - tstart[i]
    
    for (j in beginTime:endTime) {
      Q_j <- W_t - K[,,j] %*% X_t
      S_j <- t(X_t) %*% F_inv[,,j] %*% X_t +
        t(Q_j) %*% N[,,j] %*% Q_j
      s_j <- t(X_t) %*% u[,j] + t(W_t) %*% r[,j]
      S_j_inv <- try(solve(S_j), silent=TRUE)
      if(class(S_j_inv) == "try-error"){S_j_inv <- MASS::ginv(S_j)}
      delta[,j] <- S_j_inv %*% s_j
    }
  }
  
  rownames(t_value) <- paste0("t_", stateName)
  
  taste <- list(
    stats = data.frame(id=dynrModel$data$id, time=dynrModel$data$time, 
                       chi.l=latChi, chi.o=obsChi, 
                       t(t_value),
                       t(delta)),
    stateName=stateName, tstart=tstart,
    observedDim=observedDim, latentDim=latentDim)
  class(taste) <- "dynr.taste"
  return(taste)
}


computeJacobian <- function(cookDebug, jacobian, stateName, params, time){
  envList <- as.list(params)
  for(i in 1:length(stateName)) envList[[stateName[i]]] <- cookDebug$eta_smooth_final[i, time]
  
  J <- matrix(NA, length(stateName), length(stateName), dimnames=list(stateName, stateName))
  for(i in 1:(length(stateName))^2){
    cj <- as.character(jacobian[[i]])
    rc <- strsplit(cj[2], ' ~ ')[[1]]
    J[rc[1], rc[2]] <- eval(parse(text=cj[3]), envList)
  }
  return(J)
}

# example use from demo/NonlinearODE.R
#computeJacobian(res, model$dynamics$jacobianOriginal[[1]], model$measurement$state.names, coef(res), 1)



