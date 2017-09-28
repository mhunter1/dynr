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
##' @param cookDebug a dynrCook class fitted with `debug_flag=TRUE'.
##' @param dynrModel a dynrModel class.
##' @param conf.level a numeric of confidence level. The default is 0.99.
##' @return a list of shock points.
##' 
##'  @references 
##'  Chow, S.-M., Hamaker, E. L., & Allaire, J. C. (2009). 
##'  Using innovative outliers to detectdiscrete shifts in dynamics in group-based state-space models. _Multivariate BehavioralResearch_, 44, 465–496. 
dynr.taste <- function(cookDebug, dynrModel, conf.level=0.99) {
  if ( !("residual_cov" %in% slotNames(cookDebug)) ) {
    stop("Please cook again with 'debug_flag=TRUE'.")
  }
  # check for non-regime switching
  if (dynrModel$num_regime > 2) {
    stop("This test is for non-regime switching models.")
  }
  
  stateName <- dynrModel$measurement$state.names
  latentDim <- length(stateName)
  observedDim <- length(dynrModel$measurement$obs.names)
  timeDim <- dim(cookDebug$residual_cov)[3]
  personID <- unique(dynrModel$data$id)
  
  # fitted parameters
  coefx <- coef(cookDebug)
  
  
  # DONGJUN, Are the next two 4-line blocks of code there just to put the estimated
  #  parameters back into the model?  If so, this is much easier and more
  #  reliable with
  #  coef(dynrModel) <- coefx # modifies all values throughout dynrModel
  #  Lambda <- dynrModel$measurement$values.load[[1]]
  #  B <- dynrModel$dynamics$values.dyn[[1]]
  
  # Measurement matrix
  # Lambda in Chow, Hamaker, and Allaire
  # number matches the parameters in coefx
  # e.g., 0 0
  #       3 0
  #       4 0
  # 3, 4th elements in coefx match
  Lambda_p <- dynrModel$measurement$params.load[[1]]
  Lambda <- dynrModel$measurement$values.load[[1]]
  Lambda_idx <- which(Lambda_p !=0, arr.ind=TRUE)
  Lambda[Lambda_idx] <- coefx[ Lambda_p[Lambda_idx] ]
  
  # Dynamics matrix
  # B in Chow, Hamaker, and Allaire
  B_p <- dynrModel$dynamics$params.dyn[[1]]
  B <- dynrModel$dynamics$values.dyn[[1]]
  B_idx <- which(B_p !=0, arr.ind=TRUE)
  B[B_idx] <- coefx[ B_p[B_idx] ]
  
  # Array of inverse covariance matrices (i.e. information matrices) for the observed variables
  # F^-1  in Chow, Hamaker, and Allaire
  F_inv <- array(apply(cookDebug$residual_cov, 3, solve), 
                 c(observedDim, observedDim, timeDim))
  
  # Compute Kalman gain for every person/time combination
  P_pred <- cookDebug$error_cov_predicted
  t_Lambda <- t(Lambda)
  v <- cookDebug$innov_vec
  K <- array(NA, c(observedDim, latentDim, timeDim))
  obsChi <- numeric(timeDim)
  
  for(i in 1:timeDim){
    F_inv_i <- F_inv[,,i]
    K[,,i] <- P_pred[,,i] %*% t_Lambda %*% F_inv_i
    v_i <- v[,i] # matrix(v[,i], ncol=1)
    obsChi[i] <- t(v_i) %*% F_inv_i %*% v_i
  }
  
  #qqplot(qchisq(ppoints(timeDim), df=observedDim), obsChi)
  #qqline(obsChi, distribution = function(p) qchisq(p, df = observedDim), prob = c(0.1, 0.6))
  
  # For each person, loop backward from their final time to their first time
  # computing r and N
  r <- matrix(NA, latentDim, timeDim)
  N <- array(NA, c(latentDim, latentDim, timeDim))
  model_tstart <- dynrModel$data$tstart
  latChi <- numeric(timeDim)
  for(j in 1:length(personID)){
    beginTime <- model_tstart[j] + 1
    endTime <- model_tstart[j+1]
    # set endTime r and N to 0
    r[,endTime] <- 0
    N[,,endTime] <- 0
    for(i in (endTime-1):beginTime){
      # Probably could/should compute shocks in this same loop!!!
      ri <- matrix(r[,i+1], latentDim, 1)
      Ni <- matrix(N[,,i+1], latentDim, latentDim)
      Ki <- matrix(K[,,i+1], latentDim, observedDim)
      Li <- B - Ki %*% Lambda
      obsInfI <- matrix(F_inv[,,i+1], observedDim, observedDim) #TODO check for off by one index
      vi <- matrix(v[,i+1], nrow=observedDim, ncol=1) #TODO check for off by one index
      ui <- obsInfI %*% vi - t(Ki) %*% ri # TODO save ui? It is used in s_ij Eq 26.
      # Could also just compute shocks in this loop
      rnew <- t_Lambda %*% ui + t(B) %*% ri
      r[,i] <- rnew
      Nnew <- t_Lambda %*% obsInfI %*% Lambda + t(Li) %*% Ni %*% Li
      N[,,i] <- Nnew
      Ninv <- try(solve(Nnew), silent=TRUE)
      if(class(Ninv) == "try-error"){Ninv <- MASS::ginv(Nnew)}
      latChi[i] <- t(rnew) %*% Ninv %*% rnew
    }
  }
  #qqplot(qchisq(ppoints(timeDim), df=latentDim), latChi)
  #qqline(latChi, distribution = function(p) qchisq(p, df = latentDim), prob = c(0.1, 0.6))
  
  #qqplot(qchisq(ppoints(timeDim), df=observedDim+latentDim), obsChi + latChi)
  #qqline(obsChi + latChi, distribution = function(p) qchisq(p, df = observedDim+latentDim), prob = c(0.1, 0.6))
  
  # start time point of each subject
  t_start <- (model_tstart + 1)[-length(model_tstart)]
  # end time point of each subject
  t_end <- model_tstart[-1]
  # number of time points of each subject, which is df for t-test
  t_df <- t_end - model_tstart[-length(model_tstart)]
  
  ### chi-square test
  # p-values of chi-sqaure
  chi_obs_pval <- pchisq(obsChi, df=observedDim, lower.tail=FALSE)
  chi_lat_pval <- pchisq(latChi, df=latentDim, lower.tail=FALSE)
  # split chi-square p-values for each subject
  chi_obs_pval_spl <- lapply(1:length(personID), function(i) {
    chi_obs_pval[ t_start[i]:t_end[i] ]
  })
  chi_lat_pval_spl <- lapply(1:length(personID), function(i) {
    chi_obs_pval[ t_start[i]:t_end[i] ]
  })
  # locate shock points for each subject.
  chi_obs_shock_spl <- lapply(chi_obs_pval_spl, function(chi_obs) {
    which( chi_obs < 1 - conf.level )# TRUE for significant chi
  })
  chi_lat_shock_spl <- lapply(chi_lat_pval_spl, function(chi_lat) {
    which( chi_lat < 1 - conf.level )# TRUE for significant chi
  })
  
  list(obsShock=chi_obs_shock_spl,
       latShock=chi_obs_shock_spl)
  # ### t test
  # t_value <- matrix(NA, latentDim, timeDim)
  # for (t in 1:timeDim) {# CHOW, HAMAKER, ALLAIRE (2009) eq. (27)
  #   t_value[,t] <- r[,t] / sqrt(diag(N[,,t]))  
  # }
  # # split t-values for each subject
  # t_value_spl <- lapply(1:length(personID), function(i) {
  #   t_value[, t_start[i]:t_end[i]]
  # })
  # # calculate p-values of t-values
  # alternative <- match.arg(alternative)
  # if (alternative=="two.sided") {
  #   t_pval_spl <- lapply(1:length(personID), function(i) {
  #     2 * pt(-abs(t_value_spl[[i]]), df=(t_df[i]-latentDim))
  #   })
  # }	else if (alternative=="greater") {
  #   t_pval_spl <- lapply(1:length(personID), function(i) {
  #     pt(t_value_spl[[i]], df=(t_df[i]-latentDim), lower.tail=FALSE)
  #   })
  # }	else {
  #   t_pval_spl <- lapply(1:length(personID), function(i) {
  #     pt(t_value_spl[[i]], df=(t_df[i]-latentDim), lower.tail=TRUE)
  #   })
  # }
  # 
  # #qqplot(qt(ppoints(timeDim), df=t_df-latentDim), t_values)
  # #qqline(t_values, distribution = function(p) qt(p, df = t_df-latentDim), prob = c(0.1, 0.6))
  # 
  # # locate shocks from t-test. [observed_variable, t] for each subject
  # t_shock <- lapply(t_pval_spl, function(i) {
  #   i < 1 - conf.level# TRUE for significant t
  # })
  # 
  # res <- mapply(function(chishock, tshock, tvalue, tpval) {
  #   # locate time points of t-test that pass chi-square test
  #   loc <- sweep(tshock, 2, chishock, FUN="&")# TRUE for pass
  #   loc_lt <- which(loc, arr.ind=TRUE)# row: latent, col: time
  #   data.frame(
  #     latent=stateName[loc_lt[, c("row")]],
  #     time=loc_lt[, c("col")],
  #     t_value=tvalue[loc],
  #     p=tpval[loc] )
  # },
  # chi_shock, t_shock, t_value_spl, t_pval_spl, SIMPLIFY=FALSE)
  
  # names(res) <- personID
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



