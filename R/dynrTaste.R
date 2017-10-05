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
##' @param file a character of file name. A pdf file with the name will
##' be saved.
##' @param numSub number of subject to print with the largest Chi-square
##' @param idtoPlot a vector of ids.
##' @param names.state names of state to print out
##' @return a list of shock points.
##' 
##' @references 
##' Chow, S.-M., Hamaker, E. L., & Allaire, J. C. (2009). 
##' Using innovative outliers to detectdiscrete shifts in dynamics in group-based state-space models. _Multivariate BehavioralResearch_, 44, 465–496. 
dynr.taste <- function(cookDebug, dynrModel, conf.level=0.99, file=NULL,
                       numSubj=NULL, idtoPlot=NULL, names.state=NULL) {
  if ( !("residual_cov" %in% slotNames(cookDebug)) ) {
    stop("Please cook again with 'debug_flag=TRUE'.")
  }
  # N.B. maybe we could re-cook (re-heat!) the model with debug_flag=TRUE and optimization_flag=FALSE.  That is, send the model to the backend, but do not optimization.  Instead, just re-compute the Kalman scores etc.
  
  # check for non-regime switching
  if (dynrModel$num_regime > 2) {
    stop("This test is for non-regime switching models.")
  }
  
  stateName <- dynrModel$measurement$state.names
  if ( is.null(names.state) ) {
    names.state <- stateName
  } else 
    if ( !all(sn <- names.state %in% stateName) ) {
      stop(paste("cannot recognize the state name(s):", 
                 paste(names.state[!sn], collapse=", ")))
    }
  latentDim <- length(stateName)
  observedDim <- length(dynrModel$measurement$obs.names)
  timeDim <- dim(cookDebug$residual_cov)[3]
  personID <- unique(dynrModel$data$id)
  nID <- length(personID)
  
  # fitted parameters
  coefx <- coef(cookDebug)
  coef(dynrModel) <- coefx # modifies all values throughout dynrModel
  Lambda <- dynrModel$measurement$values.load[[1]]
  B <- dynrModel$dynamics$values.dyn[[1]]
  
  # Array of inverse covariance matrices (i.e. information matrices) for the observed variables
  # F^-1  in Chow, Hamaker, and Allaire
  F_inv <- array(apply(cookDebug$residual_cov, 3, solve),
                 c(observedDim, observedDim, timeDim))
  
  # Compute Kalman gain for every person/time combination
  P_pred <- cookDebug$error_cov_predicted
  t_Lambda <- t(Lambda)
  v <- cookDebug$innov_vec
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
  
  #qqplot(qchisq(ppoints(timeDim), df=observedDim), obsChi)
  #qqline(obsChi, distribution = function(p) qchisq(p, df = observedDim), prob = c(0.1, 0.6))
  
  # For each person, loop backward from their final time to their first time
  # computing r and N
  r <- matrix(NA, latentDim, timeDim)
  N <- array(NA, c(latentDim, latentDim, timeDim))
  u <- matrix(0, observedDim, timeDim)
  model_tstart <- dynrModel$data$tstart
  latChi <- numeric(timeDim)
  
  for(j in 1:length(personID)){
    beginTime <- model_tstart[j] + 1 # Add 1 because model_tstart indexes from 0 rather than 1 for interface to C
    endTime <- model_tstart[j+1] # Use the 0-indexed beginning position of the next person as the 1-indexed end of the current person
    
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
  chi_obs_pval_spl <- lapply(1:nID, function(i) {
    chi_obs_pval[ t_start[i]:t_end[i] ]
  })
  chi_lat_pval_spl <- lapply(1:nID, function(i) {
    chi_lat_pval[ t_start[i]:t_end[i] ]
  })
  # locate shock points for each subject.
  chi_obs_shk_spl <- lapply(chi_obs_pval_spl, function(chi_obs) {
    which(chi_obs < 1 - conf.level)# TRUE for significant chi
  })
  chi_lat_shk_spl <- lapply(chi_lat_pval_spl, function(chi_lat) {
    which(chi_lat < 1 - conf.level)# TRUE for significant chi
  })
  
  ##### delta estimate
  delta <- matrix(NA, timeDim, latentDim)
  colnames(delta) <- stateName
  rownames(delta) <- 1:timeDim
  # assume independent shock indicator, TODO: function arguments??
  # assume delta [latentDim,1]
  # p 469. delta = [q, 1]. how to choose q??
  W_t <- diag(1, latentDim)
  # no measurement shock assumed, TODO: function arguments??
  X_t <- matrix(0, observedDim, latentDim)
  
  for(i in 1:nID){
    beginTime <- model_tstart[i] + 1
    endTime <- model_tstart[i+1]
    time_i <- endTime - model_tstart[i]
    
    for (j in beginTime:endTime) {
      Q_j <- W_t - K[,,j] %*% X_t
      S_j <- t(X_t) %*% F_inv[,,j] %*% X_t +
        t(Q_j) %*% N[,,j] %*% Q_j
      s_j <- t(X_t) %*% u[,j] + t(W_t) %*% r[,j]
      S_j_inv <- try(solve(S_j), silent=TRUE)
      if(class(S_j_inv) == "try-error"){S_j_inv <- MASS::ginv(S_j)}
      delta[j,] <- S_j_inv %*% s_j
    }
  }
  
  delta_spl <- lapply(1:nID, function(i) {
    delta[ t_start[i]:t_end[i], ]
  })
  
  ########### plots ###################
  # split eta_smooth for each subject
  etaSmooth <- cookDebug$eta_smooth_final
  etaSmooth_spl <- lapply(1:nID, function(i) {
    etaSmooth[, t_start[i]:t_end[i] ]
  })
  # split error_cov_smooth for each subject
  errSmooth <- cookDebug$error_cov_smooth_final
  errSmooth_spl <- lapply(1:nID, function(i) {
    errSmooth[,, t_start[i]:t_end[i] ]
  })
  # split latent Chi-square for each subject
  latChi_spl <- lapply(1:nID, function(i) {
    latChi[ t_start[i]:t_end[i] ]
  })
  
  if (is.null(numSubj)) {
    numSubj_idx <- !vector(length=nID) # all TRUE
  } else {
    if (!is.null(idtoPlot)) {
      stop("'numSubj' and 'idotoPlot' cannot be both activatied.")
    }
    # max chi for each subject. save as vector
    chiMax <- sapply(latChi_spl, function(i) {
      max(i)
    })
    numSubj_idx <- rank(-chiMax) %in% 1:numSubj
  }
  
  if (is.null(idtoPlot)) {
    idtoPlot_idx <- !vector(length=nID) # all TRUE
  } else {
    if ( !all(idp <- idtoPlot %in% personID) ) {
      stop(paste("cannot recognize the id(s):", 
                 paste(idtoPlot[!idp], collapse=", ")))
    }
    idtoPlot_idx <- personID %in% idtoPlot
  }
  
  plots <- mapply(function(etaSm_i, errSm_i, latChi_i, lat_shk_i, id_i, numSubj_idx_i, idtoPlot_idx_i) {
    if (numSubj_idx_i && idtoPlot_idx_i) {
      plots_i <- lapply(1:length(names.state), function(s) {
        sName <- names.state[s]
        sProcess <- etaSm_i[s, ]
        errProcess <- errSm_i[s,s, ]
        shked_t <- lat_shk_i + 1 # shocked time points
        shked_s <- sProcess[shked_t] # state values at shocked time points
        shked_e <- 1.96 * sqrt(errProcess[shked_t]) # err s.d. at shocked time points
        
        sDF <- data.frame(t=1:length(sProcess), state=sProcess)
        names(sDF)[names(sDF)=="state"] <- sName
        shkDF <- data.frame(t=shked_t, shock=shked_s)
        errbarDF <- data.frame(t=shked_t,
                               lo=shked_s - shked_e,
                               up=shked_s + shked_e)
        
        ggplot2::ggplot(sDF, aes_string(x="t", y=sName)) +
          ggplot2::geom_line(alpha=0.5, colour="blue") +
          ggplot2::geom_point(data=shkDF, aes_string(x="t", y="shock"),
                              colour="red", size=2.5, alpha=0.5) +
          ggplot2::xlab(NULL)
        #ggplot2::geom_errorbar(data=errbarDF, aes(ymin=lo, ymax=up), width=0.1)
      } )
      chiDF <- data.frame(t=1:length(latChi_i), Chi=latChi_i)
      chiPlot <- list(
        ggplot2::ggplot(chiDF, aes_string(x="t", y="Chi")) +
          ggplot2::geom_line(alpha=0.5, colour="black") +
          ggplot2::geom_hline(yintercept=qchisq(conf.level, length(stateName)), colour="red", alpha=0.5)
      )
      plotArgs <- list(ncol=1, nrow=length(stateName) + 1, align="v")
      plots_ii <- do.call(ggpubr::ggarrange, c(plots_i, chiPlot, plotArgs))
      ggpubr::annotate_figure(plots_ii,
                              top=ggpubr::text_grob(id_i, color="black", face="bold", size=15))
    }
  },
  etaSmooth_spl, errSmooth_spl, latChi_spl, chi_lat_shk_spl, personID, numSubj_idx, idtoPlot_idx,
  SIMPLIFY=FALSE)
  
  # N.B. out pdf file name
  if ( is.null(file) ) file <- "state_shock"
  ggpubr::ggexport(plots, filename=paste0(file, ".pdf"))
  #####################################
  
  list(lat_shock=chi_lat_shk_spl,
       obs_shock=chi_obs_shk_spl,
       delta=delta_spl)
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



