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
##' @param conf.level a numeric of confidence level. The default is 0.99.
##' @param alternative a character string specifying the alternative hypothesis of t-test,
##' must be one of ``two.sided'' (default), ``greater''  or ``less''
##' @return an object of `dynrTaste' class.
##' The function summary is used to obtain a shock detection summary and optional plots.
##' 
##' @references 
##' Chow, S.-M., Hamaker, E. L., & Allaire, J. C. (2009). 
##' Using innovative outliers to detectdiscrete shifts in dynamics in group-based state-space models. _Multivariate BehavioralResearch_, 44, 465–496. 
dynr.taste <- function(dynrModel, dynrCook=NULL, conf.level=0.99,
                       alternative=c("two.sided", "less", "greater")) {
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
  dimLat <- length(stateName)
  dimObs <- length(dynrModel$measurement$obs.names)
  dimTime <- length(dynrModel$data$time)
  IDs <- unique(dynrModel$data$id)
  nID <- length(IDs)
  
  # fitted parameters
  coefx <- coef(dynrCook)
  coef(dynrModel) <- coefx # modifies all values throughout dynrModel
  Lambda <- dynrModel$measurement$values.load[[1]]
  B <- dynrModel$dynamics$values.dyn[[1]]
  
  # Array of inverse covariance matrices (i.e. information matrices) for the observed variables
  # F^-1  in Chow, Hamaker, and Allaire
  F_inv <- array(apply(dynrCook$residual_cov, 3, solve),
                 c(dimObs, dimObs, dimTime))
  
  # Compute Kalman gain for every person/time combination
  P_pred <- dynrCook$error_cov_predicted
  t_Lambda <- t(Lambda)
  v <- dynrCook$innov_vec
  K <- array(NA, c(dimLat, dimObs, dimTime))
  chiObs <- numeric(dimTime)
  
  for(i in 1:dimTime){
    F_inv_i <- F_inv[,,i]
    K[,,i] <- P_pred[,,i] %*% t_Lambda %*% F_inv_i # [lat, obs]
    v_i <- v[,i] # [obs, 1]
    chiObs[i] <- t(v_i) %*% F_inv_i %*% v_i
  }
  # N.B. the first element in P_pred is the initial, predicted latent covariance matrix, i.e. from dynrInitial
  #  The second element is the predicted cov for time=2 given the updated cov from time=1.
  
  # For each person, loop backward from their final time to their first time
  # computing r and N
  r <- matrix(NA, dimLat, dimTime)
  N <- array(NA, c(dimLat, dimLat, dimTime))
  u <- matrix(0, dimObs, dimTime)
  tstart <- dynrModel$data$tstart
  chiLat <- numeric(dimTime)
  
  for(j in 1:nID){
    beginTime <- tstart[j] + 1 # Add 1 because model_tstart indexes from 0 rather than 1 for interface to C
    endTime <- tstart[j+1] # Use the 0-indexed beginning position of the next person as the 1-indexed end of the current person
    
    # set endTime r and N to 0
    r[,endTime] <- 0
    N[,,endTime] <- 0
    
    for(i in endTime:(beginTime+1)){
      # Probably could/should compute shocks in this same loop!!!
      ri <- matrix(r[,i], dimLat, 1)
      Ni <- matrix(N[,,i], dimLat, dimLat)
      Ki <- matrix(K[,,i], dimLat, dimObs)
      Li <- B - Ki %*% Lambda
      obsInfI <- matrix(F_inv[,,i], dimObs, dimObs) #TODO check for off by one index
      vi <- matrix(v[,i], nrow=dimObs, ncol=1) #TODO check for off by one index
      ui <- obsInfI %*% vi - t(Ki) %*% ri 
      u[, i] <- ui
      # Could also just compute shocks in this loop
      rnew <- t_Lambda %*% ui + t(B) %*% ri
      r[,i-1] <- rnew
      Nnew <- t_Lambda %*% obsInfI %*% Lambda + t(Li) %*% Ni %*% Li
      N[,,i-1] <- Nnew
      Ninv <- try(solve(Nnew), silent=TRUE)
      if(class(Ninv) == "try-error"){Ninv <- MASS::ginv(Nnew)}
      chiLat[i-1] <- t(rnew) %*% Ninv %*% rnew
    }
  }

  # t-values
  t_value <- matrix(NA, dimLat, dimTime)
  rownames(t_value) <- stateName
  for (t in 1:dimTime) {# CHOW, HAMAKER, ALLAIRE (2009) eq. (27)
    t_value[,t] <- r[,t] / sqrt(diag(N[,,t]))  
  }
  
  ################ delta estimate #############################
  delta <- matrix(NA, dimLat, dimTime)
  rownames(delta) <- stateName
  # assume independent shock indicator, TODO: function arguments??
  # assume delta [latentDim,1]
  # p 469. delta = [q, 1]. how to choose q??
  W_t <- diag(1, dimLat)
  # no measurement shock assumed, TODO: function arguments??
  X_t <- matrix(0, dimObs, dimLat)
  
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
  
  # N.B. 'data.table' can solve combersome conversions and splits
  ############ chi-square test ############################
  id <- as.factor(dynrModel$data$id)
  chiLat_sp <- split(chiLat, id)
  chiObs_sp <- split(chiObs, id)
  # p-values of chi-sqaure
  chiLat_pval <- pchisq(chiLat, df=dimLat, lower.tail=FALSE)
  chiLat_pval_sp <- split(chiLat_pval, id)
  chiObs_pval <- pchisq(chiObs, df=dimObs, lower.tail=FALSE)
  chiObs_pval_sp <- split(chiObs_pval, id)
  # locate shock points, TRUE for significance
  chiLat_shk <- chiLat_pval < (1 - conf.level)
  chiLat_shk_sp <- split(chiLat_shk, id)
  chiObs_shk <- chiObs_pval < (1 - conf.level)
  chiObs_shk_sp <- split(chiObs_shk, id)
  
  ################### t-test ###############################
  # [time, dimLat] matrix, for each subject
  t_sp <- split(as.data.frame(t(t_value)), id)
  # assign function calculating p-values
  alternative <- match.arg(alternative)
  t_calcPval <- if (alternative=="two.sided") {
    function(subj) {
      2 * pt( -abs(as.matrix(subj)), df=(nrow(subj)-dimLat) )
    }
  }	else if (alternative=="greater") {
    function(subj) {
      pt( as.matrix(subj), df=(nrow(subj)-dimLat), lower.tail=FALSE )
    }
  }	else {
    function(subj) {
      pt( as.matrix(subj), df=(nrow(subj)-dimLat), lower.tail=TRUE )
    }
  }
  # [time, dimLat] matrix, for each subject
  t_pval_sp <- lapply(t_sp, FUN=t_calcPval)
  # locate shocks from t-test
  # [time, dimLat] matrix, for each subject
  t_shk_sp <- lapply(t_pval_sp, function(subj) {
    subj < (1 - conf.level)
  })
  # split delta 
  # [time, dimLat] data.frame, for each subject
  rownames(delta) <- paste0("del_", rownames(delta))
  delta_sp <- split(as.data.frame(t(delta)), id)
  
  # output data.frame for each subject
  time_sp <- split(dynrModel$data$time, id)
  res <- mapply(FUN=function(id_i, time, 
                             chiLat, chiLat_pval, chiLat_shk, 
                             chiObs, chiObs_pval, chiObs_shk, 
                             tval, t_pval, t_shk, delta) {
    colnames(tval) <- paste0("t_", colnames(tval))
    colnames(t_pval) <- paste0("t.p_", colnames(t_pval))
    colnames(t_shk) <- paste0("t.shk_", colnames(t_shk))
    # t shock that pass chi shock
    chi_t_shk <- sweep(t_shk, 1, chiLat_shk, FUN="&")
    # delta that will be input to 'dynr.detox'
    delta_dtx <- delta
    delta_dtx[!chi_t_shk] <- 0
    row.names(delta_dtx) <- 1:nrow(delta_dtx)
    
    list(
      taste=data.frame(
        id=id_i, time=time,
        chi.L=chiLat, chi.L.p=chiLat_pval, chi.L.shk=chiLat_shk,
        chi.O=chiObs, chi.O.p=chiObs_pval, chi.O.shk=chiObs_shk,
        tval, t_pval, t_shk, chi.by=chi_t_shk, delta
      ),
      delta.dtx=delta_dtx
    )
  }, SIMPLIFY=FALSE,
  unique(id), time_sp, chiLat_sp, chiLat_pval_sp, chiLat_shk_sp, 
  chiObs_sp, chiObs_pval_sp, chiObs_shk_sp, 
  t_sp, t_pval_sp, t_shk_sp, delta_sp)
  
  # TODO. display output for users
  # res <- list(res=res1, cookTaste=dynrCook)
  class(res) <- "dynrTaste"
  invisible(res)
}

##' @param dynrModel an object of dynrModel class.
##' @param dynrTaste an object of dynrTaste class
dynr.taste2 <- function(dynrModel, dynrTaste) {
  # combine delta through subjects
  delta <- do.call("rbind",
                   lapply(dynrTaste, function(taste_i) {
                     delta_i <- taste_i$delta.dtx
                     # apply delta to 'shock.time + 1', so called 
                     # 'the time the shock appears'
                     rbind( rep(0, ncol(delta_i)), delta_i[-nrow(delta_i),]  )
                   }) )
  # all parameter names + "fixed", to be used for params.xxx
  parNames <- c(names(dynrModel), "fixed")
  # to substitute 'fixed'
  numForFixed <- length(parNames)
  
  deltaName <- names(delta)
  nState <- ncol(delta)
  # build dynr.matrixDynamics
  padyn <- dynrModel@dynamics@params.dyn[[1]]
  padyn[padyn==0] <- numForFixed
  paramsDyn <- parNames[padyn]# vector
  dim(paramsDyn) <- dim(padyn)# to matrix
  new_dynamics <- prep.matrixDynamics(
    values.dyn=dynrModel@dynamics@values.dyn[[1]],
    params.dyn=paramsDyn,
    values.exo=diag(1, nrow=nState, ncol=nState),
    params.exo=matrix("fixed", nrow=nState, ncol=nState),
    covariates=deltaName,
    isContinuousTime=FALSE)
  
  # modify dynrModel@data
  dynrModel@data$covariate.names <- deltaName
  names(delta) <- paste0("covar", 1:nState)
  dynrModel@data$covariates <- delta
  
  # build measurement
  measParLoad <- dynrModel@measurement@params.load[[1]]
  measParLoad[measParLoad==0] <- numForFixed
  paramsLoad <- parNames[measParLoad]# vector
  dim(paramsLoad) <- dim(measParLoad)# to matrix
  new_measurement <- prep.measurement(
    values.load=dynrModel@measurement@values.load[[1]],
    params.load=paramsLoad, 
    state.names=dynrModel@measurement@state.names,
    obs.names=dynrModel@measurement@obs.names
  )
  
  # build noise
  noParLat <- dynrModel@noise@params.latent[[1]]
  noParLat[noParLat==0] <- numForFixed
  paramsLatent <- parNames[noParLat]# vector
  dim(paramsLatent) <- dim(noParLat)# to matrix
  
  noParObs <- dynrModel@noise@params.observed[[1]]
  noParObs[noParObs==0] <- numForFixed
  paramsObserved <- parNames[noParObs]# vector
  dim(paramsObserved) <- dim(noParObs)# to matrix
  new_noise <- prep.noise(
    values.latent=dynrModel@noise@values.latent[[1]], 
    params.latent=paramsLatent, 
    values.observed=dynrModel@noise@values.observed[[1]], 
    params.observed=paramsObserved
  )
  
  new_dynrModel <- dynr.model(
    dynamics = new_dynamics,
    measurement = new_measurement,
    noise = new_noise,
    initial = dynrModel@initial,
    data = dynrModel@data,
    outfile = "new_taste.c")
  # cook!
  dynrDetox <- dynr.cook(new_dynrModel, verbose=FALSE, debug_flag=TRUE)
  return(dynrDetox)
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



