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
##' Detect outliers in state space models.
##'
##' @param dynrModel an object of `dynrModel' class.
##' @param dynrCook the `dynrCook' object fitted with `debug_flag=TRUE' for the `dynrModel' object. The default is NULL.
##' If the dynrCook object were not provided, or the object were cooked
##' with `debug_flag=FALSE',
##' \code{dynr.taste} will fit the dynrModel object with `debug_flag=TRUE' internally.
##' @param conf.level a numeric of confidence level that is used for
##' outliers detection tests (chi-square test and t-test). The default is 0.99.
##' @param alternative a character string specifying the alternative hypothesis of t-test,
##' must be one of ``two.sided'' (default), ``greater''  or ``less''.
##' @param debug_flag a logical. 'TRUE' for output of by-products related to t-value calculation
##'
##' @return an object of `dynrTaste' class
##' that is a list containing results lists of results from the outlier detection process.
##' Vectors of ID and measured time points are included for later use,
##' such as in dynr.taste2.
##' The values, p-values, and shock points related to 
##' `joint' chi-square, `independent' chi-square, and t statistic
##' for innovative and additive outliers are following in that order.
##' The estimated delta for innovative and additive components are in the last.
##' If \code{debug_flag} is \code{TRUE}, 
##' Q, S, s, F_inv, N, u, r would be added at the end.
##' See the reference for definition of the notations.
##'
##' The `delta_chi' list comprises magnitude of innovative (Latent) and additive (Observed) outliers, `delta.L' and `delta.O',
##' when chi-square statitics is used to detect outliers.
##' The `delta_t' list comprises magnitude of innovative (Latent) and additive (Observed) outliers, `delta.L' and `delta.O',
##' when t statitics is used to detect outliers.
##'
##' @references
##' Chow, S.-M., Hamaker, E. L., & Allaire, J. C. (2009).
##' Using innovative outliers to detectdiscrete shifts in dynamics in group-based state-space models. _Multivariate Behavioral Research_, 44, 465–496.
dynr.taste <- function(dynrModel, dynrCook=NULL, conf.level=0.99,
                       alternative=c("two.sided", "less", "greater"),
                       #outliers=c("joint", "innovative", "additive"),
                       debug_flag=FALSE) {
  if ( !is(dynrModel, 'dynrModel') ) {
    stop("dynrModel object is required.") }
  
  # check for non-regime switching
  if (dynrModel$num_regime > 2) {
    stop("This test is for non-regime switching models.")
  }
  if ( is.null(dynrCook) ||
       !("residual_cov" %in% slotNames(dynrCook)) ) {# re-cooking
    # N.B. 'optimization_flag=FALSE' produces different logLik with 'TRUE'
    dynrCook <- dynr.cook(dynrModel, verbose=FALSE, debug_flag=TRUE)
  } else {
    if ( !is(dynrCook, 'dynrCook') ) {
      stop("dynrCook object is required.") }
  }
  latName <- dynrModel$measurement$state.names
  obsName <- dynrModel$measurement$obs.names
  dimLat <- length(latName)
  dimObs <- length(obsName)
  dimTime <- length(dynrModel$data$time)
  id <- dynrModel$data$id
  IDs <- unique(id)
  nID <- length(IDs)
  
  # replace values.xx in dynrModel with coef(dynrCook)
  coefCook <- coef(dynrCook)
  tryCatch(coef(dynrModel) <- coefCook,
           error=function(e) {
             stop(paste0("The estimated initial covariance matrix
                         are not positive definite. Please re-fit (re-cook)
                         the model with a modified initial condition using 'prep.initial' function, and then re-run 'dynr.taste'.",
                         "\n", e))
           })
  
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
  # save Ninv to calculate W_t (for delta)
  Ninv <- array(NA, c(dimLat, dimLat, dimTime))
  u <- matrix(NA, dimObs, dimTime)
  tstart <- dynrModel$data$tstart
  chiLat <- numeric(dimTime)
  tpoint <- dynrModel$data$time
  
  for(j in 1:nID){
    beginTime <- tstart[j] + 1 # Add 1 because model_tstart indexes from 0 rather than 1 for interface to C
    endTime <- tstart[j+1] # Use the 0-indexed beginning position of the next person as the 1-indexed end of the current person
    
    # set endTime r and N to 0
    #  --> set r and N 0 at first when allocate
    r[,endTime] <- 0
    N[,,endTime] <- 0
    
    for(i in endTime:(beginTime+1)){
      ri <- matrix(r[,i], dimLat, 1)
      Ni <- matrix(N[,,i], dimLat, dimLat)
      Ki <- matrix(K[,,i], dimLat, dimObs)
      Li <- B - Ki %*% Lambda
      obsInfI <- matrix(F_inv[,,i], dimObs, dimObs)
      vi <- matrix(v[,i], nrow=dimObs, ncol=1)
      ui <- obsInfI %*% vi - t(Ki) %*% ri
      u[, i] <- ui
      rnew <- t_Lambda %*% ui + t(B) %*% ri
      r[,i-1] <- rnew
      Nnew <- t_Lambda %*% obsInfI %*% Lambda + t(Li) %*% Ni %*% Li
      N[,,i-1] <- Nnew
      Ninv_i <- try(solve(Nnew), silent=TRUE)
      if(class(Ninv_i) == "try-error"){Ninv_i <- MASS::ginv(Nnew)}
      Ninv[,,i-1] <- Ninv_i
      chiLat[i-1] <- t(rnew) %*% Ninv_i %*% rnew
    }
    ri <- matrix(r[,beginTime], dimLat, 1)
    Ki <- matrix(K[,,beginTime], dimLat, dimObs)
    Li <- B - Ki %*% Lambda
    obsInfI <- matrix(F_inv[,,beginTime], dimObs, dimObs)
    vi <- matrix(v[,beginTime], nrow=dimObs, ncol=1)
    ui <- obsInfI %*% vi - t(Ki) %*% ri
    u[, beginTime] <- ui
  }
  
  ################ delta estimate #############################
  # De Jong and Penzer (1988) 80p. below eq 14 about X, W and delta
  delta <- matrix(NA, dimObs+dimLat, dimTime)
  rownames(delta) <- c(obsName, latName)
  # (X',W') = I. (De Jong and Penzer, 1988)
  XW <- diag(1, dimObs+dimLat)
  X_t <- XW[1:dimObs, ]
  W_t <- XW[(dimObs+1):(dimObs+dimLat), ]
  # t-values
  t_value <- matrix(NA, dimObs+dimLat, dimTime)
  rownames(t_value) <- c(obsName, latName)
  
  ###### save by-products related to t ############
  if (debug_flag) {
    Q <- array(0, c(nrow(W_t), ncol(W_t), dimTime))
    S <- array(0, c(ncol(W_t), ncol(W_t), dimTime))
    s <- matrix(0, ncol(W_t), dimTime)
    
    for(i in 1:nID){
      beginTime <- tstart[i] + 1
      endTime <- tstart[i+1]
      
      for (j in beginTime:endTime) {
        Q_j <- W_t - K[,,j] %*% X_t
        Q[,,j] <- Q_j##
        S_j <- t(X_t) %*% F_inv[,,j] %*% X_t +
          t(Q_j) %*% N[,,j] %*% Q_j
        S[,,j] <- S_j##
        s_j <- t(X_t) %*% u[,j] + t(W_t) %*% r[,j]
        s[,j] <- s_j##
        S_j_inv <- try(solve(S_j), silent=TRUE)
        if(class(S_j_inv) == "try-error"){S_j_inv <- MASS::ginv(S_j)}
        delta[,j] <- S_j_inv %*% s_j
        t_value[,j] <- s_j / sqrt(diag(S_j))
      }
      # 0/0 at endTime
      t_value[,endTime] <- 0
    }
    
  } else {
    for(i in 1:nID){
      beginTime <- tstart[i] + 1
      endTime <- tstart[i+1]
      
      for (j in beginTime:endTime) {
        Q_j <- W_t - K[,,j] %*% X_t
        S_j <- t(X_t) %*% F_inv[,,j] %*% X_t +
          t(Q_j) %*% N[,,j] %*% Q_j
        s_j <- t(X_t) %*% u[,j] + t(W_t) %*% r[,j]
        S_j_inv <- try(solve(S_j), silent=TRUE)
        if(class(S_j_inv) == "try-error"){S_j_inv <- MASS::ginv(S_j)}
        delta[,j] <- S_j_inv %*% s_j
        t_value[,j] <- s_j / sqrt(diag(S_j))
      }
      # 0/0 at endTime
      t_value[,endTime] <- 0
    }
  }
  
  ############ chi-square test ############################
  # de Jong's method.
  chiJnt <- chiLat + chiObs
  
  chiJnt_pval <- pchisq(chiJnt, df=dimLat+dimObs, lower.tail=FALSE)
  # p-values of chi-sqaure
  chiLat_pval <- pchisq(chiLat, df=dimLat, lower.tail=FALSE)
  chiObs_pval <- pchisq(chiObs, df=dimObs, lower.tail=FALSE)
  
  # locate shock points, TRUE for significance
  chiJnt_shk <- chiJnt_pval < (1 - conf.level)
  chiLat_shk <- chiLat_pval < (1 - conf.level)
  chiObs_shk <- chiObs_pval < (1 - conf.level)
  
  ################### t-test ###############################
  alternative <- match.arg(alternative)
  t_O_calcp <- if (alternative=="two.sided") {
    function(subj) {
      2 * pt( -abs(subj), df=(ncol(subj)-dimObs) )
    }
  }	else if (alternative=="greater") {
    function(subj) {
      pt( subj, df=(ncol(subj)-dimObs), lower.tail=FALSE )
    }
  }	else {
    function(subj) {
      pt( subj, df=(ncol(subj)-dimObs), lower.tail=TRUE )
    }
  }
  
  t_L_calcp <- if (alternative=="two.sided") {
    function(subj) {
      2 * pt( -abs(subj), df=(ncol(subj)-dimLat) )
    }
  }	else if (alternative=="greater") {
    function(subj) {
      pt( subj, df=(ncol(subj)-dimLat), lower.tail=FALSE )
    }
  }	else {
    function(subj) {
      pt( subj, df=(ncol(subj)-dimLat), lower.tail=TRUE )
    }
  }
  
  # t_value for observed
  t_O <- t_value[1:dimObs, ]
  # t_value for latent
  t_L <- t_value[(dimObs+1):(dimObs+dimLat), ]
  
  t_O_pval <- matrix(NA, nrow(t_O), ncol(t_O))
  rownames(t_O_pval) <- obsName
  t_L_pval <- matrix(NA, nrow(t_L), ncol(t_L))
  rownames(t_L_pval) <- latName
  for(i in 1:nID){
    begT <- tstart[i] + 1
    endT <- tstart[i+1]
    t_O_pval[, begT:endT] <- t_O_calcp(t_O[, begT:endT])
    t_L_pval[, begT:endT] <- t_L_calcp(t_L[, begT:endT])
  }
  
  # locate shocks from t-test
  t_O_shk <- t_O_pval < (1 - conf.level)
  # rownames(t_O_shk) <- obsName
  t_L_shk <- t_L_pval < (1 - conf.level)
  # rownames(t_L_shk) <- latName
  
  ####### delta structure ######################
  rownames(delta) <- paste0("d_", c(obsName, latName))
  # delta for observed
  delta_O <- delta[1:dimObs, ]
  # delta for latent
  delta_L <- delta[(dimObs+1):(dimObs+dimLat), ]
  
  res <- list(
    tstart=tstart, id=id, time=tpoint,
    chi.jnt=chiJnt, chi.jnt.pval=chiJnt_pval, chi.jnt.outlier=chiJnt_shk,
    chi.inn=chiLat, chi.inn.pval=chiLat_pval, chi.inn.outlier=chiLat_shk,
    chi.add=chiObs, chi.add.pval=chiObs_pval, chi.add.outlier=chiObs_shk,
    t.inn=t_L, t.inn.pval=t_L_pval, t.inn.outlier=t_L_shk,
    t.add=t_O, t.add.pval=t_O_pval, t.add.outlier=t_O_shk,
    delta.inn=delta_L, delta.add=delta_O)
  if (debug_flag) {
    res <- c(res, 
             list(Q=Q, S=S, s=s, F_inv=F_inv, N=N, u=u, r=r))
  }
  
  class(res) <- "dynrTaste"
  invisible(res)
}

if (outliers=="joint") {
  # delta that will be input to 'dynr.taste2'
  row.names(delta_W) <- 1:nrow(delta_W)
  row.names(delta_X) <- 1:nrow(delta_X)
  delta_W_t <- delta_W_chi <- delta_W
  delta_X_t <- delta_X_chi <- delta_X
  delta_W_chi[!chiLat_shk,] <- 0
  delta_X_chi[!chiObs_shk,] <- 0
  delta_W_t[!t_W_shk] <- 0
  delta_X_t[!t_X_shk] <- 0
}   else if (outliers=="innovative") {
  # t shock that pass chi shock
  row.names(delta_X) <- 1:nrow(delta_X)
  row.names(delta_W) <- 1:nrow(delta_W)
  delta_W_t <- delta_W_chi <- delta_W
  delta_X_t <- delta_X_chi <- delta_X
  delta_W_chi[!chiLat_shk,] <- 0
  delta_X_chi[] <- 0
  delta_W_t[!t_W_shk] <- 0
  delta_X_t[] <- 0
} else {#outliers=="additive"
  row.names(delta_X) <- 1:nrow(delta_X)
  row.names(delta_W) <- 1:nrow(delta_W)
  # delta that will be input to 'dynr.taste2'
  delta_W_t <- delta_W_chi <- delta_W
  delta_X_t <- delta_X_chi <- delta_X
  delta_W_chi[] <- 0
  delta_X_chi[!chiObs_shk,] <- 0
  delta_W_t[] <- 0
  delta_X_t[!t_X_shk] <- 0
}

##' @description
##' The function \code{dynr.taste2{}} update the \code{dynrModel}
##' object applying outliers from the \code{dynrTaste} object,
##' or outliers from users. The function then re-cook the model.
##'
##' @param dynrModel an object of dynrModel class.
##' @param dynrCook an object of dynrCook class.
##' @param dynrTaste an object of dynrTaste class. The default is NULL.
##' @param outlierTest a character string specifying the test
##'   used to detect outliers (Chi-squared or t test),
##'   must be one of ``t'' (default), ``chi''.
##' @param newOutfile a character string for \code{outfile}
##'  argument of \code{dynr.model} function
##'  to create new \code{dynrModel} object.
##'  The default is "new_taste.c".
##'  @param cook a logical specifying whether the newly built model
##'   would be cooked by 'dynr.cook' function.
##'   The default is TRUE. When 'cook=FALSE', only the newly built model will be saved for the output.
##' @param verbose a logical specifying the verbose argument
##'  of the new cook object. The default is FALSE.
##' @param delta_L a data.frame containing user-specified latent outliers.
##' The number of rows should equal to the total time points, and the number of columns should equal to the number of latent variables.
##' @param delta_O a data.frame containing user-specified observed outliers.
##' The number of rows should equal to the total time points, and the number of columns should equal to the number of observed variables.
##'
##' @details
##' The argument \code{dynrTaste} should be the dynrTaste object
##' that is output of the \code{dynr.taste} function the argument \code{dynrModel} is applied.
##'
##' The argument \code{dynrTaste} can be \code{NULL},
##' if user-specified outliers are offered by the arguments
##' \code{delta_L} and \code{delta_O}.
##'
##' @return a list with the two arguments;
##' a new \code{dynrModel} object the outliers are applied,
##' and a \code{dynrCook} object the new \code{dynrModel} object is cooked.
dynr.taste2 <- function(dynrModel, dynrCook, dynrTaste=NULL,
                        outlierTest=c("t", "chi"),
                        newOutfile="new_taste.c",
                        cook=TRUE,
                        verbose=TRUE,
                        delta_L=NULL, delta_O=NULL) {
  if ( !is(dynrModel, 'dynrModel') ) {
    stop("dynrModel object is required.") }
  if ( !is(dynrCook, 'dynrCook') ) {
    stop("dynrCook object is required.") }

  coefx <- coef(dynrCook)
  tryCatch(coef(dynrModel) <- coefx,
           error=function(e) {
             # do nothing, and re-fit with the initial of dynrModel
             coefx <- coefx
           })

  outlierTest <- match.arg(outlierTest)
  deltaTest <- paste0("delta_", outlierTest)
  # all parameter names + "fixed", to be used for params.xxx
  parNames <- c(names(dynrModel), "fixed")
  # to substitute 'fixed'
  numForFixed <- length(parNames)

  latName <- dynrModel$measurement$state.names
  obsName <- dynrModel$measurement$obs.names
  nDeltaLat <- length(latName)
  nDeltaObs <- length(obsName)

  if ( !is.null(delta_L) ) {
    if (nDeltaLat != ncol(delta_L) || length(dynrModel@data$time) != nrow(delta_L)) stop("The number of column differs with the number of the latent variables.")
    deltaLat <- delta_L
  } else
  if ( !is.null(dynrTaste) ) {
    if ( !is(dynrTaste, 'dynrTaste') ) {
      stop("dynrTaste object is required.") }
  # combine delta through subjects
  deltaLat <- do.call("rbind",
                      lapply(dynrTaste, function(taste_i) {
                        deltaL_i <- taste_i[[deltaTest]]$delta.L
                        # apply delta to 'shock.time + 1', so called
                        # 'the time the shock appears'
                        rbind( rep(0, ncol(deltaL_i)), deltaL_i[-nrow(deltaL_i),]  )
                        #rbind( deltaL_i[-1,], rep(0, ncol(deltaL_i)) )
                        #deltaL_i
                      }) )
  } else {
    stop("Both 'dynrTaste' and 'delta_L' are NULL.")
  }

  if ( !is.null(delta_O) ) {
    if (nDeltaObs != ncol(delta_O) || length(dynrModel@data$time) != nrow(delta_L)) stop("The number of column differs with the number of the observed variables.")
    deltaObs <- delta_O
  } else
  if ( !is.null(dynrTaste) ) {
    if ( !is(dynrTaste, 'dynrTaste') ) {
      stop("dynrTaste object is required.") }
    deltaObs <- do.call("rbind",
                        lapply(dynrTaste, function(taste_i) {
                          taste_i[[deltaTest]]$delta.O
                        }) )
  } else {
    stop("Both 'dynrTast' and 'delta_O' are NULL.")
  }

  deltaLatName <- names(deltaLat)
  deltaObsName <- names(deltaObs)

  # build dynr.data
  data_org <- as.data.frame(dynrModel@data$original.data)
  id_org <- dynrModel@data$id
  time_org <- dynrModel@data$time
  obsName_org <- dynrModel@data$observed.names
  covName_org <- dynrModel@data$covariate.names
  data_all <- data.frame(id=id_org, time=time_org,
                         data_org[, c(obsName_org, covName_org), drop=FALSE],
                         deltaLat, deltaObs)
  new_data <- dynr.data(data_all, observed=obsName_org,
                        covariates=c(covName_org, deltaLatName, deltaObsName))

  # build prep.initial
  paInis <- dynrModel@initial@params.inistate[[1]]
  paInis[paInis==0] <- numForFixed
  paramInis <- parNames[paInis]# vector
  dim(paramInis) <- dim(paInis)# to matrix
  paInic <- dynrModel@initial@params.inicov[[1]]
  paInic[paInic==0] <- numForFixed
  paramInic <- parNames[paInic]# vector
  dim(paramInic) <- dim(paInic)# to matrix
  new_initial <- prep.initial(
    values.inistate=dynrModel@initial@values.inistate[[1]],
    params.inistate=paramInis,
    values.inicov=dynrModel@initial@values.inicov[[1]],
    params.inicov=paramInic)

  # build dynr.matrixDynamics
  padyn <- dynrModel@dynamics@params.dyn[[1]]
  padyn[padyn==0] <- numForFixed
  paramsDyn <- parNames[padyn]# vector
  dim(paramsDyn) <- dim(padyn)# to matrix
  if ( length(dynrModel@dynamics@values.exo) ) {
    dynValExo <- cbind( dynrModel@dynamics@values.exo[[1]],
                         diag(1, nrow=nDeltaLat, ncol=nDeltaLat) )
  } else {
    dynValExo <- diag(1, nrow=nDeltaLat, ncol=nDeltaLat)
  }
  if ( length(dynrModel@dynamics@params.exo) ) {
    dynParExo1 <- dynrModel@dynamics@params.exo[[1]]
    dynParExo1[dynParExo1==0] <- numForFixed
    dynParExo2 <- parNames[dynParExo1]# vector
    dim(dynParExo2) <- dim(dynParExo1)# to matrix
    dynParExo <- cbind(dynParExo2,
                       matrix("fixed", nrow=nDeltaLat, ncol=nDeltaLat) )
  } else {
    dynParExo <- matrix("fixed", nrow=nDeltaLat, ncol=nDeltaLat)
  }
  new_dynamics <- prep.matrixDynamics(
    values.dyn=dynrModel@dynamics@values.dyn[[1]],
    params.dyn=paramsDyn,
    values.exo=dynValExo,
    params.exo=dynParExo,
    covariates=c(dynrModel@dynamics@covariates, deltaLatName),
    isContinuousTime=FALSE)

  # # modify dynrModel@data
  # if ( is.null(dynrModel@data$covariate.names) ) {#no orginal covariates
  #   dynrModel@data$covariate.names <- c(deltaLatName, deltaObsName)
  #   deltaLatCopy <- deltaLat
  #   names(deltaLatCopy) <- paste0("covar", 1:nDeltaLat)
  #   deltaObsCopy <- deltaObs
  #   names(deltaObsCopy) <- paste0("covar", (1:nDeltaObs) + nDeltaLat)
  #   dynrModel@data$covariates <- cbind(deltaLatCopy, deltaObsCopy)
  # } else {# original covariates exist
  # nPreCovariate <- length(dynrModel@data$covariate.names)
  # dynrModel@data$covariate.names <- c(dynrModel@data$covariate.names,
  #                                     deltaLatName, deltaObsName)
  # names(deltaLat) <- paste0("covar",
  #                           (1:nDeltaLat) + nPreCovariate)
  # names(deltaObs) <- paste0("covar",
  #                           (1:nDeltaObs) + nPreCovariate + nDeltaLat)
  # dynrModel@data$covariates <- cbind(dynrModel@data$covariates,
  #                                    deltaLat, deltaObs)
  # }

  # build measurement
  measParLoad <- dynrModel@measurement@params.load[[1]]
  measParLoad[measParLoad==0] <- numForFixed
  paramsLoad <- parNames[measParLoad]# vector
  dim(paramsLoad) <- dim(measParLoad)# to matrix
  if ( length(dynrModel@measurement@values.exo) ) {
    measValExo <- cbind( dynrModel@measurement@values.exo[[1]],
                         diag(1, nrow=nDeltaObs, ncol=nDeltaObs) )
  } else {
    measValExo <- diag(1, nrow=nDeltaObs, ncol=nDeltaObs)
  }
  if ( length(dynrModel@measurement@params.exo) ) {
    measParExo1 <- dynrModel@measurement@params.exo[[1]]
    measParExo1[measParExo1==0] <- numForFixed
    measParExo2 <- parNames[measParExo1]# vector
    dim(measParExo2) <- dim(measParExo1)# to matrix
    measParExo <- cbind(measParExo2,
                        matrix("fixed", nrow=nDeltaObs, ncol=nDeltaObs) )
  } else {
    measParExo <- matrix("fixed", nrow=nDeltaObs, ncol=nDeltaObs)
  }
  new_measurement <- prep.measurement(
    values.load=dynrModel@measurement@values.load[[1]],
    params.load=paramsLoad,
    values.exo=measValExo,
    params.exo=measParExo,
    exo.names=c(dynrModel@measurement@exo.names, deltaObsName),
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
    initial = new_initial,
    data = new_data,#dynrModel@data,
    outfile = newOutfile)
  # cook?
  if (cook) {
    new_dynrCook <- dynr.cook(new_dynrModel, verbose=verbose, debug_flag=TRUE)
    list(new_dynrModel=new_dynrModel,
         new_dynrCook=new_dynrCook)
  } else {
    list(new_dynrModel=new_dynrModel)
  }
}

computeJacobian <- function(cookDebug, jacobian, latName, params, time){
  envList <- as.list(params)
  for(i in 1:length(latName)) envList[[latName[i]]] <- cookDebug$eta_smooth_final[i, time]

  J <- matrix(NA, length(latName), length(latName), dimnames=list(latName, latName))
  for(i in 1:(length(latName))^2){
    cj <- as.character(jacobian[[i]])
    rc <- strsplit(cj[2], ' ~ ')[[1]]
    J[rc[1], rc[2]] <- eval(parse(text=cj[3]), envList)
  }
  return(J)
}

#example use from demo/NonlinearODE.R
#computeJacobian(res, dynm$jacobian[[1]], model$measurement$state.names, coef(res), 50)
#cf t=1 vs t=50

##' @param numSubjDemo The number of subjects to be randomly selected for plotting.
##' @param idtoPlot Values of the ID variable to plot.
##' @param names.state (optional) The names of the states to be plotted, which should be a subset of the state.names slot of the measurement slot of dynrModel.
##' @param names.observed (optional) The names of the observed variables to be plotted, which should be a subset of the obs.names slot of the measurement slot of dynrModel.
autoplot.dynrTaste(dynrtaste, numSubjDemo=2, idtoPlot=c(),
                   names.state, names.observed, ...) {

}

rbind(id=taste_shk$`1`$taste[["id"]],
           time=taste_shk$`1`$taste[["time"]],
           taste_shk$`1`$delta_chi$delta.L)

delta_chi_L <- function(tasteOut) {
  shocks <- lapply(tasteOut, function(tasteout_i) {
    data.frame(id=tasteout_i$taste[["id"]],
               time=tasteout_i$taste[["time"]],
               tasteout_i$delta_chi$delta.L)
  } )
  do.call(rbind.data.frame, )
}

shockSignature <- function(dynrModel, dynrCook,
                           T=20, shockTime=5, W=NULL) {
  if ( T <= shockTime ) {
    stop("'shockTime' must be smaller than 'T'.")
    }
  coef(dynrModel) <- coef(dynrCook)
  B <- dynrModel$dynamics@values.dyn[[1]]
  L <- dynrModel$measurement@values.load[[1]]
  nameLat <- dynrModel$measurement@state.names
  nameObs <- dynrModel$measurement@obs.names
  dimL <- dim(L)
  if ( is.null(W) ) {
    W <- diag(1, dimL[2])
  } else {
    if ( !all.equal(rep(dimL[2], 2) , dim(W)) ) {
      stop("Please check the dimension of W.")
    }
  }
  # if ( is.null(X) ) {
  #   X <- diag(1, dimL[1])
  # } else {
  #   if ( !all.equal(rep(dimL[1], 2) , dim(X)) ) {
  #     stop("Please check the dimension of X.")
  #   }
  # }
  # D_A <- array(NA, c(dimL, T),
  #              dimnames=list(nameObs, nameLat, 1:T))
  # D <- L %*% B
  # D_A[,,1] <- D
  # for (i in 2:time) {
  #   D <- D %*% B
  #   D_A[,,i] <- D
  # }
  dimB <- dim(B)
  B_all <- array(NA, c(dimB, (T-shockTime)))
  B_all[,,1] <- diag(1, dimB[1], dimB[2])
  BB <- B
  B_all[,,2] <- BB
  for ( t in 3:(T-shockTime) ) {
    BB <- BB %*% B
    B_all[,,t] <- BB
  }
  D_all <- array(NA, c(dimL, (T-shockTime)))
  for ( t in 1:(T-shockTime) ) {
    D_all[,,t] <- L %*% B_all[,,t] %*% W
  }
  invisible(D_all)
}

# panaModelc <- panaModel
# coef(panaModelc) <- coef(panaCook)
# B <- panaModelc$dynamics@values.dyn[[1]]
# L <- panaModelc$measurement@values.load[[1]]
# B
# L
# L %*% matrix(c(1, 0, 0, 1), ncol=2)
# nameLat <- panaModel$measurement@state.names
# aa <- shockSignature(panaModel, panaCook)
# stime <- 5
# ttime <- 20
# bb <- data.frame(
# rbind(
# cbind(matrix(0, nrow=stime, ncol=2), 1:stime),
# cbind(t(aa[1,,]), (stime+1):ttime)
# ) )
#  before <- data.frame(value=rep(0, stime), T=1:stime)
# # bb1 <- data.frame(t(aa[4,,]), T=(stime+1):ttime)
# names(bb) <- c(nameLat, "T")
# ttext <- min(t(aa[4,,])) + max(t(aa[4,,])) / 2
# # aa[1,3,]
# aaa <- reshape2::melt(bb, id='T')
# aaa
# ggplot(aaa, aes(x=T, y=value)) +
#   geom_line(aes(colour=variable)) +
#   theme_bw() + theme_minimal() +
#   geom_line(data=before, aes(x=T, y=value)) +
#   labs(y="Shock Coefficient", colour="States") +
#   geom_vline(xintercept=stime, linetype="dotted") +
#   theme(panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         axis.line=element_line(colour="grey60"),
#         axis.ticks=element_blank()) +
#   geom_text(aes(x=5, y=ttext, label="Shock"), colour="grey40", angle=90, vjust = -0.5) +
#   scale_x_continuous(breaks=c(5, 30), limits=c(0, 30))
# plotArgs <- list(ncol=1, nrow=4, align="v")
# latchi <- do.call(ggpubr::ggarrange, c(list(chipana, chipana_I,
#                                             tPA, tNA), plotArgs))
# ggpubr::annotate_figure(latchi,
#                         bottom=ggpubr::text_grob("time", color="black", size=11))
