# Author: Michael D. Hunter
# Date: 2017-09-14
# Filename: dynrTaste.R
# Purpose: Compute shocks and chi-squared diagnostics

##' Detect outliers in state space models.
##'
##' Compute shocks and chi-squared diagnostics following
##' Chow, Hamaker, and Allaire (2009).  Using Innovative Outliers to
##' Detect Discrete Shifts in Dynamics in Group-Based State-Space Models
##' 
##' @param dynrModel an object of `dynrModel' class.
##' @param dynrCook the `dynrCook' object fitted with `debug_flag=TRUE' for the `dynrModel' object. The default is NULL.
##' If the dynrCook object were not provided, or the object were cooked
##' with `debug_flag=FALSE',
##' \code{dynr.taste} will fit the dynrModel object with `debug_flag=TRUE' internally.
##' @param which.state a character vector of the names of latent variables. 
##' The outlier detection process will be applied only to the chosen variable. 
##' If the argument is NA, all the latent variables will be excluded in the outlier detection process. 
##' If the argument is missing (defalut), all the latent variables will be chosen.
##' @param which.obs a character vector of the names of measured or observed variables.
##' The outlier detection process will be applied only to the chosen variable. 
##' If the argument is NA, all the measured variables will be excluded in the outlier detection process. 
##' If the argument is missing (defalut), all the measured variables will be chosen.
##' @param conf.level a numeric of confidence level that is used for
##' outliers detection tests (chi-square test and t-test). The default is 0.99.
##' @param alternative a character string specifying the alternative hypothesis of t-test,
##' must be one of ``two.sided'' (default), ``greater''  or ``less''.
##' @param debug_flag a logical. 'TRUE' for output of by-products related to t-value calculation
##'
##' @return an object of `dynrTaste' class
##' that is a list containing lists of results from the outlier detection process.
##' Vectors of ID and measured time points are included for later use,
##' such as in dynr.taste2.
##' The values, p-values, and shock points related to 
##' `joint' chi-square, `independent' chi-square, and t statistic
##' for innovative and additive outliers are following in that order.
##' The estimated delta for innovative and additive components are in the last.
##' If \code{debug_flag} is \code{TRUE}, 
##' The by-products of the Kalman filter and smoother (Q, S, s, F_inv, N, u, r) would be added at the end.
##' See the reference for definition of the notations.

##' The t statistic (estimate of an outlier divided by standard error of the outlier) of the last time point is NA, 
##' because the Kalman smoothing process starts with setting r and N to zero for the last time point 
##' (core elements of calculating estimates and the standard errors of outliers) 
##' that lead to 0/0 of the t statistic of the last time point.
##' For the time-varing models, more NAs would appear at the end of times because the Kalman smoother needs more time points to obtain all elements of r nad N from limited number of observed variables in the model. 
##'
##' The `delta_chi' list comprises magnitude of innovative (Latent) and additive (Observed) outliers, `delta.L' and `delta.O',
##' when chi-square statitics is used to detect outliers.
##' The `delta_t' list comprises magnitude of innovative (Latent) and additive (Observed) outliers, `delta.L' and `delta.O',
##' when t statitics is used to detect outliers.
##'
##' @references
##' Chow, S.-M., Hamaker, E. L., & Allaire, J. C. (2009).
##' Using innovative outliers to detect discrete shifts in dynamics in group-based state-space models. _Multivariate Behavioral Research_, 44, 465-496.
##' 
##' @examples
##' \dontrun{
##' # See the demo for outlier detection, OutlierDetection.R
##' dynrCook <- dynr.cook(dynrModel)
##' dynrTaste <- dynr.taste(dynrModel, dynrCook)
##' 
##' # Detect outliers related to 'eta1' out of, say, three latent
##' # variables c("eta1", "eta2", "eta3"), and all measured variables.
##' dynrTaste <- dynr.taste(dynrModel, dynrCook, which.state=c("eta1"))
##' }
dynr.taste <- function(dynrModel, dynrCook=NULL,
                       which.state, which.obs,
                       conf.level=0.99,
                       alternative=c("two.sided", "less", "greater"),
                       debug_flag=FALSE) {
  if ( !inherits(dynrModel, 'dynrModel') ) {
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
    if ( !inherits(dynrCook, 'dynrCook') ) {
      stop("dynrCook object is required.") }
  }
  lat_name <- dynrModel$measurement$state.names
  obs_name <- dynrModel$measurement$obs.names
  dimLat <- length(lat_name)
  dimObs <- length(obs_name)
  dimTime <- length(dynrModel$data$time)
  id <- dynrModel$data$id
  IDs <- unique(id)
  nID <- length(IDs)
  
  # choose latent variables
  if( missing(which.state) ) {
    W_t_exclude <- rep(FALSE, dimLat)
    
  } else if ( any(is.na(which.state)) ) {
    W_t_exclude <- rep(TRUE, dimLat)
    
  } else {
    which.state <- unique(which.state)
    lcv <- which.state %in% lat_name
    if ( !all(lcv) ) {
      stop( paste("Cannot recognize the latent variable(s):",
                  paste(which.state[!lcv], collapse=", ")),
            call.=FALSE)
    }
    W_t_exclude <- !(lat_name %in% which.state)
  }
  
  # choose measured variables
  if( missing(which.obs) ) {
    X_t_exclude <- rep(FALSE, dimObs)
    
  } else if ( any(is.na(which.obs)) ) {
    X_t_exclude <- rep(TRUE, dimObs)
    
  } else {
    which.obs <- unique(which.obs)
    ocv <- which.obs %in% obs_name
    if ( !all(ocv) ) {
      stop( paste("Cannot recognize the measured variable(s):",
                  paste(which.obs[!ocv], collapse=", ")),
            call.=FALSE)
    }
    X_t_exclude <- !(obs_name %in% which.obs)
  }
  
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
    F_inv_i <- matrix(F_inv[,,i], dimObs, dimObs)
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
  time <- dynrModel$data$time
  
  if ( inherits(dynrModel$dynamics, 'dynrDynamicsMatrix') ) {
    B <- dynrModel$dynamics$values.dyn[[1]]
    
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
        if("try-error" %in% class(Ninv_i)){Ninv_i <- MASS::ginv(Nnew)}
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
  }
  if ( inherits(dynrModel$dynamics, 'dynrDynamicsFormula') ) {
    # if ( is.null(dynrDynamics) ) {
    #   stop("The 'dynrDynamicsFormula' object MUST be provided, which was created by 'prep.formulaDynamics'.")
    # }
    # if ( !inherits(dynrDynamics, 'dynrDynamicsFormula') ) {
    #   stop("dynrDynamics should be a 'dynrDynamicsFormula' object.") }
    
    eta <- dynrCook$eta_smooth_final
    rownames(eta) <- lat_name
    ## replace 'param[]' to parnames in jacobian formula
    # formula to character
    to_character <- function(x) {
      Reduce( paste, deparse(x) )
    }
    mdj <- dynrModel@dynamics@jacobian[[1]]
    parnames <- dynrModel@dynamics@paramnames
    n_params <- length(parnames)
    pch <- paste0("param\\[", 0:(n_params-1), "\\]")
    jacobian <- lapply(mdj, function(ja) {
      jac <- to_character(ja)
      for (i in 1:n_params) {
        jac <- gsub(pch[i], parnames[i], jac)
      }
      as.formula(jac)
    })
    ####
    #jacobian <- dynrDynamics$jacobian[[1]]
    params <- coefCook[parnames]
    pp <- matrix(params, nrow=n_params, ncol=dimTime)
    rownames(pp) <- parnames
    # check the order of differential equations
    jac_n1 <- lapply(jacobian, function(jj) {
      jj[[2]]
    })
    jac_n2 <- vector("list", length=dimLat*dimLat)
    ii <- 1
    for (rr in lat_name) {
      for (cc in lat_name) {
        jac_n2[[ii]] <- as.formula(paste(c(rr, cc), collapse=" ~ "))
        ii <- ii + 1
      }}
    stopifnot(all.equal(jac_n1, jac_n2))
    
    ddenv <- as.data.frame( t(rbind(pp, eta)) )
    jac_l <- lapply(jacobian, function(jac) {
      eval(jac[[3]], ddenv)
    })
    jac_m <- do.call(rbind, jac_l)
    Bjac <- array(NA, dim=c(dimLat, dimLat, dimTime))
    for (rr in 1:dimTime) {
      Bjac[,,rr] <- matrix(jac_m[,rr],
                           nrow=dimLat, ncol=dimLat, byrow=TRUE)
    } 
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
        B <- Bjac[,,i]
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
        if("try-error" %in% class(Ninv_i)){Ninv_i <- MASS::ginv(Nnew)}
        Ninv[,,i-1] <- Ninv_i
        chiLat[i-1] <- t(rnew) %*% Ninv_i %*% rnew
      }
      ri <- matrix(r[,beginTime], dimLat, 1)
      Ki <- matrix(K[,,beginTime], dimLat, dimObs)
      Li <- Bjac[,,beginTime] - Ki %*% Lambda
      obsInfI <- matrix(F_inv[,,beginTime], dimObs, dimObs)
      vi <- matrix(v[,beginTime], nrow=dimObs, ncol=1)
      ui <- obsInfI %*% vi - t(Ki) %*% ri
      u[, beginTime] <- ui
    }
  }
  
  ################ delta estimate #############################
  # De Jong and Penzer (1988) 80p. below eq 14 about X, W and delta
  delta <- matrix(NA, dimObs+dimLat, dimTime)
  rownames(delta) <- c(obs_name, lat_name)
  # (X',W') = I. (De Jong and Penzer, 1988)
  XW <- diag(1, dimObs+dimLat)
  X_t <- XW[1:dimObs, , drop=FALSE]
  W_t <- XW[(dimObs+1):(dimObs+dimLat), , drop=FALSE]
  W_t[W_t_exclude, ] <- 0
  X_t[X_t_exclude, ] <- 0
  
  # t-values
  t_value <- matrix(NA, dimObs+dimLat, dimTime)
  rownames(t_value) <- c(obs_name, lat_name)
  
  ###### save by-products related to t ############
  if (debug_flag) {
    Q <- array(0, c(nrow(W_t), ncol(W_t), dimTime))
    S <- array(0, c(ncol(W_t), ncol(W_t), dimTime))
    s <- matrix(0, ncol(W_t), dimTime)
    
    for(i in 1:nID){
      beginTime <- tstart[i] + 1
      endTime <- tstart[i+1]
      
      for (j in beginTime:endTime) {
        Q_j <- W_t - matrix(K[,,j], dimLat, dimObs) %*% X_t
        Q[,,j] <- Q_j##
        S_j <- t(X_t) %*% matrix(F_inv[,,j], dimObs, dimObs) %*% X_t +
          t(Q_j) %*% matrix(N[,,j], dimLat, dimLat) %*% Q_j
        S[,,j] <- S_j##
        s_j <- t(X_t) %*% u[,j, drop=FALSE] + t(W_t) %*% r[,j, drop=FALSE]
        s[,j] <- s_j##
        S_j_inv <- try(solve(S_j), silent=TRUE)
        if("try-error" %in% class(S_j_inv)){S_j_inv <- MASS::ginv(S_j)}
        delta[,j] <- S_j_inv %*% s_j
        t_value[,j] <- s_j / sqrt(diag(S_j))
      }
      # 0/0 at endTime
      # t_value[(dimObs+1):(dimObs+dimLat), endTime] <- 0
    }
    
  } else {
    for(i in 1:nID){
      beginTime <- tstart[i] + 1
      endTime <- tstart[i+1]
      
      for (j in beginTime:endTime) {
        Q_j <- W_t - matrix(K[,,j], dimLat, dimObs) %*% X_t
        S_j <- t(X_t) %*% matrix(F_inv[,,j], dimObs, dimObs) %*% X_t +
          t(Q_j) %*% matrix(N[,,j], dimLat, dimLat) %*% Q_j
        s_j <- t(X_t) %*% u[,j, drop=FALSE] + t(W_t) %*% r[,j, drop=FALSE]
        S_j_inv <- try(solve(S_j), silent=TRUE)
        if("try-error" %in% class(S_j_inv)){S_j_inv <- MASS::ginv(S_j)}
        delta[,j] <- S_j_inv %*% s_j
        t_value[,j] <- s_j / sqrt(diag(S_j))
      }
      # 0/0 at endTime
      # t_value[(dimObs+1):(dimObs+dimLat), endTime] <- 0
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
  t_O <- t_value[1:dimObs, , drop=FALSE]
  # t_value for latent
  t_L <- t_value[(dimObs+1):(dimObs+dimLat), , drop=FALSE]
  
  t_O_pval <- matrix(NA, nrow(t_O), ncol(t_O))
  rownames(t_O_pval) <- obs_name
  t_L_pval <- matrix(NA, nrow(t_L), ncol(t_L))
  rownames(t_L_pval) <- lat_name
  for(i in 1:nID){
    begT <- tstart[i] + 1
    endT <- tstart[i+1]
    t_O_pval[, begT:endT] <- t_O_calcp(t_O[, begT:endT, drop=FALSE])
    t_L_pval[, begT:endT] <- t_L_calcp(t_L[, begT:endT, drop=FALSE])
  }
  
  # locate shocks from t-test
  t_O_shk <- t_O_pval < (1 - conf.level)
  t_L_shk <- t_L_pval < (1 - conf.level)
  
  ####### delta structure ######################
  #rownames(delta) <- paste0("d_", c(obs_name, lat_name))
  rownames(delta) <- c(obs_name, lat_name)
  # delta for observed
  delta_O <- delta[1:dimObs, , drop=FALSE]
  # delta for latent
  delta_L <- delta[(dimObs+1):(dimObs+dimLat), , drop=FALSE]
  
  res <- list(
    conf.level=conf.level, tstart=tstart, id=id, time=time,
    chi.jnt=chiJnt, chi.jnt.pval=chiJnt_pval, chi.jnt.shock=chiJnt_shk,
    chi.inn=chiLat, chi.inn.pval=chiLat_pval, chi.inn.shock=chiLat_shk,
    chi.add=chiObs, chi.add.pval=chiObs_pval, chi.add.shock=chiObs_shk,
    t.inn=t_L, t.inn.pval=t_L_pval, t.inn.shock=t_L_shk,
    t.add=t_O, t.add.pval=t_O_pval, t.add.shock=t_O_shk,
    delta.inn=delta_L, delta.add=delta_O)
  #if ( !is.null(dynrDynamics) ) {
  #  res <- c(res, dynrDynamics=dynrDynamics)
  #}
  if (debug_flag) {
    res <- c(res, 
             list(Q=Q, S=S, s=s, F_inv=F_inv, N=N, u=u, r=r))
  }
  
  class(res) <- "dynrTaste"
  invisible(res)
}


##' Re-fit state-space model using the estimated outliers.
##' 
##' The function \code{dynr.taste2{}} update the \code{dynrModel}
##' object applying outliers from the \code{dynrTaste} object,
##' or outliers from users. The function then re-cook the model.
##'
##' @param dynrModel an object of dynrModel class.
##' @param dynrCook an object of dynrCook class.
##' @param dynrTaste an object of dynrTaste class. The default is NULL.
##' @param delta_inn a character string for a method detecting `inn'ovative outliers,
##' which must be one of ``t'' (default), ``ind'', ``jnt'' or ``null''.
##' According to the method,
##' corresponding delta estimates (magnitude of estimated outliers) 
##' will be included in the new dynrModel in output.
##' `t' represents the t statistic, `ind' represents the independent chi-square statistic,
##' `jnt' represents the joint chi-square statistic.
##' If no outliers are assumed, ``null'' can be used. 
##' @param delta_add a character string for a method detecting `add'itive outliers,
##' which must be one of ``t'' (default), ``ind'', ``jnt'' or ``null''.
##' According to the method,
##' corresponding delta estimates will be included in the new dynrModel.
##' @param delta_L a data.frame containing user-specified latent outliers.
##' The delta estimates from \code{dynrTaste} will be ignored.
##' The number of rows should equal to the total time points, and the number of columns should equal to the number of latent variables.
##' @param delta_O a data.frame containing user-specified observed outliers.
##' The delta estimates from \code{dynrTaste}, and 
##' arguments of \code{delta_inn} and \code{delta_add} will be ignored.
##' The number of rows should equal to the total time points, and the number of columns should equal to the number of observed variables.
##' @param cook a logical specifying whether the newly built model
##'   would be cooked by 'dynr.cook' function.
##'   The default is TRUE. When 'cook=FALSE', only the newly built model will be saved for the output.
##' @param verbose a logical specifying the verbose argument
##'  of the new cook object. The default is FALSE.
##' @param newOutfile a character string for \code{outfile}
##'  argument of \code{dynr.model} function
##'  to create new \code{dynrModel} object.
##'  The default is "new_taste.c".
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
##' 
##' @examples
##' \dontrun{
##' # See the demo for outlier detection, OutlierDetection.R
##' dynrCook <- dynr.cook(dynrModel)
##' dynrTaste <- dynr.taste(dynrModel, dynrCook)
##' 
##' # Detect outliers related to 'eta1' out of, say, three latent
##' # variables c("eta1", "eta2", "eta3"), and all measured variables.
##' taste2 <- dynr.taste2(dynrModel, dynrCook, dynrTaste)
##' }
dynr.taste2 <- function(dynrModel, dynrCook, dynrTaste,
                        delta_inn=c("t", "ind", "jnt", "null"),
                        delta_add=c("t", "ind", "jnt", "null"),
                        delta_L=NULL, delta_O=NULL,
                        cook=TRUE, verbose=FALSE,
                        newOutfile="new_taste.c") {
  if ( !inherits(dynrModel, 'dynrModel') ) {
    stop("dynrModel object is required.")
  }
  if ( !inherits(dynrCook, 'dynrCook') ) {
    stop("dynrCook object is required.") 
  }
  if ( !inherits(dynrTaste, "dynrTaste") ) {
    stop("dynrTaste object is required.") 
  }
  
  id <- dynrModel$data$id
  id_unq <- unique(id) 
  id_n <- length(id_unq)
  tstart <- dynrModel$data$tstart
  time <- dynrModel$data$time
  lat_name <- dynrModel$measurement$state.names
  obs_name <- dynrModel$measurement$obs.names
  lat_n <- length(lat_name)
  obs_n <- length(obs_name)
  
  # coefx <- coef(dynrCook)
  # tryCatch(coef(dynrModel) <- coefx,
  #          error=function(e) {
  #            # do nothing, and re-fit with the initial of dynrModel
  #            coefx <- coefx
  #          })

  deltaLat <- dynrTaste$delta.inn
  delta_inn <- match.arg(delta_inn)
  if (delta_inn == "t") {
    deltaLat[!dynrTaste$t.inn.shock] <- 0
    
  } else if (delta_inn == "ind") {
    deltaLat[!dynrTaste$chi.inn.shock] <- 0
        
  } else if (delta_inn == "jnt") {
    deltaLat[!dynrTaste$chi.jnt.shock] <- 0
  
  } else {# "null"
    deltaLat[] <- 0
  }
  # apply deltaLat to 'shock.time + 1', so called
  # 'the time the shock appears'
  for(j in 1:id_n) {
    begT <- tstart[j] + 1
    endT <- tstart[j+1]
    deltaLat[, begT:endT] <- 
      cbind( rep(0, lat_n), deltaLat[, begT:(endT-1)] )
  }
  deltaLat <- as.data.frame( t(deltaLat) )
  
  deltaObs <- dynrTaste$delta.add
  delta_add <- match.arg(delta_add)
  if (delta_add == "t") {
    deltaObs[!dynrTaste$t.add.shock] <- 0
    
  } else if (delta_inn == "ind") {
    deltaObs[!dynrTaste$chi.add.shock] <- 0
    
  } else if (delta_inn == "jnt") {
    deltaObs[!dynrTaste$chi.jnt.shock] <- 0
    
  } else {# "null"
    deltaObs[] <- 0
  }
  deltaObs <- as.data.frame( t(deltaObs) )
  
  if ( !is.null(delta_L) ) {
    if ( !is.data.frame(delta_L) ) {
      stop("'delta_L' must be a data.frame.") }
    if ( lat_n != ncol(delta_L) ) {
      stop("The number of columns of 'delta_L' differs with the number of the latent variables.") }
    if ( length(time) != nrow(delta_L) ) {
      stop("The number of rows of 'delta_L' differs with the number of the total time points.") }
    deltaLat <- delta_L
  }
  dlname <- names(deltaLat)
  names(deltaLat) <- paste0("d_", dlname)
  
  if ( !is.null(delta_O) ) {
    if ( !is.data.frame(delta_O) ) {
      stop("'delta_O' must be a data.frame. ") }
    if ( obs_n != ncol(delta_O) ) {
      stop("The number of columns of 'delta_O' differs with the number of the measured variables.") }
    if ( length(time) != nrow(delta_O) ) {
      stop("The number of rows of 'delta_L' differs with the number of the total time points.") }
    deltaObs <- delta_O
  }
  doname <- names(deltaObs)
  names(deltaObs) <- paste0("d_", doname)

  delta_lat_name <- names(deltaLat)
  delta_obs_name <- names(deltaObs)
  
  # all parameter names + "fixed", to be used for params.xxx
  par_name <- c(names(dynrModel), "fixed")
  # to substitute 'fixed'
  fixed_pos <- length(par_name)

  if ( inherits(dynrModel$dynamics, 'dynrDynamicsMatrix') ) {  
  # build dynr.matrixDynamics
  padyn <- dynrModel@dynamics@params.dyn[[1]]
  padyn[padyn==0] <- fixed_pos
  paramsDyn <- par_name[padyn]# vector
  dim(paramsDyn) <- dim(padyn)# to matrix
  if ( length(dynrModel@dynamics@values.exo) ) {
    dynValExo <- cbind( dynrModel@dynamics@values.exo[[1]],
                         diag(1, nrow=lat_n, ncol=lat_n) )
  } else {
    dynValExo <- diag(1, nrow=lat_n, ncol=lat_n)
  }
  if ( length(dynrModel@dynamics@params.exo) ) {
    dynParExo1 <- dynrModel@dynamics@params.exo[[1]]
    dynParExo1[dynParExo1==0] <- fixed_pos
    dynParExo2 <- par_name[dynParExo1]# vector
    dim(dynParExo2) <- dim(dynParExo1)# to matrix
    dynParExo <- cbind(dynParExo2,
                       matrix("fixed", nrow=lat_n, ncol=lat_n) )
  } else {
    dynParExo <- matrix("fixed", nrow=lat_n, ncol=lat_n)
  }
  new_dynamics <- prep.matrixDynamics(
    values.dyn=dynrModel@dynamics@values.dyn[[1]],
    params.dyn=paramsDyn,
    values.exo=dynValExo,
    params.exo=dynParExo,
    covariates=c(dynrModel@dynamics@covariates, delta_lat_name),
    isContinuousTime=FALSE)
  }
  
  if ( inherits(dynrModel$dynamics, 'dynrDynamicsFormula') ) {
    ## copied the function from 'statnet.common' package
    append_rhs <- function(object, newterms, keep.onesided = FALSE) {
      #if (inherits(newterms, "formula")) 
      #  newterms <- list_rhs.formula(newterms)
      for (i in seq_along(newterms)) {
        newterm <- newterms[[i]]
        #termsign <- if (NVL(attr(newterms, "sign")[i], +1) > 0) 
        #  "+"
        #else "-"
        termsign <- "+"
        if (length(object) == 3) 
          object[[3]] <- call(termsign, object[[3]], newterm)
        else if (keep.onesided) 
          object[[2]] <- call(termsign, object[[2]], newterm)
        else object[[3]] <- if (termsign == "+") 
          newterm
        else call(termsign, newterm)
      }
      object
    }
    
    #new_dynamics <- dynrTaste$dynrDynamics
    #old_dynamics <- dynrModel@dynamics
    to_character <- function(x) {
      Reduce( paste, deparse(x) )
    }
    parnames <- dynrModel@dynamics@paramnames
    n_params <- length(parnames)
    pch <- paste0("param\\[", 0:(n_params-1), "\\]")
    
    # mdj <- dynrModel@dynamics@jacobian[[1]]
    # new_jacobian <- lapply(mdj, function(ja) {
    #   jac <- to_character(ja)
    #   for (i in 1:n_params) {
    #     jac <- sub(pch[i], parnames[i], jac)
    #   }
    #   as.formula(jac)
    # })
    
    mdf <- dynrModel@dynamics@formula[[1]]
    new_formula <- lapply(mdf, function(fo) {
      foc <- to_character(fo)
      for (i in 1:n_params) {
        foc <- gsub(pch[i], parnames[i], foc)
      }
      as.formula(foc)
    })
    new_formula <- list(
      lapply(new_formula, function(ff) {
      del <- paste0("d_", ff[[2]])
      append_rhs(ff, list(as.name(del)))
    }) )
    new_dynamics <- prep.formulaDynamics(
      formula=new_formula, 
      startval=dynrModel@dynamics@startval,
      isContinuousTime=dynrModel@dynamics@isContinuousTime
    )
  }
  
  #coefx <- coef(dynrCook)
  #tryCatch(coef(dynrModel) <- coefx,
  #         error=function(e) {
  #           # do nothing, and re-fit with the initial of dynrModel
  #           coefx <- coefx
  #         })
  
  # build dynr.data
  data_org <- as.data.frame(dynrModel@data$original.data)
  obs_name_org <- dynrModel@data$observed.names
  cov_name_org <- dynrModel@data$covariate.names
  data_all <- data.frame(id=id, time=time,
                         data_org[, c(obs_name_org, cov_name_org), drop=FALSE],
                         deltaLat, deltaObs)
  new_data <- dynr.data(data_all, observed=obs_name_org,
                        covariates=c(cov_name_org, delta_lat_name, delta_obs_name))
  
  # build prep.initial
  paInis <- dynrModel@initial@params.inistate[[1]]
  paInis[paInis==0] <- fixed_pos
  paramInis <- par_name[paInis]# vector
  dim(paramInis) <- dim(paInis)# to matrix
  paInic <- dynrModel@initial@params.inicov[[1]]
  paInic[paInic==0] <- fixed_pos
  paramInic <- par_name[paInic]# vector
  dim(paramInic) <- dim(paInic)# to matrix
  new_initial <- prep.initial(
    values.inistate=dynrModel@initial@values.inistate[[1]],
    params.inistate=paramInis,
    values.inicov=dynrModel@initial@values.inicov[[1]],
    params.inicov=paramInic)
  
  # build measurement
  measParLoad <- dynrModel@measurement@params.load[[1]]
  measParLoad[measParLoad==0] <- fixed_pos
  paramsLoad <- par_name[measParLoad]# vector
  dim(paramsLoad) <- dim(measParLoad)# to matrix
  if ( length(dynrModel@measurement@values.exo) ) {
    measValExo <- cbind( dynrModel@measurement@values.exo[[1]],
                         diag(1, nrow=obs_n, ncol=obs_n) )
  } else {
    measValExo <- diag(1, nrow=obs_n, ncol=obs_n)
  }
  if ( length(dynrModel@measurement@params.exo) ) {
    measParExo1 <- dynrModel@measurement@params.exo[[1]]
    measParExo1[measParExo1==0] <- fixed_pos
    measParExo2 <- par_name[measParExo1]# vector
    dim(measParExo2) <- dim(measParExo1)# to matrix
    measParExo <- cbind(measParExo2,
                        matrix("fixed", nrow=obs_n, ncol=obs_n) )
  } else {
    measParExo <- matrix("fixed", nrow=obs_n, ncol=obs_n)
  }
  new_measurement <- prep.measurement(
    values.load=dynrModel@measurement@values.load[[1]],
    params.load=paramsLoad,
    values.exo=measValExo,
    params.exo=measParExo,
    exo.names=c(dynrModel@measurement@exo.names, delta_obs_name),
    state.names=dynrModel@measurement@state.names,
    obs.names=dynrModel@measurement@obs.names
  )

  # build noise
  noParLat <- dynrModel@noise@params.latent[[1]]
  noParLat[noParLat==0] <- fixed_pos
  paramsLatent <- par_name[noParLat]# vector
  dim(paramsLatent) <- dim(noParLat)# to matrix

  noParObs <- dynrModel@noise@params.observed[[1]]
  noParObs[noParObs==0] <- fixed_pos
  paramsObserved <- par_name[noParObs]# vector
  dim(paramsObserved) <- dim(noParObs)# to matrix
  new_noise <- prep.noise(
    values.latent=dynrModel@noise@values.latent[[1]],
    params.latent=paramsLatent,
    values.observed=dynrModel@noise@values.observed[[1]],
    params.observed=paramsObserved
  )

  dynrModel_new <- dynr.model(
    dynamics = new_dynamics,
    measurement = new_measurement,
    noise = new_noise,
    initial = new_initial,
    data = new_data,#dynrModel@data,
    outfile = newOutfile)
  # cook?
  if (cook) {
    dynrCook_new <- dynr.cook(dynrModel_new, verbose=verbose, debug_flag=TRUE)
    list(dynrModel_new=dynrModel_new,
         dynrCook_new=dynrCook_new)
  } else {
    list(dynrModel_new=dynrModel_new)
  }
}

computeJacobian <- function(cookDebug, jacobian, lat_name, params, time){
  envList <- as.list(params)
  for(i in 1:length(lat_name)) envList[[lat_name[i]]] <- cookDebug$eta_smooth_final[i, time]

  J <- matrix(NA, length(lat_name), length(lat_name), dimnames=list(lat_name, lat_name))
  for(i in 1:(length(lat_name))^2){
    cj <- as.character(jacobian[[i]])
    rc <- strsplit(cj[2], ' ~ ')[[1]]
    J[rc[1], rc[2]] <- eval(parse(text=cj[3]), envList)
  }
  return(J)
}

#example use from demo/NonlinearODE.R
#computeJacobian(res, dynm$jacobian[[1]], model$measurement$state.names, coef(res), 50)
#cf t=1 vs t=50

cj <- function(cookDebug, jacobian, lat_name, params){
  eta <- cookDebug$eta_smooth_final
  rownames(eta) <- lat_name
  tt <- ncol(eta)
  lat_n <- length(lat_name)
  pp <- matrix(params, nrow=length(params), ncol=tt)
  rownames(pp) <- names(params)
  ddenv <- as.data.frame( t(rbind(pp, eta)) )
  jac_n <- lapply(jacobian, function(jj) {
    jj[[2]]
  })
  #paste(lat_name, collapse=" ~ ")
  ff <- vector("list", length=4)
  ii <- 1
  for (rr in lat_name) {
    for (cc in lat_name) {
      ff[[ii]] <- as.formula(paste(c(rr, cc), collapse=" ~ "))
      ii <- ii + 1
    }
  }
  all.equal(jac_n, ff)
  
  # jjl <- lapply(jacobian, function(jj) {
  #   eval(jj[[3]], ddenv)
  # })
  # jjm <- do.call(rbind, jjl)
  # jja <- array(NA, dim=c(lat_n, lat_n, tt))
  # for (rr in 1:tt) {
  #   jja[,,rr] <- matrix(jjm[,rr], nrow=lat_n, ncol=lat_n, byrow=TRUE)
  # }
  # jja
}

#example use from demo/NonlinearODE.R
#computeJacobian(res, dynm$jacobian[[1]], model$measurement$state.names, coef(res), 3)
#str(aa)
#aa[,,3]
#cf t=1 vs t=50

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
