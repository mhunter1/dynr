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
##' must be one of ``two.sided'' (default), ``greater''  or ``less''
##' @param outliers a character string specifying the outlier detection
##' method. must be one of ``joint'' (default), ``innovative'' or ``additive''.
##' When ``joint'' is selected, dynr.taste detects innovative and additive outliers together.
##' @param debug_flag a logical. 'TRUE' for output of by-products related to t-value calculation
##' 
##' @return an object of `dynrTaste' class
##' that is a list containing results lists for each participant.
##' The result list for a paticipant includes
##' a data.frame called `taste', 
##' a list called `delta_chi' for detected magnitude outliers when the chi-square test was used,
##' a list called `delta_t' for detected magnitude outliers when the t-test was used,
##' and a list called `debug' for additional information,
##' which are Q, S, s, F_inv, N, u, r. 
##' See the reference for definition of the notations.
##' 
##' The `debug' will be saved only when the argument \code{debug_flag=TRUE}.
##' The `delta_chi' list comprises magnitude of innovative (Latent) and additive (Observed) outliers, `delta.L' and `delta.O',
##' when chi-square statitics is used to detect outliers.
##' The `delta_t' list comprises magnitude of innovative (Latent) and additive (Observed) outliers, `delta.L' and `delta.O',
##' when t statitics is used to detect outliers.
##' 
##' @references 
##' Chow, S.-M., Hamaker, E. L., & Allaire, J. C. (2009). 
##' Using innovative outliers to detectdiscrete shifts in dynamics in group-based state-space models. _Multivariate BehavioralResearch_, 44, 465–496. 
dynr.taste <- function(dynrModel, dynrCook=NULL, conf.level=0.99,
                       alternative=c("two.sided", "less", "greater"),
                       outliers=c("joint", "innovative", "additive"),
                       debug_flag=FALSE) {
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
  obsName <- dynrModel$measurement$obs.names
  dimLat <- length(stateName)
  dimObs <- length(obsName)
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
  r <- matrix(0, dimLat, dimTime)
  N <- array(0, c(dimLat, dimLat, dimTime))
  # save Ninv to calculate W_t (for delta)
  Ninv <- array(0, c(dimLat, dimLat, dimTime))
  u <- matrix(0, dimObs, dimTime)
  tstart <- dynrModel$data$tstart
  chiLat <- numeric(dimTime)
  
  for(j in 1:nID){
    beginTime <- tstart[j] + 1 # Add 1 because model_tstart indexes from 0 rather than 1 for interface to C
    endTime <- tstart[j+1] # Use the 0-indexed beginning position of the next person as the 1-indexed end of the current person
    
    # set endTime r and N to 0
    #  --> set r and N 0 at first when allocate
    #r[,endTime] <- 0
    #N[,,endTime] <- 0
    
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
  }
  
  ################ delta estimate #############################
  outliers <- match.arg(outliers)
  if (outliers=="joint") {
    # De Jong and Penzer (1988) 80p. below eq 14 about X, W and delta
    delta <- matrix(NA, dimObs+dimLat, dimTime)
    rownames(delta) <- c(obsName, stateName)
    # (X',W') = I. (De Jong and Penzer, 1988)
    XW <- diag(1, dimObs+dimLat)
    X_t <- XW[1:dimObs, ]
    W_t <- XW[(dimObs+1):(dimObs+dimLat), ]
    # t-values
    t_value <- matrix(NA, dimObs+dimLat, dimTime)
    rownames(t_value) <- c(obsName, stateName)
    
  } else if (outliers=="innovative") {
    delta <- matrix(NA, dimLat, dimTime)
    rownames(delta) <- c(stateName)
    # (X',W') = I. (De Jong and Penzer, 1988)
    # care for dimensions
    # X: [dimObs, dimObs + dimLat], W: [dimLat, dimObs + dimLat]
    X_t <- matrix(0, dimObs, dimLat)
    W_t <- diag(1, dimLat)
    # t-values
    t_value <- matrix(NA, dimLat, dimTime)
    rownames(t_value) <- c(stateName)
  
  } else {#outliers=="additive"
    delta <- matrix(NA, dimObs, dimTime)
    rownames(delta) <- c(obsName)
    X_t <- diag(1, dimObs)
    W_t <- matrix(0, dimLat, dimObs)
    # t-values
    t_value <- matrix(NA, dimObs, dimTime)
    rownames(t_value) <- c(obsName)
  }
  ###### save by-products related to t ############
  if (debug_flag) {
    Q <- array(0, c(dim(W_t)[1], dim(W_t)[2], dimTime))
    S <- array(0, c(dim(W_t)[2], dim(W_t)[2], dimTime))
    s <- matrix(0, dim(W_t)[2], dimTime)
    # F_inv, N, u, r
    
    for(i in 1:nID){
      beginTime <- tstart[i] + 1
      endTime <- tstart[i+1]
      time_i <- endTime - tstart[i]
      
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
    
    ## TEMPORAL: T-REALTED!
    Q_sp <- vector("list", length=nID)
    S_sp <- vector("list", length=nID)
    s_sp <- vector("list", length=nID)
    F_inv_sp <- vector("list", length=nID)
    N_sp <- vector("list", length=nID)
    u_sp <- vector("list", length=nID)
    r_sp <- vector("list", length=nID)
    for(i in 1:nID){
      beginTime <- tstart[i] + 1
      endTime <- tstart[i+1]
      
      Q_sp[[i]] <- Q[,,beginTime:endTime]
      S_sp[[i]] <- S[,,beginTime:endTime]
      s_sp[[i]] <- s[,beginTime:endTime]
      F_inv_sp[[i]] <- F_inv[,,beginTime:endTime]
      N_sp[[i]] <- N[,,beginTime:endTime]
      u_sp[[i]] <- u[,beginTime:endTime]
      r_sp[[i]] <- r[,beginTime:endTime]
    }
  } else {
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
        t_value[,j] <- s_j / sqrt(diag(S_j))
      }
      # 0/0 at endTime
      t_value[,endTime] <- 0
    }
  }
  
  ############ chi-square test ############################
  id <- as.factor(dynrModel$data$id)
  if (outliers=="joint") {
    # de Jong's method.
    chiBoth <- chiLat + chiObs
    #chiBoth_sp <- split(chiBoth, id)
    chiLat_sp <- chiObs_sp <- split(chiBoth, id)
    #chiBoth_pval <- pchisq(chiBoth, df=dimLat+dimObs, lower.tail=FALSE)
    chiLat_pval <- chiObs_pval <- 
      pchisq(chiBoth, df=dimLat+dimObs, lower.tail=FALSE)
    #chiBoth_pval_sp <- split(chiBoth_pval, id)
    #chiBoth_shk <- chiBoth_pval < (1 - conf.level)
    #chiBoth_shk_sp <- split(chiBoth_shk, id) 
  } else {
    chiLat_sp <- split(chiLat, id)
    chiObs_sp <- split(chiObs, id)
    # p-values of chi-sqaure
    chiLat_pval <- pchisq(chiLat, df=dimLat, lower.tail=FALSE)
    chiObs_pval <- pchisq(chiObs, df=dimObs, lower.tail=FALSE)
  }
  chiLat_pval_sp <- split(chiLat_pval, id)
  chiObs_pval_sp <- split(chiObs_pval, id)
  # locate shock points, TRUE for significance
  chiLat_shk <- chiLat_pval < (1 - conf.level)
  chiLat_shk_sp <- split(chiLat_shk, id)
  chiObs_shk <- chiObs_pval < (1 - conf.level)
  chiObs_shk_sp <- split(chiObs_shk, id)
  
  ################### t-test ###############################
  # assign function calculating p-values
  alternative <- match.arg(alternative)
  tLat_calcPval <- if (alternative=="two.sided") {
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
  
  tObs_calcPval <- if (alternative=="two.sided") {
    function(subj) {
      2 * pt( -abs(as.matrix(subj)), df=(nrow(subj)-dimObs) )
    }
  }	else if (alternative=="greater") {
    function(subj) {
      pt( as.matrix(subj), df=(nrow(subj)-dimObs), lower.tail=FALSE )
    }
  }	else {
    function(subj) {
      pt( as.matrix(subj), df=(nrow(subj)-dimObs), lower.tail=TRUE )
    }
  }
  # t-values for each 'outliers' option
  if (outliers=="joint") {
    # t_value for observed
    t_X <- t_value[1:dimObs, ]
    # t_value for latent
    t_W <- t_value[(dimObs+1):(dimObs+dimLat), ]
    
    # [time_i, dimObs] matrix, for each subject
    t_X_sp <- split(as.data.frame(t(t_X)), id)
    # [time_i, dimLat] matrix, for each subject
    t_W_sp <- split(as.data.frame(t(t_W)), id)
    
    # [time_i, dimObs] matrix, for each subject
    t_X_pval_sp <- lapply(t_X_sp, FUN=tObs_calcPval)
    # [time_i, dimLat] matrix, for each subject
    t_W_pval_sp <- lapply(t_W_sp, FUN=tLat_calcPval)
    
    # locate shocks from t-test
    # [time_i, dimObs] matrix, for each subject
    t_X_shk_sp <- lapply(t_X_pval_sp, function(subj) {
      subj < (1 - conf.level)
    })
    t_W_shk_sp <- lapply(t_W_pval_sp, function(subj) {
      subj < (1 - conf.level)
    })
  } else if (outliers=="innovative") {
    # t_value for observed
    t_X <- matrix(NA, 1, dimTime)
    rownames(t_X) <- "t.L"
    # t_value for latent
    t_W <- t_value
    
    # [time_i, dimObs] matrix, for each subject
    t_X_sp <- split(as.data.frame(t(t_X)), id)
    # [time_i, dimLat] matrix, for each subject
    t_W_sp <- split(as.data.frame(t(t_W)), id)
    
    # [time_i, dimObs] matrix, for each subject
    t_X_pval_sp <- t_X_sp
    # [time_i, dimLat] matrix, for each subject
    t_W_pval_sp <- lapply(t_W_sp, FUN=tLat_calcPval)
    
    # locate shocks from t-test
    # [time_i, dimObs] matrix, for each subject
    t_X_shk_sp <- t_X_sp
    t_W_shk_sp <- lapply(t_W_pval_sp, function(subj) {
      subj < (1 - conf.level)
    })
  } else {#outliers=="additive"
    # t_value for observed
    t_X <- t_value
    # t_value for latent
    t_W <- matrix(NA, 1, dimTime)
    rownames(t_W) <- "t.O"
    
    # [time_i, dimObs] matrix, for each subject
    t_X_sp <- split(as.data.frame(t(t_X)), id)
    # [time_i, dimLat] matrix, for each subject
    t_W_sp <- split(as.data.frame(t(t_W)), id)
    
    # [time_i, dimObs] matrix, for each subject
    t_X_pval_sp <- lapply(t_X_sp, FUN=tObs_calcPval)
    # [time_i, dimLat] matrix, for each subject
    t_W_pval_sp <- t_W_sp
    
    # locate shocks from t-test
    # [time_i, dimObs] matrix, for each subject
    t_X_shk_sp <- lapply(t_X_pval_sp, function(subj) {
      subj < (1 - conf.level)
    })
    t_W_shk_sp <- t_W_sp
  }

  if (outliers=="joint") {
    rownames(delta) <- paste0("d_", c(obsName, stateName))
    # delta for observed
    delta_X <- delta[1:dimObs, ]
    # delta for latent
    delta_W <- delta[(dimObs+1):(dimObs+dimLat), ]
    # split delta  
    # [time_i, dimObs] data.frame, for each subject
    delta_X_sp <- split(as.data.frame(t(delta_X)), id)
    # [time_i, dimLat] data.frame, for each subject
    delta_W_sp <- split(as.data.frame(t(delta_W)), id)
  
  } else if (outliers=="innovative") {
    rownames(delta) <- paste0("d_", stateName)
    # delta for observed
    delta_X <- matrix(0, dimObs, dimTime)
    rownames(delta_X) <- obsName
    # delta for latent
    delta_W <- delta
    # split delta  
    # [time_i, dimObs] data.frame, for each subject
    delta_X_sp <- split(as.data.frame(t(delta_X)), id)
    # [time_i, dimLat] data.frame, for each subject
    delta_W_sp <- split(as.data.frame(t(delta_W)), id)
  
  } else {#outliers=="additive"
    rownames(delta) <- paste0("d_", obsName)
    # delta for observed
    delta_X <- delta
    # delta for latent
    delta_W <- matrix(0, dimLat, dimTime)
    rownames(delta_W) <- stateName
    # split delta  
    # [time_i, dimObs] data.frame, for each subject
    delta_X_sp <- split(as.data.frame(t(delta_X)), id)
    # [time_i, dimLat] data.frame, for each subject
    delta_W_sp <- split(as.data.frame(t(delta_W)), id)
  }
  # id vector
  idv <- unique(id)
  
  # output data.frame for each subject
  time_sp <- split(dynrModel$data$time, id)
  if (debug_flag) {
    res <- mapply(FUN=function(id_i, time,
                               #chiBoth, chiBoth_pval, chiBoth_shk,
                               chiLat, chiLat_pval, chiLat_shk, 
                               chiObs, chiObs_pval, chiObs_shk,
                               t_X, t_W, t_X_pval, t_W_pval, 
                               t_X_shk, t_W_shk, delta_X, delta_W,
                               Q, S, s, F_inv, N, u, r) 
    {
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
        
        list(
          taste=data.frame(
            id=id_i, time=time,
            chi.L=chiLat, chi.L.p=chiLat_pval, chi.L.shk=chiLat_shk,
            chi.O=chiObs, chi.O.p=chiObs_pval, chi.O.shk=chiObs_shk,
            t.L=t_W, t.O=t_X, t.L.p=t_W_pval, t.O.p=t_X_pval,
            t.L.shk=t_W_shk, t.O.shk=t_X_shk
          ),
          delta_chi=list(
            delta.L=delta_W_chi,
            delta.O=delta_X_chi),
          delta_t=list(
            delta.L=delta_W_t,
            delta.O=delta_X_t),
          debug=list(
            Q=Q, S=S, s=s, F_inv=F_inv, N=N, u=u, r=r
          )
        )
        
      } else if (outliers=="innovative") {
        # t shock that pass chi shock
        row.names(delta_X) <- 1:nrow(delta_X)
        row.names(delta_W) <- 1:nrow(delta_W)
        delta_W_t <- delta_W_chi <- delta_W
        delta_X_t <- delta_X_chi <- delta_X
        delta_W_chi[!chiLat_shk,] <- 0
        delta_X_chi[] <- 0
        delta_W_t[!t_W_shk] <- 0
        delta_X_t[] <- 0
        
        list(
          taste=data.frame(
            id=id_i, time=time,
            chi.L=chiLat, chi.L.p=chiLat_pval, chi.L.shk=chiLat_shk,
            t.L=t_W, t.L.p=t_W_pval,
            t.L.shk=t_W_shk
          ),
          delta_chi=list(
            delta.L=delta_W_chi,
            delta.O=delta_X_chi),
          delta_t=list(
            delta.L=delta_W_t,
            delta.O=delta_X_t),
          debug=list(
            Q=Q, S=S, s=s, F_inv=F_inv, N=N, u=u, r=r
          )
        )
        
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
        
        list(
          taste=data.frame(
            id=id_i, time=time,
            chi.O=chiObs, chi.O.p=chiObs_pval, chi.O.shk=chiObs_shk,
            t.O=t_X, t.O.p=t_X_pval,
            t.O.shk=t_X_shk
          ),
          delta_chi=list(
            delta.L=delta_W_chi,
            delta.O=delta_X_chi),
          delta_t=list(
            delta.L=delta_W_t,
            delta.O=delta_X_t),
          debug=list(
            Q=Q, S=S, s=s, F_inv=F_inv, N=N, u=u, r=r
          )
        )
      }
    }, SIMPLIFY=FALSE,
    idv, time_sp,
    chiLat_sp, chiLat_pval_sp, chiLat_shk_sp, 
    chiObs_sp, chiObs_pval_sp, chiObs_shk_sp, 
    t_X_sp, t_W_sp, t_X_pval_sp, t_W_pval_sp, t_X_shk_sp, t_W_shk_sp,
    delta_X_sp, delta_W_sp,
    Q_sp, S_sp, s_sp, F_inv_sp, N_sp, u_sp, r_sp)
    
  } else {
    res <- mapply(FUN=function(id_i, time,
                               chiLat, chiLat_pval, chiLat_shk, 
                               chiObs, chiObs_pval, chiObs_shk,
                               t_X, t_W, t_X_pval, t_W_pval, 
                               t_X_shk, t_W_shk, delta_X, delta_W) 
    {
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
        
        list(
          taste=data.frame(
            id=id_i, time=time,
            chi.L=chiLat, chi.L.p=chiLat_pval, chi.L.shk=chiLat_shk,
            chi.O=chiObs, chi.O.p=chiObs_pval, chi.O.shk=chiObs_shk,
            t.L=t_W, t.O=t_X, t.L.p=t_W_pval, t.O.p=t_X_pval,
            t.L.shk=t_W_shk, t.O.shk=t_X_shk
          ),
          delta_chi=list(
            delta.L=delta_W_chi,
            delta.O=delta_X_chi),
          delta_t=list(
            delta.L=delta_W_t,
            delta.O=delta_X_t)
        )
        
      } else if (outliers=="innovative") {
        # t shock that pass chi shock
        row.names(delta_X) <- 1:nrow(delta_X)
        row.names(delta_W) <- 1:nrow(delta_W)
        delta_W_t <- delta_W_chi <- delta_W
        delta_X_t <- delta_X_chi <- delta_X
        delta_W_chi[!chiLat_shk,] <- 0
        delta_X_chi[] <- 0
        delta_W_t[!t_W_shk] <- 0
        delta_X_t[] <- 0
        
        list(
          taste=data.frame(
            id=id_i, time=time,
            chi.L=chiLat, chi.L.p=chiLat_pval, chi.L.shk=chiLat_shk,
            t.L=t_W, t.L.p=t_W_pval,
            t.L.shk=t_W_shk
          ),
          delta_chi=list(
            delta.L=delta_W_chi,
            delta.O=delta_X_chi),
          delta_t=list(
            delta.L=delta_W_t,
            delta.O=delta_X_t)
        )
        
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
        
        list(
          taste=data.frame(
            id=id_i, time=time,
            chi.O=chiObs, chi.O.p=chiObs_pval, chi.O.shk=chiObs_shk,
            t.O=t_X, t.O.p=t_X_pval,
            t.O.shk=t_X_shk
          ),
          delta_chi=list(
            delta.L=delta_W_chi,
            delta.O=delta_X_chi),
          delta_t=list(
            delta.L=delta_W_t,
            delta.O=delta_X_t)
        )
      }
    }, SIMPLIFY=FALSE,
    idv, time_sp,
    chiLat_sp, chiLat_pval_sp, chiLat_shk_sp, 
    chiObs_sp, chiObs_pval_sp, chiObs_shk_sp, 
    t_X_sp, t_W_sp, t_X_pval_sp, t_W_pval_sp, t_X_shk_sp, t_W_shk_sp,
    delta_X_sp, delta_W_sp)
  }
  
  names(res) <- idv
  class(res) <- "dynrTaste"
  invisible(res)
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
                        verbose=FALSE,
                        delta_L=NULL, delta_O=NULL) {
  if ( length(dynrModel@dynamics@values.exo) != 0 ||
       length(dynrModel@measurement@values.exo) != 0) {
    stop("Currently, a model without covariates can be used.")
  }
  coefx <- coef(dynrCook)
  coef(dynrModel) <- coefx # modifies all values throughout dynrModel
  outlierTest <- match.arg(outlierTest)
  deltaTest <- paste0("delta_", outlierTest)
  # all parameter names + "fixed", to be used for params.xxx
  parNames <- c(names(dynrModel), "fixed")
  # to substitute 'fixed'
  numForFixed <- length(parNames)
  
  stateName <- dynrModel$measurement$state.names
  obsName <- dynrModel$measurement$obs.names
  nDeltaLat <- length(stateName)
  nDeltaObs <- length(obsName)
  
  if ( !is.null(dynrTaste) ) {
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
    if ( is.null(delta_L) ) stop("Both 'dynrTaste' and 'delta_L' are NULL.")
    if (nDeltaLat != ncol(delta_L) || length(dynrModel@data$time) != nrow(delta_L)) stop("The number of column differs with the number of the latent variables.")
    deltaLat <- delta_L
  }
  
  if ( !is.null(dynrTaste) ) {
    deltaObs <- do.call("rbind",
                        lapply(dynrTaste, function(taste_i) {
                          taste_i[[deltaTest]]$delta.O
                        }) )
  } else {
    if ( is.null(delta_O) ) stop("Both 'dynrTast' and 'delta_O' are NULL.")
    if (nDeltaObs != ncol(delta_O) || length(dynrModel@data$time) != nrow(delta_L)) stop("The number of column differs with the number of the observed variables.")
    deltaObs <- delta_O
  }
  
  deltaLatName <- names(deltaLat)
  deltaObsName <- names(deltaObs)
  
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
  new_dynamics <- prep.matrixDynamics(
    values.dyn=dynrModel@dynamics@values.dyn[[1]],
    params.dyn=paramsDyn,
    values.exo=diag(1, nrow=nDeltaLat, ncol=nDeltaLat),
    params.exo=matrix("fixed", nrow=nDeltaLat, ncol=nDeltaLat),
    covariates=deltaLatName,
    isContinuousTime=FALSE)
  
  # modify dynrModel@data
  if ( is.null(dynrModel@data$covariate.names) ) {#no orginal covariates
    dynrModel@data$covariate.names <- c(deltaLatName, deltaObsName)
    deltaLatCopy <- deltaLat
    names(deltaLatCopy) <- paste0("covar", 1:nDeltaLat)
    deltaObsCopy <- deltaObs
    names(deltaObsCopy) <- paste0("covar", (1:nDeltaObs) + nDeltaLat)
    dynrModel@data$covariates <- cbind(deltaLatCopy, deltaObsCopy)
  } else {# original covariates exist
  nPreCovariate <- length(dynrModel@data$covariate.names)
  dynrModel@data$covariate.names <- c(dynrModel@data$covariate.names,
                                      deltaLatName, deltaObsName)
  names(deltaLat) <- paste0("covar", 
                            (1:nDeltaLat) + nPreCovariate)
  names(deltaObs) <- paste0("covar", 
                            (1:nDeltaObs) + nPreCovariate + nDeltaLat)
  dynrModel@data$covariates <- cbind(dynrModel@data$covariates,
                                     deltaLat, deltaObs)
  }
  
  # build measurement
  measParLoad <- dynrModel@measurement@params.load[[1]]
  measParLoad[measParLoad==0] <- numForFixed
  paramsLoad <- parNames[measParLoad]# vector
  dim(paramsLoad) <- dim(measParLoad)# to matrix
  new_measurement <- prep.measurement(
    values.load=dynrModel@measurement@values.load[[1]],
    params.load=paramsLoad,
    values.exo=diag(1, nrow=nDeltaObs, ncol=nDeltaObs),
    params.exo=matrix("fixed", nrow=nDeltaObs, ncol=nDeltaObs),
    exo.names=deltaObsName,
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
    data = dynrModel@data,
    outfile = newOutfile)
  # cook!
  new_dynrCook <- dynr.cook(new_dynrModel, verbose=verbose, debug_flag=TRUE)
  list(new_dynrModel=new_dynrModel,
       new_dynrCook=new_dynrCook)
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

#example use from demo/NonlinearODE.R
#computeJacobian(res, dynm$jacobian[[1]], model$measurement$state.names, coef(res), 50)
#cf t=1 vs t=50


