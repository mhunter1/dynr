##' Simulated single-subject time series to capture features of facial electromyography data
##' 
##' A dataset simulated using an autoregressive model of order (AR(1)) with
##' regime-specific AR weight, intercept, and slope for a covariate. This model
##' is a special case of Model 1 in Yang and Chow (2010) in which the moving average
##' coefficient is set to zero.
##'  
##'  Reference:
##'  Yang, M-S. & Chow, S-M. (2010). Using state-space models with regime switching to 
##'  represent the dynamics of facial electromyography (EMG) data. Psychometrika, 74(4), 744-771
##'
##' The variables are as follows:
##' 
##' \itemize{
##'   \item id. ID of the participant (= 1 in this case, over 500 time points)
##'   \item EMG. Hypothetical observed facial electromyograhy data
##'   \item self. Covariate - the individual's concurrent self-reports
##'   \item truestate. The true score of the individual's EMG at each time point
##'   \item trueregime. The true underlying regime for the individual at each time point
##' }
##' 
##' @docType data
##' @keywords datasets
##' @name EMGsim
##' @usage data(EMGsim)
##' @format A data frame with 500 rows and 6 variables
NULL

##' Single-subject time series of facial electromyography data
##'
##' A dataset obtained and analyzed in Yang and Chow (2010).
##'
##'  Reference:
##'  Yang, M-S. & Chow, S-M. (2010). Using state-space models with regime switching to
##'  represent the dynamics of facial electromyography (EMG) data. Psychometrika, 74(4), 744-771
##'
##' The variables are as follows:
##'
##' \itemize{
##'   \item id. ID of the participant (= 1 in this case, over 695 time points)
##'   \item time Time in seconds
##'   \item iEMG. Observed integrated facial electromyograhy data
##'   \item SelfReport. Covariate - the individual's concurrent self-reports
##' }
##'
##' @docType data
##' @keywords datasets
##' @name EMG
##' @usage data(EMG)
##' @format A data frame with 695 rows and 4 variables
NULL

##' Simulated multi-subject time series based on a dynamic factor analysis model with nonlinear relations at the latent level
##' 
##' A dataset simulated using a discrete-time nonlinear dynamic factor analysis model
##' with 6 observed indicators for identifying two latent factors: individuals'
##' positive and negative emotions. Proposed by Chow and Zhang (2013), the model was inspired 
##' by models of affect and it posits that the two latent factors follow a vector autoregressive
##' process of order 1 (VAR(1)) with parameters that vary between two possible regimes:
##' (1) an "independent" regime in which the lagged influences between positive and negative
##' emotions are zero; (2) a "high-activation" regime to capture instances
##' on which the lagged influences between PA and NA intensify when an individual's previous 
##' levels of positive and negative emotions were unusually high or low (see Model 2 in Chow & Zhang).
##'  
##'  Reference:
##'  Chow, S-M, & Zhang, G. (2013). Regime-switching nonlinear dynamic factor analysis 
##'  models. Psychometrika, 78(4), 740-768.
##' 
##' \itemize{
##'   \item id. ID of the participant (1 to 10)
##'   \item time. Time index (300 time points from each subject)
##'   \item y1-y3. Observed indicators for positive emotion
##'   \item y4-y6. Observed indicators for negative emotion
##' }
##' 
##' @docType data
##' @keywords datasets
##' @name NonlinearDFAsim
##' @usage data(NonlinearDFAsim)
##' @format A data frame with 3000 rows and 8 variables
NULL

##' Simulated time series data for multiple eco-systems based on a predator-and-prey model
##' 
##' A dataset simulated using a continuous-time nonlinear predator-and-prey model
##' with 2 observed indicators for identifying two latent factors. The variables are as follows:
##' 
##' \itemize{
##'   \item id. ID of the systems (1 to 20)
##'   \item time. Time index (50 time points for each system)
##'   \item prey. The true population of the prey species
##'   \item predator. The true population of the predator species
##'   \item x. Observed indicator for the population of the prey species
##'   \item y. Observed indicator for the population of the predator species
##' }
##' 
##' @docType data
##' @keywords datasets
##' @name PPsim
##' @usage data(PPsim)
##' @format A data frame with 1000 rows and 6 variables
NULL

##' Simulated time series data for multiple eco-systems based on a regime-switching predator-and-prey model
##' 
##' A dataset simulated using a regime-switching continuous-time nonlinear predator-and-prey model
##' with 2 observed indicators for identifying two latent factors. The variables are as follows:
##' 
##' \itemize{
##'   \item id. ID of the systems (1 to 20)
##'   \item time. Time index (300 time points for each system)
##'   \item prey. The true population of the prey species
##'   \item predator. The true population of the predator species
##'   \item x. Observed indicator for the population of the prey species
##'   \item y. Observed indicator for the population of the predator species
##'   \item cond. A time-varying covariate indicating the conditions of the respective eco-system across time which 
##'   affects the regime-switching transition matrix
##'   \item regime. The true regime indicators across time (1 and 2).
##' }
##' 
##' @docType data
##' @keywords datasets
##' @name RSPPsim
##' @usage data(RSPPsim)
##' @format A data frame with 6000 rows and 8 variables
NULL

##' Simulated time series data of a damped linear oscillator
##' 
##' A dataset simulated using a damped linear oscillator model in continuous time
##' with 1 observed indicator for identifying two latent factors (position and velocity).
##' The variables are as follows:
##' 
##' \itemize{
##'   \item id. ID of the systems (1 to 1 because this is a single person)
##'   \item y1. Noisy observed position
##'   \item times. Time index (1000 time points) spaced at one unit intervals
##'   \item x1. True latent position
##'   \item x2. True latent velocity
##' }
##' 
##' @docType data
##' @keywords datasets
##' @name Oscillator
##' @usage data(Oscillator)
##' @format A data frame with 1000 rows and 5 variables
##' 
##' @examples
##' # The following was used to generate the data
##' #--------------------------------------
##' # Data Generation
##' #require(mvtnorm)
##' #require(Matrix)
##' #
##' #xdim <- 2
##' #udim <- 1
##' #ydim <- 1
##' #tdim <- 1000
##' #set.seed(315)
##' #tA <- matrix(c(0, -.3, 1, -.7), xdim, xdim)
##' #tB <- matrix(c(0), xdim, udim)
##' #tC <- matrix(c(1, 0), ydim, xdim)
##' #tD <- matrix(c(0), ydim, udim)
##' #tQ <- matrix(c(0), xdim, xdim); diag(tQ) <- c(0, 2.2)
##' #tR <- matrix(c(0), ydim, ydim); diag(tR) <- c(1.5)
##' #
##' #x0 <- matrix(c(0, 1), xdim, 1)
##' #P0 <- diag(c(1), xdim)
##' #tdx <- matrix(0, xdim, tdim+1)
##' #tx <- matrix(0, xdim, tdim+1)
##' #tu <- matrix(0, udim, tdim)
##' #ty <- matrix(0, ydim, tdim)
##' #
##' #tT <- matrix(0:tdim, nrow=1, ncol=tdim+1)
##' #
##' #tI <- diag(1, nrow=xdim)
##' #
##' #tx[,1] <- x0
##' #for(i in 2:(tdim+1)){
##' #	q <- t(rmvnorm(1, rep(0, xdim), tQ))
##' #	tdx[,i] <- tA %*% tx[,i-1] + tB %*% tu[,i-1] + q
##' #	expA <- as.matrix(expm(tA * (tT[,i]-tT[,i-1])))
##' #	intA <- solve(tA) %*% (expA - tI)
##' #	tx[,i] <- expA %*% tx[, i-1] + intA %*% tB %*% tu[,i-1] + intA %*% q
##' #	ty[,i-1] <- tC %*% tx[,i] + tD %*% tu[,i-1] + t(rmvnorm(1, rep(0, ydim), tR))
##' #}
##' #
##' #
##' #
##' #rownames(ty) <- paste('y', 1:ydim, sep='')
##' #rownames(tx) <- paste('x', 1:xdim, sep='')
##' #simdata <- cbind(id=rep(1, tdim), t(ty), times=tT[,-1], t(tx)[-1,])
##' # write.table(simdata, file='Oscillator.txt', row.names=FALSE, col.names=TRUE)
##' #
##' #plot(tx[1,], type='l')
##' #plot(tT[,-1], ty[1,], type='l')
NULL

##' Simulated time series data for a stochastic linear damped oscillator model with logistic time-varying setpoints
##' 
##' A dataset simulated using a continuous-time stochastic linear damped oscillator model.
##' The variables are as follows:
##' 
##' \itemize{
##'   \item id. ID of the systems (1 to 10)
##'   \item times. Time index (241 time points for each system)
##'   \item x. Latent level variable
##'   \item y. Latent first derivative variable
##'   \item z. True values of time-varying setpoints
##'   \item obsy. Observed level
##' }
##' 
##' @docType data
##' @keywords datasets
##' @name LogisticSetPointSDE
##' @usage data(LogisticSetPointSDE)
##' @format A data frame with 2410 rows and 6 variables
##' @examples
##' # The following was used to generate the data
##' #--------------------------------------
##' #require(Sim.DiffProc)
##' #freq <- -1
##' #damp <- -.1
##' #mu <- -2
##' #r <- .5
##' #b <- .1
##' #sigma1 <- 0.1
##' #sigma2 <- 0.1
##' #fx <- expression(y, freq*(x-z) + damp*y, r*z*(1-b*z))
##' #gx <- expression(0, sigma1, 0)
##' #r3dall <- c()
##' #for (j in 1:10){
##' #  r3dtemp <- c(-5,0,.1)
##' #  r3d <- r3dtemp
##' #  for (i in seq(0.125, 30, by=0.125)){
##' #    mod3dtemp <- snssde3d(drift=fx, diffusion=gx, M=1, t0=i-0.125,
##' #        x0=as.numeric(r3dtemp), T=i, N=500, type="str",
##' #        method="smilstein")
##' #    r3dtemp <- rsde3d(mod3dtemp,at=i)
##' #    r3d <-rbind(r3d,r3dtemp)
##' #  }
##' #  r3dall <- rbind(r3dall, cbind(r3d, id=j))
##' #}
##' #
##' #r3dall$obsy <- r3dall$x+rnorm(length(r3dall$x),0,1)
##' #write.table(r3dall, file="LogisticSetPointSDE.txt")
NULL

##' Simulated time series data of a multisubject process factor analysis
##' 
##' A multiple subject dataset simulated using a two factor process factor analysis model in discrete time
##' with 6 observed indicators for identifying two latent factors.
##' The variables are as follows:
##' 
##' \itemize{
##'   \item ID. Person ID variable (1 to 50 because there are 50 simulated people)
##'   \item Time. Time ID variable (1 to 50 because there are 50 time points)
##'   \item V1. Noisy observed variable 1
##'   \item V2. Noisy observed variable 2
##'   \item V3. Noisy observed variable 3
##'   \item V4. Noisy observed variable 4
##'   \item V5. Noisy observed variable 5
##'   \item V6. Noisy observed variable 6
##'   \item F1. True latent variable 1 scores
##'   \item F2. True latent variable 2 scores
##' }
##' 
##' Variables V1, V2, and V3 load on F1, whereas variables V4, V5, V6 load on F2.  The true values of the factor loadings are 1, 2, 1, 1, 2, and 1, respectively.  The true measurement error variance is 0.5 for all variables.  The true dynamic noise covariance has F1 with a variance of 2.77, F2 with a variance of 8.40, and their covariance is 2.47.  The across-time dynamics have autoregressive effects of 0.5 for both F1 and F2 with a cross-lagged effect from F1 to F2 at 0.4.  The cross-lagged effect from F2 to F1 is zero.  The true initial latent state distribution as mean zero and a diagonal covariance matrix with var(F1) = 2 and var(F2) = 1.  The generating model is the same for all individuals.
##' 
##' @docType data
##' @keywords datasets
##' @name PFAsim
##' @usage data(PFAsim)
##' @format A data frame with 2,500 rows and 10 variables
##' 
##' @examples
##' # The following was used to generate the data
##' #set.seed(12345678)
##' #library(mvtnorm)
##' ## setting up matrices
##' #time      <- 50
##' ## Occasions to throw out to wash away the effects of initial condition
##' #npad      <- 0
##' #np        <- 50
##' #ne        <- 2 #Number of latent variables
##' #ny        <- 6 #Number of manifest variables
##' ## Residual variance-covariance matrix
##' #psi       <- matrix(c(2.77, 2.47,
##' #                      2.47, 8.40),
##' #                    ncol = ne, byrow = T)
##' ## Lambda matrix containing contemporaneous relations among
##' ## observed variables and 2 latent variables. 
##' #lambda    <- matrix(c(1, 0,
##' #                      2, 0,
##' #                      1, 0,
##' #                      0, 1,
##' #                      0, 2,
##' #                      0, 1),
##' #                    ncol = ne, byrow = TRUE)
##' ## Measurement error variances
##' #theta     <- diag(.5, ncol = ny, nrow = ny)
##' ## Lagged directed relations among variables
##' #beta      <- matrix(c(0.5, 0,
##' #                      0.4, 0.5), 
##' #                    ncol = ne, byrow = TRUE)
##' #a0        <- mvtnorm::rmvnorm(1, mean = c(0, 0),
##' #                                 sigma = matrix(c(2,0,0,1),ncol=ne))
##' #yall <- matrix(0,nrow = time*np, ncol = ny)
##' #eall <- matrix(0,nrow = time*np, ncol = ne)
##' #for (p in 1:np){
##' #  # Latent variable residuals
##' #  zeta      <- mvtnorm::rmvnorm(time+npad, mean = c(0, 0), sigma = psi)
##' #  # Measurement errors
##' #  epsilon   <- rmvnorm(time, mean = c(0, 0, 0, 0, 0, 0), sigma = theta)
##' #  # Set up matrix for contemporaneous variables
##' #  etaC      <- matrix(0, nrow = ne, ncol = time + npad)
##' #  # Set up matrix for lagged variables
##' #  etaL      <- matrix(0, nrow = ne, ncol = time + npad + 1)
##' #  
##' #  etaL[,1]  <- a0
##' #  etaC[,1] <- a0
##' #  # generate factors
##' #  for (i in 2:(time+npad)){
##' #    etaL[ ,i] <- etaC[,i-1]
##' #    etaC[ ,i]   <- beta %*% etaL[ ,i] + zeta[i, ]
##' #  }
##' #  etaC <- etaC[,(npad+1):(npad+time)]
##' #  eta <- t(etaC)
##' #  
##' #  # generate observed series
##' #  y   <- matrix(0, nrow = time, ncol = ny)
##' #  for (i in 1:nrow(y)){
##' #    y[i, ] <- lambda %*% eta[i, ] + epsilon[i, ]
##' #  }
##' #  yall[(1+(p-1)*time):(p*time),] <- y
##' #  eall[(1+(p-1)*time):(p*time),] <- eta
##' #}
##' #yall <- cbind(rep(1:np,each=time),rep(1:time,np),yall)
##' #yeall <- cbind(yall,eall)
##' #write.table(yeall,'PFAsim.txt',row.names=FALSE,
##' #  col.names=c("ID", "Time", paste0("V", 1:ny), paste0("F", 1:ne)))
NULL

##' Simulated time series data for detecting outliers.
##' 
##' This is a list object containing true outliers, the dataset, and the saved result from running dynr.taste.
##' 
##' The true outliers for observed variables are saved in `Outliers$generated$shockO'.
##' \itemize{
##'   \item id. Six outliers were added for each ID. 
##'   \item time_O. Time points where the outliers were added.
##'   \item obs. Variable indices where the outliers were added.
##'   \item shock.O. The magnitude of outliers.
##' }
##' 
##' The true outliers for state variables are saved in `Outliers$generated$shockL'.
##' \itemize{
##'   \item id. Three outliers were added for each ID. 
##'   \item time_L. Time points where the outliers were added.
##'   \item lat. Variable indices where the outliers were added.
##'   \item shock.L. The magnitude of outliers.
##' }
##' 
##' A dataset simulated based on state-space model including the outliers. The data is saved in `Outliers$generated$y'.
##' The variables are as follows:
##' 
##' \itemize{
##'   \item id. ID of the systems (1 to 100)
##'   \item times. Time indices (100 time points for each participant)
##'   \item V1 - V6. observed variables
##' }
##' 
##' The detected innovative outliers from dynr.taste for this dataset, which is used for testing whether the dynr.taste replicate the same result. The data is saved in `Outliers$detect_O'.
##' The variables are as follows:
##' 
##' \itemize{
##'   \item id. IDs
##'   \item time_L. Time points where the outliers were detected
##'   \item obs. Variable indices for observed variables where the outliers were detected
##' }
##' 
##' The detected additive outliers from dynr.taste for this dataset, which is used for testing whether the dynr.taste replicate the same result. The data is saved in `Outliers$detect_L'.
##' The variables are as follows:
##' 
##' \itemize{
##'   \item id. IDs
##'   \item time_L. Time points where the outliers were detected
##'   \item obs. Variable indices for latent variables where the outliers were detected
##' }
##' 
##' @docType data
##' @keywords datasets
##' @name Outliers
##' @usage data(Outliers)
##' @format A data frame with 6000 rows and 6 variables
##' @examples
##' # The following was used to generate the data
##' #--------------------------------------
##' # lambda <- matrix(c(1.0, 0.0,
##' # 0.9, 0.0,
##' # 0.8, 0.0,
##' # 0.0, 1.0,
##' # 0.0, 0.9,
##' # 0.0, 0.8), ncol=2, byrow=TRUE)
##' # psi <- matrix(c(0.3, -0.1,
##' #                 -0.1, 0.3), ncol=2, byrow=TRUE)
##' # beta <- matrix(c(0.8, -0.2,
##' #                  -0.2,  0.7), ncol=2, byrow=TRUE)
##' # theta <- diag(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2), ncol=6, nrow=6)
##' # nlat <- 2; nobs <- 6
##' # mean_0 <- rep(0, nlat)
##' # psi_inf <- diag(1, 2*2) - kronecker(beta, beta)
##' # psi_inf_inv <- try(solve(psi_inf), silent=TRUE)
##' # if("try-error" %in% class(psi_inf_inv)) {
##' #   psi_inf_inv <- MASS::ginv(psi_inf)}
##' # psi_0 <- psi_inf_inv %*% as.vector(psi)
##' # dim(psi_0) <- c(2, 2)
##' # # measurement error covariance matrix
##' # mea_cov <- lambda %*% psi_0 %*% t(lambda) + theta
##' # resL <- lapply(1:100, function(subj) {
##' #   # initial state
##' #   eta_0 <- mvtnorm::rmvnorm(1, mean=mean_0, sigma=psi_0)#[1,nlat]
##' #   zeta_0 <- mvtnorm::rmvnorm(1, mean=rep(0, nlat), sigma=psi)
##' #   eta <- matrix(0, nrow=time, ncol=nlat)
##' #   eta[1, ] <- beta %*% t(eta_0) + t(zeta_0) 
##' #   zeta <- mvtnorm::rmvnorm(time, mean=rep(0, nlat), sigma=psi)
##' #   # random shock generation
##' #   # to avoid shock appearing too early or late (first and last 3)
##' #   shkLat_time <- sample(4:(time-3), nshockLat)
##' #   shk_lat <- sample(1:nlat, nshockLat, replace=TRUE)
##' #   shockLatIdx <- matrix(c(shkLat_time, shk_lat), ncol=2)
##' #   shockSignL <- sample(c(1,-1), nshockLat, replace=TRUE)
##' #   colnames(shockLatIdx) <- c("time_L","lat")
##' #   shockLatV <- shockSignL*( shockMag*sqrt(diag(shockPsi)))[shockLatIdx[,"lat"]]
##' #   shockLatM <- matrix(0, time, nlat)
##' #   shockLatM[shockLatIdx] <- shockLatV
##' #   shkObs_time <- sample(4:(time-3), nshockObs)
##' #   shk_obs <- sample(1:nobs, nshockObs, replace=TRUE)
##' #   shockObsIdx <- matrix(c(shkObs_time, shk_obs), ncol=2)
##' #   shockSignO <- sample(c(1,-1), nshockObs, replace=TRUE)
##' #   colnames(shockObsIdx) <- c("time_O","obs")
##' #   shockObsV <- shockSignO*( shockMag*sqrt(diag(mea_cov)) )[shockObsIdx[,"obs"]]
##' #   shockObsM <- matrix(0, time, nobs)
##' #   shockObsM[shockObsIdx] <- shockObsV
##' #   # generate state process WITH shock
##' #   for (t in 1:(time-1)) {
##' #     eta[t+1, ] <- shockLatM[t, ] + beta %*% eta[t, ] + zeta[t, ]
##' #   }
##' #   # generate observed process
##' #   y <- shockObsM + eta %*% t(lambda) +
##' #     mvtnorm::rmvnorm(time, mean=rep(0, nobs), sigma=theta)# epsilon
##' # }
NULL

##' Simulated time series data for multiple imputation in dynamic modeling.
##' 
##' A dataset simulated using a vector autoregressive (VAR) model of order 1 with 
##' two observed variables and two covariates. Data are generated following 
##' the simulation design illustrated by Ji and colleagues (2018). Specifically, 
##' missing data are generated following the missing at random (MAR) condition under which 
##' the probability of missingness in both dependent variables and covariates is conditioned
##' on two completely observed auxiliary variables.   
##' 
##' The variables are as follows:
##' \itemize{
##'   \item ID. ID of the participant (1 to 100)
##'   \item Time. Time index (100 time points from each subject)
##'   \item ca. Covariate 1
##'   \item cn. Covariate 2
##'   \item wp. Dependent variable 1
##'   \item hp. Dependent variable 2
##'   \item x1. Auxiliary variable 1
##'   \item x2. Auxiliary variable 2
##' }
##' 
##' @references
##' Ji, L., Chow, S-M., Schermerhorn, A.C., Jacobson, N.C., & Cummings, E.M. (2018). Handling 
##' Missing Data in the Modeling of Intensive Longitudinal Data. Structural Equation Modeling: 
##' A Multidisciplinary Journal, 1-22.
##' 
##' @docType data
##' @keywords datasets
##' @name VARsim
##' @usage data(VARsim)
##' @format A data frame with 10000 rows and 8 variables
NULL

##' Simulated multilevel multi-subject time series of a Van der Pol Oscillator
##' 
##' A dataset simulated using methods described in the reference below.
##'  
##'  Reference:
##'  Chow, S., Lu, Z., Sherwood, A., and Zhu, H. (2016). Fitting Nonlinear Ordinary
##'  Differential Equation Models with Random Effects and Unknown Initial Conditions
##'  Using the Stochastic Approximation Expectation-Maximization (SAEM) Algorithm.
##'  Psychometrika, 81(1), 102-134.
##'
##' The variables are as follows:
##' 
##' \itemize{
##'   \item batch. Batch number from simulation
##'   \item kk. Unclear
##'   \item trueInit. True initial condition
##'   \item id. Person ID
##'   \item time. Continuous time of measurement
##'   \item y1. Observed score 1
##'   \item y2. Observed score 2
##'   \item y3. Observed score 3
##'   \item co1. Covariate 1
##'   \item co2. Covariate 2
##' }
##' 
##' @docType data
##' @keywords datasets
##' @name TrueInit_Y14
##' @usage data(TrueInit_Y14)
##' @format A data frame with 60,000 rows and 10 variables
NULL

##' Simulated time series data for a deterministic linear damped oscillator model 
##' 
##' The variables are as follows:
##' 
##' \itemize{
##'   \item ID. ID of the systems (1 to 10)
##'   \item x. Latent level variable
##'   \item theTimes. Measured time Points
##' }
##' 
##' @docType data
##' @keywords datasets
##' @name LinearOsc
##' @usage data(LinearOsc)
##' @format A data frame with 1000 rows and 3 variables
##' @examples
##' # The following was used to generate the data
##' #--------------------------------------
##' #Osc<- function(t, prevState, parms) {
##' #  x1 <- prevState[1] # x1[t]
##' #  x2 <- prevState[2] # x2[t]
##' #  eta1 = parms[1]
##' #  zeta1 = parms[2]
##' #  with(as.list(parms), {
##' #   dx1 <- x2
##' #    dx2 <- eta1*x1 + zeta1*x2 
##' #    res<-c(dx1,dx2)
##' #    list(res)
##' #  }
##' #  )
##' #}
##' #n = 10 #Number of subjects
##' #T = 100 #Number of time points
##' #deltaT = .1 #dt
##' #lastT = deltaT*T #Value of t_{i,T}
##' #theTimes  = seq(0, lastT, length=T)  #A list of time values
##' #
##' #eta = -.8
##' #zeta = -.1
##' #out1 = matrix(NA,T*n,1)
##' #trueOut = matrix(NA,T*n,1)
##' #parms = c(eta, zeta)
##' #  for (i in 1:n){
##' #  xstart = c(rnorm(1,0,2),rnorm(1,0,.5))
##' #  out <- lsoda(as.numeric(xstart), theTimes, Osc, parms)
##' #  trueOut[(1+(i-1)*T):(i*T)] = out[,2]
##' #  out1[(1+(i-1)*T):(i*T)] = out[,2]+rnorm(T,0,1)
##' #  }
##' #
##' #LinearOsc= data.frame(ID=rep(1:n,each=T),x=out1[,1],
##' #                  theTimes=rep(theTimes,n))
##' #save(LinearOsc,file="LinearOsc.rda")
NULL

##' Another simulated multilevel multi-subject time series of a Van der Pol Oscillator
##' 
##' A dataset simulated using methods described in the reference below.
##'  
##'  Reference:
##'  Chow, S., Lu, Z., Sherwood, A., and Zhu, H. (2016). Fitting Nonlinear Ordinary
##'  Differential Equation Models with Random Effects and Unknown Initial Conditions
##'  Using the Stochastic Approximation Expectation-Maximization (SAEM) Algorithm.
##'  Psychometrika, 81(1), 102-134.
##'
##' The variables are as follows:
##' 
##' \itemize{
##'   \item batch. Batch number from simulation
##'   \item kk. Unclear
##'   \item trueInit. True initial condition
##'   \item id. Person ID
##'   \item time. Continuous time of measurement
##'   \item y1. Observed score 1
##'   \item y2. Observed score 2
##'   \item y3. Observed score 3
##'   \item u1. Covariate 1
##'   \item u2. Covariate 2
##'   \item trueb. True value of person-specific random effect	
##' }
##' 
##' @docType data
##' @keywords datasets
##' @name vdpData
##' @usage data(vdpData)
##' @format A data frame with 10,000 rows and 11 variables
NULL

##' Another simulated multilevel multi-subject time series of a damped oscillator model
##'
##' The variables are as follows:
##' 
##' \itemize{
##'   \item id. Person ID
##'   \item times. Continuous time of measurement
##'   \item y1. Observed score 1
##'   \item u1. Covariate 1
##'   \item u2. Covariate 2
##'   \item trueb. True value of person-specific random effect	
##' }
##' 
##' @docType data
##' @keywords datasets
##' @name oscData
##' @usage data(oscData)
##' @format A data frame with 1,800 rows and 6 variables
NULL
