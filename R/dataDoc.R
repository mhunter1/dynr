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
##' @name TrueInit_Y1
##' @usage data(TrueInit_Y1)
##' @format A data frame with 60,000 rows and 10 variables
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
