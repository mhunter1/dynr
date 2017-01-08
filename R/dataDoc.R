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
