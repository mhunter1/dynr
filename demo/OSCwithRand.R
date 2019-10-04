#------------------------------------------------------------------------------
# Author: Sy-Miin Chow & Hui-Ju Hung
# Date: 2019-09-28
# Filename: OSCwithRand.R
# Purpose: An illustrative example for using dynr to estimate random effects
#          for damped oscillator model with random effects in eta, 
#          and random initial conditions (ICs). Fixed and random effects on
#          eta, a frequency-related parameter, is specified via the function
#          prep.thetaFormula
#------------------------------------------------------------------------------

library('dynr')
options(scipen=999)

# ---- Read in the data ----
data(oscData)
data <- dynr.data(oscData, id="id", time="times",
                 observed=c('y1'),
                 covariates=c("u1","u2"))

# ---- Prepare the recipes (i.e., specifies modeling functions) ----

# Measurement (factor loadings)				 
meas <- prep.measurement(
	values.load=matrix(c(1,0), 1, 2),
    params.load=matrix(c('fixed','fixed'), 1, 2),
    obs.names=c('y1'),
    state.names=c('x', 'dx'))

	
# Initial conditions on the latent state and covariance
initial <- prep.initial(
    values.inistate=c(10, 0),
    params.inistate=c("mu_x0", "mu_dx0"),
    values.inicov=matrix(c(20, 0,
                            0, 5), ncol=2, byrow=TRUE), 
    params.inicov=matrix(c( 'v_x0','fixed',
                           'fixed','v_dx0'), ncol=2, byrow=T))

# Noises
mdcov <- prep.noise(
    values.latent=diag(0, 2),
    params.latent=diag(c("fixed","fixed"), 2),
    values.observed=diag(.5, 1),
    params.observed=diag("var_e", 1))


# dynamics	
formula = list(x ~ dx,
               dx ~ eta_i * x + zeta*dx)

#Specifies the structure of the fixed and random effects of the person-specific
#parameter, eta_i in the ODE (dynamic) formula
theta.formula  = list (eta_i ~ 1 * eta0  + u1 * eta1 + u2 * eta2 + 1 * b_eta)

dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(eta0=-1,
                                      eta1=.1,
                                      eta2=-.1,
                                      zeta=-.02),
                            isContinuousTime=TRUE,
						    theta.formula=theta.formula,
							random.names=c('b_eta'),
							random.params.inicov=matrix(c('sigma2_b_eta'), ncol=1,byrow=TRUE),
							random.values.inicov=matrix(c(0.1), ncol=1,byrow=TRUE))

							
# -----Cooking materials ----

# Put all the recipes together in a Model Specification
model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data,
                    outfile="osc.cpp")


# Estimate free parameters		
fitted_model <- dynr.cook(model, verbose=FALSE)


# ----Examine results----

# True values:
# eta0 = -.8, eta1 = .1, eta2 = -.3
# zeta = -.01
# sigma2_b_eta = .1^2 = .01
summary(fitted_model)

# Check the correlation between estimated b and true b:
# estimated b in fitted_model@eta_smooth_final
# true b in the in the 6th column of the input data set
cor(fitted_model@eta_smooth_final[3,oscData$times==30],
    oscData$trueb[oscData$times==30])

	
	
# ----Other useful fucntions----

# Get the estimated parameters from a cooked model
coef(fitted_model)

# Get the log likelihood, AIC, and BIC from a cooked model
logLik(fitted_model)
AIC(fitted_model)
BIC(fitted_model)




