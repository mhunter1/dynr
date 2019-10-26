#------------------------------------------------------------------------------
# Author: Hui-Ju Hung
# Date: 2019-10-03
# Filename: VDPwithRand.R
# Purpose: An illustrative example for using dynr to estimate random effects  
#          for Van Der Pol model.
#------------------------------------------------------------------------------

library('dynr')

# ---- Read in the data ----

data(vdpData)
data <- dynr.data(vdpData, id="id", time="time",
                 observed=c('y1','y2','y3'),
                 covariates=c("u1","u2"))

# ---- Prepare the recipes (i.e., specifies modeling functions) ----

# Measurement (factor loadings)	
meas <- prep.measurement(
	values.load=matrix(c(1, 1, 1, 0, 0, 0), 3, 2),
    params.load=matrix(c('fixed', 'lambda_21', 'lambda_31', 'fixed', 'fixed', 'fixed'), 3, 2),
    obs.names = c('y1', 'y2', 'y3'),
    state.names=c('x1', 'x2')) 

# Initial conditions on the latent state and covariance
initial <- prep.initial(
    values.inistate=c(3, 1),
    params.inistate=c("mu_x1", "mu_x2"),
    values.inicov=matrix(c(.5,.2,
                           .2,.6),ncol=2,byrow=TRUE), 
    params.inicov=matrix(c('sigma2_bx1','sigma_bx1x2',
                           'sigma_bx1x2','sigma2_bx2'),ncol=2,byrow=TRUE))

# Noises
mdcov <- prep.noise(
    values.latent=diag(0, 2),
    params.latent=diag(c("fixed","fixed"), 2),
    values.observed=diag(rep(0.3,3)),
    params.observed=diag(c("var_1","var_2","var_3"),3)
)


# Dynamics	
formula = list(x1 ~ x2,
               x2 ~ -61.68503 * x1 + zeta_i * (1 - x1^2) * x2)

#To have random effect in the parameter, zeta_i, specified via theta.formula
theta.formula  = list (zeta_i ~ 1 * zeta0  + u1 * zeta1 + u2 * zeta2 + 1 * b_zeta)

dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(zeta0=-1,
                                      zeta1=.5,
                                      zeta2=.2),
                           isContinuousTime=TRUE,
						   theta.formula=theta.formula,
						   random.names=c('b_zeta'),
						   random.params.inicov=matrix(c('sigma2_b_zeta'), ncol=1,byrow=TRUE),
						   random.values.inicov = matrix(c(0.9), ncol=1,byrow=TRUE))

						   
# -----Cooking materials ----

# Put all the recipes together in a Model Specification								
model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data, #saem=TRUE,
                    outfile="VanDerPol.cpp")


# Estimate free parameters	
fitted_model <- dynr.cook(model, verbose=FALSE)



# ----Examine results----

# True values:
# zeta0 = 3, zeta1 = .5, zeta2 = .5     
# lambda_21 = .7, lambda_31 = 1.2
# var_1 = .5, var_2 = .5, var_3 = .5
# mu_x1 = 0, mu_x2 = 0
# sigma2_bx1 = 1, sigma_bx1x2 = .3, sigma2_bx2 = 1, sigma2_b_zeta = .5

summary(fitted_model)

# Check the correlation between estimated b and true b:
# estimated b in fitted_model@eta_smooth_final
# true b in the in the 10th column of the input data set
cor(fitted_model@eta_smooth_final[3,data$tstart[2:51]],
    vdpData$trueb[data$tstart[2:51]])



# ----Other useful fucntions----

# Get the estimated parameters from a cooked model
coef(fitted_model)

# Get the log likelihood, AIC, and BIC from a cooked model
logLik(fitted_model)
AIC(fitted_model)
BIC(fitted_model)

