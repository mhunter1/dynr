#------------------------------------------------------------------------------
# Author: Hui-Ju Hung
# Date: 2020-11-01
# Filename: errorCheckSAEM.R
# Purpose: Check that errors related to EstimateRandomAsLV.
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

require(dynr)

#------------------------------------------------------------------------------
# example model: OSC 



formula <- list(x ~ dx,
               dx ~ eta_i * x + zeta*dx)
theta.formula <- list (eta_i ~ 1 * eta0  + u1 * eta1 + u2 * eta2 + 1 * b_eta)

data(oscData)
data <- dynr.data(oscData, id="id", time="times",
                 observed=c('y1'),
                 covariates=c("u1","u2"))

meas <- prep.measurement(
    values.load=matrix(c(1,0), 1, 2),
    params.load=matrix(c('fixed','fixed'), 1, 2),
    obs.names=c('y1'),
    state.names=c('x', 'dx'))

initial <- prep.initial(
    values.inistate=c(10, 0),
    params.inistate=c("mu_x0", "mu_dx0"),
    values.inicov=matrix(c(20, 0,
                            0, 5), ncol=2, byrow=TRUE), 
    params.inicov=matrix(c( 'v_x0','fixed',
                           'fixed','v_dx0'), ncol=2, byrow=T))

mdcov <- prep.noise(
    values.latent=diag(0, 2),
    params.latent=diag(c("fixed","fixed"), 2),
    values.observed=diag(.5, 1),
    params.observed=diag("var_e", 1))

dynm <- prep.formulaDynamics(formula=formula, 
    startval=c(eta0=-1, eta1=.1, eta2=-.1, zeta=-.01),
    isContinuousTime=TRUE,
    theta.formula=theta.formula,
    random.names=c('b_eta'),
    random.params.inicov=matrix(c('sigma2_b_eta'), ncol=1,byrow=TRUE),
    random.values.inicov=matrix(c(0.1), ncol=1,byrow=TRUE))

model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data)

#------------------------------------------------------------------------------

# Examination of prep.formulaDynamics: slot names should be random.names instead of random.name
testthat::expect_error(dynm<-prep.formulaDynamics(formula=formula, theta.formula = theta.formula, startval=c(eta0=-1, eta1=.1, eta2=-.1, zeta=-.02), isContinuousTime=TRUE, random.name=c('b_eta')), regexp="You passed some invalid names to the ... argument. Check with US Customs or the ?prep.formulaDynamics help page.", fixed=TRUE)

#------------------------------------------------------------------------------

# Examination of dynr.cook: if random.names are not specified in model@dynamics@random.names
dynm_w <- prep.formulaDynamics(formula=formula, 
    startval=c(eta0=-1, eta1=.1, eta2=-.1, zeta=-.02),
    isContinuousTime=TRUE,
    theta.formula=theta.formula,
    random.params.inicov=matrix(c('sigma2_b_eta'), ncol=1,byrow=TRUE),
    random.values.inicov=matrix(c(0.1), ncol=1,byrow=TRUE))



testthat::expect_error(model_w <- dynr.model(dynamics=dynm_w, measurement=meas,
                                             noise=mdcov, initial=initial, data=data), regexp="In theta.formula, there are variables that are not specified in startvals: b_eta", fixed=TRUE)

#------------------------------------------------------------------------------

# Examination of dynr.model: Check that number of assigned values from coef<- is the same as number available
testthat::expect_error(
    coef(model) <- rep(.1, 5),
    regexp="Number of model coeficients (10) does not match number assigned (5).", fixed=TRUE)

#------------------------------------------------------------------------------
#fix seed
set.seed(42)
	
# Examination of dynr.cook: the results
fitted_model <- dynr.cook(model, verbose=FALSE)


# True values:
# eta0 = -.8, eta1 = .1, eta2 = -.3
# zeta = -.01
# sigma2_b_eta = .1^2 = .01

# compare true parameters to estimated ones
trueParams <- c(eta0 = -.8, eta1 = .1, eta2 = -.3, zeta = -.01 , sigma2_b_eta= .01)
data.frame(name=c('eta0', 'eta1', 'eta2', 'zeta', 'sigma2_b_eta'), true=trueParams, estim=coef(fitted_model)[c('eta0', 'eta1', 'eta2', 'zeta', 'sigma2_b_eta')])

(CI <- confint(fitted_model, level = 0.99)[c('eta0', 'eta1', 'eta2', 'zeta', 'sigma2_b_eta'),])

# Check that all true parameters are within the confidence intervals of the estimated params
withinIntervals <- CI[,1] < trueParams & trueParams < CI[,2]
withinIntervals
testthat::expect_true(all(withinIntervals))

#------------------------------------------------------------------------------
# test on VDP model
data(vdpData)
#vdpData = vdpData[vdpData$id<=30,]
data <- dynr.data(vdpData, id="id", time="time",
                 observed=c('y1','y2','y3'),
                 covariates=c("u1","u2"))

meas <- prep.measurement(
    values.load=matrix(c(1, 1, 1, 0, 0, 0), 3, 2),
    params.load=matrix(c('fixed', 'lambda_21', 'lambda_31', 'fixed', 'fixed', 'fixed'), 3, 2),
    obs.names = c('y1', 'y2', 'y3'),
    state.names=c('x1', 'x2')) 

initial <- prep.initial(
    values.inistate=c(0, 0),
    params.inistate=c("mu_x1", "mu_x2"),
    values.inicov=matrix(c(.8,.3,
                           .3,.7),ncol=2,byrow=TRUE), 
    params.inicov=matrix(c('sigma2_bx1','sigma_bx1x2',
                           'sigma_bx1x2','sigma2_bx2'),ncol=2,byrow=TRUE))

mdcov <- prep.noise(
    values.latent=diag(0, 2),
    params.latent=diag(c("fixed","fixed"), 2),
    values.observed=diag(rep(0.6,3)),
    params.observed=diag(c("var_1","var_2","var_3"),3)
)


formula <- list(x1 ~ x2,
               x2 ~ -61.68503 * x1 + zeta_i * (1 - x1^2) * x2)

theta.formula <- list (zeta_i ~ 1 * zeta0  + u1 * zeta1 + u2 * zeta2 + 1 * b_zeta)

dynm <- prep.formulaDynamics(formula=formula,
    startval=c(zeta0=-3,
               zeta1=.8,
               zeta2=.2),
    isContinuousTime=TRUE,
    theta.formula=theta.formula,
    random.names=c('b_zeta'),
    random.params.inicov=matrix(c('sigma2_b_zeta'), ncol=1, byrow=TRUE),
    random.values.inicov = matrix(c(0.9), ncol=1, byrow=TRUE))


model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data)
msg <- paste0("Number of model coeficients (", length(coef(model)),
    ") does not match number assigned (7).")
testthat::expect_error(
    coef(model) <- rep(.1, 7),
    regexp=msg, fixed=TRUE)



# Examination of dynr.cook: the results
fitted_model <- dynr.cook(model, verbose=FALSE)


# True values:
# zeta0 = 3, zeta1 = .5, zeta2 = .5     
# lambda_21 = .7, lambda_31 = 1.2
# var_1 = .5, var_2 = .5, var_3 = .5
# mu_x1 = 0, mu_x2 = 0
# sigma2_bx1 = 1, sigma_bx1x2 = .3, sigma2_bx2 = 1, sigma2_b_zeta = .5

# compare true parameters to estimated ones
trueParams <- c(zeta0 = 3, zeta1 = .5, zeta2 = .5, lambda_21 = .7, lambda_31 = 1.2, var_1 = .5, var_2 = .5, var_3 = .5, mu_x1 = 0, mu_x2 = 0, sigma2_bx1 = 1, sigma_bx1x2 = .3, sigma2_bx2 = 1, sigma2_b_zeta = .5)
#trueParams <- c(zeta0 = 3, zeta1 = .5, zeta2 = .5, lambda_21 = .7, lambda_31 = 1.2, var_1 = .5, var_2 = .5, var_3 = .5, sigma2_b_zeta = .5)


data.frame(name=c('zeta0', 'zeta1', 'zeta2', 'lambda_21', 'lambda_31', 'var_1', 'var_2', 'var_3', 'mu_x1', 'mu_x2','sigma_bx1', 'sigma_bx1x2', 'sigma_bx2','sigma2_b_zeta'), true=trueParams, estim=coef(fitted_model))
#data.frame(name=c('zeta0', 'zeta1', 'zeta2', 'lambda_21', 'lambda_31','var_1', 'var_2', 'var_3','sigma2_b_zeta'), true=trueParams, estim=coef(fitted_model))


(CI <- confint(fitted_model, level=0.99))

# Check that all true parameters are within the confidence intervals of the estimated params
withinIntervals <- CI[,1] < trueParams & trueParams < CI[,2]
withinIntervals
testthat::expect_true(all(withinIntervals))

#------------------------------------------------------------------------------
