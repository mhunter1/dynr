#------------------------------------------------------------------------------
# Author: Hui-Ju Hung
# Date: 2021-07-22
# Filename: errorCheckSAEMDifferentiation_VDP.R
# Purpose: Check that errors related to SAEM differntiation matrices
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

require(dynr)

#------------------------------------------------------------------------------
# example model: Damped OSC and its extension



data(vdpData)
data <- dynr.data(vdpData, id="id", time="time",
                 observed=c('y1', 'y2', 'y3'),
                 covariates=c("u1", "u2"))

meas <- prep.measurement(
    #values.load=matrix(c(1, -0.0087, -0.0126, 0, 0, 0), 3, 2), #original setting
    values.load=matrix(c(1, 0.6994, 1.1974, 0, 0, 0), 3, 2), # in SAEMTesting
    params.load=matrix(c('fixed', 'lambda_21', 'lambda_31', 'fixed', 'fixed', 'fixed'), 3, 2),
    obs.names = c('y1', 'y2', 'y3'),
    state.names=c('x1', 'x2'),
    #values.int=matrix(c(0, 0, 0), ncol=1), #original setting
	values.int=matrix(c(-0.0018, -0.0051, -0.0007), ncol=1), #SAEMTesting setting
    params.int=matrix(c('mu1', 'mu2', 'mu3'), ncol=1))


initial <- prep.initial(
    values.inistate=c(1, 1),
    #params.inistate=c("mu_x1", "mu_x2"),
    params.inistate=c("fixed", "fixed"),
    #values.inicov=matrix(c( 0, -1.204,
    #                        -1.204,0), ncol=2, byrow=TRUE), #enter values in unconstrainted scale(exp(0) = 1, exp(-1.204)=0.3)
    values.inicov=matrix(c( 1, 0.3,
                            0.3,1), ncol=2, byrow=TRUE), #enter values in constrainted scale(exp(0) = 1, exp(-1.204)=0.3)
    #values.inicov=matrix(c(1.14, .26,
    #                        .26,1.15), ncol=2, byrow=TRUE), #original setting in SAEM
    #params.inicov=matrix(c('sigma2_bx1','sigma_bx1x2',
    #                       'sigma_bx1x2','sigma2_bx2'), ncol=2, byrow=TRUE)
    
	params.inicov=matrix(c('fixed','fixed',
                           'fixed','fixed'), ncol=2, byrow=TRUE) 
)

mdcov <- prep.noise(
    values.latent=diag(0, 2), 
    params.latent=diag(c("fixed","fixed"), 2),
    #values.observed=diag(rep(-0.693,3)), # enter values in unconstrained scale (exp(-0.693) = 0.5)
    #values.observed=diag(rep(0.5,3)), # enter values in unconstrained scale (exp(-0.693) = 0.5)
	values.observed=diag(c(0.5076, 0.5027, 0.5140)),
    params.observed=diag(c("var_1","var_2","var_3"),3)
)

formula=
    list(x1 ~ x2,
         x2 ~ -61.68503 * x1 + zeta_i * (1 - x1^2) * x2
    )

theta.formula  = list (zeta_i ~ 1 * zeta0 + u1 * zeta1 + u2 * zeta2 + 1 * b_zeta)





dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(zeta0= 3,
                                      zeta1=.5,
                                      zeta2=.5),
                           isContinuousTime=TRUE,
                           theta.formula=theta.formula,
                           random.names=c('b_zeta'),
                           random.params.inicov = matrix(c('sigma2_b_zeta'), ncol=1,byrow=TRUE),
                           random.values.inicov = matrix(c(0.5), ncol=1,byrow=TRUE),
						   #random.values.inicov = matrix(c(1.0227), ncol=1,byrow=TRUE), original setting
                           random.lb = -5, 
                           random.ub = 5,
                           saem=TRUE
                           )



model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data, #saem=TRUE, 
                    #outfile=paste0(tempfile(),'.cpp'),
                    outfile="vdp_.cpp")

#------------------------------------------------------------------------------

					
# examine the dimension of differentiation of dynamic formulas
testthat::expect_true(all(dim(model@dynamics@jacobian[[1]]) - c(2,2)==0))
testthat::expect_true(all(dim(model@dynamics@jacobianOriginal[[1]]) - c(2,2)==0))
testthat::expect_true(all(dim(model@dynamics@dfdtheta[[1]]) - c(1,2)==0))
testthat::expect_true(all(dim(model@dynamics@dfdxdtheta[[1]]) - c(4,1)==0))
testthat::expect_true(all(dim(model@dynamics@dfdthetadx[[1]]) - c(2,2)==0))
testthat::expect_true(all(dim(model@dynamics@dfdx2[[1]]) - c(4,2)==0))
testthat::expect_true(all(dim(model@dynamics@dfdtheta2[[1]]) - c(2,1)==0))

# examine the dimension of differentiation variance matrices
testthat::expect_true(all(dim(model@dLambdparLamb) - c(2,6)==0))
testthat::expect_true(all(dim(model@dLambdparLamb2) - c(12,2)==0))
testthat::expect_true(all(dim(model@dmudparMu) - c(3,3)==0))
testthat::expect_true(all(dim(model@dmudparMu2) - c(9,3)==0))
testthat::expect_true(all(dim(model@dSigmaede) - c(3,9)==0))
testthat::expect_true(all(dim(model@dSigmaede2) - c(27,3)==0))
testthat::expect_true(all(dim(model@dSigmabdb) - c(1,1)==0))
testthat::expect_true(all(dim(model@dSigmabdb2) - c(1,1)==0))


#------------------------------------------------------------------------------
{"mode":"full","isActive":false}