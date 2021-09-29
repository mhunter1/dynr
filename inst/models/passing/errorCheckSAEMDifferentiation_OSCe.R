#------------------------------------------------------------------------------
# Author: Hui-Ju Hung
# Date: 2021-07-22
# Filename: errorCheckSAEMDifferentiation_OSCe.R
# Purpose: Check that errors related to SAEM differntiation matrices
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

require(dynr)

#------------------------------------------------------------------------------
# example model: Damped OSC and its extension



data(Oscillator)
Oscillator$y2 = 2 + 1.2*Oscillator$y1 + rnorm(length(Oscillator$y1),0,2)
Oscillator$u1 <- 1
Oscillator$u2 <- 0
Oscillator$u3 <- 0

data <- dynr.data(Oscillator, id="id", time="times", observed=c("y1","y2"),
                  covariates=c("u1", "u2","u3"))


initial <- prep.initial(
  values.inistate=c(0, 1,0, 1),
  params.inistate=c('inipos1', 'fixed','inipos2', 'inidx2'), #initial position is free parameter 5, initial slope is fixed at 1
  values.inicov=diag(1, 4),
  params.inicov=diag('fixed', 4))

mdcov <- prep.noise(
  values.latent=diag(c(0, 1,0,1), 4), params.latent=diag(c('fixed', 'dnoise','fixed', 'dnoise2'), 4), # uses free parameter 3
  values.observed=diag(1.5, 2), params.observed=diag(c('mnoise1','mnoise2'), 2)) # uses free parameter 4


formula=
  list(x1 ~ x2,
       x2 ~ -k1_i * (x1-base1_i) -zeta1_i* x2 + c1*cos(x1),
       x3 ~ x4,
       x4 ~ -k2_i * (x3-base2_i) -zeta2_i* x1
  )

theta.formula  = list (zeta1_i ~ 1 * zeta0 + u1 * zeta1 +  1 * b_zeta,
                       k1_i~1*k10+1*b_k1,
                       k2_i~1*k20+ u3*k21 + 1*b_k2,
                       base1_i~1*base10+ u3*base11 + 1*b_base1,
                       base2_i~1*base20+ u2*base21 + 1*b_base2)


dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(k10=1,base10=0, base11 = 1.2,
                                      base20=0, base21 = .4,
                                      k20=1,k21 = .5, 
                                      zeta0= 3,
                                      zeta1=.5, c1=2,zeta2_i=0.2),
                           isContinuousTime=TRUE,
                           theta.formula=theta.formula,
                           random.names=c('b_zeta',"b_k1","b_k2",
                                          'b_base1','b_base2'),
                           random.params.inicov = 
                             matrix(c('sigma2_b_zeta',0,'c13',0,0,
                                      0,'sigma2_b_k1',0,'c24',0,
                                      'c13',0,'sigma2_b_k2',0,0,
                                      0,'c24',0,'sigma2_b_base1','c45',
                                      0,0,0,'c45','sigma2_b_base1'),ncol=5,byrow=TRUE),
                           random.values.inicov = 
                              matrix(c(2,0,.5,0,0,
                                        0,.8,0,.001,0,
                                        .5,0,1,0,0,
                                        0,.001,0,5,.8,
                                        0,0,0,.8,1.2),ncol=5,byrow=TRUE),
                           random.lb = -10, 
                           random.ub = 2,
                           saem=TRUE)
						   
#------------------------------------------------------------------------------

meas <- prep.measurement(
  values.load=matrix(c(1, 0), 1, 2), # starting values and fixed values
  params.load=matrix(c('fixed', 'fixed'), 1, 2),
  state.names=c("x1","x2"))

					
# Examination of prep.measurement: state names are in prep.measurement and prep.formulaDynamics are not matched
testthat::expect_error(model_w <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data, 
                    outfile="do.cpp"), regexp="Found (4) latent states in dynamics formula, but expected (2) latent states from measurement model.", fixed=TRUE)

#------------------------------------------------------------------------------
meas <- prep.measurement(
  values.load=matrix(c(1, 0, 1, 0, 1, 0, 1, 0), 2, 4), # starting values and fixed values
  params.load=matrix(c('a', 'fixed', 0, 0, 'fixed', 'fixed', 0, 0), 2, 4),
  state.names=c("x1","x2", 'x3', 'x4'),
  obs.names=c("y1", "y2"),
  values.int=matrix(c(0,1),ncol=1),
  params.int=matrix(c('fixed','mu_x2'), ncol=1))

testthat::expect_warning( model <- dynr.model(dynamics=dynm, measurement=meas,
							       noise=mdcov, initial=initial, data=data,  
                                   outfile="do.cpp"),
						  regexp="Currently you can use the SAEM with ordinary differential equations (i.e., null matrix for the \"values.latent\" in dynrNoise)", fixed=TRUE) 
					
# examine the dimension of differentiation of dynamic formulas
testthat::expect_true(all(dim(model@dynamics@jacobian[[1]]) - c(4,4)==0))
testthat::expect_true(all(dim(model@dynamics@jacobianOriginal[[1]]) - c(4,4)==0))
testthat::expect_true(all(dim(model@dynamics@dfdtheta[[1]]) - c(5,4)==0))
testthat::expect_true(all(dim(model@dynamics@dfdxdtheta[[1]]) - c(16,5)==0))
testthat::expect_true(all(dim(model@dynamics@dfdthetadx[[1]]) - c(20,4)==0))
testthat::expect_true(all(dim(model@dynamics@dfdx2[[1]]) - c(16,4)==0))
testthat::expect_true(all(dim(model@dynamics@dfdtheta2[[1]]) - c(20,5)==0))

# examine the dimension of differentiation variance matrices
testthat::expect_true(all(dim(model@dLambdparLamb) - c(1,8)==0))
testthat::expect_true(all(dim(model@dLambdparLamb2) - c(8,1)==0))
testthat::expect_true(all(dim(model@dmudparMu) - c(1,2)==0))
testthat::expect_true(all(dim(model@dmudparMu2) - c(2,1)==0))
testthat::expect_true(all(dim(model@dmudparMu) - c(1,2)==0))
testthat::expect_true(all(dim(model@dmudparMu2) - c(2,1)==0))
testthat::expect_true(all(dim(model@dSigmaede) - c(2,4)==0))
testthat::expect_true(all(dim(model@dSigmaede2) - c(8,2)==0))
testthat::expect_true(all(dim(model@dSigmabdb) - c(8,81)==0))
testthat::expect_true(all(dim(model@dSigmabdb2) - c(648,8)==0))


#------------------------------------------------------------------------------
