#------------------------------------------------------------------------------
# Author: Hui-Ju Hung
# Date: 2019-05-23
# Filename: SAEMVDPFree.R
# Purpose: Model script for dynr/SAEM gateway function, Van Der Pol fucntion 
#------------------------------------------------------------------------------

library('dynr')
library('plyr')


# nPeople = 200
# nTimes = 300
# vdpData <- read.csv("../data/TrueInitY1.txt", header=FALSE)
# colnames(vdpData) <- c('batch', 'kk', 'trueInit', 'time', 'y1','y2','y3', 'u1', 'u2')
# vdpData$id <- rep(1:nPeople, each=nTimes)


data(vdpData)
data <- dynr.data(vdpData, id="id", time="time",
                 observed=c('y1', 'y2', 'y3'),
                 covariates=c("u1", "u2"))

meas <- prep.measurement(
    values.load=matrix(c(1, 0.7, 1.2, 0, 0, 0), 3, 2),
    params.load=matrix(c('fixed', 'lambda_21', 'lambda_31', 'fixed', 'fixed', 'fixed'), 3, 2),
    obs.names = c('y1', 'y2', 'y3'),
    state.names=c('x1', 'x2'),
    values.int=matrix(c(0, 0, 0), ncol=1),
    params.int=matrix(c('mu1', 'mu2', 'mu3'), ncol=1))


initial <- prep.initial(
    values.inistate=c(3, 1),
    params.inistate=c("mu_x1", "mu_x2"),
    #params.inistate=c("fixed", "fixed"),
    values.inicov=matrix(c(1.14, .26,
                            .26,1.15), ncol=2, byrow=TRUE), 
    params.inicov=matrix(c('sigma2_bx1','sigma_bx1x2',
                           'sigma_bx1x2','sigma2_bx2'), ncol=2, byrow=TRUE)
    
	#params.inicov=matrix(c('fixed','fixed',
    #                       'fixed','fixed'), ncol=2, byrow=TRUE) 
)

mdcov <- prep.noise(
    values.latent=diag(0, 2), 
    params.latent=diag(c("fixed","fixed"), 2),
    values.observed=diag(rep(0.5,3)),
    params.observed=diag(c("var_1","var_2","var_3"),3)
)

formula=
    list(x1 ~ x2,
         x2 ~ -61.68503 * x1 + zeta_i * (1 - x1^2) * x2
    )

theta.formula  = list (zeta_i ~ 1 * zeta0  + u1 * zeta1 + u2 * zeta2 + 1 * b_zeta)





dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(zeta0= 3,
                                      zeta1=.5,
                                      zeta2=.5),
                           isContinuousTime=TRUE,
                           theta.formula=theta.formula,
                           random.names=c('b_zeta'),
                           random.params.inicov = matrix(c('sigma2_b_zeta'), ncol=1,byrow=TRUE),
                           random.values.inicov = matrix(c(0.9), ncol=1,byrow=TRUE),
                           random.lb = -5, 
                           random.ub = 5,
                           saem=TRUE
                           )



model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data, #saem=TRUE, 
                    outfile="VanDerPol.cpp")

print('here')
print(model@freeIC)
#-10
#exp(-10) 
#InfDS.lowBound(11:16) = log(10e-8);


#Using $ with lb and ub assumes that you are adjusting the parameters
#model@lb[names(model@lb) %in% c("var_1","var_2","var_3")] = log(10e-8)
model$ub[names(model@ub)] = 10
print(model@random.params.inicov)

saemp <- prep.saemParameter(MAXGIB = 10, 
                            MAXITER = 100, 
                            maxIterStage1 = 1005, 
                            gainpara = 0.600000, 
                            gainparb = 3.000000, 
                            gainpara1 = 0.900000, 
                            gainparb1 = 1.000000, 
                            bAdaptParams = c(0.5, 2.5, 0.5),
							KKO=30)

timestart<-Sys.time()

fitted_model <- dynr.cook(model, optimization_flag = TRUE, hessian_flag = TRUE, verbose=TRUE, debug_flag=TRUE, saemp = saemp)

timeend<-Sys.time()

print(timeend-timestart)


