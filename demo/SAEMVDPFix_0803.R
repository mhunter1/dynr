#------------------------------------------------------------------------------
# Author: Hui-Ju Hung
# Date: 2019-05-23
# Filename: SAEMVDPFree.R
# Purpose: Model script for dynr/SAEM gateway function, Van Der Pol fucntion 
#------------------------------------------------------------------------------

library('dynr')
library('plyr')

#setwd("C:/Users/Cynthia/Documents/gits/dynr/")

nPeople = 200
#nTimes = 300
#vdpData <- read.csv("./data/VDPFix.rda", header=FALSE)
#colnames(vdpData) <- c('batch', 'kk', 'trueInit', 'time', 'y1','y2','y3', 'u1', 'u2')
#vdpData$id <- rep(1:nPeople, each=nTimes)
data(VDPFix)


#data(vdpData)
data <- dynr.data(vdpData, id="id", time="time",
                  observed=c('y1', 'y2', 'y3'),
                  covariates=c("u1", "u2"))

trueb <- data$trueb[data$tstart[1:nPeople+1], ]

#truebData <- read.csv("./data/truebFile220803.txt", header=FALSE)
#trueb <- truebData[,4]

meas <- prep.measurement(
  #values.load=matrix(c(1, -0.0087, -0.0126, 0, 0, 0), 3, 2), #original setting
  values.load=matrix(c(1, 0.7, 1.2, 0, 0, 0), 3, 2), # in SAEMTesting
  params.load=matrix(c('fixed', 'fixed', 'fixed', 'fixed', 'fixed', 'fixed'), 3, 2),
  obs.names = c('y1', 'y2', 'y3'),
  state.names=c('x1', 'x2'),
  #values.int=matrix(c(0, 0, 0), ncol=1), #original setting
  values.int=matrix(c(0,0,0), ncol=1), #SAEMTesting setting
  params.int=matrix(c('fixed', 'fixed', 'fixed'), ncol=1))


initial <- prep.initial(
  values.inistate=c(1,1),
  params.inistate=c("mu_x1", "mu_x2"),
  #params.inistate=c("fixed", "fixed"),
  values.inicov=matrix(c( 1, .3,
                          .3, 1), ncol=2, byrow=TRUE), #enter values in unconstrainted scale(exp(0) = 1, exp(-1.204)=0.3)
  #values.inicov=matrix(c( 0, 0,
  #                        0,0), ncol=2, byrow=TRUE), #enter values in constrainted scale(exp(0) = 1, exp(-1.204)=0.3)
  #values.inicov=matrix(c(1.14, .26,
  #                        .26,1.15), ncol=2, byrow=TRUE), #original setting in SAEM
  #params.inicov=matrix(c('sigma2_bx1','sigma_bx1x2',
  #                       'sigma_bx1x2','sigma2_bx2'), ncol=2, byrow=TRUE)
  params.inicov=matrix(c('fixed','fixed',
                         'fixed','fixed'), ncol=2, byrow=TRUE) 
)


sampleCovformula=
  list(Sigma11 ~ par1*u1^3,
       Sigma12 ~ par2*u1, 
       Sigma22 ~ par3*u1^3)

# Todo: mdcov of startval needs to transformed to unconstrained scale
mdcov <- prep.noise(
  values.latent=diag(0, 2), 
  params.latent=matrix(c('Sigma11', 'Sigma12', 'Sigma12', 'Sigma22'), nrow = 2, byrow = TRUE),
  #values.observed=diag(rep(-0.693,3)), # enter values in unconstrained scale (exp(-0.693) = 0.5)
  #values.observed=diag(rep(0.5,3)), # enter values in unconstrained scale (exp(-0.693) = 0.5)
  values.observed=diag(c(0.5, 0.5, 0.5)),
  params.observed=diag(c("var1","var2","var3"),3),
  latent.formula = sampleCovformula,
  covariates = c("u1"),
  latent.startval = c(par1=-.6, par2=-.7, par3 = -.3)
)

#ldl.transformed = matrix(c(0.3726659, 0.3119856, 0.3119856,1.480466),2,2)
#ldl.transformed (subj 1)
#          [,1]      [,2]
#[1,] 0.3726659 0.3119856
#[2,] 0.3119856 1.4804660
#ldl.transformed (subj 50)
#          [,1]      [,2]
#[1,] 0.0523459 0.0211216
#[2,] 0.0211216 0.3783712

#UI problem
#1. Without specifying params.latent=matrix(c('Sigma11', 'Sigma12', 'Sigma12', 'Sigma22'), nrow = 2, byrow = TRUE),
#   is it okay to assume that a user will know that it would be a matrix form in column/row major?


formula=
  list(x1 ~ x2,
       x2 ~ -61.68503 * x1 + (1 * 3 + u1 * 0.5 + u2 * 0.5 + 0.2) * (1 - x1^2) * x2
  )

theta.formula  = list (zeta_i ~  1 * b_zeta)





dynm<-prep.formulaDynamics(formula=formula,
                           #startval=c(zeta0= 3,
                           #           zeta1=.5,
                           #           zeta2=.5),
                           isContinuousTime=TRUE,
                           #theta.formula=theta.formula,
                           #random.names=c('b_zeta'),
                           #random.params.inicov = matrix(c('fixed'), ncol=1,byrow=TRUE),
                           #random.values.inicov = matrix(c(0.5), ncol=1,byrow=TRUE),
                           #random.values.inicov = matrix(c(1.0227), ncol=1,byrow=TRUE), original setting
                           #random.lb = -10, 
                           #random.ub = 10,
                           covariate.formula = sampleCovformula,
                           covariate.names = c('delta_t'),
                           saem=FALSE
)



model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data, #saem=TRUE, 
                    #outfile=paste0(tempfile(),'.cpp'))
                    outfile="vdp7.cpp"
                    #ub =c(zeta0=10,zeta1=10,lambda_21 =10, lambda_31=10, mu1=10, mu2=10, mu3=10, var_1=10, var_2=10, var_3=10)
)

#model@ub = c(zeta0=10,zeta1=10,zeta2=10,lambda_21 =10, lambda_31=10, mu1=10, mu2=10, mu3=10, var_1=10, var_2=10, var_3=10, sigma2_b_zeta=10)
#model@ub[10:12] = log(5)
#model@lb[1:9] = -5
#model@lb[10:12] = -log(5)
#setwd("C:/Users/Cynthia/Documents/gits/dynr/temp")
#model@outfile = "vdp4.cpp"
print(model@freeIC)
print(model@xstart)
print(model$xstart)
#-10
#exp(-10) 
#InfDS.lowBound(11:16) = log(10e-8);


#Using $ with lb and ub assumes that you are adjusting the parameters
#model@lb[names(model@lb) %in% c("var_1","var_2","var_3")] = log(10e-8)
#model$ub[names(model@ub)] = 10
#model$ub[names(model@lb)] = 10
#print(model@random.params.inicov)

saemp <- prep.saemParameter(MAXGIB = 1, 
                            MAXITER = 1, 
                            seed = 9,
                            setScaleb = 0, trueb =as.matrix(trueb)
                            
                            #setAccept= 0.7
                            #scaleb = 10,
                            #maxIterStage1 = 100, 
                            #gainpara = 0.600000, 
                            #gainparb = 3.000000, 
                            #gainpara1 = 0.900000, 
                            #gainparb1 = 1.000000, 
                            #bAdaptParams = c(5, 2.5, 0.5),
                            #KKO=5
)

timestart<-Sys.time()
#setwd("C:/Users/Cynthia/Documents/gits/dynr/temp")
#fitted_model <- dynr.cook(model, optimization_flag = FALSE, hessian_flag = FALSE, verbose=TRUE, debug_flag=TRUE, saemp = saemp)
fitted_model <- dynr.cook(model, optimization_flag = FALSE, hessian_flag = FALSE, verbose=TRUE, debug_flag=TRUE, maxeval = 1)
#print(fitted_model)

timeend<-Sys.time()

print(timeend-timestart)

#save.image('~/Dropbox/Brekfis/SAEM/Armadillo/SymiinTestExamples/VDP_test.Rdata')


#vdpTrueX <- read.csv("C:\\Users\\Cynthia\\Dropbox\\Brekfis\\SAEM\\Armadillo\\SymiinTestExamples\\OutputforDebugging\\NewTrueInit_trueXG1.txt", header=FALSE)
#colnames(vdpTrueX) <- c('batch', 'kk', 'trueInit', 'time', 'id','truex','truedx')
#View(vdpTrueX)

#fitted_model$dXtildAll[t][[1]][theta,nx, nsubj]
