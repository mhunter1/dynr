#------------------------------------------------------------------------------
# Author: Hui-Ju Hung
# Date: 2019-05-23
# Filename: vanderpol_3.R
# Purpose: Model script for dynr/SAEM gateway function, Van Der Pol fucntion 
# Note: Workable for developer dynr on arma branch
#------------------------------------------------------------------------------

#setwd('C:\\Users\\Cynthia\\Documents\\gits\\dynr\\demo')
library('dynr')
library('plyr')

state.names = c('x1', 'x2')
#beta.names = c('zeta0', 'zeta1', 'zeta2', 'mu_x1', 'mu_x2')
beta.names = c('zeta0', 'zeta1', 'zeta2')
covariate.names = c('u1', 'u2')
theta.names = c('zeta_i') 
#[todo] Z*b's b
random.names = c('b_zeta')
intercept.names = c('mu1', 'mu2', 'mu3')


#random data
# N = 200
# T = 300
# vdpData <- data.frame(id=rep(1:N,each=T), time=rep(seq(0.005,1.5,by=0.005),N),
                      # y1=rnorm(100*300),  y2=rnorm(100*300),  y3=rnorm(100*300),
                      # u1 = rnorm(100*300), u2 = rnorm(100*300))
# colnames(vdpData) <- c("id","time","y1","y2", "y3","u1", "u2") 
# data <- dynr.data(vdpData, id="id", time="time",
                   # observed=c('y1', 'y2', 'y3'),
                   # covariates=c('u1','u2'))

nPeople = 200
nTimes = 300
#vdpData <- read.csv("../data/TrueInitY1.txt", header=FALSE)
#vdpData <- read.csv("C:\\Users\\Cynthia\\Documents\\gits\\dynr\\data\\TrueInitY1.txt", header=FALSE)
vdpData <- read.csv("~/gits/dynr/data/TrueInit_Y1.txt", header=FALSE)

colnames(vdpData) <- c('batch', 'kk', 'trueInit', 'time', 'y1','y2','y3', 'u1', 'u2')
vdpData$id <- rep(1:nPeople, each=nTimes)
data <- dynr.data(vdpData, id="id", time="time",
                 observed=c('y1','y2','y3'),
                 covariates=c("u1","u2"))

meas <- prep.measurement(
	values.load=matrix(c(1, 1, 1, 0, 0, 0), 3, 2),
    params.load=matrix(c('fixed', 'lambda_21', 'lambda_31', 'fixed', 'fixed', 'fixed'), 3, 2),
    obs.names = c('y1', 'y2', 'y3'),
    state.names=state.names) #,
	#values.int=c(3, 1, 0),
	#params.int=intercept.names) #intercept.names = c('mu1', 'mu2', 'mu3')


initial <- prep.initial(
    values.inistate=c(3, 1),
    params.inistate=c("mu_x1", "mu_x2"),
    values.inicov=matrix(c(.5,.2,
                           .2,.6),ncol=2,byrow=TRUE), 
    params.inicov=matrix(c('sigma2_bx1','sigma_bx1x2',
                           'sigma_bx1x2','sigma2_bx2'),ncol=2,byrow=T)
)

mdcov <- prep.noise(
    values.latent=diag(0, 2),
    params.latent=diag(c("fixed","fixed"), 2),
    values.observed=diag(rep(0.3,3)),
    params.observed=diag(c("var_1","var_2","var_3"),3)
)

formula=
    list(x1 ~ x2,
         x2 ~ -61.68503 * x1 + zeta_i * (1 - x1^2) * x2 #,
#         zeta0 ~0,
#         zeta1 ~0,
#         zeta2 ~0,
#         mu_x1 ~0,
#         mu_x2 ~0
    )

theta.formula  = list (zeta_i ~ 1 * zeta0  + u1 * zeta1 + u2 * zeta2 + 1 * b_zeta)


#theta.formula = list( zeta_i ~ 1 * zeta0 + u1 * zeta1 + u2 * zeta2 + 1 * b_zeta,
#                      x1_0 ~ 1 * mu_x1 + 1 * b_x1,
#                      x2_0 ~ 1 * mu_x2 + 1 * b_x2)




dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(zeta0=-1,
                                           zeta1=.5,
                                           zeta2=.2),
                                isContinuousTime=TRUE,
								#state.names=state.names,
								theta.formula=theta.formula,
								#theta.names=theta.names,
								#beta.names=beta.names,
								#intercept.names=intercept.names, 
								random.names=random.names,
								random.params.inicov = matrix(c('sigma2_b_zeta'), ncol=1,byrow=TRUE),
								random.values.inicov = matrix(c(0.9), ncol=1,byrow=TRUE),
							    random.lb = -5, 
				                random.ub = 5#,
								#saem=TRUE
								)

								
model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data, #saem=TRUE,
                    outfile="VanDerPol.cpp")


#-10
#exp(-10) 
#InfDS.lowBound(11:16) = log(10e-8);


#Using $ with lb and ub assumes that you are adjusting the parameters
#model@lb[names(model@lb) %in% c("var_1","var_2","var_3")] = log(10e-8)

#model$ub[names(model@ub)] = 10

#print(model@random.params.inicov)

saemp <- prep.saemParameter(MAXGIB = 10, MAXITER = 10, maxIterStage1 = 1005, gainpara = 0.600000, gainparb = 3.000000, gainpara1 = 0.900000, gainparb1 = 1.000000, bAdaptParams = c(0.5, 2.5, 0.5))
								
timestart<-Sys.time()

fitted_model <- dynr.cook(model, optimization_flag = TRUE, hessian_flag = TRUE, verbose=TRUE, debug_flag=TRUE, saemp = saemp)

timeend<-Sys.time()

print(timeend-timestart)
#print(fitted_model)



# -------
# previous perp.random interface
# ran <- prep.random(random.names=random.names, 
				   # num.subj = N, 
				   # random.lb = c(-.1, -.1, -.1), 
				   # random.ub = c(.1, .1, .1),
				   # params.inicov=matrix(c('b_zeta', 'c01', 'c02',
										     # 'c01','b_x1', 'c12',
                                             # 'c02', 'c12','b_x2'),ncol=3,byrow=T),
				   # values.inicov=matrix(c(1, 0, 0,
										  # 0,.5,.6,
                                          # 0,.6,.2),ncol=3,byrow=T))
#print(ran)

#model <- dynr.model(dynamics=dynm, measurement=meas,
#                    noise=mdcov, initial=initial, data=data, random=ran, armadillo=TRUE,
#                    outfile="VanDerPol.c")

# -------
# previous implemented function: write a matrix in R to Armadillo Code format
# matrix2ArmadilloCode <- function(variable.name, matrix.input){
#     str = paste0(variable.name, " = \"")
#     dim(matrix.input)
#     for(i in c(1:dim(matrix.input)[1])){
#         for(j in c(1:dim(matrix.input)[2])){
#             str = paste0(str, " ", as.character(matrix.input[i,j]))
#         }
#         str = paste0(str, ";")
#     }
#     str = paste0(str, "\";")
#     return(str)
# }

