#------------------------------------------------------------------------------
# Author: Hui-Ju Hung
# Date: 2019-05-23
# Filename: vanderpol_3.R
# Purpose: Model script for dynr/SAEM gateway function, Van Der Pol fucntion 
# Note: Workable for developer dynr on arma branch
#------------------------------------------------------------------------------

library('dynr')
state.names = c('x1', 'x2')
beta.names = c('zeta0', 'zeta1', 'zeta2', 'mu1', 'mu2')
covariate.names = c('u1', 'u2')
theta.names = c('zeta_i', 'zeta_x1', 'zeta_x2') 
#Do we need the user to set up the names of 'zeta_x1', 'zeta_x2'?
#[todo] Z*b's b
random.names = c('b_zeta')
intercept.names = c('mu1', 'mu2', 'mu3')


NxState = length(state.names)
Nbeta = length(beta.names)
Nx = NxState + Nbeta
Ntheta = length(theta.names)



#random data
N = 100
T = 300
vdpData <- data.frame(id=rep(1:N,each=T), time=rep(seq(0.005,1.5,by=0.005),N),
                      y1=rnorm(100*300),  y2=rnorm(100*300),  y3=rnorm(100*300),
                      u1 = rnorm(100*300), u2 = rnorm(100*300))
colnames(vdpData) <- c("id","time","y1","y2", "y3","u1", "u2") 
data <- dynr.data(vdpData, id="id", time="time",
                  observed=c('y1', 'y2', 'y3'),
                  covariates=c('u1','u2'))


meas <- prep.measurement(
    #values.load=matrix(c(1,0), 1, 2),
    #params.load=matrix(c('fixed'), 1, 2),
    values.load=matrix(c(1,0, 1,0,1,0), 3, 2),
    params.load=matrix(c('fixed'), 3, 2),
    obs.names = c('y1', 'y2', 'y3'),
    state.names=state.names,
	values.int=c(3, 1, 0),
	params.int=intercept.names) #intercept.names = c('mu1', 'mu2', 'mu3')


initial <- prep.initial(
    values.inistate=c(3, 1),
    params.inistate=c("mu_x1", "mu_x2"),
    values.inicov=matrix(c(.5,.2,
                           .2,.6),ncol=2,byrow=T), 
    params.inicov=matrix(c('sigma2_x10','c12',
                           'c12','sigma2_x20'),ncol=2,byrow=T)
)

mdcov <- prep.noise(
    values.latent=diag(0, 2),
    params.latent=diag(c("fixed","fixed"), 2),
    values.observed=diag(rep(0.3,2)),
    params.observed=diag(c("var_1","var_2"),2) #sigma_e
)


formula=
    list(x1 ~ x2,
         x2 ~ -61.68503 * x1 + zeta_i * (1 - x1^2) * x2,
         zeta0 ~0,
         zeta1 ~0,
         zeta2 ~0,
         mu_x1 ~0,
         mu_x2 ~0
    )
theta.formula  = list (zeta_i ~ 1 * zeta0  + u1 * zeta1 + u2 * zeta2 + 1 * b_zeta,
x1_0 ~ 1 * 0,
x2_0 ~ 1 * 0)

#theta.formula = list( zeta_i ~ 1 * zeta0 + u1 * zeta1 + u2 * zeta2 + 1 * b_zeta,
#                      x1_0 ~ 1 * mu_x1 + 1 * b_x1,
#                      x2_0 ~ 1 * mu_x2 + 1 * b_x2)

#theta.formula2 = prep.thetaFormula(theta.formula, intercept.names, random.names)
#print(theta.formula2)
dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(zeta0=-1,
                                           zeta1=.5,
                                           zeta2=.2),
                                isContinuousTime=FALSE,
								state.names=state.names,
								theta.formula=theta.formula,
								theta.names=theta.names,
								beta.names=beta.names,
								intercept.names=intercept.names, 
								random.names=random.names,
							    random.lb = -5, 
				                random.ub = 5,
								saem=TRUE)
print(dynm$random.lb)
								
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

model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data, armadillo=TRUE,
                    outfile="VanDerPol.c")
print(model@random.params.inicov)
print(model@random.values.inicov)
#print(model@random)
# to do consist the formula in Line 71 and here
#model@dynamics@theta.formula = list( zeta_i ~ 1*zeta_0 + u1*zeta_1 + u2*zeta_2 + 1*b_zeta,
#                     zeta_i_2 ~ 1*mu1 + 1*b_x1,
#                      zeta_i_3 ~ 1*mu2 + 1*b_x2)

fitted_model <- dynr.cook(model, saem=TRUE, optimization_flag = TRUE, hessian_flag = TRUE, verbose=TRUE, debug_flag=TRUE)
print(fitted_model)

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

