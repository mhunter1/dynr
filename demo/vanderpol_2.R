#------------------------------------------------------------------------------
# Author: Hui-Ju Hung
# Date: 2018-05-15
# Filename: vanderpol_2.R
# Purpose: Fit van der Pol oscillator model with SAEM
#------------------------------------------------------------------------------


library('dynr')
#vdPdata <- read.table('directory')
state.names = c('x1', 'x2')
#[todo] beta=> fix
beta.names = c('zeta0', 'zeta1', 'zeta2', 'mu1', 'mu2')
covariate.names = c('u1', 'u2')
theta.names = c('zeta_i', 'zeta_i_2', 'zeta_i_3')
#[todo] Z*b's b
random.names = c('b1')
DV.names = c('y1', 'y2', 'y3')


NxState = length(state.names)
Nbeta = length(beta.names)
Nx = NxState + Nbeta
Ntheta = length(theta.names)

#random data
#N = 100
#T = 300
#vdpData <- data.frame(id=rep(1:N,each=T), time=rep(0:299,N), 
#                      y=rnorm(100*300),
#                      u1 = rnorm(100*300), u2 = rnorm(100*300))
#colnames(vdpData) <- c("id","time","y","u1", "u2") # try

nPeople = 200
nTimes = 300
table = read.csv("TrueInitY1.txt", header=FALSE)
colnames(table) <- union(union(c('batch', 'kk', 'trueInit', 'time'), DV.names), covariate.names)
vdpData <- data.frame(id=rep(1:N, each=nTimes), time=table['time'],
                      y1=table['y1'], y2=table['y2'], y3=table['y3'],
                      u1 = table['u1'], u2 = table['u2'])
data <- dynr.data(vdpData, id="id", time="time", 
                  observed=DV.names,
                  covariates=covariate.names)

#????
meas <- prep.measurement(
    values.load = matrix(c(1, 0), 3, 2),
    params.load = matrix(c('fixed'), 3, 2),
    obs.names = DV.names,
    state.names = state.names)


#TODO adjust initial condition
# Initial conditions on the latent state and covariance
# initial <- prep.initial(
	# values.inistate=c(3, 1),
	# params.inistate=c("fixed", "fixed"),
	# values.inicov=diag(c(0.01,0.01)), 
	# params.inicov=diag("fixed",2)
# )
initial <- prep.initial(
    values.inistate=c(3, 1),
    params.inistate=c("mu1", "mu2"),
    values.inicov=matrix(c(.5,.2,
                           .2,.6),ncol=2,byrow=T), 
    params.inicov=matrix(c('v10','c120',
                           'c120','v20'),ncol=2,byrow=T)
)
if (length(unlist(initial$params.inistate[!initial$params.inistate==
                                   "fixed"]))>0)
initial$params.inistate

# TODO adjust noise
#measurement and dynamics covariances
mdcov <- prep.noise(
	values.latent=diag(0, 2),
	params.latent=diag(c("fixed","fixed"), 2),
	values.observed=diag(rep(0.3,2)),
	params.observed=diag(c("var_1","var_2"),2)
)
		 

formula=
    list(x1 ~ x2,
         x2 ~ -61.68503 * x1 + zeta_i * (1 - x1^2) * x2,
         zeta0 ~0,
         zeta1 ~0,
         zeta2 ~0,
         mu1 ~0,
         mu2 ~0
    )
theta.formula  = list (zeta_i ~ zeta0  + u1 * zeta1 + u2 * zeta2)


#beta.names = c('param[0]', 'param[1]', 'param[2]', 'mu1', 'mu2')
dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(zeta0=-1,
                                      zeta1=.5,
                                      zeta2=.2),
								isContinuousTime=FALSE,
								state.names=state.names,
								theta.formula=theta.formula,
								theta.names=theta.names,
								beta.names=beta.names,
								saem=TRUE)



meas@state.names = c('x1', 'x2', 'zeta0', 'zeta1', 'zeta2', 'mu1', 'mu2')
data$covariate.names = covariate.names

model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data,armadillo=TRUE,
                    outfile="VanDerPol.c")


#covariate.names=covariate.names,beta.names=beta.names,

