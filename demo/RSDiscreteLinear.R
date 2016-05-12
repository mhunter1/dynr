#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Author: Michael D. Hunter
# Date: 2016-04-22
# Filename: RSDiscreteLinear.R
# Purpose: Show regime-switching 
# measurement and dynamics in linear discrete time model
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#---- (1) Load packages ----
require(dynr)
#---- (2) Read in data ----
# create dynr data object

data(EMG)
ds <- EMG
ds$ID <- rep(1, nrow(EMG))
ds$t <- 1:nrow(EMG)
dd <- dynr.data(ds, id='ID', time='t', observed='EMG', covariates='self')


#---- (3) Specify recipes for all model pieces ----

#---- (3a) Measurement ----
recMeas <- prep.measurement(
	values.load=rep(list(matrix(1, 1, 1)), 2),
	values.int=list(matrix(0, 1, 1), matrix(1, 1, 1)),
	params.int=list(matrix('mu0', 1, 1), matrix('mu1', 1, 1)),
	values.exo=list(matrix(0, 1, 1), matrix(1, 1, 1)),
	params.exo=list(matrix('beta0', 1, 1), matrix('beta1', 1, 1)),
	obs.names = c('EMG'),
	state.names=c('lEMG'),
	exo.names=c("self"))

p = printex(recMeas)

#---- (3b) Dynamic and measurement noise cov structures----

recNoise <- prep.noise(
	values.latent=matrix(1, 1, 1),
	params.latent=matrix('dynNoise', 1, 1),
	values.observed=matrix(0, 1, 1),
	params.observed=matrix('fixed', 1, 1))

# ---- (3c) Regimes-switching model ----

recReg <- prep.regimes(
	values=matrix(0, 2, 2),
	params=matrix(c('p00', 'p10', 'fixed', 'fixed'), 2, 2))

#---- (3d) Initial condition specification ----

recIni <- prep.initial(
	values.inistate=matrix(0, 1, 1),
	params.inistate=matrix('fixed', 1, 1),
	values.inicov=matrix(1, 1, 1),
	params.inicov=matrix('fixed', 1, 1),
	values.regimep=c(1, 0),
	params.regimep=c('fixed', 'fixed'))


#---- (3e) Dynamic model ----

recDyn <- prep.matrixDynamics(
	values.dyn=list(matrix(.1, 1, 1), matrix(.8, 1, 1)),
	params.dyn=list(matrix('phi0', 1, 1), matrix('phi1', 1, 1)),
	isContinuousTime=FALSE)

#---- (4) Create model and cook it all up ----

rsmod <- dynr.model(dynamics=recDyn, measurement=recMeas, noise=recNoise, initial=recIni, regimes=recReg, outfile="cooked")

printex(rsmod,show=FALSE)

yum <- dynr.cook(rsmod, dd, debug_flag=TRUE)

#---- (5) Serve it! ----

summary(yum)
p1 = dynr.ggplot(yum, data.dynr=dd, states=1, 
                 names.regime=c("Deactivated","Activated"),
                 names.state=c("EMG"),
                 title="Results from RS-AR model", numSubjDemo=1,
                 shape.values = c(1),
                 text=element_text(size=16))

print(p1)
#true parameters
truep <- c(phi0=.3, phi1=.9, beta0=0, beta1=.5, mu0=3, mu1=4, measnoise=0, dynnoise=.5^2, p00=.99, p10=.01)

rsmod@xstart <- coef(yum)
r1 <- c(rsmod$xstart[which(rsmod$param.names=="p00")],0)
exp(r1)/sum(exp(r1)) #first row of transition probability matrix
r2 <- c(rsmod$xstart[which(rsmod$param.names=="p10")],0)
exp(r2)/sum(exp(r2)) #second row of transition probability matrix

#---- End of demo ---- 
