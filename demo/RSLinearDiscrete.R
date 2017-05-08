#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2016-05-24
# Filename: RSLinearDiscrete.R
# Purpose: Show regime-switching 
# measurement and dynamics in linear discrete time model
#------------------------------------------------------------------------------


#---- (1) Load packages ----
require(dynr)
#---- (2) Read in data ----
# create dynr data object

data(EMGsim)
dd <- dynr.data(EMGsim, id='id', time='time', observed='EMG', covariates='self')


#---- (3) Specify recipes for all model pieces ----

#---- (3a) Measurement ----
recMeas <- prep.measurement(
	values.load=rep(list(matrix(1, 1, 1)), 2),
	values.int=list(matrix(0, 1, 1), matrix(1, 1, 1)),
	params.int=list(matrix('mu_0', 1, 1), matrix('mu_1', 1, 1)),
	values.exo=list(matrix(0, 1, 1), matrix(1, 1, 1)),
	params.exo=list(matrix('beta_0', 1, 1), matrix('beta_1', 1, 1)),
	obs.names = c('EMG'),
	state.names=c('lEMG'),
	exo.names=c("self"))

#---- (3b) Dynamic and measurement noise cov structures----

recNoise <- prep.noise(
	values.latent=matrix(1, 1, 1),
	params.latent=matrix('dynNoise', 1, 1),
	values.observed=matrix(0, 1, 1),
	params.observed=matrix('fixed', 1, 1))

# ---- (3c) Regimes-switching model ----

recReg <- prep.regimes(
	values=matrix(0, 2, 2),
	params=matrix(c('c11', 'c21', 'fixed', 'fixed'), 2, 2))

recReg2 <- prep.regimes(
	values=matrix(0, 2, 2),
	params=matrix(c('dev1', 'base1', 'fixed', 'fixed'), 2, 2),
	deviation=TRUE) # refRow gets set to refCol=2

recReg1 <- prep.regimes(
	values=matrix(0, 2, 2),
	params=matrix(c('base1', 'dev2', 'fixed', 'fixed'), 2, 2),
	deviation=TRUE, refRow=1) #refRow get set to non-default 1

#recReg <- prep.regimes(
#	values=matrix(0, 2, 2),
#	params=matrix(c('p00', 'p10', 'fixed', 'fixed'), 2, 2),
#	deviation=TRUE, refRow=3) ##rightly causes error, refRow out of bounds




#---- (3d) Initial condition specification ----

recIni <- prep.initial(
	values.inistate=matrix(0, 1, 1),
	params.inistate=matrix('fixed', 1, 1),
	values.inicov=matrix(1, 1, 1),
	params.inicov=matrix('fixed', 1, 1),
	values.regimep=c(10, 0),
	params.regimep=c('fixed', 'fixed'))


#---- (3e) Dynamic model ----

recDyn <- prep.matrixDynamics(
	values.dyn=list(matrix(.1, 1, 1), matrix(.8, 1, 1)),
	params.dyn=list(matrix('phi_0', 1, 1), matrix('phi_1', 1, 1)),
	isContinuousTime=FALSE)

#---- (4) Create model and cook it all up ----

rsmod <- dynr.model(dynamics=recDyn, measurement=recMeas, noise=recNoise, initial=recIni, regimes=recReg, data=dd, outfile="RSLinearDiscrete.c")

rsmod$lb['phi_0'] <- -0.01

rsmod1 <- dynr.model(dynamics=recDyn, measurement=recMeas, noise=recNoise, initial=recIni, regimes=recReg1, data=dd, outfile="RSLinearDiscreteDev1.c")

rsmod1$lb['phi_0'] <- -0.01

rsmod2 <- dynr.model(dynamics=recDyn, measurement=recMeas, noise=recNoise, initial=recIni, regimes=recReg2, data=dd, outfile="RSLinearDiscreteDev2.c")

rsmod2$lb['phi_0'] <- -0.01
rsmod2$ub['beta_0'] <- 1.0


# Inspect three versions of the same model
#  The only difference is how the regimes are created
printex(rsmod, ParameterAs=rsmod$param.names, printInit=TRUE,printRS=TRUE,
        outFile="RSLinearDiscrete.tex")
#tools::texi2pdf("RSLinearDiscrete.tex")
#system(paste(getOption("pdfviewer"), "RSLinearDiscrete.pdf"))
printex(rsmod1, ParameterAs=rsmod1$param.names, printInit=TRUE,printRS=TRUE,
        outFile="RSLinearDiscrete1.tex")
#tools::texi2pdf("RSLinearDiscrete1.tex")
#system(paste(getOption("pdfviewer"), "RSLinearDiscrete1.pdf"))
printex(rsmod2, ParameterAs=rsmod2$param.names, printInit=TRUE,printRS=TRUE,
        outFile="RSLinearDiscrete2.tex")
#tools::texi2pdf("RSLinearDiscrete2.tex")
#system(paste(getOption("pdfviewer"), "RSLinearDiscrete2.pdf"))


yum <- dynr.cook(rsmod, debug_flag=TRUE)
yum1 <- dynr.cook(rsmod1)
yum2 <- dynr.cook(rsmod2)



#---- (5) Serve it! ----

summary(yum)

dynr.ggplot(yum, dynrModel = rsmod, style = 1, 
            names.regime=c("Deactivated","Activated"),
            title="Results from RS-AR model", numSubjDemo=1,
            shape.values = c(1),
            text=element_text(size=16))
#ggsave("RSLinearDiscreteggPlot1.pdf")
dynr.ggplot(yum, dynrModel = rsmod, style = 2,
            names.regime=c("Deactivated","Activated"),
            title="Results from RS-AR model", numSubjDemo=1,
            text=element_text(size=16))
#ggsave("RSLinearDiscreteggPlot2.pdf")
plot(yum, dynrModel = rsmod, style = 1, textsize = 5)
plot(yum, dynrModel = rsmod, style = 2, textsize = 5)
#---- End of demo ---- 
#save(rsmod,yum,file="RSLinearDiscrete.RData")


#---- (6) Compare true and estimated params ----

#true parameters
truep <- c(phi0=.3, phi1=.9, beta0=0, beta1=.5, mu0=3, mu1=4, dynnoise=.5^2, p00=.99, p10=.01)
estp <- coef(yum)

r1 <- c(coef(yum)[which(rsmod$param.names=="c11")],0)
(r1 <- exp(r1)/sum(exp(r1))) #first row of transition probability matrix
r2 <- c(coef(yum)[which(rsmod$param.names=="c21")],0)
(r2 <- exp(r2)/sum(exp(r2))) #second row of transition probability matrix

estp[8:9] <- c(r1[1], r2[1])

plot(estp, truep)
abline(a=0, b=1)

rms <- function(x, y){sqrt(mean((x-y)^2))}

testthat::expect_true(rms(estp, truep) < 0.07)

# Check that all true parameters are within the confidence intervals of the estimated params
# other than the regime-switching params
withinIntervals <- yum@conf.intervals[1:7,1] < truep[1:7] & truep[1:7] < yum@conf.intervals[1:7,2]
testthat::expect_true(all(withinIntervals))


#---- (7) Compare true and estimated states ----

#---- (7a) latent states ----

#pdf('plotFilterSmoothEMG.pdf')

smoothCor <- cor(yum@eta_smooth_final[1,], EMGsim$truestate)
smoothRMS <- rms(yum@eta_smooth_final[1,], EMGsim$truestate)

plot(yum@eta_smooth_final, EMGsim$truestate, xlab='Smoothed Latent State', ylab='True Latent State', main='Smoothed Estimates')
text(x=1, y=-2, labels=paste0('r = ', round(smoothCor,3)), adj=c(1,1))
text(x=1, y=-2.5, labels=paste0('RMSE = ', round(smoothRMS,3)), adj=c(1,1))

estRegime <- apply(yum@pr_t_given_T, 2, which.max) - 1


fstate2 <- yum$eta_filtered[1,]

updateCor <- cor(fstate2, EMGsim$truestate)
updateRMS <- rms(fstate2, EMGsim$truestate)

plot(fstate2, EMGsim$truestate, xlab='Updated Latent State', ylab='True Latent State', main='Updated Estimates')
text(x=1, y=-2, labels=paste0('r = ', round(updateCor,3)), adj=c(1,1))
text(x=1, y=-2.5, labels=paste0('RMSE = ', round(updateRMS,3)), adj=c(1,1))


#---- (7b) latent regimes ----


(rtab <- table(trueRegime=EMGsim$trueregime, estRegime=estRegime))

plot(x=jitter(apply(yum@pr_t_given_T, 2, which.max) - 1), y=jitter(EMGsim$trueregime), xlab='Estimated Smoothed Regime', ylab='True Regime', xaxt='n', yaxt='n', main='Smoothed Regimes')
axis(side=1, at=c(0,1))
axis(side=2, at=c(0,1))
text(x=c(0,0,1,1), y=c(0+.25,1-.25,0+.25,1-.25), labels=c(rtab))

#dev.off()


#---- (8) Compare deviation and non-deviation estimates ----

# Check that regime-switching probabilities match
testthat::expect_equal(coef(yum)['c11'], coef(yum1)['base1'], tolerance=0.001, check.names=FALSE)
testthat::expect_equal(coef(yum)['c21'], coef(yum2)['base1'], tolerance=0.001, check.names=FALSE)
testthat::expect_equal(coef(yum)['c21'], sum(coef(yum1)[c('base1', 'dev2')]), tolerance=0.001, check.names=FALSE)
testthat::expect_equal(coef(yum)['c11'], sum(coef(yum2)[c('base1', 'dev1')]), tolerance=0.001, check.names=FALSE)

# Check that other parameters match
estNonRegimeParams <- cbind(coef(yum), coef(yum1), coef(yum2))[1:7,]
for(i in 1:nrow(estNonRegimeParams)){
	testthat::expect_equal(estNonRegimeParams[i, 1], estNonRegimeParams[i, 2], tolerance=0.001, check.names=FALSE)
	testthat::expect_equal(estNonRegimeParams[i, 2], estNonRegimeParams[i, 3], tolerance=0.001, check.names=FALSE)
}

