#------------------------------------------------------------------------------
# Author: Hui-Ju Hung
# Date: 2018-05-15
# Filename: vanderpol_dynr.R
# Purpose: Fit van der Pol oscillator model (cook workable)
# Note: Workable for developer dynr on master branch 0.1.14-17
#------------------------------------------------------------------------------

library('dynr')

nPeople = 200
nTimes = 300
vdpData <- read.csv("TrueInit_Y14.txt", header=FALSE)
colnames(vdpData) <- c('batch', 'kk', 'trueInit', 'time', 'y1','y2','y3', 'co1', 'co2')
vdpData$id <- rep(1:nPeople, each=nTimes)
data <- dynr.data(vdpData, id="id", time="time",
                  observed=c('y1','y2','y3'),
                  covariates=c("co1","co2"))


mdcov <- prep.noise(
    values.latent=diag(0, 2),
    params.latent=diag(c("fixed","fixed"), 2),
    values.observed=diag(rep(0.3,3)),
    params.observed=diag(c("var_1","var_2","var_3"),3)
)

meas <- prep.measurement(
    values.load = matrix(c(1, 0), 3, 2),
    params.load = matrix(c('fixed'), 3, 2),
    obs.names = c('y1', 'y2', 'y3'),
    state.names = c('x1', 'x2'))


initial <- prep.initial(
    values.inistate=c(3, 1),
    params.inistate=c("fixed", "fixed"),
    values.inicov=matrix(c(.5,.2,
                           .2,.6),ncol=2,byrow=T), 
    params.inicov=diag("fixed",2)
)
if (length(unlist(initial$params.inistate[!initial$params.inistate==
                                          "fixed"]))>0)
    initial$params.inistate





formula=
    list(x1 ~ x2,
         x2 ~ -61.68503 * x1 + (zeta0  + co1 * zeta1 + co2 * zeta2) * (1 - x1^2) * x2
    )


dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(zeta0=-1,
                                      zeta1=.5,
                                      zeta2=.2, mu1= 3, mu2 = 1),
                           isContinuousTime=TRUE)

model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data, armadillo=FALSE,
                    outfile="VanDerPol_.c")
model$xstart<-c(3,0.5,0.3,0,0,0.5,0.5,0.5)
fitted_model <- dynr.cook(model, optimization_flag = TRUE, hessian_flag = TRUE, verbose=TRUE, debug_flag=TRUE)
