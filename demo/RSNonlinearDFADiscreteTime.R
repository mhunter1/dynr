#rm(list=ls(all=TRUE))
# Example for fitting a nonlinear dynamic factor analysis model with
# nonlinear vector autoregressive (VAR) relation at the factor level

require(dynr)
setwd("~/Dropbox/Symiin_Lu/dynr")

#Reading in simulated data
thedata <- read.table(paste0("./data/NonlinearVARsimT300n10.txt"),na.strings = "NaN",sep=",")
thedata2 <- t(thedata)
n = 10; T = 300
thedata2 <- cbind(rep(1:n,each=T),rep(1:T,n),thedata2)

colnames(thedata2)<-c("id", "Time", "y1", "y2", "y3", "y4", "y5", "y6")
data <- dynr.data(thedata2, id="id", time="Time",observed=colnames(thedata2)[c(3:8)])

formula=list(x1~param[6]*x1+param[8]*(exp(abs(x2))/(1+exp(abs(x2))))*x2,
             x2~param[7]*x2+param[9]*(exp(abs(x1))/(1+exp(abs(x1))))*x1)
jacob=list(x1~x1~param[6],
           x1~x2~param[8]*(exp(abs(x2))/(exp(abs(x2))+1)+x2*sign(x2)*exp(abs(x2))/pow(1+exp(abs(x2)),2)),
           x2~x2~param[7],
           x2~x1~param[9]*(exp(abs(x1))/(exp(abs(x1))+1)+x1*sign(x1)*exp(abs(x1))/pow(1+exp(abs(x1)),2)))
dynm<-dynr.nonlindynamics(formula,jacob,isContinuosTime=FALSE)
writeLines(dynm)

pstart <- log(.9/(1-.9))

#True values = 
#c(1.2, 1.2, 1.1, .95, .98, .85, .2, .25, -.6, -.8,
# .28, .10, .12, .13, .12, .11,
# .35, .3)

truepar <- c(1.2, 1.2, 1.1, .95, 
             log(.98/(1-.98)), log(.85/(1-.85)), 
             .2, .25, -.6, -.8,
            log(c(.28, .10, .12, .13, .12, .11)),
            log(c(.35, .3)))

model <- dynr.model(
  num_regime=2,
  dim_latent_var=2,
  xstart=c(1,1,1,1,pstart,pstart,.3,.4,-.5,-.5,
           rep(log(.1),6),
           rep(log(.3),2)),
  #ub=rep(9999,18),
  #lb=rep(9999,18),
  ub=c(rep(3,4),5,5,1.5,1.5,1.5,1.5,rep(log(5),6),rep(log(5),2)),
  lb=c(rep(0,4),-5,-5,-1.5,-1.5,-1.5,-1.5,rep(log(.001),6),rep(log(.001),2)),
  options=list(maxtime=60*60, maxeval=1000,ftol_rel=as.numeric(1e-10),
               xtol_rel=as.numeric(1e-7)),
  isContinuousTime=FALSE,
  infile="./demo/RSNonlinearDiscrete.c", 
  outfile="./demo/RSNonlinearDiscrete2", 
  verbose=TRUE,
  compileLib=TRUE
)



tfun <- function(x){c(x[1:4],
                      exp(x[5])/(1+exp(x[5])), exp(x[6])/(1+exp(x[6])),
                      x[7:10],
                      exp(x[11:18]))}
res <- dynr.run(model, data,tfun)


#True values should be
#c(log(.2), log(.1), log(.3), log(.2),  100, log(9.0), log(9.0), -4, 8.5, -1, 1,-2, -1)


summary(res)
#save.image(file="./demo/RSLinearODETestBrief.RData")
#plot(res, data=data, graphingPar=list(cex.main=1, cex.axis=1, cex.lab=1.2), numSubjDemo=2)


#dynr.ggplot(res, data.dynr=data, states=c(1,2), title="Smoothed State Values", numSubjDemo=2)

