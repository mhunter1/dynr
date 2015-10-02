require(dynr)

thedata = read.table('~/Dropbox/SyMiin_Lu/dynr/data/New2CovModel1T1000n20batch1ODEsimData.txt')
thedata$V6<-as.numeric(thedata$V6)
data <- dynr.data(thedata, id="V1", time="V2",observed=paste0('V', 3:4), 
                  covariates=paste0('V', 5:6))
model <- list(num_sbj=20,
              dim_latent_var=2,
              dim_obs_var=2,
              dim_co_variate=2, 
              num_regime=2,
              xstart=c(log(.1), log(.1), log(.1), log(.1), log(10.0), log(10.0), -3.0, 9.0, -1.5, -0.5, 95.0,-.3,-.3),
              num_func_param=13,
              ub=c(10.0,10.0,10.0,10.0,10.0,10.0,20.0,20.0,20.0,20.0,1000.0,20.0,20.0),
              lb=c(-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-20.0,-20.0,-20.0,-20.0,0.0,-20.0,-20.0)
)


x <- dynr.run(model, data)
summary(x)
plot(x, data=data, graphingPar=list(cex.main=1, cex.axis=1, cex.lab=1.2), numSubjDemo=2)

require(ggplot2)
require(reshape2)
require(plyr)
dynr.ggplot(x,data.dynr=data,states=c(1,2),names.state=paste0("state",states),title="Smoothed State Values",numSubjDemo=2)

