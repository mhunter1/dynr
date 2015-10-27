
require(dynr)

thedata = read.table('./data/New2CovModel1T1000n20batch1ODEsimData.txt')
thedata$V6 <- as.numeric(thedata$V6)
data <- dynr.data(thedata, id="V1", time="V2",observed=paste0('V', 3:4), 
                  covariates=paste0('V', 5:6))
#options=list(maxtime=1*60, maxeval=20) #may add this to the "model" line below

model <- dynr.model(
              num_regime=2,
              dim_latent_var=2,
              xstart=c(rep(log(.1), 4), log(10.0), log(10.0), -3.0, 9.0, -1.5, -0.5, 95.0,-.3,-.3),
              ub=c(rep(10, 6), rep(10, 4), 200, 10, 10),
              lb=c(rep(-10, 6), rep(-10, 4), 0, -10, -10),
	     options=list(maxtime=60*60, maxeval=1000)
)

tfun <- function(x){c(exp(x[1:6]), x[7:13])}
res <- dynr.run(model, data, tfun)
#True values should be
#c(log(.2), log(.1), log(.3), log(.2), log(9.0), log(9.0), -4, 8.5, -1, 1, 100, -2, -1)

summary(res)
save.image(file="./demo/RSLinearODE.RData")
plot(res, data=data, graphingPar=list(cex.main=1, cex.axis=1, cex.lab=1.2), numSubjDemo=2)


dynr.ggplot(res, data.dynr=data, states=c(1,2), title="Smoothed State Values", numSubjDemo=2)

