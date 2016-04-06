
require(dynr)

nT = 500; n = 10; batch = 1
thedata = read.table(paste0("./data/New2CovT", nT,"n", n, "batch", batch, "ODEsimData.txt"))
thedata$V6 <- as.numeric(thedata$V6)
data <- dynr.data(thedata, id="V1", time="V2",observed=paste0('V', 3:4), 
                  covariates=paste0('V', 5:6))


model <- dynr.model(
              num_regime=2,
              dim_latent_var=2,
              xstart=c(rep(log(.1), 4), 95.0, log(10.0), log(10.0), -3.0, 9.0, -1.5, 0.5, -1,-.3),
              ub=rep(9999,13),lb=rep(9999,13),
              options=list(maxtime=60*60, 
                           maxeval=500,
                           ftol_rel=as.numeric(1e-8),
                           xtol_rel=as.numeric(1e-5)),
              isContinuousTime=TRUE,
              infile="./RSlinearODE.c", 
              outfile="./RSODEmodel2", 
              verbose=TRUE,
              compileLib=FALSE
)




tfun <- function(x){c(exp(x[1:4]),x[5],exp(x[6:7]), x[8:13])}
res <- dynr.cook(model, data, tfun)


#True values should be
trueValues <- tfun(c(log(.2), log(.1), log(.3), log(.2),  100, log(9.0), log(9.0), -4, 8.5, -1, 1,-2, -1))


summary(res)
#save.image(file="./demo/RSLinearODETestBrief.RData")

#testthat::expect_equal(coef(res), trueValues, tolerance=.01)

withinIntervals <- res@conf.intervals[,1] < trueValues & trueValues < res@conf.intervals[,2]
testthat::expect_true(all(withinIntervals))

#plot(res, data=data, graphingPar=list(cex.main=1, cex.axis=1, cex.lab=1.2), numSubjDemo=2)


#dynr.ggplot(res, data.dynr=data, states=c(1,2), title="Smoothed State Values", numSubjDemo=2)

