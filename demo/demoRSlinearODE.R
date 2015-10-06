require(dynr)

data(New2CovModel1T1000n20batch1ODEsimData)
thedata = New2CovModel1T1000n20batch1ODEsimData
thedata$V6 <- as.numeric(thedata$V6)
data <- dynr.data(thedata, id="V1", time="V2",observed=paste0('V', 3:4), 
                  covariates=paste0('V', 5:6))
model <- dynr.model(
              num_regime=2,
              dim_latent_var=2,
              xstart=c(rep(log(.1), 4), log(10.0), log(10.0), -3.0, 9.0, -1.5, -0.5, 95.0,-.3,-.3),
              ub=c(rep(10, 6), rep(20, 4), 1000, 20, 20),
              lb=c(rep(-10, 6), rep(-20, 4), 0, -20, -20),
              options=list(maxtime=30*60, maxeval=1)
)


res <- dynr.run(model, data)
summary(res)
save.image(file="RSLinearODE.RData")
plot(res, data=data, graphingPar=list(cex.main=1, cex.axis=1, cex.lab=1.2), numSubjDemo=2)


dynr.ggplot(res, data.dynr=data, states=c(1,2),
            #mancolorPalette=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
            #manfillPalette=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
            title="Smoothed State Values", numSubjDemo=2)

