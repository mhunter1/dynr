if(!require(numDeriv)){
	install.packages('numDeriv')
}
if(!require(ggplot2)){
    install.packages('ggplot2')
}
if(!require(reshape2)){
    install.packages('reshape2')
}
if(!require(plyr)){
    install.packages('plyr')
}
require(dynr)
require(numDeriv)



data(dataPANAsim)
data <- dynr.data(dataPANAsim, id="V1", time="V2",observed=paste0('V', 3:4), covariates=paste0('V', 5))

# the model
model <- dynr.model(num_regime=1,
                    dim_latent_var=4,
                    xstart=c(log(1),log(2),0,0,-10,-10),
                    ub=rep(5, 6),
                    lb=c(-5,-5,-5,-5,-15, -15)
                  #, options=list(maxtime=30, maxeval=1)
                    )

# initial values and bounds


tfun <- function(x){c(exp(x[1]), exp(x[2]), x[3:6])}
x <- dynr.run(model, data, tfun)
summary(x)
plot(x, data=data, graphingPar=list(cex.main=1, cex.axis=1, cex.lab=1.2), numSubjDemo=2)


dynr.ggplot(x, data.dynr=data, states=c(1,3), numSubjDemo=2)

