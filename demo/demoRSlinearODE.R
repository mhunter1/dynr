#rm(list=ls(all=TRUE))
require(dynr)
options(scipen=999)
nT = 500; n = 10; batch = 1
thedata = read.table(paste0("./data/New2CovT",nT,"n",n,"batch",batch,"ODEsimData.txt"))
thedata$V6 <- as.numeric(thedata$V6)
data <- dynr.data(thedata, id="V1", time="V2",observed=paste0('V', 3:4), 
                  covariates=paste0('V', 5:6))
#options=list(maxtime=1*60, maxeval=20) #may add this to the "model" line below

model <- dynr.model(
              num_regime=2,
              dim_latent_var=2,
              xstart=c(rep(log(.1), 4), 95.0, log(10.0), log(10.0), -3.0, 9.0, -1.5, 0.5, -1,-.3),
              ub=rep(9999,13),lb=rep(9999,13),
              #ub=c(rep(1.5, 4), 0, 5, 5, rep(10, 6)),
              #lb=c(rep(-10, 4), 200, -10, -10, rep(-10, 6)),
              options=list(maxtime=60*60, 
                           maxeval=500,
                           ftol_rel=as.numeric(1e-8),
                           xtol_rel=as.numeric(1e-5)),
              isContinuousTime=TRUE,
              infile="./demo/RSODEmodel.c", 
              outfile="./demo/RSODEmodel2", 
              verbose=TRUE,
              compileLib=FALSE
)


#func_noise_cov_txt="void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov){
#size_t i;
#for (i=0;i<eta_noise_cov->size1;i++){
#gsl_matrix_set(eta_noise_cov,i,i,-20);
#}
#gsl_matrix_set(y_noise_cov,0,0, param[5]);
#gsl_matrix_set(y_noise_cov,1,1, param[6]);
#
#}
#" 

#func_address=dynr.funcaddress(file="./demo/RSODEmodel.c",verbose=FALSE,model=model)
tfun <- function(x){c(exp(x[1:4]),x[5],exp(x[6:7]), x[8:13])}
res <- dynr.cook(model, data,tfun)


#True values should be
#c(log(.2), log(.1), log(.3), log(.2),  100, log(9.0), log(9.0), -4, 8.5, -1, 1,-2, -1)


summary(res)
#save.image(file="./demo/RSLinearODETestBrief.RData")
#plot(res, data=data, graphingPar=list(cex.main=1, cex.axis=1, cex.lab=1.2), numSubjDemo=2)


#dynr.ggplot(res, data.dynr=data, states=c(1,2), title="Smoothed State Values", numSubjDemo=2)

