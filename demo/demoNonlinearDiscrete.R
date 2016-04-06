#rm(list=ls(all=TRUE))
require(dynr)

thedata = read.table(paste0("./data/HSDParcels.txt"),na.strings = "-999")
colnames(thedata)<-c("id", "Day", "Time", "pos1", "pos2", "pos3", "neg1", "neg2", "neg3", "isSocCompany",
"posEvents", "negEvents", "pss", "extrav", "agree", "consc", "emoStab", "open")
data <- dynr.data(thedata, id="id", time="Time",observed=colnames(thedata)[c(6,5,4,9,8,7)])

formula=list(x1~param[4]*x1+param[6]*(exp(abs(x2))/(1+exp(abs(x2))))*x2,
             x2~param[5]*x2+param[7]*(exp(abs(x1))/(1+exp(abs(x1))))*x1)
jacob=list(x1~x1~param[4],
           x1~x2~param[6]*(exp(abs(x2))/(exp(abs(x2))+1)+x2*sign(x2)*exp(abs(x2))/pow(1+exp(abs(x2)),2)),
           x2~x2~param[5],
           x2~x1~param[7]*(exp(abs(x1))/(exp(abs(x1))+1)+x1*sign(x1)*exp(abs(x1))/pow(1+exp(abs(x1)),2)))
dynm<-dynr.nonlindynamics(formula,jacob,isContinuosTime=FALSE)
writeLines(dynm)


model <- dynr.model(
              num_regime=1,
              dim_latent_var=2,
              xstart=c(1,1,1,1,.5,.5,-.3,-.3,
                       rep(log(.5),6),
                       rep(log(.3),2)),
              ub=rep(9999,16),
              lb=rep(9999,16),
              options=list(maxtime=1*60, 
                           maxeval=1,
                           ftol_rel=as.numeric(1e-8),
                           xtol_rel=as.numeric(1e-5)),
              isContinuousTime=FALSE,
              infile="./demo/NonlinearDiscrete.c", 
              outfile="./demo/NonlinearDiscrete2", 
              verbose=TRUE,
              compileLib=TRUE
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
tfun <- function(x){c(x[1:8],exp(x[9:16]))}
res <- dynr.cook(model, data,tfun)


#True values should be
#c(log(.2), log(.1), log(.3), log(.2),  100, log(9.0), log(9.0), -4, 8.5, -1, 1,-2, -1)


summary(res)
#save.image(file="./demo/RSLinearODETestBrief.RData")
#plot(res, data=data, graphingPar=list(cex.main=1, cex.axis=1, cex.lab=1.2), numSubjDemo=2)


#dynr.ggplot(res, data.dynr=data, states=c(1,2), title="Smoothed State Values", numSubjDemo=2)

