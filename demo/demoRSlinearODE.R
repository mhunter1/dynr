require(dynr)

thedata = read.table('~/Dropbox/Brekfis/dynr/data/New2CovModel1T1000n20batch1ODEsimData.txt')
thedata$V6 <- as.numeric(thedata$V6)
data <- dynr.data(thedata, id="V1", time="V2",observed=paste0('V', 3:4), 
                  covariates=paste0('V', 5:6))
model <- dynr.model(
              num_regime=2,
              dim_latent_var=2,
              xstart=c(rep(log(.1), 4), log(10.0), log(10.0), -3.0, 9.0, -1.5, -0.5, 95.0,-.3,-.3),
              ub=c(rep(10, 6), rep(20, 4), 1000, 20, 20),
              lb=c(rep(-10, 6), rep(-20, 4), 0, -20, -20),
              options=list(maxtime=1*60, maxeval=20)
)

func_noise_cov_txt="void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov){
size_t i;
for (i=0;i<eta_noise_cov->size1;i++){
gsl_matrix_set(eta_noise_cov,i,i,-20);
}
gsl_matrix_set(y_noise_cov,0,0, param[4]);
gsl_matrix_set(y_noise_cov,1,1, param[5]);

}
" 
func_address=dynr.funcaddresses(file="~/Dropbox/Brekfis/dynr/RSODEmodel.c",verbose=FALSE)
tfun <- function(x){c(exp(x[1:4]), x[5:13])}
res <- dynr.run(model, data,func_address,tfun)
summary(res)


save.image(file="~/Dropbox/Symiin_Lu/dynr/RSLinearODE.RData")
plot(res, data=data, graphingPar=list(cex.main=1, cex.axis=1, cex.lab=1.2), numSubjDemo=2)
dynr.ggplot(res, data.dynr=data, states=c(1,2), title="Smoothed State Values", numSubjDemo=2)

