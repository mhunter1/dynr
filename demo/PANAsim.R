
require(dynr)

data(dataPANAsim)
data <- dynr.data(dataPANAsim, id="V1", time="V2",observed=paste0('V', 3:4), covariates=paste0('V', 5))

# the model
model <- list(num_sbj=217,
            dim_latent_var=4,
            dim_obs_var=2,
            dim_co_variate=1, 
            num_regime=1,
            xstart=c(log(1),log(2),0,0,-10,-10),
            num_func_param=6,
            ub=c(5, 5, 5, 5, 5, 5),
            lb=c(-5,-5,-5,-5,-15, -15)
)
# initial values and bounds

starttime <- proc.time()
x <- dynr.run(model, data)
(time <- proc.time()-starttime)
x
