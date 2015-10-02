# TODO
# create dynrModel objects
# They should output a list like this
#model <- list(num_sbj=20,
#              dim_latent_var=2,
#              dim_obs_var=2,
#              dim_co_variate=2, 
#              num_regime=2,
#              xstart=c(rep(log(.1), 4), log(10.0), log(10.0), -3.0, 9.0, -1.5, -0.5, 95.0,-.3,-.3),
#              num_func_param=13,
#              ub=c(rep(10, 6), rep(20, 4), 1000, 20, 20),
#              lb=c(rep(-10, 6), rep(-20, 4), 0, -20, -20)
#)

# Really all the model needs right now is
# num_regime, xstart, ub, lb
# Everything else can be gathered from the data.
# Importantly, the thing handed to the backend must
# remain a list exactly like the above.

