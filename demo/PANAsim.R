Tutorial: Using *dynr* to fit a regime-switching linear ODE model
--------------------------------------------------------------------
  This example illustrates how to use functions in the dynr package
to fit a regime-switching linear ODE model of the form:
  
  $\frac{dx_1(t)}{dt} = -r_1*x_1(t) + a_{12}*[x_2(t)-x_1(t)]$
  
  $\frac{dx_2(t)}{dt} = -r_2*(x_2(t)-x_{20}) + a_{21}*[x_2(t)-x_1(t)]$ 
  
  - if $S_t=1$: $r_1$ and $r_2$ are freely estimated; $a_{12}$ $=$ $a_{21}$ $=$ 0;
- if $S_t=2$: $r_1$ $=$ $r_2$ $=$ 0; $a_{12}$ and $a_{21}$ $=$ are freely estimated.

```{r test-child}
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

tfun <- function(x){c(exp(x[1]), exp(x[2]), x[3:6])}
starttime <- proc.time()
x <- dynr.run(model, data, tfun)
(time <- proc.time()-starttime)
x
data.frame(x$transformed.parameters, x$standard.errors, x$CI)

```
