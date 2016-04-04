# 
# dynr model CLASS
# 

setClass("dynrOptOptions",
         representation(maxtime="numeric", 
                        maxeval="integer",
                        ftol_rel="numeric",
                        xtol_rel="numeric"
         )
)
setClass("dynrModel",
         representation(response="list", # response models
                        prior="ANY", # the prior model (multinomial)
                        dens="array", # response densities (B)
                        init="array", # usually called pi 
                        nstates="numeric",
                        nresp="numeric",
                        ntimes="numeric",
                        npars="numeric", # number of parameters
                        num_regime="integer",
                        dim_latent_var="integer",
                        xstart="vector",
                        ub="vector",
                        lb="vector",
                        options="list",
                        isContinuousTime="logical",
                        infile="character",
                        outfile="character",
                        verbose="logical",
                        compileLib="logical"
         )
)

