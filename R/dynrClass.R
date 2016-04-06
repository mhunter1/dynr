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
setClass(Class =  "dynrModel",
         representation = representation(
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
                        compileLib="logical",
			            dynamics =  "dynrDynamics",
			            measurement = "dynrMeasurement",
			            initial = "dynrInitial",
						regimes= "dynrRegimes",
						Class = "dynrNoise"
         )
)
