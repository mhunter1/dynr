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
setClass(Class =  "dynrControl",
         representation = representation(
           xstart="vector",
           ub="vector",
           lb="vector",
           options="dynrOptOptions",
           isContinuousTime="logical",
           infile="character",
           outfile="character",
           verbose="logical",
           compileLib="logical",
         ),
         contains="dynrOptOptions"
)
setClass(Class =  "dynrModel",
         representation = representation(
                     	  num_regime="integer",
                        dim_latent_var="integer",
                        opt.control="dynrControl",
                        dynamics =  "dynrDynamics",
			                  measurement = "dynrMeasurement",
			                  initial = "dynrInitial",
						            regimes= "dynrRegimes",
						            noise = "dynrNoise"
         ),
         contains="dynrControl"
)

PrepFuncNametochange<-function(dynrModel){
  
}
