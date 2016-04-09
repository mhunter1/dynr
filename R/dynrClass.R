# 
# dynr model CLASS
# 
setClass(Class =  "dynrModel",
         representation = representation(
           dynamics =  "dynrDynamics",
           measurement = "dynrMeasurement",
           noise = "dynrNoise",
           initial = "dynrInitial",
           regimes= "dynrRegimes",
           num_regime="integer",
           dim_latent_var="integer",
           infile="character",
           outfile="character",
           isContinuousTime="logical",
           verbose="logical",
           compileLib="logical",
           xstart="vector",
           ub="vector",
           lb="vector",
           options="list"
         ),
         prototype(
           num_regime=as.integer(1),
           isContinuousTime=TRUE,
           verbose=TRUE,
           compileLib=TRUE,
           options=default.model.options
         )
)

setMethod("initialize", "dynrModel",
          function(.Object, x){
            for(i in names(x)){
              slot(.Object, name=i, check = TRUE) <- x[[i]]
            }
            return(.Object)
          }
)

setMethod("$", "dynrModel",
          function(x, name){slot(x, name)}
)
