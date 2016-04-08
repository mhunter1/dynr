# 
# dynr model CLASS
# 
setClass("dynrOptOptions",
         representation(maxtime="integer", 
                        maxeval="integer",
                        ftol_rel="numeric",
                        xtol_rel="numeric"
         )
)
setMethod("initialize", "dynrOptOptions",
          function(.Object, opt.list){
            .Object@maxtime=as.integer(opt.list$maxtime)
            .Object@maxeval=as.integer(opt.list$maxeval)
            .Object@ftol_rel=as.numeric(opt.list$ftol_rel)
            .Object@xtol_rel=as.numeric(opt.list$xtol_rel)
            return(.Object)
          }
)
setClass(Class =  "dynrControl",
         representation = representation(
           xstart="vector",
           ub="vector",
           lb="vector",
           options="dynrOptOptions"
         )
)
setMethod("initialize", "dynrControl",
          function(.Object,num_func_param){
            .Object@xstart=rep(0,num_func_param)
            .Object@ub=rep(9999,num_func_param)
            .Object@lb=rep(9999,num_func_param)
            .Object@options=new("dynrOptOptions",opt.list=list(maxtime=30*60, 
                                                               maxeval=5000,
                                                               ftol_rel=as.numeric(1e-8),
                                                               xtol_rel=as.numeric(1e-8)))
            return(.Object)
          }
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
						            noise = "dynrNoise",
						            infile="character",
						            outfile="character",
						            isContinuousTime="logical",
						            verbose="logical",
						            compileLib="logical"
         )
)
setMethod("initialize", "dynrModel",
          function(.Object, meas, noise, initial, dynamics){
            .Object@num_regime="integer",
            .Object@dim_latent_var="integer",
            .Object@opt.control="dynrControl",
            .Object@dynamics =  "dynrDynamics",
            .Object@measurement = meas,
            .Object@initial = initial,
            .Object@regimes= "dynrRegimes",
            .Object@noise = noise
            return(.Object)
          }
)
dynr.prep.Nametochange<-function(dynrModel){
  #take in a dynrModel object
  
  #modify the object slot, including starting values, etc.
}

dynr.cook.Nametochange<-function(dynrModel,data){
  #1. dynr.model convert dynrModel to a model list
  model<-dynr.model(
    dynrModel@num_regime,
    dynrModel@dim_latent_var,
    dynrModel@opt.control@xstart,
    dynrModel@opt.control@ub,
    dynrModel@opt.control@lb,
    list(dynrModel@opt.control@optionsmaxtime, 
         dynrModel@opt.control@optionsmaxeval,
         dynrModel@opt.control@optionsftol_rel,
         dynrModel@opt.control@optionsxtol_rel),
    dynrModel@isContinuousTime,
    dynrModel@infile,
    dynrModel@outfile,
    dynrModel@verbose,
    dynrModel@compileLib
  )
  #Others keep unchanged
  #model$func_addresss = dynr.funcaddress()
}