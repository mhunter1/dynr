# 
# dynr model CLASS
# 

setClass(Class =  "dynrModel",
         representation = representation(
           num_regime="integer",
           dim_latent_var="integer",
           dynamics =  "dynrDynamics",
           measurement = "dynrMeasurement",
           initial = "dynrInitial",
           regimes= "dynrRegimes",
           noise = "dynrNoise",
           infile="character",
           outfile="character",
           isContinuousTime="logical",
           verbose="logical",
           compileLib="logical",
           xstart="vector",
           ub="vector",
           lb="vector",
           options="list"
         )
)
setMethod("initialize", "dynrModel",
          function(.Object, meas, noise, initial, dynamics){
            .Object@num_regime="integer",
            .Object@dim_latent_var="integer",
            .Object@dynamics =  "dynrDynamics",
            .Object@measurement = meas,
            .Object@initial = initial,
            .Object@regimes= "dynrRegimes",
            .Object@noise = noise,
            .Object@infile="character",
            .Object@outfile="character",
            .Object@isContinuousTime=FALSE,
            .Object@verbose=FALSE,
            .Object@compileLib=TRUE,
            .Object@xstart="vector",
            .Object@ub="vector",
            .Object@lb="vector",
            .Object@options=default.model.options
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
    num_regime=dynrModel@num_regime,
    dim_latent_var=dynrModel@dim_latent_var,
    xstart=dynrModel@opt.control@xstart,
    ub=dynrModel@opt.control@ub,
    lb=dynrModel@opt.control@lb,
    options=dynrModel@opt.control@options,
    isContinuousTime=dynrModel@isContinuousTime,
    infile=dynrModel@infile,
    outfile=dynrModel@outfile,
    compileLib=dynrModel@compileLib,
    verbose=dynrModel@verbose
  )
  #Others keep unchanged
  #model$func_addresss = dynr.funcaddress()
}