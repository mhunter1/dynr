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

dynr.prep.Nametochange<-function(dynamics,measurement,noise,regimes,initial){
  #take in a dynrModel object
  obj.dynrModel=new("dynrModel",dynamics=dynamics,measurement=measurement,noise=noise,regimes=regimes,initial=initial)
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