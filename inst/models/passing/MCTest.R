###########################################################
# A Monte Carlo simulation example
#
# Author: Sy-Miin Chow
# Last modified: 11/9/2015
# Two sample size configurations:
# (1) n = 20, nT = 1000 (mirroring the empirical data)
# (2) n = 40, nT = 500 (same total sample size)
##########################################################

require('deSolve') #Library to needed for ODE data generation
require('dynr')
require('matrixcalc')
library(pryr)#check memory usuage

datadir = "./data"


#These are the arguments provided in the shell script
nconfig = 1# sample size configurations
rStart = 1# first MC run
noMC = 2 #number of MC replications
r1 = 1 #batch number
Model = 2
fitModel=2
options(scipen=999)
#onServer = 0
#isHome = 1
#if (onServer==0){
#nconfig = 1# sample size configurations
#rStart = 1# first MC run
#noMC = 2 #number of MC replications
#r1 = 1 #batch number
#Model = 3 #Model number (1 = without covariate in transition LO; 2 = with covariate in transition LO)
#dynrdir = "~/Dropbox/Symiin_Lu/dynr"
#fitModel = 3
#setwd(dynrdir)
#MCdir = "~/Dropbox/NLMEode/RegimeSwitching/MonteCarlo"
#}

# -------------------------------------------------------------
#(1) First we define some functions that are helpful for 
#computing outcome measures for the simulation
# -------------------------------------------------------------
isSimulation=1


if (nconfig==1){
  n = 20; nT = 1000
}else{
  n = 40; nT = 500  
}


set.seed(219*n/nT*r1)
#options(width=110)

ne = 2 # size of latent variable vector or "prevState" in the ODE model
ny = 2
isSimulation = 1
#-------------------------------------------------------------
#(2) Set up control info that is invariant across MC
# -------------------------------------------------------------
#lastT = 50
SDerror1 = 3
SDerror2 = 3
switch(Model,
       {#Model 1 - 2R model with no covariate in transition LO function
       nx=0
       LOpar = c(-4,4+4.5,0,0,0,0)
       parms <- c(log(.2), log(.1),log(.3), log(.2), 
                  100,log(SDerror1^2),log(SDerror2^2),LOpar[1:2]) 
       truePar = c(exp(parms[1:4]),parms[5],
                   SDerror1^2,SDerror2^2,LOpar[1:2])
       },
       {#Model 2 - 2R model with covariates in transition LO function
       nx=2
       LOpar = c(-4,4+4.5,-1,1,-2,-1)   
       parms <- c(log(.2), log(.1),log(.3), log(.2), 
                   100,log(SDerror1^2),log(SDerror2^2),LOpar) 
       truePar = c(exp(parms[1:4]),parms[5],
                   SDerror1^2,SDerror2^2,LOpar)
        },
       {#Model 3 - 1R model
         nx=0
         LOpar = NULL
         parms <- c(log(.2), log(.1),log(.3), log(.2), 
                    100,log(SDerror1^2),log(SDerror2^2)) 
         truePar = c(exp(parms[1:4]),parms[5],
                     SDerror1^2,SDerror2^2)
       }
)

switch(fitModel,
       {
         noDynPar = 4
         noBase = 1
         noLOPar = 2
         noPar = noDynPar + noLOPar + noBase + ny
         xstart0=c(rep(log(.1), 4), 95.0, log(10.0), log(10.0), -3.0, 9.0)
         ub=c(rep(1.5, 4), 200, 5, 5, rep(30, 2))
         lb=c(rep(-20, 4), 50, -10, -10, rep(-30, 2))
         tfun <- function(x){c(exp(x[1:4]),x[5],exp(x[6:7]), x[8:9])}
       },
       {
         noDynPar = 4
         noBase = 1
         noLOPar = 6
         noPar = noDynPar + noLOPar + noBase + ny
         xstart0=c(rep(log(.1), 4), 95.0, log(10.0), log(10.0), -3.0, 9.0, -1.5, 0.5, -1,-.3)
         ub=c(rep(1.5, 4), 200, 5, 5, rep(30, 6))
         lb=c(rep(-20, 4), 50, -10, -10, rep(-30, 6))
         tfun <- function(x){c(exp(x[1:4]),x[5],exp(x[6:7]), x[8:13])}
       },
       {         
         noDynPar = 4
         noBase = 1
         noLOPar = 0
         noPar = noDynPar + noLOPar + noBase + ny
         xstart0=c(rep(log(.1), 4), 95.0, log(10.0), log(10.0))
         ub=c(rep(1.5, 4), 200, 5, 5)
         lb=c(rep(-20, 4), 50, -10, -10)
         tfun <- function(x){c(exp(x[1:4]),x[5],exp(x[6:7]))}
       })


allExitFlag = matrix(NA,noMC,2)
allPar = matrix(NA, noMC, noPar+1)
allSE = matrix(NA, noMC, noPar+1)
allInt = matrix(NA,noMC,noPar*2+1)
include0Matrix = matrix(NA,noMC,noPar+1)
if (Model==fitModel) includeTrueParMatrix = matrix(NA,noMC,noPar+1)
allLVRMSE = matrix(NA,noMC,ne+1)
allLik = matrix(NA,noMC, 4)
print(fitModel)



switch(fitModel,
{
  infile = paste0("./RSODEmodelNoCov.c")
  outfile = paste0("./RSODEmodelNoCov")
},
{
  infile=paste0("./RSODEmodel.c")
  outfile = paste0("./RSODEmodel")
},
{
  infile=paste0("./RSODEmodel1R.c")
  outfile = paste0("./RSODEmodel1R")  
}
)
themodel <- dynr.model(
  num_regime=ifelse(fitModel < 3, 2, 1),
  dim_latent_var=2,
  xstart=xstart0,
  ub=ub,lb=lb,
  options=list(maxtime=60*60, maxeval=1000,ftol_rel=as.numeric(1e-10),
               xtol_rel=as.numeric(1e-7)),
  isContinuousTime=TRUE,
  infile=infile, 
  outfile=outfile, 
  verbose=TRUE,
  compileLib=FALSE
)
sessionInfo()
#Perform a bunch of Monte Carlo replications
#---------------------------------------------

for (r in rStart:noMC){
  run = r+(r1-1)*noMC
  print(paste(c("MC run = ", r)))
  allExitFlag[r,1] = run
  allPar[r,1] = run
  allSE[r,1] = run
  allInt[r,1] = run
  include0Matrix[r,1] = run
  if (fitModel==Model)includeTrueParMatrix[r,1] = run
  allLVRMSE[r,1] = run
  allLik[r,1] = run
  
  if (nx > 0){
  thedata = read.table(paste0(datadir,"/Model",fitModel,"T",nT,"n",n,"batch",run,"ODEsimData.txt"),col.names=c("ID","time",paste0("DV",1:ne),paste0("TV",1:nx))) 
  thedata[,"TV2"] = as.numeric(thedata[,"TV2"])
  data <- dynr.data(thedata, id="ID", time="time",observed=paste0('DV', 1:2),covariates=paste0('TV', 1:2))
  }else{
    thedata = read.table(paste0(datadir,"/Model",Model,"T",nT,"n",n,"batch",run,"ODEsimData.txt"),col.names=c("ID","time",paste0("DV",1:ne))) 
    data <- dynr.data(thedata, id="ID", time="time",observed=paste0('DV', 1:2))  
}

  print(pryr::mem_used())#function from the pryr library  
  if (exists("res")) {
   #rm(res)
   print(mem_change(rm(res)))
  }
  
  #print(getLoadedDLLs())
  #print(themodel$func_address)
  res <- dynr.run(themodel, data,tfun)
  print(pryr::mem_used())#function from the pryr library   
  summary(res)
         
  #------Compile results from dynr run ---------
  if (!is.null(res)){
    if (length(res@standard.errors[is.na(res@standard.errors)]) == 0 && res@exitflag >= 6 && res@exitflag < 10){
      
      parhat = coef(res)
      sehat = res@standard.errors
      include0Flag = rep(NA,noPar)
      includeTrueParFlag = rep(NA,noPar)
      for (k in 1:noPar){
        toCompare = res@conf.intervals[k,]
        #See if interval includes 0. If so = 1, else  = 0
        include0Flag[k] = ifelse(toCompare[1] < 0.000 && toCompare[2] > 0.000,1,0)
        if (fitModel==Model){
          includeTrueParFlag[k] = ifelse(toCompare[1] <= truePar[k] && toCompare[2] >= truePar[k],1,0)}
      }
      
      include0Matrix[r,2:(noPar+1)] = include0Flag
      if (fitModel==Model){includeTrueParMatrix[r,2:(noPar+1)] = includeTrueParFlag}
      allPar[r,2:(noPar+1)] = parhat
      allSE[r,2:(noPar+1)] = sehat
      allInt[r,2:(noPar*2 + 1)] = t(matrix(unlist(t(res@conf.intervals)),ncol=1,byrow=F))
      
      allExitFlag[r,2] = res@exitflag
      allLik[r,2] = 2*res@neg.log.likelihood
      allLik[r,3] = AIC(res)
      allLik[r,4] = BIC(res) 
    }
  }#Check if dynr object exists
  
  #End of 1 replication
  #---------------------------------------------

}#End of noMC loop

#q(save="no")
