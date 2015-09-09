#dyn.load(file.path("~/Dropbox/Brekfis/dynr/src",paste("wrappernegloglike",.Platform$dynlib.ext,sep="")))

#require(dynr)

dynr.run<- function(model,data) {
  tmp=.Call('main_R',model,data, PACKAGE='dynr')
  return(tmp)
}

dynr.data<-function(dataframe,id,time,observed,covariates){
  ids=unique(dataframe[,id])
  tstart=c(sapply(1:length(ids),function(i){min(which(dataframe[,id]%in%ids[i]))})-1,dim(dataframe)[1])
  data.object<-list(tstart=tstart,time=dataframe[,time],observed=data.frame(dataframe[,observed]),covariates=data.frame(dataframe[,covariates]))
  names(data.object$observed)<-paste0("obs",1:length(observed))
  names(data.object$covariates)<-paste0("covar",1:length(covariates))
  return(data.object)
}


#-------------------------------------------------------------------------

#mydata<-read.table("~/Dropbox/Brekfis/dynr/data/dataPANAsim.txt")#Missing data NA
#data<-dynr.data(mydata, id="V1", time="V2",observed=paste0('V', 3:4), covariates=paste0('V', 5))
#model.defaults <- list(blah)#starting values, ub,lb

# the model
model<-list(num_sbj=217,
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

starttime=proc.time()
x<-dynr.run(model,data)
(time=proc.time()-starttime)
x
