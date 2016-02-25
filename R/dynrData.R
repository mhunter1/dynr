
#setClass(Class = "dynrData",
#	representation = representation(
#	data
#	)
#)

#setMethod("initialize", "dynrData",
#	function(.Object){
#		return(.Object)
#	}
#)

#setMethod("summary", "dynrData",
#	function(.Object){
#		summary(.Object@data)
#	}
#)

#setMethod("print", "dynrData", function(x,...) { 
#	displayDynrData(x) 
#})

#setMethod("show", "dynrData", function(x,...) { 
#	displayDynrData(x) 
#})

dynr.data <- function(dataframe, id, time, observed, covariates){
  ids <- unique(dataframe[,id])
  tstart <- c(sapply(1:length(ids),function(i){min(which(dataframe[,id]%in%ids[i]))})-1,dim(dataframe)[1])
   if (!missing(covariates)){
     data.object <- list(id=dataframe[,id],tstart=as.integer(tstart),time=as.double(dataframe[,time]),observed=data.frame(as.double(dataframe[,observed])),covariates=data.frame(as.double(dataframe[,covariates])))
    names(data.object$covariates) <- paste0("covar",1:length(covariates))
   }
  else{
     data.object <- list(id=dataframe[,id],tstart=as.integer(tstart),time=as.double(dataframe[,time]),observed=data.frame(as.double(dataframe[,observed])))
     }
  names(data.object$observed) <- paste0("obs",1:length(observed))
    return(data.object)
}



