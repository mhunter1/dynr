
dynr.data <- function(dataframe, id, time, observed, covariates){
  ids <- unique(dataframe[,id])
  tstart <- c(sapply(1:length(ids),function(i){min(which(dataframe[,id]%in%ids[i]))})-1,dim(dataframe)[1])
  data.object <- list(tstart=tstart,time=dataframe[,time],observed=data.frame(dataframe[,observed]),covariates=data.frame(dataframe[,covariates]))
  names(data.object$observed) <- paste0("obs",1:length(observed))
  names(data.object$covariates) <- paste0("covar",1:length(covariates))
  return(data.object)
}



