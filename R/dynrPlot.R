#------------------------------------------------------------------------------
# Plotting methods

findRegime = function(prop){
	RegimeIndices = 1:length(prop)
	MostLikelyRegime = RegimeIndices[prop==max(prop)]
	return(MostLikelyRegime[1])
	#If the probabilities are equal, return the first regime
}

Mode <- function(y) {
	uy <- unique(y)
	uy[which.max(tabulate(match(y, uy)))]
}

setMethod("plot", "dynrRun",
	function(x, y=NULL, data, graphingPar=par(no.readonly = TRUE), ylab = "Smoothed state values", xlab = "Time", numSubjDemo=2, legend.cex=1.2){
		opar = par(no.readonly = TRUE)
		par(graphingPar)
		thesmooth = data.frame(t(x@eta_smooth_final))  
		colnames(thesmooth) = paste0("state",1:model$dim_latent_var)
		ID = data[["id"]]
		rowIndex = 1:length(ID)
		uniID <- sort(unique(ID))
		thes = sort(sample(1:length(uniID), numSubjDemo))
		par(mfrow=c(ifelse(numSubjDemo%%2 > 0,
			numSubjDemo,numSubjDemo/2),
			ifelse(numSubjDemo%%2 > 0,
			1,2)))
		if (model$num_regime > 1){
			thePr = t(x@pr_t_given_T)
			mostLikelyRegime = apply(thePr,1,findRegime)  
			theR = Mode(mostLikelyRegime)
		}
		for (s in 1:numSubjDemo){
			T <- length(ID[ID %in% uniID[thes[s]]])
			therow = rowIndex[ID %in% uniID[thes[s]]]
	 		plot(
				c(1,T),
				c(min(thesmooth)-1, max(thesmooth)+quantile(unlist(thesmooth),.1)),
				ylab=ifelse(exists("ylab"),ylab,"State values"), 
				xlab=ifelse(exists("xlab"),xlab,"Time"),
				main=ifelse(exists("main"),main,""),
				type='n')
			if (model$num_regime > 1){
				times=1:T
				thepri = mostLikelyRegime[therow]
				rect(times[thepri==theR]-.5,quantile(unlist(thesmooth),.01),times[thepri==theR]+.5,quantile(unlist(thesmooth),.99),col="yellow",density=30)
				#legend('topleft',paste0("Regime",theR), 	
				#bty="n",cex=1.4,col=c(NA),fill=c("yellow"),
				#density=c(100))
			}
			time2 = if(T>500){
				seq(1,T,9)
			}else {1:T}
			
			for (j in 1:model$dim_latent_var){
				lines(1:T,thesmooth[therow,j], lty=1, lwd=1, col=j)
				points(time2, thesmooth[therow,j][time2], pch=as.character(j), col=j)
			}
			stateNames = NULL
			for(j in 1:model$dim_latent_var){
				stateNames = c(stateNames, paste0('State ',j))
			}
			legend('topright',
			paste0("State",as.character(1:model$dim_latent_var)),
				lty=1:model$dim_latent_var,
				lwd=2, col=1:model$dim_latent_var, bty="n",
				pch=as.character(1:model$dim_latent_var), cex=legend.cex)
		}
		on.exit(par(opar))
	}
)


dynr.ggplot <- function(x, data.dynr, states, names.state=paste0("state", states), title="Smoothed State Values", numSubjDemo=2){
	dims=dim(x@eta_regime_regime_t_pred)
	dim_latent_var=dims[1]
	num_regime=dims[2]
	num_sbj=length(unique(data.dynr$id))
	randid=sample(unique(data.dynr$id),numSubjDemo)
	
	data.plot<-data.frame(id=as.factor(data.dynr$id),time=data.dynr$time)
	#data.plot<-cbind(data.plot,data.dynr$observed)
	
	addendtime <- function(data.frame){
		data.frame$endtime <- c(data.frame$time[2:nrow(data.frame)], data.frame$time[nrow(data.frame)]+min(diff(data.frame$time)))
		return(data.frame)
	}
	if(num_regime==1){
		names.id.vars <- c("id","time")
		names.measure.vars <- names.state
		data.plot <- cbind(data.plot, t(x@eta_smooth_final[states,]))
		names(data.plot) <- c(names.id.vars, names.measure.vars)
		
		data_long<- melt(data.plot[data.plot$id%in%randid,],
			id.vars=names.id.vars,
			measure.vars=names.measure.vars,
			value.name="value")
		ggplot(data_long, aes(x=time, y=value, colour = variable)) + geom_line(size=1, alpha=0.8) + facet_wrap(~id) + 
			ggtitle(title) + 
			theme(plot.title = element_text(lineheight=.8, face="bold"))
	}else{
		
		data.plot <- ddply(data.plot,"id",addendtime)
		data.plot$regime<-as.factor(apply(x@pr_t_given_T,2,findRegime))
		names.id.vars<-c("id","time","endtime","regime")
		names.measure.vars<-paste0("states",states)
		data.plot<-cbind(data.plot,t(x@eta_smooth_final[states,]))
		names(data.plot)<-c(names.id.vars,names.measure.vars)
		
		data_long<- melt(data.plot[data.plot$id%in%randid,],
			id.vars=names.id.vars,
			measure.vars=names.measure.vars,
			value.name="value")
		
		ggplot(data_long,aes(x=time, y=value, group=variable)) +
			geom_point(size=2, aes(color=variable)) +
			geom_line(size=1, aes(color=variable)) +
			geom_rect(aes(xmin=time, xmax=endtime, ymin=-Inf, ymax=Inf, fill=regime), alpha=.15) +
			facet_wrap(~id)+
			ggtitle(title) +
			theme(plot.title=element_text(lineheight=.8, face="bold"), legend.position="bottom", legend.text=element_text(lineheight=0.8, face="bold"))
		
	}
}






