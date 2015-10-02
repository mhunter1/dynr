#------------------------------------------------------------------------------
# Plotting methods

findR = function(y){
  Rindex2 = 1:model$num_regime
  Rtoret = Rindex2[y==max(y)]
  return(Rtoret)
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
		thes = sort(sample(1:model$"num_sbj", numSubjDemo))
		par(mfrow=c(ifelse(numSubjDemo%%2 > 0,
			numSubjDemo,numSubjDemo/2),
			ifelse(numSubjDemo%%2 > 0,
			1,2)))
		if (model$num_regime > 1){
			thePr = t(x@pr_t_given_T)
			mostLikelyRegime = apply(thePr,1,findR)  
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
				lines(1:T,thesmooth[therow,j],lty=1,lwd=1,col=j)
				points(time2,thesmooth[therow,j][time2],pch=as.character(j),col=j)
			}
			stateNames = NULL
			for(j in 1:model$dim_latent_var){
				stateNames = c(stateNames, paste0('State ',j))
			}
			legend('topright',
			paste0("State",as.character(1:model$dim_latent_var)),
				lty=1:model$dim_latent_var,
				lwd=2,col=1:model$dim_latent_var,bty="n",
				pch=as.character(1:model$dim_latent_var),cex=legend.cex)
		}
		on.exit(par(opar))
	}
)

dynr.ggplot <- function(){
	#library(ggplot2)
	#library(reshape2)
	
	## ----plotting the 6-state solution---------------------------------------
	economics6FitProbs <- posterior(economics6Fit)
	pDats <- economics
	
	
	pDats$EndDate <- c(economics$date[2:nrow(economics)], max(economics$date)+31)
	#---------Add this End time variable 
	pDats$pHigh <- economics6FitProbs$S6
	
	pDats$State <- as.factor(economics6FitProbs$state)
	#---------State is factor
	
	#-----Melt the data: Define id and measure variables-----------
	pMelted <- melt(pDats, id=c("date", "EndDate", "State"), measure=c("median.unemployment.duration.in.week", "pHigh"))
	
	#——Facet plots :Probabilities vs. values
	ggplot(pMelted, aes(date, value), 
		main = "median unemployment duration and probability of High state") + 
		geom_line() + 
		facet_grid(variable~., scales="free_y")
	
	#-----rectangular state background, Likely State By Color
	ggplot(pMelted, aes(date, value), 
		main = "median unemployment duration and probability of High state") + 
		geom_line() + 
		geom_rect(aes(xmin=date, xmax=EndDate, ymin=-Inf, ymax=Inf, fill=State), alpha=.25)
}


