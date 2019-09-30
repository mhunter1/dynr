#------------------------------------------------------------------------------
# Plotting methods

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

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

#Old plot function for dynrCook object
#setMethod("plot", "dynrCook",
#          function(x, y=NULL, data.dynr, graphingPar=par(no.readonly = TRUE), ylab = #"Smoothed state values", xlab = "Time", numSubjDemo=2, legend.cex=1.2){
#            opar = par(no.readonly = TRUE)
#            par(graphingPar)
#            thesmooth = data.frame(t(res@eta_smooth_final))  
#            
#            dim_latent_var=dim(res@eta_smooth_final)[1]
#            num_regime=dim(res@pr_t_given_T)[1]
#            
#            colnames(thesmooth) = paste0("state",1:dim_latent_var)
#            
#            ID = data.dynr[["id"]]
#            rowIndex = 1:length(ID)
#            uniID <- sort(unique(ID))
#            thes = sort(sample(1:length(uniID), numSubjDemo))
#            
#            par(mfrow=c(ifelse(numSubjDemo%%2 > 0,
#                               numSubjDemo,numSubjDemo/2),
#                        ifelse(numSubjDemo%%2 > 0,
#                               1,2)))
#            if (num_regime > 1){
#              thePr = t(x@pr_t_given_T)
#              mostLikelyRegime = apply(thePr,1,findRegime)  
#              theR = Mode(mostLikelyRegime)
#            }
#            for (s in 1:numSubjDemo){
#              T <- length(ID[ID %in% uniID[thes[s]]])
#              therow = rowIndex[ID %in% uniID[thes[s]]]
#              plot(
#                c(1,T),
#                c(min(thesmooth)-1, max(thesmooth)+quantile(unlist(thesmooth),.1)),
#                ylab=ifelse(exists("ylab"),ylab,"State values"), 
#                xlab=ifelse(exists("xlab"),xlab,"Time"),
#                main=ifelse(exists("main"),main,""),
#                type='n')
#              if (num_regime > 1){
#                times=1:T
#                thepri = mostLikelyRegime[therow]
#                rect(times[thepri==theR]-.5,quantile(unlist(thesmooth),.01),times[th#epri==theR]+.5,quantile(unlist(thesmooth),.99),col="yellow",density=30)
#                #legend('topleft',paste0("Regime",theR), 	
#                #bty="n",cex=1.4,col=c(NA),fill=c("yellow"),
#                #density=c(100))
#              }
#              time2 = if(T>500){
#                seq(1,T,9)
#              }else {1:T}
#              
#              for (j in 1:dim_latent_var){
#                lines(1:T,thesmooth[therow,j], lty=1, lwd=1, col=j)
#                points(time2, thesmooth[therow,j][time2], pch=as.character(j), col=j#)
#            }
#            stateNames = NULL
#            for(j in 1:dim_latent_var){
#              stateNames = c(stateNames, paste0('State ',j))
#            }
#            legend('topright',
#                   paste0("State",as.character(1:dim_latent_var)),
#                   lty=1:dim_latent_var,
#                   lwd=2, col=1:dim_latent_var, bty="n",
#                   pch=as.character(1:dim_latent_var), cex=legend.cex)
#          }
#          on.exit(par(opar))
#        }
#)

##' Plot method for dynrCook objects
##' 
##' @param x dynrCook object
##' @param dynrModel model object
##' @param style The style of the plot in the first panel. If style is 1 (default), user-selected smoothed state variables are plotted. If style is 2, user-selected observed-versus-predicted values are plotted.
##' @param names.state (optional) The names of the states to be plotted, which should be a subset of the state.names slot of the measurement slot of dynrModel.
##' @param names.observed (optional) The names of the observed variables to be plotted, which should be a subset of the obs.names slot of the measurement slot of dynrModel.
##' @param printDyn A logical value indicating whether or not to plot the formulas for the dynamic model
##' @param printMeas A logical value indicating whether or not to plot the formulas for the measurement model
##' @param textsize numeric. Font size used in the plot.
##' @param ... Further named arguments
##' 
##' @details
##' This is a wrapper around \code{\link{dynr.ggplot}}.  A great benefit of it is that it shows the model equations in a plot.
plot.dynrCook <- function(x, dynrModel, style = 1, names.state, names.observed, printDyn=TRUE, printMeas=TRUE, textsize=4, ...) {
  #The first panel is the ggplot
  if(missing(names.observed)){names.observed=dynrModel@measurement@obs.names}
  if(missing(names.state)){names.state=dynrModel@measurement@state.names}
  p1 <- dynr.ggplot(x, dynrModel, style, numSubjDemo=1, names.state=names.state, names.observed=names.observed, ...)

  #If there are more than 2 regimes, the second panel shows the histogram of the most probable regimes across time and subjects.
  if (dynrModel$num_regime>1){
    regime <- NULL
    highProbR <- data.frame(regime=apply(x@pr_t_given_T,2,which.max))
    p2 <- ggplot2::ggplot(data = highProbR, ggplot2::aes(factor(regime))) +
      ggplot2::geom_bar() +
      ggplot2::scale_fill_brewer(palette = 3) +
      ggplot2::xlab('Regime') + ggplot2::ylab('Counts') + ggplot2::labs(fill = '') + 
      ggplot2::ggtitle('Counts by most probable regime') +
      ggplot2::theme(axis.title = ggplot2::element_text(size=14),
                     axis.text=ggplot2::element_text(size=12),
                     plot.title = ggplot2::element_text(size = 12, colour = "black", face = "bold", hjust = 0.5, vjust=0.01))
  }
  
  #The third panel plots the model formulae
  p3 <- plotFormula(dynrModel, ParameterAs=sprintf("%.2f", x@transformed.parameters), printDyn=printDyn, printMeas=printMeas, textsize=textsize)

  #Organize the panels using the multiplot function
  if (dynrModel$num_regime > 1){
    multiplot(p1, p2, p3, cols = 1, layout=matrix(c(1,1,2,3), nrow=2, byrow=TRUE))
  }else{
    multiplot(p1, p3, layout=matrix(c(1,1,2,2), nrow=2, byrow=TRUE))
  }

}

plotdf <- function(vec_tex){
  dataframe=data.frame(text=sapply(paste0("$",vec_tex,"$"),function(x){as.character(latex2exp::TeX(x))}))
  dataframe$x<-0
  return(dataframe)
}

##' Plot of the estimated frequencies of the regimes across all individuals and time points
##' based on their smoothed regime probabilities
##' 
##' @param res The dynr object returned by dynr.cook().
##' @param dynrModel The model object to plot.
##' @param names.regime (optional) Names of the regimes (must match the length of the number of regimes)
##' @param title (optional) Title of the plot.
##' @param xlab (optional) Label of the x-axis.
##' @param ylab (optional) Label of the y-axis.
##' @param textsize (default = 12) Text size for the axis labels and title (= textsize + 2).
##' @param print (default = TRUE) A flag for whether the plot should be printed.
dynr.plotFreq <- function(res, dynrModel, names.regime, title, xlab, ylab, textsize=12,print=TRUE) {
  if(missing(names.regime)){names.regime<-paste0("Regime",1:dynrModel@num_regime)}
  if(missing(title)){title<-'Counts by most probable regime'}
  if(missing(xlab)){xlab<-'Regime'}
  if(missing(ylab)){ylab<-'Counts'}

    regime <- NULL
    highProbR <- data.frame(regime=factor(apply(res@pr_t_given_T,2,which.max),
                                   levels=1:dynrModel@num_regime,
                                   labels=names.regime))
    p2 <- ggplot2::ggplot(data = highProbR, ggplot2::aes(factor(regime))) +
      ggplot2::geom_bar() +
      ggplot2::scale_fill_brewer(palette = 3) +
      ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::labs(fill = '') + 
      ggplot2::ggtitle(title) +
      ggplot2::theme(axis.title = ggplot2::element_text(size=textsize+2),
                     axis.text=ggplot2::element_text(size=textsize),
                     plot.title = ggplot2::element_text(size = textsize, colour = "black", face = "bold", hjust = 0.5, vjust=0.01))
    if (print) print(p2)
    return(p2)
  }
  

##' Plot the formula from a model
##' 
##' @param dynrModel The model object to plot.
##' @param ParameterAs The parameter values or names to plot. The underscores in parameter names are saved for 
##' use of subscripts.  Greek letters can be specified as corresponding LaTeX symbols without backslashes (e.g., "lambda") 
##' and printed as greek letters.
##' @param printDyn A logical value indicating whether or not to plot the formulas for the dynamic model.
##' @param printMeas A logical value indicating whether or not to plot the formulas for the measurement model
##' @param printRS logical. Whether or not to print the regime-switching model. The default is FALSE.
##' @param textsize The text size use in the plot.
##' 
##' @details
##' This function typesets a set of formulas that represent the model.  Typical inputs to the \code{ParameterAs} argument are (1) the starting values for a model, (2) the final estimated values for a model, and (3) the parameter names.  These are accessible with (1) \code{model$xstart}, (2) \code{coef(cook)}, and (3) \code{model$param.names} or \code{names(coef(cook))}, respectively.
plotFormula <- function(dynrModel, ParameterAs, printDyn=TRUE, printMeas=TRUE,
                        printRS=FALSE, textsize=4){
  
  dynrModel <- PopBackModel(dynrModel, LaTeXnames(ParameterAs, latex = FALSE))

  #Dynamic model
  if (printDyn) {  
    if (class(dynrModel$dynamics) == "dynrDynamicsMatrix"){
      state.names <- (dynrModel$measurement)$state.names
      exo.names <- (dynrModel$dynamics)$covariates
      dynrModel@dynamics@values.dyn <- lapply((dynrModel$dynamics)$values.dyn, preProcessNames,state.names,state.names)
      dynrModel@dynamics@values.exo <- lapply((dynrModel$dynamics)$values.exo, preProcessNames,state.names,exo.names)
      dynrModel@dynamics@values.int <- lapply((dynrModel$dynamics)$values.int, preProcessNames,state.names)
    }

    dyn.df <- data.frame(text="bold('Dynamic Model')",x=0)
    nRegime=ifelse(class(dynrModel@dynamics)=="dynrDynamicsFormula",
                   length(dynrModel@dynamics@formula),
                   length(dynrModel@dynamics@values.dyn))
    dyn_tex=printex(dynrModel@dynamics,AsMatrix=FALSE)
    nEq <- length(dyn_tex[[1]])
    # noise_tex <- paste0(" + ", ifelse(dynrModel@dynamics@isContinuousTime, "d", ""), "w_{",1:nEq, "}(t)")
    noise_tex <- paste0(ifelse(dynrModel@dynamics@isContinuousTime, "d", ""), 
                        "w_{", 1:nEq, "}(t)")
    for (i in 1:nRegime){
      if (nRegime>1){
        dyn.df <- rbind(dyn.df,data.frame(text=paste0("'Regime ",i,":'"),x=0))
      }
      if (length(dynrModel@noise@values.latent)==1){
        values.latent.mat <- dynrModel@noise@values.latent[[1]] 
      }else if (length(dynrModel@noise@values.latent)==nRegime){
        values.latent.mat <- dynrModel@noise@values.latent[[i]]
      }else{stop("The number of regimes implied by the dynamic noise structure does not match the number of regimes in the dynamic model.")}
      pos.zero <- diag(values.latent.mat) == 0
      # noise_tex[diag(values.latent.mat)==0] <- ""
      # dyn.df <- rbind(dyn.df,plotdf(LaTeXnames(paste0(dyn_tex[[i]], noise_tex))))
      dyn.df <- rbind(dyn.df,
                      plotdf(LaTeXnames(
                        paste0(dyn_tex[[i]], 
                               ifelse(pos.zero, "", 
                                      paste0(" + ", noise_tex))))),
                      plotdf(LaTeXnames(
                        paste0(ifelse(pos.zero, "", 
                                      paste0(noise_tex, " \\sim ", "N(0,\\,", 
                                             diag(values.latent.mat), ")"))))) )
    }
  } else {
    dyn.df <- plotdf("")
  }
  
	#Measurement model
	if (printMeas) {
		meas.df <- data.frame(text="bold('Measurement Model')", x=0)
		nRegime <- length(dynrModel$measurement$values.load)
		meas_tex <- printex(dynrModel$measurement, AsMatrix=FALSE)
		nEq <- length(meas_tex[[1]])
		# noise_tex <- paste0(" + ", "epsilon_{", 1:nEq, "}")
		noise_tex <- paste0("epsilon_{", 1:nEq, "}")
		for (i in 1:nRegime){
			if (nRegime > 1){
				meas.df <- rbind(meas.df, data.frame(text=paste0("'Regime ",i,":'"), x=0))
			}
			if (length(dynrModel$noise$values.observed) == 1){
				values.observed.mat <- dynrModel$noise$values.observed[[1]] 
			}else if (length(dynrModel$noise$values.observed) == nRegime){
				values.observed.mat <- dynrModel$noise$values.observed[[i]]
			}else{stop("The number of regimes implied by the measurement noise structure does not match the number of regimes in the measurement model.")}
			pos.zero <- diag(values.observed.mat) == 0
			# noise_tex[diag(values.observed.mat)==0] <- ""
			# meas.df<-rbind(meas.df,plotdf(LaTeXnames(paste0(meas_tex[[i]], noise_tex))))
			meas.df <- rbind(meas.df,
							plotdf(LaTeXnames(
								paste0(meas_tex[[i]],
										ifelse(pos.zero, "", paste0(" + ", noise_tex, ",")),
										"\\;\\,",
										ifelse(pos.zero,
											"",
											paste0(noise_tex, " \\sim ", "N(0,\\,", diag(values.observed.mat), ")"))))))
		}
	} else {
		meas.df <- plotdf("")
	}
	
  # Regime-switching model
  if (printRS) {
    rs.df <- data.frame(text="bold('Regime-switching Model')", x=0)
    if (any(dim(dynrModel@regimes@params) == 0)) {
      warning("No Regime-switching Model. Please turn off 'printRS'.")
    } else {
      # helper function to make sub. e.g., a_111 to a_{111}
      to_sub <- function(parnames) {
        idx_sub <- grepl('_', parnames)
        if ( sum(idx_sub)==0 ) { return(parnames) }
        subs <- parnames[idx_sub]
        subs_spl <- strsplit(subs, '_')
        subbed <- sapply(subs_spl, function(sub) {
          paste0(sub[1], '_', paste0('{', sub[2], '}') )
        })
        parnames[idx_sub] <- subbed
        parnames 
      }
      covars <- dynrModel@regimes@covariates
      if ( length(covars) == 0 ) {
        covars <- c("")
      } else {
        covars <- c("", paste0("*", covars) )
      }
      covars_n <- length(covars)
      # TODO: check nregime to fit with others
      nregime <- nrow(dynrModel@regimes@params)
      # column index for the intercepts of regimes
      inter_col <- cumsum(c(1, rep(covars_n, nregime-2)))
      # i:(i+ncov-1)
      parnames <- to_sub( dynrModel@regimes@paramnames )
      param_rs <- matrix(parnames, nregime)
      param_inter <- param_rs[, inter_col, drop=FALSE]# intercepts
      refrow <- dynrModel@regimes@refRow
      if (length(refrow) != 0) {
        for ( co in inter_col ) {
          param_rs[-refrow, co] <- 
            paste0(param_rs[refrow, co], ' + ', param_rs[-refrow, co])
        }
      }
      param_rs <- matrix(sweep(param_rs, 2, covars, FUN=paste0), 
                            nrow=nregime)
      ll <- vector('list', length=length(inter_col))
      for ( i in 1:length(inter_col) ) {
        ll[[i]] <- 
          param_rs[, inter_col[i]:(inter_col[i]+covars_n-1), 
                   drop=FALSE]
      }
      sl <- sapply(ll, function(l) {
        apply(l, 1, paste, collapse=" + ")
      })
      # transition probability subs
      transub <- array(c(row(sl), col(sl)), c(nrow(sl), ncol(sl), 2))
      transubed <- apply(transub, c(1,2), paste, collapse="")
      rs.df <- rbind(rs.df,
                     plotdf(paste0("Log-odds(p_{",
                                   transubed, "}) = ",
                                   sl)))
  } 
    } else {
    rs.df <- plotdf("")
  }
  
  plot.df <- rbind(meas.df, dyn.df, rs.df)

  plot.df$y=seq((dim(plot.df)[1]-1)*5+1, 1, by=-5)
  
  x <- NULL
  y <- NULL
  fig<-ggplot2::ggplot(plot.df, ggplot2::aes(x=x, y=y, label=text))+
    ggplot2::geom_text(parse=TRUE,size=textsize)+
    ggplot2::theme_void() + ggplot2::ylim(min(plot.df$y)-5, max(plot.df$y)+5)
  
  return(fig)
}

##' The ggplot of the smoothed state estimates and the most likely regimes
##'
##' @aliases autoplot.dynrCook
##' 
##' @param res The dynr object returned by \code{dynr.cook()}.
##' @param dynrModel The model object to plot.
##' @param style The style of the plot. If style is 1 (default), user-selected smoothed state variables are plotted. If style is 2, user-selected observed-versus-predicted values are plotted.
##' @param numSubjDemo The number of subjects to be randomly selected for plotting.
##' @param idtoPlot Values of the ID variable to plot.
##' @param names.state (optional) The names of the states to be plotted, which should be a subset of the state.names slot of the measurement slot of dynrModel.
##' @param names.observed (optional) The names of the observed variables to be plotted, which should be a subset of the obs.names slot of the measurement slot of dynrModel.
##' @param names.regime (optional) The names of the regimes to be plotted, which can be missing.
##' @param shape.values (optional) A vector of values that correspond to the shapes of the points, which can be missing. See the R documentation on pch for details on possible shapes.
##' @param title (optional) A title of the plot.
##' @param ylab (optional) The label of the y axis.
##' @param is.bw Is plot in black and white? The default is FALSE.
##' @param colorPalette A color palette for lines and dots. It is a value passed to the palette argument of the \code{ggplot2::scale_colour_brewer()} function. These palettes are in the R package \pkg{RColorBrewer}. One can find them by attaching the package with \code{library(RColorBrewer)} and run \code{display.brewer.all()}.
##' @param fillPalette A color palette for blocks. It is a value passed to the palette argument of the \code{ggplot2::scale_fill_brewer()} function. These palettes are in the package \pkg{RColorBrewer}. One can find them by attaching the package with \code{library(RColorBrewer)} and run \code{display.brewer.all()}.  
##' @param mancolorPalette (optional) A color palette for manually scaling the colors of lines and dots. It is a vector passed to the values argument of the \code{ggplot2::scale_colour_manual} function.
##' @param manfillPalette (optional) A color palette for manually scaling the colors of filled blocks. It is a vector passed to the values argument of the \code{ggplot2::scale_fill_manual} function.
##' @param ... A list of elements that modify the existing ggplot theme. Consult the \code{ggplot2::theme()} function in the R package \pkg{ggplot2} for more options.
##' 
##' @details
##' This function outputs a ggplot layer that can be modified using functions in the package \pkg{ggplot2}. That is, one can add layers, scales, coords and facets with the "+" sign. In an example below, the \code{ggplot2::ylim()} function is used to modify the limits of the y axis of the graph. More details can be found on \url{http://ggplot2.tidyverse.org} and \url{http://ggplot2.tidyverse.org/reference/}.
##'
##' The two functions \code{dynr.ggplot()} and \code{autoplot()} as identical aliases of one another.  The \code{autoplot()} function is an S3 method from the package \pkg{ggplot2} that allows many objects to be plotted and works like the base \code{plot()} function.
##'
##' @examples
##' # The following code is part of a demo example in dynr 
##' # One can obtain the yum and rsmod objects needed below by running demo(RSLinearDiscreteYang).
##' # p <- dynr.ggplot(yum, dynrModel = rsmod, style = 1,
##' # 	names.regime = c("Deactivated", "Activated"),
##' # 	title = "(B) Results from RS-AR model", numSubjDemo = 1,
##' # 	shape.values = c(1),
##' # 	text = element_text(size = 16),
##' # 	is.bw = TRUE)
##' # One can modify the limits on the y axis by using '+'
##' # p + ggplot2::ylim(-2, 4)
##'
##' # autoplot(yum, dynrModel = rsmod, style = 1,
##' #	names.regime = c("Deactivated", "Activated"),
##' #	title = "(B) Results from RS-AR model", numSubjDemo = 1,
##' #	shape.values = c(1),
##' #	text = element_text(size = 16),
##' #	is.bw = TRUE)
dynr.ggplot <- function(res, dynrModel, style = 1,
                        numSubjDemo=2, idtoPlot=c(),
                        names.state, 
                        names.observed,
                        names.regime,
                        shape.values,
                        title, 
                        ylab, 
                        is.bw=FALSE, 
                        colorPalette="Set2", 
                        fillPalette="Set2", 
                        mancolorPalette, manfillPalette, ...){
 
  data.dynr=dynrModel@data
  dim_latent_var=dim(res@eta_smooth_final)[1]
  num_regime=dim(res@pr_t_given_T)[1]
  if(missing(names.observed)){names.observed=dynrModel@measurement@obs.names}
  observed=which(dynrModel@measurement@obs.names%in%names.observed)
  if ((length(observed)>8)&(missing(mancolorPalette))&(style==2)){stop("You provided too many variables than the default color palette can handle.\nPlease consider specify the mancolorPalette argument.")}
  if(missing(names.state)){names.state=dynrModel@measurement@state.names}
  states=which(dynrModel@measurement@state.names%in%names.state)
  if(missing(shape.values)){
    if (style==1) shape.values=states#48+states
    if (style==2) shape.values=c(rep(32,length(observed)),observed)#pre-obs
  }else{
  	if (style==2) shape.values=c(rep(32,length(shape.values)), shape.values)#pre-obs
  }
  if (style==1) line.values=rep(1,length(states))
  if (style==2) line.values=rep(c(1,0),each=length(observed))#pre-obs
  num_sbj=length(unique(data.dynr$id))
  if(length(idtoPlot)<1){
    randid =sample(unique(data.dynr$id),numSubjDemo)
  }else {randid=idtoPlot}
  if(missing(title)){
    if (style==1) title=" "
    if (style==2) title="Observed vs Predicted"
  }
  if(missing(ylab)){
    if (style==1) ylab="Smoothed State Values"
    if (style==2) ylab="Values"
  }
  
  data.plot<-data.frame(id=as.factor(data.dynr$id),time=data.dynr$time)
  
  if(num_regime==1){
    names.id.vars <- c("id","time")
    if (style==1){
      names.measure.vars <- names.state
      data.plot <- cbind(data.plot, matrix(res@eta_smooth_final[states,],ncol=length(states),byrow=TRUE))
      lines.var <- names.state
      points.var <- names.state
    }#end of style 1
    if (style==2){
      dynrModel=PopBackModel(dynrModel, LaTeXnames(res@transformed.parameters, latex = FALSE))
      dynrModel@measurement@values.load=lapply(dynrModel@measurement@values.load, function(x){matrix(as.numeric(x),ncol=ncol(x))})
      num_regime_meas <- length(dynrModel@measurement@values.load)
      if (length(dynrModel@measurement@values.exo)!=0){
        dynrModel@measurement@values.exo=lapply(dynrModel@measurement@values.exo, function(x){matrix(as.numeric(x),ncol=ncol(x))})
      }
      if (length(dynrModel@measurement@values.int)!=0){
        dynrModel@measurement@values.int=lapply(dynrModel@measurement@values.int, function(x){matrix(as.numeric(x),ncol=ncol(x))})
      }
      #predicted scores
      if (num_regime_meas==1){#one-regime measurement model
        predicted=dynrModel@measurement@values.load[[1]]%*%res@eta_smooth_final
        if (length(dynrModel@measurement@values.int)!=0){
          predicted=predicted+dynrModel@measurement@values.int[[1]]%*%matrix(rep(1,ncol(res@eta_smooth_final)), nrow=1)
        }
        if (length(dynrModel@measurement@values.exo)!=0){
          exo.mat=t(as.matrix(dynrModel@data$covariates[,paste0("covar",sapply(dynrModel@measurement@exo.names, function(x){which(dynrModel@data$covariate.names==x)}))]))
          predicted=predicted+dynrModel@measurement@values.exo[[1]]%*%exo.mat
        }
      }#end of one-regime meas
      names.measure.vars <- c(paste0(names.observed, ".predicted"), paste0(names.observed, ".observed"))
      data.plot <- cbind(data.plot, matrix(predicted[observed,],ncol=length(observed),byrow=TRUE, dimnames=list(NULL, paste0(names.observed, ".predicted"))))
      data.plot <- cbind(data.plot, data.dynr$observed[,observed])
      lines.var <- paste0(names.observed, ".predicted")
      points.var <- paste0(names.observed, ".observed")
      if (length(observed)<8){
        fancy8 <- c("#edaac1","#a7dab0","#dfb3e2","#8bc5ed","#bcb8ec","#e6b296","#d2d39d","#82d9d6")
        mancolorPalette=rep(fancy8[observed],2)
      }
    }#end of style 2
    names(data.plot) <- c(names.id.vars, names.measure.vars)
    
    data_long <- reshape2::melt(data.plot[data.plot$id%in%randid,],
                                id.vars=names.id.vars,
                                measure.vars=names.measure.vars,
                                value.name="value")
	data_long$variable<-factor(data_long$variable, levels=sort(names.measure.vars)) 
    value <- NULL
    variable <- NULL
    partial.plot<-ggplot2::ggplot(data_long, ggplot2::aes(x=time, y=value, colour = variable, shape = variable, linetype = variable)) +
      ggplot2::geom_line(data=data_long[data_long$variable%in%lines.var,], size=1) +
      ggplot2::geom_point(data=data_long[data_long$variable%in%points.var,], size=4) +
      ggplot2::scale_linetype_manual(labels=sort(names.measure.vars), values=line.values[order(names.measure.vars)])+
      ggplot2::scale_shape_manual(labels=sort(names.measure.vars), values=shape.values[order(names.measure.vars)])+
      ggplot2::facet_wrap(~id)+
      ggplot2::labs(title = title, y=ylab, ggtitle="")+
      ggplot2::theme(plot.title=ggplot2::element_text(lineheight=.8, face="bold", hjust = 0.5, vjust=0.01), legend.position="bottom",legend.text=ggplot2::element_text(lineheight=0.8, face="bold"),...)
    
  }else{#more than two regimes
    if(missing(names.regime)){names.regime = 1:num_regime}

    findRegime = function(prop){
      MostLikelyRegime = names.regime[prop==max(prop)]
      return(MostLikelyRegime[1])
      #If the probabilities are equal, return the first regime
    }
    addendtime <- function(data_frame){
      data_frame$endtime <- c(data_frame$time[2:nrow(data_frame)], data_frame$time[nrow(data_frame)]+min(diff(data_frame$time)))
      return(data_frame)
    }
    data.plot <- plyr::ddply(data.plot,"id",addendtime)
    data.plot$regime<-as.factor(apply(res@pr_t_given_T,2,findRegime))
    names.id.vars<-c("id","time","endtime","regime")
    if (style==1){
      names.measure.vars<-names.state
      data.plot<-cbind(data.plot,matrix((res@eta_smooth_final[states,]),ncol=length(states),byrow=TRUE))
      lines.var <- names.state
      points.var <- names.state
    }#end of style 1
    if (style==2){
      dynrModel=PopBackModel(dynrModel, LaTeXnames(res@transformed.parameters, latex = FALSE))
      dynrModel@measurement@values.load=lapply(dynrModel@measurement@values.load, function(x){matrix(as.numeric(x),ncol=ncol(x))})
      num_regime_meas <- length(dynrModel@measurement@values.load)
      if (length(dynrModel@measurement@values.exo)!=0){
        dynrModel@measurement@values.exo=lapply(dynrModel@measurement@values.exo, function(x){matrix(as.numeric(x),ncol=ncol(x))})
      }
      if (length(dynrModel@measurement@values.int)!=0){
        dynrModel@measurement@values.int=lapply(dynrModel@measurement@values.int, function(x){matrix(as.numeric(x),ncol=ncol(x))})
      }
      #predicted scores
      if (num_regime_meas==1){#one-regime measurement model
        predicted=matrix(dynrModel@measurement@values.load[[1]][observed,],nrow=length(observed))%*%res@eta_smooth_final
        if (length(dynrModel@measurement@values.int)!=0){
          predicted=predicted+matrix(dynrModel@measurement@values.int[[1]][observed,],nrow=length(observed))%*%matrix(rep(1,ncol(res@eta_smooth_final)), nrow=1)
        }
        if (length(dynrModel@measurement@values.exo)!=0){
          exo.mat=t(as.matrix(dynrModel@data$covariates[,paste0("covar",sapply(dynrModel@measurement@exo.names, function(x){which(dynrModel@data$covariate.names==x)}))]))
          predicted=predicted+matrix(dynrModel@measurement@values.exo[[1]][observed,],nrow=length(observed))%*%exo.mat
        }
        #end of one-regime meas
      }else if (num_regime_meas==num_regime){
        ntotaltime=ncol(res@eta_smooth_final)
        predicted=matrix(rep(0,length(observed)*ntotaltime),ncol=ntotaltime)
        for (index in 1:ntotaltime){
          predicted[,index]=matrix(dynrModel@measurement@values.load[[data.plot$regime[index]]][observed,],nrow=length(observed))%*%matrix(res@eta_smooth_final[,index], ncol=1)
          if (length(dynrModel@measurement@values.int)!=0){
            predicted[,index]=matrix(predicted[,index]+dynrModel@measurement@values.int[[data.plot$regime[index]]][observed,],nrow=length(observed))
          }
          if (length(dynrModel@measurement@values.exo)!=0){
            exo.mat=t(as.matrix(dynrModel@data$covariates[index,paste0("covar",sapply(dynrModel@measurement@exo.names, function(x){which(dynrModel@data$covariate.names==x)}))]))
            predicted[,index]=matrix(predicted[,index]+dynrModel@measurement@values.exo[[data.plot$regime[index]]][observed,],nrow=length(observed))%*%exo.mat
          }
        }
        #end of regime specific measurement
      }
      
      names.measure.vars <- c(paste0(names.observed, ".predicted"), paste0(names.observed, ".observed"))
      data.plot <- cbind(data.plot, matrix(predicted,ncol=length(observed), byrow=TRUE, dimnames=list(NULL, paste0(names.observed, ".predicted"))))
      data.plot <- cbind(data.plot, data.dynr$observed[,observed])
      lines.var <- paste0(names.observed, ".predicted")
      points.var <- paste0(names.observed, ".observed")
      if (length(observed)<=8){
        fancy8 <- c("#edaac1","#a7dab0","#dfb3e2","#8bc5ed","#bcb8ec","#e6b296","#d2d39d","#82d9d6")
        manfillPalette=fancy8[9-(1:num_regime)]
        mancolorPalette=rep(fancy8[observed],2)
      }
    }#end of style 2
    
    names(data.plot)<-c(names.id.vars,names.measure.vars)
    
    data_long<- reshape2::melt(data.plot[data.plot$id%in%randid,],
                               id.vars=names.id.vars,
                               measure.vars=names.measure.vars,
                               value.name="value")
    data_long$variable<-factor(data_long$variable, levels=sort(names.measure.vars)) 
    #data_long$statenumber<-as.factor(sub("state","",data_long$variable))
    endtime <- NULL
    regime <- NULL
    partial.plot<-ggplot2::ggplot(data_long,ggplot2::aes(x=time, y=value, group=variable)) +
      ggplot2::geom_rect(ggplot2::aes(xmin=time, xmax=endtime, ymin=-Inf, ymax=Inf, fill=regime), alpha=.15) +
      ggplot2::geom_line(data=data_long[data_long$variable%in%lines.var,], size=1, ggplot2::aes(color=variable, linetype=variable)) +
      ggplot2::geom_point(data=data_long[data_long$variable%in%points.var,], size=3, ggplot2::aes(color=variable, shape=variable)) +
      ggplot2::scale_linetype_manual(labels=sort(names.measure.vars), values=line.values[order(names.measure.vars)])+
      ggplot2::scale_shape_manual(labels=sort(names.measure.vars), values=shape.values[order(names.measure.vars)])+
      #geom_text(size=1, ggplot2::aes(label=statenumber,color=variable))+
      ggplot2::facet_wrap(~id)+ 
      ggplot2::labs(title = title,y=ylab,ggtitle="")+
      ggplot2::theme(plot.title=ggplot2::element_text(lineheight=.8, face="bold", hjust = 0.5, vjust=0.01), legend.position="bottom",legend.text=ggplot2::element_text(lineheight=0.8, face="bold"),...)
    
  }
  
  if(is.bw){
    bw.plot<-partial.plot+
      ggplot2::scale_fill_grey(start = 0.1, end = 0.9)+
      ggplot2::scale_color_grey(labels=sort(names.measure.vars), start = 0.2, end = 0.7)
    return(bw.plot)
  }else{
    default.plot<-partial.plot+
      ggplot2::scale_colour_brewer(labels=sort(names.measure.vars), palette=colorPalette)+
      ggplot2::scale_fill_brewer(palette=fillPalette)		  
    if(!missing(mancolorPalette)){
      manual.plot<-partial.plot+ggplot2::scale_colour_manual(labels=sort(names.measure.vars), values=mancolorPalette[order(names.measure.vars)])
    }
    if(!missing(manfillPalette)){
      if(!exists("manual.plot")){manual.plot<-partial.plot}
      manual.plot<-manual.plot+ggplot2::scale_fill_manual(values=manfillPalette)
    }
    if (exists("manual.plot")){
      return(manual.plot) #print(manual.plot)
    }else{
      return(default.plot)#print(default.plot)
    }
  }
}

##' @rdname dynr.ggplot
##' 
##' @param object The same as res. The dynr object returned by dynr.cook().
autoplot.dynrCook <- function(object, dynrModel, style = 1,
                        numSubjDemo=2, idtoPlot=c(),
                        names.state, 
                        names.observed,
                        names.regime,
                        shape.values,
                        title, 
                        ylab, 
                        is.bw=FALSE, 
                        colorPalette="Set2", 
                        fillPalette="Set2", 
                        mancolorPalette, manfillPalette, ...){
	dynr.ggplot(object, dynrModel=dynrModel, style = style,
                        numSubjDemo=numSubjDemo, idtoPlot=idtoPlot,
                        names.state=names.state, 
                        names.observed=names.observed,
                        names.regime=names.regime,
                        shape.values=shape.values,
                        title=title, 
                        ylab=ylab, 
                        is.bw=is.bw, 
                        colorPalette=colorPalette, 
                        fillPalette=fillPalette, 
                        mancolorPalette=mancolorPalette, manfillPalette=manfillPalette, ...)
}

##' The ggplot of the outliers estimates.
##' 
##' @param object A dynrTaste object.
##' @param numSubjDemo The number of subjects, who have 
##' largest joint chi-square statistic, to be selected  for plotting.
##' @param idtoPlot Values of the ID variable to plot.
##' @param names.state (optional) The names of the states to be plotted, which should be a subset of the state.names slot of the measurement slot of dynrModel. If NULL, the t statistic plots for all state variables will be included. 
##' @param names.observed (optional) The names of the observed variables to be plotted, which should be a subset of the obs.names slot of the measurement slot of dynrModel. If NULL, the t statistic plots for all observed variables will be included.
##' @param ... Place holder for other arguments. Please do not use.
##' 
##' @return a list of ggplot objects for each ID. 
##' The plots of chi-square statistics (joint and independent),
##' and the plots of t statistic for \code{names.state} and \code{names.observed} will be included.
##' Users can modify the ggplot objects using ggplot grammar.
##' If a \code{filename} is provided, a pdf of plots will be saved additionally.
autoplot.dynrTaste <- function(object, 
                               numSubjDemo=2, idtoPlot=NULL,
                               names.state=NULL, names.observed=NULL, ...) {
  if ( !inherits(object, "dynrTaste") ) {
    stop("dynrTaste object is required.") 
  }
  stopifnot(numSubjDemo >= 1)
  numSubjDemo <- as.integer(numSubjDemo)
  
  conf_level <- object$conf.level
  tstart <- object$tstart
  chi_jnt <- object$chi.jnt
  id <- object$id
  id_unq <- unique(id)
  id_n <- length(id_unq)
  lat_name <- rownames(object$t.inn)
  obs_name <- rownames(object$t.add)
  lat_n <- length(lat_name)
  obs_n <- length(obs_name)
  
  if( is.null(names.state) ) {
    lat_toplot <- lat_name
  } else {
    if ( !all(names.state %in% lat_name) ) {
      stop("'names.state' should be a subset of the latent variables.")
    }
    lat_toplot <- names.state
  }
  
  if( is.null(names.observed) ) {
    obs_toplot <- obs_name
  } else {
    if ( !all(names.observed %in% obs_name) ) {
      stop("'names.observed' should be a subset of the observed variables.")
    }
    obs_toplot <- names.observed
  }
  
  if ( is.null(idtoPlot) ) {
    chi_max <- vector("numeric", id_n)
    for(i in 1:id_n){
      begT <- tstart[i] + 1
      endT <- tstart[i+1]
      chi_max[i] <- max( chi_jnt[begT:endT] )
    }
    id_toplot <- id_unq[ order(chi_max, decreasing=TRUE) ][
      1:numSubjDemo ]
  } else {
    id_toplot <- idtoPlot
    if ( !all(id_toplot %in% id_unq) ) {
      stop("Not all ID are in the data.")
    }
  }
  
  ### create ggplot objects
  ggadd_hline <- function(yintercept) {
    geom_hline(yintercept=yintercept, 
               colour="black", linetype="dashed", alpha=0.5)
  }
  ggadd_point <- function(data, aes_y) {
    geom_point(data=data, aes_string(x="time", y=aes_y),
               colour="red", size=2, alpha=0.5)
  }
  ggadd_theme <- function() {
    theme(plot.title=element_text(size = rel(0.85)),
          axis.title.x=element_blank(),
          #axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          plot.margin = margin(0, 0, 0, 0),
          axis.line=element_line(colour = "black"))
  }
  
  ## Joint chi-square
  id_toplot_n <- length(id_toplot)
  time <- object$time
  chi_jnt_shk <- NULL
  chi_df_jnt <- data.frame(id = id, time = time,
                           chi_jnt = object$chi.jnt,
                           chi_jnt_shk = object$chi.jnt.shock)
  # subset data for id to plot
  chi_df_jnt_plot <- subset(chi_df_jnt, is.element(id, id_toplot))
  # subset data only for shocks
  chi_df_jnt_shk <- subset(chi_df_jnt_plot, chi_jnt_shk)
  qchi_jnt <- qchisq(conf_level, lat_n + obs_n)
  
  plots_jnt <- lapply(id_toplot, function(id_i) {
    chi_df_jnt_shk_i <- subset(chi_df_jnt_shk, id==id_i)
    ggplot2::ggplot(subset(chi_df_jnt_plot, id==id_i), 
                    aes(x=time, y=chi_jnt)) +
      geom_line(alpha=0.5) +
      labs(title="") +
      ylab(expression(paste("Joint ", Chi^2))) +
      list(ggadd_point(chi_df_jnt_shk_i, "chi_jnt"), 
           ggadd_hline(qchi_jnt), ggadd_theme())
  })
  
  #### Independend chi-square for innovative ####
  chi_inn <- chi_inn_shk <- NULL
  chi_df_inn <- data.frame(id = id, time = time,
                           chi_inn = object$chi.inn,
                           chi_inn_shk = object$chi.inn.shock)
  chi_df_inn_plot <- subset(chi_df_inn, is.element(id, id_toplot))
  chi_df_inn_shk <- subset(chi_df_inn_plot, chi_inn_shk)
  qchi_inn <- qchisq(conf_level, lat_n)
  
  plots_inn <- lapply(id_toplot, function(id_i) {
    chi_df_inn_shk_i <- subset(chi_df_inn_shk, id==id_i)
    ggplot2::ggplot(subset(chi_df_inn_plot, id==id_i),
                    aes(x=time, y=chi_inn)) +
      geom_line(alpha=0.5) +
      labs(title="") +
      ylab(expression(paste("Indep.  ", Chi^2, " - innovative"))) +
      list(ggadd_point(chi_df_inn_shk_i, "chi_inn"),
           ggadd_hline(qchi_inn), ggadd_theme())
  })
  
  #### Independent chi-square for additive ####
  chi_add <- chi_add_shk <- NULL
  chi_df_add <- data.frame(id = id, time = time,
                           chi_add = object$chi.add,
                           chi_add_shk = object$chi.add.shock)
  chi_df_add_plot <- subset(chi_df_add, is.element(id, id_toplot))
  chi_df_add_shk <- subset(chi_df_add_plot, chi_add_shk)
  
  plots_add <- lapply(id_toplot, function(id_i) {
    chi_df_add_shk_i <- subset(chi_df_add_shk, id==id_i)
    ggplot2::ggplot(subset(chi_df_add_plot, id==id_i),
                    aes(x=time, y=chi_add)) +
      geom_line(alpha=0.5) +
      labs(title="") +
      ylab(expression(paste("Indep.  ", Chi^2, " - additive"))) +
      list(ggadd_point(chi_df_add_shk_i, "chi_add"),
           ggadd_hline(qchisq(conf_level, lat_n)), ggadd_theme())
  })
  
  #### t-statistic for innovative ####
  t_inn <- object$t.inn
  t_inn_name <- row.names(t_inn)
  row.names(t_inn) <- paste0(t_inn_name, "_t")
  t_df_inn <- data.frame(id = id, time = time,
                         t(t_inn), t(object$t.inn.shock))
  t_df_inn_plot <- subset(t_df_inn, is.element(id, id_toplot))
  
  plots_t_inn <- lapply(id_toplot, function(id_i) {
    t_df_inn_plot_i <- subset(t_df_inn_plot, id==id_i)
    lapply(lat_toplot, function(lat_i) {
      t_df_inn_shk_i <- subset(t_df_inn_plot_i,
                               t_df_inn_plot_i[[lat_i]])
      qt_i <- qt(conf_level, nrow(t_df_inn_plot_i) - obs_n)
      ggplot2::ggplot(t_df_inn_plot_i,
                      aes_string(x="time", y=paste0(lat_i, "_t"))) +
        geom_line(alpha=0.5) +
        labs(title="") +
        ylab(paste("t_[", lat_i, "]")) +
        list(ggadd_point(t_df_inn_shk_i, paste0(lat_i, "_t")),
             ggadd_hline(qt_i), ggadd_hline(-qt_i), ggadd_theme())
    })
  })
  
  #### t-statistic for additive ####
  t_add <- object$t.add
  t_add_name <- row.names(t_add)
  row.names(t_add) <- paste0(t_add_name, "_t")
  t_df_add <- data.frame(id = id, time = time,
                         t(t_add), t(object$t.add.shock))
  t_df_add_plot <- subset(t_df_add, is.element(id, id_toplot))
  
  plots_t_add <- lapply(id_toplot, function(id_i) {
    t_df_add_plot_i <- subset(t_df_add_plot, id==id_i)
    t_add_df_i <- nrow(t_df_add_plot_i)
    lapply(obs_toplot, function(obs_i) {
      t_df_add_shk_i <- subset(t_df_add_plot_i,
                               t_df_add_plot_i[[obs_i]])
      qt_i <- qt(conf_level, t_add_df_i - obs_n)
      ggplot2::ggplot(t_df_add_plot_i,
                      aes_string(x="time", y=paste0(obs_i, "_t"))) +
        geom_line(alpha=0.5) +
        labs(title="") +
        ylab(paste("t_[", obs_i, "]")) +
        list(ggadd_point(t_df_add_shk_i, paste0(obs_i, "_t")),
             ggadd_hline(qt_i), ggadd_hline(-qt_i), ggadd_theme())
    })
  })
  
  gg_objects <- mapply(function(jnt_i, inn_i, add_i, 
                                t_inn_i, t_add_i, id_toplot_i) {
    gg_i <- c(list(jnt_i, inn_i, add_i), t_inn_i, t_add_i)
    names(gg_i) <- c("chi_jnt", "chi_inn", "chi_add", 
                     paste0("t_[", lat_toplot, "]"), 
                     paste0("t_[", obs_toplot, "]"))
    gg_i
  },
  plots_jnt, plots_inn, plots_add, plots_t_inn, plots_t_add, id_toplot,
  SIMPLIFY=FALSE)
  
  names(gg_objects) <- id_toplot
  
  # if ( !is.null(filename) && is.character(filename) ) {
  #   gg_plots <- mapply(function(plots_i, id_toplot_i) {
  #     plots_arr <- ggpubr::ggarrange(
  #       plotlist=plots_i, nrow=length(plots_i), ncol=1, align="v")
  #     ggpubr::annotate_figure(
  #       plots_arr,
  #       top=ggpubr::text_grob(paste0("ID: ", id_toplot_i), 
  #                             color="black", face="bold", size=14),
  #       bottom=ggpubr::text_grob("time", color="black", size=11))
  #   }, 
  #   gg_objects, id_toplot,
  #   SIMPLIFY=FALSE)
  #   
  #   id_v <- chi_df_jnt_plot[["id"]]
  #   id_time_length <- vector("numeric", id_toplot_n)
  #   # time points for each id
  #   for (i in 1:id_toplot_n) {
  #     id_time_length[i] <- sum(is.element(id_v, id_toplot[i]))
  #   }
  #   # width: time/6, height: n of gg object * 2
  #   ggpubr::ggexport(gg_plots, filename=filename, 
  #                    width=mean(id_time_length)/6, 
  #                    height=length(gg_objects[[1]])*2)
  # }
  invisible(gg_objects)
}
