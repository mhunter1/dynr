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
##' @param printDyn A logical value indicating whether or not to plot the formulas for the dynamic model
##' @param printMeas A logical value indicating whether or not to plot the formulas for the measurement model
##' @param textsize numeric. Font size used in the plot.
##' @param ... Further named arguments
##' 
##' @details
##' This is a wrapper around \code{\link{dynr.ggplot}}.  A great benefit of it is that is shows the equations for the plot.
setMethod("plot", "dynrCook",
          function(x, dynrModel, style = 1, printDyn=TRUE, printMeas=TRUE, textsize=4, ...) {
  #The first panel is the ggplot
  p1 <- dynr.ggplot(x, dynrModel, style, numSubjDemo=1, ...)

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
                     plot.title = ggplot2::element_text(size = 12, colour = "black", face = "bold"))
  }
  
  #The third panel plots the model formulae
  p3 <- plotFormula(dynrModel, ParameterAs=sprintf("%.2f", x@transformed.parameters), printDyn=printDyn, printMeas=printMeas, textsize=textsize)

  #Organize the panels using the multiplot function
  if (dynrModel$num_regime > 1){
    multiplot(p1, p2, p3, cols = 1, layout=matrix(c(1,1,2,3), nrow=2, byrow=TRUE))
  }else{
    multiplot(p1, p3, layout=matrix(c(1,1,2,2), nrow=2, byrow=TRUE))
  }

  })

plotdf <- function(vec_tex){
  dataframe=data.frame(text=sapply(paste0("$",vec_tex,"$"),function(x){as.character(latex2exp::TeX(x))}))
  dataframe$x<-0
  return(dataframe)
}

##' Plot the formula from a model
##' 
##' @param dynrModel The model object to plot.
##' @param ParameterAs The parameter values or names to plot. The underscores in parameter names are saved for 
##' use of subscripts.  Greek letters can be specified as corresponding LaTeX symbols without backslashes (e.g., "lambda") 
##' and printed as greek letters.
##' @param printDyn A logical value indicating whether or not to plot the formulas for the dynamic model.
##' @param printMeas A logical value indicating whether or not to plot the formulas for the measurement model
##' @param textsize The text size use in the plot.
##' 
##' @details
##' This function typesets a set of formulas that represent the model.
plotFormula <- function(dynrModel, ParameterAs, printDyn=TRUE, printMeas=TRUE, textsize=4){
  
  dynrModel <- PopBackModel(dynrModel, LaTeXnames(ParameterAs, latex = FALSE))

  #Dynamic model
  if (printDyn){  
    if (class(dynrModel$dynamics) == "dynrDynamicsMatrix"){
      state.names <- (dynrModel$measurement)$state.names
      exo.names <- (dynrModel$dynamics)$covariates
      dynrModel@dynamics@values.dyn <- lapply((dynrModel$dynamics)$values.dyn, preProcessNames,state.names,state.names)
      dynrModel@dynamics@values.exo <- lapply((dynrModel$dynamics)$values.exo, preProcessNames,state.names,exo.names)
      dynrModel@dynamics@values.int <- lapply((dynrModel$dynamics)$values.int, preProcessNames,state.names)
    }

    dyn.df <- data.frame(text="'Dynamic Model'",x=0)
    nRegime=ifelse(class(dynrModel@dynamics)=="dynrDynamicsFormula",
                   length(dynrModel@dynamics@formula),
                   length(dynrModel@dynamics@values.dyn))
    dyn_tex=printex(dynrModel@dynamics,AsMatrix=FALSE)
    nEq <- length(dyn_tex[[1]])
    noise_tex <- paste0(" + ", ifelse(dynrModel@dynamics@isContinuousTime, "d", ""), "w_{",1:nEq, "}(t)")
    for (i in 1:nRegime){
      if (nRegime>1){
        dyn.df <- rbind(dyn.df,data.frame(text=paste0("'Regime ",i,":'"),x=0))
      }
      if (length(dynrModel@noise@values.latent)==1){
        values.latent.mat <- dynrModel@noise@values.latent[[1]] 
      }else if (length(dynrModel@noise@values.latent)==nRegime){
        values.latent.mat <- dynrModel@noise@values.latent[[i]]
      }else{stop("The number of regimes implied by the dynamic noise structure does not match the number of regimes in the dynamic model.")}
      noise_tex[which(diag(values.latent.mat)==0)] <- ""
      dyn.df <- rbind(dyn.df,plotdf(LaTeXnames(paste0(dyn_tex[[i]], noise_tex))))
    }
  }
  
  #Measurement model
  if (printMeas){
    meas.df<-data.frame(text="'Measurement Model'",x=0)
    nRegime=length(dynrModel@measurement@values.load)
    meas_tex=printex(dynrModel@measurement,AsMatrix=FALSE)
    nEq <- length(meas_tex[[1]])
    noise_tex <- paste0(" + ", "epsilon_{", 1:nEq, "}")
    for (i in 1:nRegime){
      if (nRegime>1){
        meas.df<-rbind(meas.df,data.frame(text=paste0("'Regime ",i,":'"),x=0))
      }
      if (length(dynrModel@noise@values.observed)==1){
        values.observed.mat <- dynrModel@noise@values.observed[[1]] 
      }else if (length(dynrModel@noise@values.observed)==nRegime){
        values.observed.mat <- dynrModel@noise@values.observed[[i]]
      }else{stop("The number of regimes implied by the measurement noise structure does not match the number of regimes in the measurement model.")}
      noise_tex[which(diag(values.observed.mat)==0)] <- ""
      meas.df<-rbind(meas.df,plotdf(LaTeXnames(paste0(meas_tex[[i]], noise_tex))))
    }
  }
  
  if (printDyn&printMeas){
    plot.df=rbind(dyn.df,meas.df)
  }else if (printDyn){
    plot.df=dyn.df
  }else if (printMeas){
    plot.df=meas.df
  }else{
    stop("There needs to be something to plot.")
  }
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
##' @param res The dynr object returned by dynr.cook().
##' @param dynrModel The model object to plot.
##' @param style The style of the plot. If style is 1 (default), user-selected smoothed state variables are plotted. If style is 2, user-selected observed-versus-predicted values are plotted.
##' @param numSubjDemo The number of subjects to be randomly selected for plotting.
##' @param idtoPlot Values of the ID variable to plot
##' @param names.state (optional) The names of the states to be plotted, which should be a subset of dynrModel@measurement@state.names.
##' @param names.observed (optional) The names of the observed variables to be plotted, which should be a subset of dynrModel@measurement@obs.names.
##' @param names.regime (optional) The names of the regimes to be plotted, which can be missing.
##' @param shape.values (optional) A vector of values that correspond to the shapes of the points, which can be missing. See the R documentation on pch for details on possible shapes.
##' @param title (optional) A title of the plot.
##' @param ylab (optional) The label of the y axis.
##' @param is.bw Is plot in black and white? The default is FALSE.
##' @param colorPalette A palette function for lines and dots that when called with a single integer argument (the number of levels in the scale) returns the values that they should take.
##' @param fillPalette A palette function for blocks that when called with a single integer argument (the number of levels in the scale) returns the values that they should take.
##' @param mancolorPalette (optional) A color palette for manually scaling the colors of lines and plots.
##' @param manfillPalette (optional) A color palette for manually scaling the colors of filled blocks.
##' @param ... A list of element name, element pairings that modify the existing ggplot2::ggplot ggplot2::theme. Consult the ggplot2::theme() function in the R package ggplot2::ggplot.
dynr.ggplot <- function(res, dynrModel, style,
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
  observed=which(names.observed%in%dynrModel@measurement@obs.names)
  if ((length(observed)>8)&(missing(mancolorPalette))&(style==2)){stop("You provided too many variables than the default color palette can handle.\nPlease consider specify the mancolorPalette argument.")}
  if(missing(shape.values)){
    if (style==1) shape.values=48+states
    if (style==2) shape.values=15+observed
  }
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
      names.measure.vars <- c(paste0(names.observed, ".observed"), paste0(names.observed, ".predicted"))
      data.plot <- cbind(data.plot, matrix(predicted[observed,],ncol=length(observed),byrow=TRUE, dimnames=list(NULL, paste0(names.observed, ".predicted"))))
      colnames(data.dynr$observed[,observed])<-paste0(names.observed, ".observed")
      data.plot <- cbind(data.plot, data.dynr$observed[,observed])
      lines.var <- paste0(names.observed, ".predicted")
      points.var <- paste0(names.observed, ".observed")
      if (length(observed)<8){
        fancy8 <- c("#edaac1","#a7dab0","#dfb3e2","#8bc5ed","#bcb8ec","#e6b296","#d2d39d","#82d9d6")
        mancolorPalette=rep(fancy8[observed],each=2)
      }
    }#end of style 2
    names(data.plot) <- c(names.id.vars, names.measure.vars)
    
    data_long <- reshape2::melt(data.plot[data.plot$id%in%randid,],
                                id.vars=names.id.vars,
                                measure.vars=names.measure.vars,
                                value.name="value")
    value <- NULL
    variable <- NULL
    partial.plot<-ggplot2::ggplot(data_long, ggplot2::aes(x=time, y=value, colour = variable)) + 
      ggplot2::geom_line(data=data_long[data_long$variable%in%lines.var,], size=1) +
      ggplot2::geom_point(data=data_long[data_long$variable%in%points.var,], size=4, ggplot2::aes(shape=variable)) +
      ggplot2::scale_shape_manual(values=shape.values)+
      ggplot2::facet_wrap(~id)+
      ggplot2::labs(title = title, y=ylab, ggtitle="")+
      ggplot2::theme(plot.title=ggplot2::element_text(lineheight=.8, face="bold"), legend.position="bottom",legend.text=ggplot2::element_text(lineheight=0.8, face="bold"),...)
    
  }else{#more than two regimes
    if(missing(names.regime)){names.regime = 1:num_regime}
    if(missing(names.state)){names.state=dynrModel@measurement@state.names}
    states=which(names.state%in%dynrModel@measurement@state.names)
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
      dynrModel=dynr:::PopBackModel(dynrModel, LaTeXnames(res@transformed.parameters, latex = FALSE))
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
        #end of one-regime meas
      }else if (num_regime_meas==num_regime){
        ntotaltime=ncol(res@eta_smooth_final)
        predicted=matrix(rep(0,length(observed)*ntotaltime),ncol=ntotaltime)
        for (index in 1:ntotaltime){
          predicted[,index]=dynrModel@measurement@values.load[[data.plot$regime[index]]]%*%matrix(res@eta_smooth_final[,index], ncol=1)
          if (length(dynrModel@measurement@values.int)!=0){
            predicted[,index]=predicted[,index]+dynrModel@measurement@values.int[[data.plot$regime[index]]]
          }
          if (length(dynrModel@measurement@values.exo)!=0){
            exo.mat=t(as.matrix(dynrModel@data$covariates[index,paste0("covar",sapply(dynrModel@measurement@exo.names, function(x){which(dynrModel@data$covariate.names==x)}))]))
            predicted[,index]=predicted[,index]+dynrModel@measurement@values.exo[[data.plot$regime[index]]]%*%exo.mat
          }
        }
        #end of regime specific measurement
      }
      
      names.measure.vars <- c(paste0(names.observed, ".observed"), paste0(names.observed, ".predicted"))
      data.plot <- cbind(data.plot, matrix(predicted[observed,],ncol=length(observed),byrow=TRUE, dimnames=list(NULL, paste0(names.observed, ".predicted"))))
      colnames(data.dynr$observed[,observed])<-paste0(names.observed, ".observed")
      data.plot <- cbind(data.plot, data.dynr$observed[,observed])
      lines.var <- paste0(names.observed, ".predicted")
      points.var <- paste0(names.observed, ".observed")
      if (length(observed)<=8){
        fancy8 <- c("#edaac1","#a7dab0","#dfb3e2","#8bc5ed","#bcb8ec","#e6b296","#d2d39d","#82d9d6")
        manfillPalette=fancy8[9-(1:num_regime)]
        mancolorPalette=rep(fancy8[observed],each=2)
      }
    }#end of style 2
    
    names(data.plot)<-c(names.id.vars,names.measure.vars)
    
    data_long<- reshape2::melt(data.plot[data.plot$id%in%randid,],
                               id.vars=names.id.vars,
                               measure.vars=names.measure.vars,
                               value.name="value")
    #data_long$statenumber<-as.factor(sub("state","",data_long$variable))
    endtime <- NULL
    regime <- NULL
    partial.plot<-ggplot2::ggplot(data_long,ggplot2::aes(x=time, y=value, group=variable)) +
      ggplot2::geom_rect(ggplot2::aes(xmin=time, xmax=endtime, ymin=-Inf, ymax=Inf, fill=regime), alpha=.15) +
      ggplot2::geom_line(data=data_long[data_long$variable%in%lines.var,], size=1, ggplot2::aes(color=variable)) +
      ggplot2::geom_point(data=data_long[data_long$variable%in%points.var,], size=4, ggplot2::aes(color=variable,shape=variable)) +
      ggplot2::scale_shape_manual(values=shape.values)+
      #geom_text(size=1, ggplot2::aes(label=statenumber,color=variable))+
      ggplot2::facet_wrap(~id)+ 
      ggplot2::labs(title = title,y=ylab,ggtitle="")+
      ggplot2::theme(plot.title=ggplot2::element_text(lineheight=.8, face="bold"), legend.position="bottom",legend.text=ggplot2::element_text(lineheight=0.8, face="bold"),...)
    
  }
  
  if(is.bw){
    bw.plot<-partial.plot+
      ggplot2::scale_fill_grey(start = 0.1, end = 0.9)+
      ggplot2::scale_color_grey(start = 0.2, end = 0.7)
    print(bw.plot)
  }else{
    default.plot<-partial.plot+
      ggplot2::scale_colour_brewer(palette=colorPalette)+
      ggplot2::scale_fill_brewer(palette=fillPalette)		  
    if(!missing(mancolorPalette)){
      manual.plot<-partial.plot+ggplot2::scale_colour_manual(values=mancolorPalette)
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




