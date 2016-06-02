# 
# dynr model CLASS
# 
##' The dynrModel Class
##' 
##' @aliases
##' $,dynrModel-method
##' $<-,dynrModel-method
##' 
##' @details
##' This is an internal class structure.  You should not use it directly.
##' Use \code{\link{dynr.model}} instead.
setClass(Class =  "dynrModel",
         representation = representation(
           dynamics =  "dynrDynamics",
           measurement = "dynrMeasurement",
           noise = "dynrNoise",
           initial = "dynrInitial",
           regimes= "dynrRegimes",
           transform="dynrTrans",
           data="list",
           num_regime="integer",
           dim_latent_var="integer",
           infile="character",
           outfile="character",
           verbose="logical",
           compileLib="logical",
           xstart="vector",
           ub="vector",
           lb="vector",
           options="list",
           param.names="character"
         ),
         prototype = prototype(
           num_regime=as.integer(1),
           verbose=TRUE,
           compileLib=TRUE,
           options=default.model.options
         )
)

setMethod("initialize", "dynrModel",
          function(.Object, x){
            for(i in names(x)){
              slot(.Object, name=i, check = TRUE) <- x[[i]]
            }
            return(.Object)
          }
)

setMethod("$", "dynrModel",
          function(x, name){slot(x, name)}
)

setReplaceMethod("$", "dynrModel",
	function(x, name, value){
		if(name %in% c('xstart', 'ub', 'lb')){
			# Check that the length is okay
			slot(object=x, name=name, check = TRUE) <- x$transform$inv.tfun.full(value)
		} else if(name %in% c('dynamics', 'measurement', 'noise', 'initial', 'regimes', 'transform')) {
			slot(object=x, name=name, check = TRUE) <- value
		} else {
			stop(paste0("You can't set the ", name, " slot of a dynrModel.", "  You're not allowed to touch me there."))
		}
		return(x)
	}
)

implode <- function(..., sep='') {
  paste(..., collapse=sep)
}


vecRegime <- function(object){
  numRegimes <- nrow(object$values)
  covariates <- object$covariates
  numCovariates <- length(covariates)
  
  Prlist <- list()
  for (j in 1:numRegimes){
  values <- object$values[j,which(object$values[j,]!=0)]
  params <- object$params[j,which(object$values[j,]!=0)]
  colIndex <- which(params==params)    #which(object$values[j,]!="0")
  colIndexSet <- ceiling((colIndex)/(numCovariates+1))
  for (q in unique(colIndexSet)){
  colIndex2 <- which(colIndexSet==q)
  a <- diag(matrix(outer(matrix(values[colIndex2],ncol=numCovariates+1), 
           matrix(c(1,object$covariates),ncol=1),
           FUN=paste,sep="*"),ncol=numCovariates+1))
  namesLO = paste0("&\\frac{Pr(p",j,q,")}{1-Pr(p",j,q,")}")
  a <- paste0(namesLO," = ", 
              implode(gsub("*1","",a,fixed=TRUE),sep=" + "))
  Prlist <- paste0(Prlist , a, "\\\\")
              }#End of loop through colIndex
  }#End of loop through regime
  return(Prlist)
}  

setMethod("printex", "dynrModel",
          function(object, ParameterAs, 
                   printDyn=TRUE, printMeas=TRUE, printInit=FALSE, printRS=FALSE,
                   outFile){
            model2<-PopBackModel(object, ParameterAs)
            
            inlist <- list(model2$dynamics, model2$measurement, model2$noise, model2$initial, model2$regimes)
            outlist <- lapply(inlist, printex, show = FALSE)
            
            if (printInit){
              #initial condition of latent variables
              if (length(outlist[[4]]$initial.state) > 1){              
                for (j in 1:length(outlist[[4]]$initial.state)){
                  space <-ifelse(j<length(outlist[[4]]$initial.state),",\\\\\n","")
                  a <- paste0("\\text{Regime ",j,":}&\\\\\n")  
                  initmu <- outlist[[4]]$initial.state[[j]]
                  initcov <- outlist[[4]]$initial.covariance[[j]]
                  initCond=paste0(paste0(a,"&",.xtableMatrix(matrix(paste0((model2$measurement)$state.names,"(0)"),ncol=1),F), "\\sim N \\Big(",
                                    initmu,", ",
                                    initcov,"\\Big)",space,"\n"))
                }#End of loop through j
              }else{
                a <- ""
                space <- ""
                initmu <- outlist[[4]]$initial.state
                initcov <- outlist[[4]]$initial.covariance
                initCond=paste0(a,.xtableMatrix(matrix(paste0((model2$measurement)$state.names,"(0)"),ncol=1),F), "\\sim N \\Big(",
                                initmu,", ",
                                initcov,"\\Big)",space,"\n")
              }#End of check if the initial condition specification is regime-specific
              initequ=paste0("\\begin{align*}\n",
                             initCond,
                             "\\end{align*}")
            }
            
            #Dynamic Model
            if (printDyn){
              #dynamic noise
              #TODO Check if this part is correct
              #Determine whether model is deterministic or stochastic
              pn <- c((model2$noise)$values.latent)
              isProcessNoise <- ifelse(length(pn[which(pn=="0")])==
                                         length((model2$measurement)$state.names)*length((model2$measurement)$state.names)
                                       ,0,1)
              lw <- length(outlist[[3]]$dynamic.noise)
              if (lw ==1){
                processNoise <- outlist[[3]]$dynamic.noise
              }
              
              # Modify state names on the LHS based on continuous- vs discrete-time models
              if ((model2$dynamics)$isContinuousTime){
                #Continuous-time dynamic model
                LHS <- paste0(.xtableMatrix(matrix(paste0("d",(model2$measurement)$state.names,"(t)"),ncol=1),F))
                state <- paste0(.xtableMatrix(matrix(paste0((model2$measurement)$state.names,"(t)"),ncol=1),F))
                
                if (isProcessNoise){
                  pnLab <- " + dw(t),"
                  pNoisePre <- paste0("&dw(t) \\sim N\\Big(",
                                      .xtableMatrix(matrix(rep(0,length((model2$measurement)$state.names)),ncol=1),F),
                                      ",")
                }else{#No Measurement Noise
                  pnLab <- ""
                  pNoisePre <- ""
                }
              }else{#Discrete-time dynamic model
                LHS <- .xtableMatrix(matrix(paste0((model2$measurement)$state.names,"(t+1)"),ncol=1),F)
                state <- .xtableMatrix(matrix(paste0((model2$measurement)$state.names,"(t)"),ncol=1),F)
                
                if (isProcessNoise){
                  pnLab <- " + w(t),"
                  pNoisePre <- paste0("&w(t) \\sim N\\Big(",
                                      .xtableMatrix(matrix(rep(0,length((model2$measurement)$state.names)),ncol=1),F),
                                      ",")
                }else{
                  pnLab <- ""
                  pNoisePre <- ""
                }
              }#End discrete-time dynamic model
              
              dynequ="\\begin{align*}\n"
              if (class(model2$dynamics) == 'dynrDynamicsFormula'){
                exp1 <- printex(model2$dynamics)
                for (j in 1:length((model2$dynamics)$formula)){
                  if(lw>1){processNoise <- outlist[[3]]$dynamic.noise[[j]]}
                  if (length((model2$dynamics)$formula)>1){
                    a <- paste0("\\text{Regime ",j,":}&\\\\\n")
                  }else{a <- ""}
                  neq <- length((model2$dynamics)$formula[[j]])
                  for (k in 1:neq){
                    space <- ifelse(((j*k)<neq*length((model2$dynamics)$formula))|isProcessNoise,",\\\\\n","")
                    if ((model2$dynamics)$isContinuousTime && isProcessNoise){
                      #TODO deal with noise properly
                      #Continuous-time dynamic model
                      pnLab <- paste0(" + dw",k,"(t)")
                    }else if(isProcessNoise){#Discrete-time dynamic model
                      pnLab <- paste0(" + w",k,"(t)")
                    }#End discrete-time dynamic model
                    dynequ <- paste0(dynequ,paste0(a,"&",exp1[[j]][k],pnLab,space))
                    a <- NULL
                  }#loop through eqs within regime j
                  if (isProcessNoise){
                    dynequ <- paste0(dynequ,paste0(pNoisePre,processNoise,"\\Big)",ifelse(j==length((model2$dynamics)$formula),"\n","\\\\\n")))}
                }#loop through regimes
                #state <- .xtableMatrix(matrix(paste0((model2$measurement)$state.names,"(t)"),ncol=1),F)
              }else{ #DynamicMatrix specification

                a <- NULL; space = ""
                for (j in 1:length(outlist[[1]]$dyn_tran)){
                  space<- ifelse((j<length(outlist[[1]]$dyn_tran))|isProcessNoise,"\\\\","")
                  if (length(outlist[[1]]$dyn_tran) > 1) {
                    if(lw>1){processNoise <- outlist[[3]]$dynamic.noise[[j]]}
                    a <- paste0("\\text{Regime ",j,":}&\\\\")
                  }#End of text edits required only for multiple-regime models
                  exo <- NULL
                  dint <- NULL
                  if (length((model2$dynamics)$covariates) > 0){
                    exo <- paste0("+",outlist[[1]]$dyn_exo[[j]],outlist[[1]]$dyn_exo.names)  
                  }#End covariate if
                  if (length((model2$dynamics)$values.int) > 0){
                    dint <- paste0(outlist[[1]]$dyn_int[[j]],"+")  
                  }#End int if
                  if ((model2$dynamics)$isContinuousTime){
                    RHSeqPre <- "("
                    RHSeqPost <- ")dt"
                  }else{
                    RHSeqPre <- ""
                    RHSeqPost <- ""
                  }
                  dynequ<-paste0(dynequ,a,"&",LHS,"=",dint,
                                               RHSeqPre,outlist[[1]]$dyn_tran[[j]],state,RHSeqPost,
                                               exo,pnLab,space)
                  if (isProcessNoise){
                    dynequ <- paste0(dynequ,
                                paste0(pNoisePre,processNoise,"\\Big)",ifelse(j==length(outlist[[1]]$dyn_tran),"\n","\\\\\n")))
                  }
                }#Loops through regimes
              }#End of formula vs. linear dynamic matrix specification
              dynequ <- paste0(dynequ,"\\end{align*}")
            }
            
            #Measurement model
            if (printMeas){
              #measurement noise
              #TODO Check if this part is correct
              #Determine whether model is deterministic or stochastic
              pn <- c((model2$noise)$values.observed)
              isMeasNoise <- ifelse(length(pn[which(pn=="0")])==
                                      length((model2$measurement)$obs.names)*length((model2$measurement)$obs.names)
                                    ,0,1)
              lmeas <- length(outlist[[3]]$measurement.noise)
              if (lmeas ==1){
                measNoise <- outlist[[3]]$measurement.noise
              }
              
              measequ="\\begin{align*}\n"
              obs <- .xtableMatrix(matrix(paste0((model2$measurement)$obs.names,"(t)"),ncol=1),F)
              
              for (j in 1:length((model2$measurement)$values.load)){
                space <-ifelse(j<length((model2$measurement)$values.load),"\\\\\n","")
                a <- NULL
                if (length((model2$measurement)$values.load) > 1) {
                  a <-  paste0("\\text{Regime ",j,":}&\\\\\n")
                  if (lmeas > 1) measNoise <- outlist[[3]]$measurement.noise[[j]]
                }
                exom <- NULL
                mint <- NULL
                if (length((model2$measurement)$exo.names) > 0){
                  exom <- paste0("+",outlist[[2]]$meas_exo[[j]],outlist[[2]]$meas_exo.names)  
                }
                if (length((model2$measurement)$values.int) > 0){
                  mint <- paste0(outlist[[2]]$meas_int[[j]]," + ")  
                }
                measequ <- paste0(measequ,paste0(a,"&",obs," = ",mint,
                                            outlist[[2]]$meas_loadings[[j]],
                                            state,exom))
                if (isMeasNoise){
                  measequ <- paste0(measequ,paste0("+ \\epsilon,",
                                              "\\\\&\\epsilon\\sim N\\Big(",
                                              .xtableMatrix(matrix(rep(0,length((model2$measurement)$obs.names)),ncol=1),F),
                                              ",",measNoise,"\\Big)"))
                }#end of (isMeasNoise check)
                measequ=paste0(measequ,space)
              }#end of loop through regimes
              measequ=paste0(measequ,"\\end{align*}")
            }
            
            #Regime-switching model
            if (printRS){
              #initial regime probabilities
              initProb <- outlist[[4]]$initial.probability
              outProb <- NULL
              if (length((model2$initial)$values.regimep)>1){
                #Only print initial regime probabilities if > 1 regime
                outProb <- paste0("&\\text{Initial regime probabilities = }",
                                  initProb,"\\\\\n")
              }
              #regime switch probability
              if (length((model2$initial)$values.regimep)>1){
                Prlist <- implode(vecRegime(model2$regimes),sep="&\\\\\n")
                #Only print initial RS probabilities if > 1 regime
                RSequ=paste0("\\begin{align*}\n",outProb,
                             Prlist,"\n\\end{align*}\n")
              }
            }
            
            #Print out the latex code
            outcode <- "\\documentclass[fleqn]{article}\n\\usepackage{amsmath}\n\\setlength{\\mathindent}{0pt}\n\n\\begin{document}\n"
            if (printMeas==TRUE){
              outcode <- paste0(outcode, "\nThe measurement model is given by:\n", measequ)
            }
            if (printDyn==TRUE){
              outcode <- paste0(outcode, "\n\nThe dynamic model is given by:\n", dynequ)
            }
            if (printInit==TRUE){
              outcode <- paste0(outcode, "\n\nThe initial condition of the dynamic model is given by:\n", initequ)
            }
            if (printRS==TRUE){
              outcode <- paste0(outcode, "\n\nThe regime-switching model is given by:\n", RSequ)
            }
            outcode <- paste0(outcode, "\\end{document}\n")
            if (missing(outFile)){
              cat(outcode)
            }else{
              cat(outcode,file=outFile)
            }
          })

# modeling is what happens to recipes.
# The alpha version of this file just takes a bunch of recipes and puts
#  them together into a C file of the user's naming.
#

.logisticCFunction <- "/**\n * This function takes a double and gives back a double\n * It computes the logistic function (i.e. the inverse of the logit link function)\n * @param x, the double value e.g. a normally distributed number\n * @return logistic(x), the double value e.g. a number between 0 and 1\n */\ndouble mathfunction_logistic(const double x){\n\tdouble value = 1.0/(1.0 + exp(-x));\n\treturn value;\n}\n"

.softmaxCFunction <- "/**\n * This function takes a gsl_vector and modifies its second argument (another gsl_vector)\n * It computes the softmax function (e.g. for multinomial logistic regression)\n * @param x, vector of double values e.g. a vector of normally distributed numbers\n * @param result, softmax(x), e.g. a vector of numbers between 0 and 1 that sum to 1\n */\nvoid mathfunction_softmax(const gsl_vector *x, gsl_vector *result){\n\t/* Elementwise exponentiation */\n\tsize_t index=0;\n\tfor(index=0; index < x->size; index++){\n\t\tgsl_vector_set(result, index, exp(gsl_vector_get(x, index)));\n\t}\n\t\n\t/* Sum for the scaling coeficient */\n\tdouble scale = 0.0;\n\tfor(index=0; index < x->size; index++){\n\t\tscale += gsl_vector_get(result, index);\n\t}\n\t\n\t/* Multiply all elements of result by 1/scale */\n\tgsl_blas_dscal(1/scale, result);\n}\n"

.cfunctions <- paste(.logisticCFunction, .softmaxCFunction, sep="\n")

##' Create a dynrModel object for parameter estimation (cooking dynr) using \code{\link{dynr.cook}}
##' 
##' @param dynamics a dynrDynamics object prepared with \code{\link{prep.formulaDynamics}} 
##' or \code{\link{prep.matrixDynamics}}
##' @param measurement a dynrMeasurement object prepared with \code{\link{prep.loadings}} 
##' or \code{\link{prep.measurement}}
##' @param noise a dynrNoise object prepared with \code{\link{prep.noise}}
##' @param initial a dynrInitial object prepared with \code{\link{prep.initial}}
##' @param data a dynrData object made with \code{\link{dynr.data}}
##' @param ... additional arguments specifying other dynrRecipe objects. Argument regimes is for 
##' a dynrRegimes object prepared with \code{\link{prep.regimes}} and argument transform is for 
##' a dynrTrans object prepared with \code{\link{prep.tfun}}.
##' @param infile a character string of the name of the input C script of model functions to be compiled 
##' for parameter estimation. The default is a temporary file name.
##' @param outfile a character string of the name of the output C script of model functions to be compiled 
##' for parameter estimation.
##' 
##' @examples
##' #rsmod <- dynr.model(dynamics=recDyn, measurement=recMeas, noise=recNoise, initial=recIni, regimes=recReg, data=dd, outfile="RSLinearDiscrete.c")
##' #For a full demo example, see:
##' #demo(RSLinearDiscrete , package="dynr")
dynr.model <- function(dynamics, measurement, noise, initial, data, ..., infile=tempfile(), outfile){
  #check the order of the names 
  if (class(dynamics) == "dynrDynamicsFormula"){
    states.dyn <- lapply(dynamics@formula, function(list){sapply(list, function(fml){as.character(as.list(fml)[[2]])})})
    if (all(sapply(states.dyn, function(x, y){all(x==y)}, y=states.dyn[[1]]))){
      states.dyn=states.dyn[[1]]
    }else{
      stop("Formulas should be specified in the same order for different regimes.")
    }
    if (!all(measurement@state.names == states.dyn)){
      stop("The state.names slot of the 'dynrMeasurement' object hould match the order of the dynamic formulas specified.")
    }
  }
  if (!all(measurement@obs.names == data$observed.names)){
    stop("The obs.names slot of the 'dynrMeasurement' object should match the observed argument passed to the dynr.data() function.")
  }
  # gather inputs
  inputs <- list(dynamics=dynamics, measurement=measurement, noise=noise, initial=initial, ...)

  # Figure out what the unique parameters are
  all.params <- unlist(sapply(inputs, slot, name='paramnames'))
  unique.params <- extractParams(all.params)
  unique.numbers <- c() #allow for model with no free parameters
  if(length(unique.params) > 0){unique.numbers <- 1L:(length(unique.params))}
  
  # Create the map between parameter values, the user-specified parameter names, and the automatically-produced parameter numbers (param.data$param.number)
  param.data <- data.frame(param.number=unique.numbers, param.name=unique.params,stringsAsFactors=FALSE)
  
  param.data$ldl.latent<-param.data$param.name%in%extractParams(inputs$noise$params.latent)
  param.data$ldl.observed<-param.data$param.name%in%extractParams(inputs$noise$params.observed)
  param.data$ldl.inicov<-param.data$param.name%in%extractParams(inputs$initial$params.inicov)
  
  #TODO write a way to extract param.data from a model object (grabs from recipes within model)
  
  #TODO write a way to assign param.data to a model object (assigns to recipes within model)
  # populate transform slots
  if(any(sapply(inputs, class) %in% 'dynrTrans')){
    inputs$transform<-createRfun(inputs$transform, param.data, 
                                 params.observed=inputs$noise$params.observed, params.latent=inputs$noise$params.latent, params.inicov=inputs$initial$params.inicov,
                                 values.observed=inputs$noise$values.observed.inv.ldl, values.latent=inputs$noise$values.latent.inv.ldl, values.inicov=inputs$initial$values.inicov.inv.ldl)
    #at this step, the paramnum slot of transform gets populated, which is needed for paramName2Number
  }else{
    inputs$transform <- createRfun(prep.tfun(), param.data, 
                                 params.observed=inputs$noise$params.observed, params.latent=inputs$noise$params.latent, params.inicov=inputs$initial$params.inicov,
                                 values.observed=inputs$noise$values.observed.inv.ldl, values.latent=inputs$noise$values.latent.inv.ldl, values.inicov=inputs$initial$values.inicov.inv.ldl)
  }
  # paramName2Number on each recipe (this changes are the params* matrices to contain parameter numbers instead of names
  inputs <- sapply(inputs, paramName2Number, names=param.data$param.name)
  
  # writeCcode on each recipe
  inputs <- sapply(inputs, writeCcode, data$covariate.names) 
  all.values <- unlist(sapply(inputs, slot, name='startval'))
  unique.values <- extractValues(all.values, all.params)
  
  if(length(inputs$transform$formula.inv)>0){
    unique.values <- inputs$transform$inv.tfun(unique.values)
  }
  param.data$param.value=unique.values
  
  #initiate a dynrModel object
  obj.dynrModel <- new("dynrModel", c(list(data=data, infile=infile, outfile=outfile, param.names=as.character(param.data$param.name)), inputs))
  obj.dynrModel@dim_latent_var <- dim(inputs$measurement$values.load[[1]])[2] #numbber of columns of the factor loadings
  
  obj.dynrModel@xstart <- param.data$param.value
  obj.dynrModel@ub <- rep(9999,length(obj.dynrModel@xstart))
  obj.dynrModel@lb <- rep(9999,length(obj.dynrModel@xstart))
  if(any(sapply(inputs, class) %in% 'dynrRegimes')){
    obj.dynrModel@num_regime<-dim(inputs$regimes$values)[1]
  }
  #write out the C script
  cparts <- unlist(sapply(inputs, slot, name='c.string'))
  includes <- "#include <math.h>\n#include <gsl/gsl_matrix.h>\n#include <gsl/gsl_blas.h>\n"
  body <- paste(cparts, collapse="\n\n")
  if( length(grep("void function_regime_switch", body)) == 0 ){ # if regime-switching function isn't provided, fill in 1 regime model
    body <- paste(body, writeCcode(prep.regimes())$c.string, sep="\n\n")
  }
  body<-gsub("NUM_PARAM",length(obj.dynrModel@xstart),body)
#   if( length(grep("void function_transform", body)) == 0 ){ # if transformation function isn't provided, fill in identity transformation
#     body <- paste(body, writeCcode(prep.tfun())$c.string, sep="\n\n")
#   }
  glom <- paste(includes, .cfunctions, body, sep="\n\n")
  if (obj.dynrModel@dynamics@isContinuousTime){
    glom <- paste(glom, dP_dt, sep="\n\n")
  }
  cat(glom, file=obj.dynrModel@outfile)
  
  return(obj.dynrModel)
  #modify the object slot, including starting values, etc.
}


