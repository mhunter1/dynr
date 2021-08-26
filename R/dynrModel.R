# 
# dynr model CLASS
# 
##' The dynrModel Class
##' 
##' @aliases
##' $,dynrModel-method
##' $<-,dynrModel-method
##' print,dynrModel-method
##' show,dynrModel-method
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
           #random = "dynrRandom",
           regimes= "dynrRegimes",
           transform="dynrTrans",
           data="list",
           num_regime="integer",
           dim_latent_var="integer",
           outfile="character",
           verbose="logical",
           compileLib="logical",
           xstart="vector",
           ub="vector",
           lb="vector",
           options="list",
           param.names="character",
           random.params.inicov = "matrix",
           random.values.inicov = "matrix",
		   dLambdparLamb2="matrix",
		   dLambdparLamb="matrix",
		   dmudparMu="matrix",
		   dmudparMu2="matrix",
		   dSigmaede="matrix",
		   dSigmaede2="matrix",
		   dSigmabdb="matrix",
		   dSigmabdb2="matrix",
		   Sigmab="matrix",
		   known.vars = "list",
		   freeIC="logical"
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


setMethod("print", "dynrModel", printRecipeOrModel)
setMethod("show", "dynrModel", function(object){printRecipeOrModel(object)})


##' Extract the free parameter names of a dynrModel object
##' 
##' @param x The dynrModel object from which the free parameter names are desired
setMethod("names", "dynrModel",
	function(x) {
		pnames <- x$param.names
		output <- c(pnames)
		output <- gsub("(\\w+\\W+.*)", "'\\1'", output)
		return(output)
	}
)

.DollarNames.dynrModel <- function(x, pattern){
	if(missing(pattern)){
		pattern <- ''
	}
	output <- slotNames(x)
	output <- gsub("(\\w+\\W+.*)", "'\\1'", output)
	return(grep(pattern, output, value=TRUE))
}

setReplaceMethod("$", "dynrModel",
	function(x, name, value){
		if(name %in% c('xstart', 'ub', 'lb')){
			# Check that the length is okay
			if(length(value) != length(x$param.names)){
				stop(paste("I'm going over my borders.", "You gave me", length(value), "things,",
					"but I need", length(x$param.names),
					"(the number of free parameters)."))
			}
			if(is.null(names(value))){
				names(value) <- x$param.names
			}
			lookup <- match(names(value), x$param.names)
			lookup <- union(na.omit(lookup), 1L:length(value))
			value[lookup] <- value
			names(value) <- x$param.names
			slot(object=x, name=name, check = TRUE) <- x$transform$inv.tfun.full(value) #suppressWarnings(expr)
			# Check that parameters are within the bounds and the bounds are in order
			if(any(na.omit(x$ub < x$lb))){
				warning("Some of the lower bounds are above the upper bounds. Bad, user!")
			}
			if(any(na.omit(x$ub < x$xstart))){
				offenders <- names(which(na.omit(x$ub < x$xstart)))
				offenders <- paste(offenders, collapse=', ')
				msg <- paste0("I spy with my little eye an upper bound that is smaller than the starting value.\n",
					"Offending parameter(s) named: ", offenders)
				warning(msg)
			}
			if(any(na.omit(x$lb > x$xstart))){
				offenders <- names(which(na.omit(x$lb > x$xstart)))
				offenders <- paste(offenders, collapse=', ')
				msg <- paste0("I spy with my little eye a lower bound that is larger than the starting value.\n",
					"Offending parameter(s) named: ", offenders)
				warning(msg)
			}
		} else if(name %in% c('dynamics', 'measurement', 'noise', 'initial', 'regimes', 'transform')) {
			slot(object=x, name=name, check = TRUE) <- value
		} else {
			stop(paste0("You can't set the ", name, " slot of a dynrModel.", "  You're not allowed to touch me there."))
		}
		return(x)
	}
)

##' Extract the number of observations for a dynrModel object
##' 
##' @param object An unfitted model object
##' @param ... Further named arguments. Ignored.
##' 
##' @details
##' We return the total number of rows of data, adding up the number of time points for each person. For some purposes, you may want the mean number of observations per person or the number of people instead.  These are not currently supported via \code{nobs}.
##' 
##' @return
##' A single number. The total number of observations across all IDs.
##' 
##' @examples
##' # Let rawModel be the output from dynr.model
##' #nobs(rawModel)
nobs.dynrModel <- function(object, ...){
	dim(object$data$observed)[1]
}


##' @rdname coef.dynrCook
coef.dynrModel <- function(object, ...){
	object$transform$tfun(object$xstart)
}

##' @rdname coef.dynrCook
##' 
##' @param value values for setting
`coef<-` <- function(object, value){
	UseMethod("coef<-")
}

##' @rdname coef.dynrCook
`coef<-.dynrModel` <- function(object, value){
	if(length(coef(object)) != length(value)){
		stop(paste0("Number of model coeficients (", length(coef(object)), ") does not match number assigned (", length(value), ")."))
	}
	object <- PopBackModel(object, value)
	return(object)
}


implode <- function(..., sep='') {
  paste(..., collapse=sep)
}


vecRegime <- function(object){
	objValues <- object$values
	objParams <- object$params
	numRegimes <- nrow(object$values)
	covariates <- object$covariates
	numCovariates <- length(covariates)
	deviation <- object$deviation
	refRow <- object$refRow
	
	if(deviation){
		if(nrow(objValues)!=0 && nrow(objParams)!=0){
			interceptSel <- seq(1, ncol(objValues), by=numCovariates+1)
			intercept.values <- matrix(objValues[refRow, interceptSel], numRegimes, 1)
			intercept.params <- matrix(objParams[refRow, interceptSel], numRegimes, 1)
			objValues[refRow, interceptSel] <- 0
			objParams[refRow, interceptSel] <- 0
		}
	} else {
		intercept.values <- matrix(0, numRegimes, 1)
		intercept.params <- matrix(0, numRegimes, 1)
	}
	
	colBeginSeq <- seq(1, ncol(objValues), by=numCovariates+1)
	colEndSeq <- colBeginSeq + numCovariates
	Prlist <- list()
	for (j in 1:numRegimes){
		for(k in 1:numRegimes){
			colSel <- colBeginSeq[k]:colEndSeq[k]
			#namesLO = paste0("&\\frac{Pr(p",j,k,")}{1-Pr(p",j,k,")}")
			namesLO = paste0("&Log Odds(p", j, k, ")")
			mat1 <- matrix(objValues[j, colSel], ncol=numCovariates+1)
			mat2 <- matrix(c(1, covariates), ncol=1)
			# drop zeros before multiplication
			mat2 <- mat2[mat1 !=0 ]
			mat1 <- mat1[mat1 !=0 ]
			a <- paste(mat1, mat2, sep="*")
			a <- gsub("*1", "", a, fixed=TRUE)
			b <- intercept.values[k, 1]
			b <- b[ b != 0]
			a <- implode(c(b, a), sep=" + ")
			if(nchar(a) > 0){
				a <- paste0(namesLO, " = ", a)
				Prlist <- paste0(Prlist , a, "\\\\")
			}
		}#End of loop through colIndex
	}#End of loop through regime
	return(Prlist)
}

LaTeXnames<-function(names, decimal = 2, latex = TRUE){
  if (latex == TRUE){
    if (class(names)=="character"){
      names<-gsub("([[:alnum:]]+[\\]*)_([[:alnum:]]+)","\\1_{\\2}",names)
      
      greek <- c("alpha", "nu", "beta", "xi", "Xi", "gamma", "Gamma", "delta", 
                 "Delta", "pi", "Pi", "epsilon", "varepsilon", "rho", "varrho", 
                 "zeta", "sigma", "Sigma", "eta", "tau", "theta", "vartheta", 
                 "Theta", "upsilon", "Upsilon", "iota", "phi", "varphi", "Phi", 
                 "kappa", "chi", "lambda", "Lambda", "psi", "Psi", "mu", "omega", "Omega")
      greekpat <- paste0("(\\<", paste0(greek, collapse = "_*\\>|\\<"), "_*\\>)")
      names<-gsub(greekpat,"\\\\\\1",names)
      
      names<-gsub("\\\\_", "_", names)
    }
  }
  
  if (class(names)=="numeric"){
    names<-sprintf(paste0("%.",decimal,"f"), names)
  }
  
  return(names)
}

setMethod("printex", "dynrModel",
          function(object, ParameterAs, 
                   printDyn=TRUE, printMeas=TRUE, printInit=FALSE, printRS=FALSE,
                   outFile){
            model2<-PopBackModel(object, LaTeXnames(ParameterAs, latex = FALSE))
            
            inlist <- list(model2$dynamics, model2$measurement, model2$noise, model2$initial, model2$regimes)
            outlist <- lapply(inlist, printex, show = FALSE)
            
            if (printInit){
              #initial condition of latent variables
              if (length(outlist[[4]]$initial.state) > 1){
                initCond=character(0)
                for (j in 1:length(outlist[[4]]$initial.state)){
                  space <-ifelse(j<length(outlist[[4]]$initial.state),",\\\\\n","")
                  a <- paste0("\\text{Regime ",j,":}&\\\\\n")  
                  initmu <- outlist[[4]]$initial.state[[j]]
                  initcov <- outlist[[4]]$initial.covariance[[j]]
                  initCond=paste0(initCond, paste0(a,"&",.xtableMatrix(matrix(paste0((model2$measurement)$state.names,"(0)"),ncol=1),F), "\\sim N \\Big(",
                                                   initmu,", ",
                                                   initcov,"\\Big)",space))
                }#End of loop through j
              }else{
                a <- ""
                space <- ""
                initmu <- outlist[[4]]$initial.state
                initcov <- outlist[[4]]$initial.covariance
                initCond=paste0(a,.xtableMatrix(matrix(paste0((model2$measurement)$state.names,"(0)"),ncol=1),F), "\\sim N \\Big(",
                                initmu,", ",
                                initcov,"\\Big)",space)
              }#End of check if the initial condition specification is regime-specific
              initequ=paste0("\\begin{align*}\n",
                             initCond,
                             "\n\\end{align*}")
            }
            
            #Dynamic Model
            if (printDyn){
              #dynamic noise
              lw <- length(outlist[[3]]$dynamic.noise)
              
              # Modify state names on the LHS based on continuous- vs discrete-time models
              if ((model2$dynamics)$isContinuousTime){
                #Continuous-time dynamic model
                LHS <- paste0(.xtableMatrix(matrix(paste0("d",(model2$measurement)$state.names,"(t)"),ncol=1),F))
                state <- paste0(.xtableMatrix(matrix(paste0((model2$measurement)$state.names,"(t)"),ncol=1),F))
              }else{#Discrete-time dynamic model
                LHS <- .xtableMatrix(matrix(paste0((model2$measurement)$state.names,"(t+1)"),ncol=1),F)
                state <- .xtableMatrix(matrix(paste0((model2$measurement)$state.names,"(t)"),ncol=1),F)
              }#End discrete-time dynamic model
              
              dynequ="\\begin{align*}\n"
              if (class(model2$dynamics) == 'dynrDynamicsFormula'){
                exp1 <- printex(model2$dynamics)
                for (j in 1:length((model2$dynamics)$formula)){
                  #Determine whether model in regime is deterministic or stochastic
                  if (lw == 1){
                    processNoise <- outlist[[3]]$dynamic.noise[[1]]
                    pn <- c(model2@noise@values.latent[[1]])
                  }else if (lw == length((model2$dynamics)$formula)){
                    processNoise <- outlist[[3]]$dynamic.noise[[j]]
                    pn <- c(model2@noise@values.latent[[j]])
                  }else{
                    stop("The number of regimes implied by the dynamic noise structure does not match the number of regimes in the dynamic model.")
                  }
                  
                  isProcessNoise <- ifelse(length(pn[which(pn=="0")])==
                                             length((model2$measurement)$state.names)*length((model2$measurement)$state.names)
                                           ,0,1)
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
                    }else {#End discrete-time dynamic model
                      pnLab <- ""}
                    dynequ <- paste0(dynequ,paste0(a,"&",exp1[[j]][k],pnLab,space))
                    a <- NULL
                  }#loop through eqs within regime j
                  if (isProcessNoise){
                    pNoisePre <- paste0(ifelse((model2$dynamics)$isContinuousTime, "&dw(t) \\sim N\\Big(", "&w(t) \\sim N\\Big("),
                                        .xtableMatrix(matrix(rep(0,length((model2$measurement)$state.names)),ncol=1),F),
                                        ",")
                    dynequ <- paste0(dynequ,paste0(pNoisePre,processNoise,"\\Big)",ifelse(j==length((model2$dynamics)$formula),"\n","\\\\\n")))}
                }#loop through regimes
                #state <- .xtableMatrix(matrix(paste0((model2$measurement)$state.names,"(t)"),ncol=1),F)
              }else{ #DynamicMatrix specification
                
                a <- NULL; space = ""
                for (j in 1:length(outlist[[1]]$dyn_tran)){
                  #Determine whether model in regime is deterministic or stochastic
                  if (lw == 1){
                    processNoise <- outlist[[3]]$dynamic.noise[[1]]
                    pn <- c(model2@noise@values.latent[[1]])
                  }else if (lw == length((model2$dynamics)$formula)){
                    processNoise <- outlist[[3]]$dynamic.noise[[j]]
                    pn <- c(model2@noise@values.latent[[j]])
                  }else{
                    stop("The number of regimes implied by the dynamic noise structure does not match the number of regimes in the dynamic model.")
                  }
                  isProcessNoise <- ifelse(length(pn[which(pn=="0")])==
                                             length((model2$measurement)$state.names)*length((model2$measurement)$state.names)
                                           ,0,1)
                  space<- ifelse((j<length(outlist[[1]]$dyn_tran))|isProcessNoise,"\\\\","")
                  if (length(outlist[[1]]$dyn_tran) > 1) {
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
                  
                  if (isProcessNoise){
                    pnLab <- ifelse((model2$dynamics)$isContinuousTime, " + dw(t),", " + w(t),")
                  }else{#No Noise
                    pnLab <- ""
                  }
                  dynequ<-paste0(dynequ,a,"&",LHS,"=",dint,
                                 RHSeqPre,outlist[[1]]$dyn_tran[[j]],state,RHSeqPost,
                                 exo,pnLab,space)
                  if (isProcessNoise){
                    pNoisePre <- paste0(ifelse((model2$dynamics)$isContinuousTime, "&dw(t) \\sim N\\Big(", "&w(t) \\sim N\\Big("),
                                        .xtableMatrix(matrix(rep(0,length((model2$measurement)$state.names)),ncol=1),F),
                                        ",")
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
              lmeas <- length(outlist[[3]]$measurement.noise)
              
              measequ="\\begin{align*}\n"
              obs <- .xtableMatrix(matrix(paste0((model2$measurement)$obs.names,"(t)"),ncol=1),F)
              
              for (j in 1:length((model2$measurement)$values.load)){
                #Determine whether model in regime is deterministic or stochastic
                if (lmeas == 1){
                  measNoise <- outlist[[3]]$measurement.noise[[1]]
                  pn <- c(model2@noise@values.observed[[1]])
                }else if (lmeas == length((model2$measurement)$values.load)){
                  measNoise <- outlist[[3]]$measurement.noise[[j]]
                  pn <- c(model2@noise@values.observed[[j]])
                }else{
                  stop("The number of regimes implied by the measurement noise structure does not match the number of regimes in the measurement model.")
                }
                isMeasNoise <- ifelse(length(pn[which(pn=="0")])==
                                        length((model2$measurement)$obs.names)*length((model2$measurement)$obs.names)
                                      ,0,1)
                
                space <-ifelse(j<length((model2$measurement)$values.load),"\\\\\n","")
                a <- NULL
                if (length((model2$measurement)$values.load) > 1) {
                  a <-  paste0("\\text{Regime ",j,":}&\\\\\n")
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
                  measequ <- paste0(measequ,paste0("+ epsilon,",
                                                   "\\\\&epsilon\\sim N\\Big(",
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
              if (model2$num_regime > 1){
                #Only print initial regime probabilities if > 1 regime
                outProb <- paste0("&\\text{Initial regime log odds = }",
                                  initProb,"\\\\\n")
              }
              #regime switch probability
              if (model2$num_regime > 1){
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
            outcode <- LaTeXnames(outcode)
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
##' @param outfile a character string of the name of the output C script of model functions to be compiled 
##' for parameter estimation. The default is the name for a potential temporary file returned by tempfile().
##' 
##' @details
##' A \code{dynrModel} is a collection of recipes.  The recipes are constructed with the functions \code{\link{prep.measurement}}, \code{\link{prep.noise}}, \code{\link{prep.formulaDynamics}}, \code{\link{prep.matrixDynamics}}, \code{\link{prep.initial}}, and in the case of regime-switching models \code{\link{prep.regimes}}.  Additionally, data must be prepared with \code{\link{dynr.data}} and added to the model.
##' 
##' Several \emph{named} arguments can be passed into the \code{...} section of the function.  These include
##' \itemize{
##' 	\item Argument \code{regimes} is for a dynrRegimes object prepared with \code{\link{prep.regimes}}
##' 	\item Argument \code{transform} is for a dynrTrans object prepared with \code{\link{prep.tfun}}.
##' 	\item Argument \code{options} a list of options. Check the NLopt website \url{https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#stopping-criteria}
##' for details. Available options for use with a dynrModel object 
##' include xtol_rel, stopval, ftol_rel, ftol_abs, maxeval, and maxtime, 
##' all of which control the termination conditions for parameter optimization. The examples below show a case where options were set.
##' }
##' 
##' There are several available methods for \code{dynrModel} objects.
##' \itemize{
##' 	\item The dollar sign ($) can be used to both get objects out of a model and to set pieces of the model.
##' 	\item \code{names} returns the names of the free parameters in a model.
##' 	\item \code{\link{printex}} prints LaTeX expressions for the equations that compose a model. The output can then be readily typeset for inclusion in presentations and papers.
##' 	\item \code{nobs} gives the total number of observations (e.g. all times across all people)
##' 	\item \code{coef} gives the free parameter starting values.  Free parameters can also be assigned with \code{coef(model) <- aNamedVectorOfCoefficients}
##' }
##' 
##' @examples
##' #rsmod <- dynr.model(dynamics=recDyn, measurement=recMeas, noise=recNoise, 
##' #    initial=recIni, regimes=recReg, data=dd, outfile="RSLinearDiscrete.c")
##'
##' #Set relative tolerance on function value via 'options':
##' #rsmod <- dynr.model(dynamics=recDyn, measurement=recMeas, noise=recNoise, 
##' #    initial=recIni, regimes=recReg, data=dd, outfile="RSLinearDiscrete.c",
##' #    options=list(ftol_rel=as.numeric(1e-6)))
##' 
##' #For a full demo example, see:
##' #demo(RSLinearDiscrete , package="dynr")
dynr.model <- function(dynamics, measurement, noise, initial, data, ..., outfile = tempfile()){
  
  if(!class(dynamics) %in% c("dynrDynamicsFormula","dynrDynamicsMatrix")){
    stop("Check to see that dynamics argument is of the correct class. Hint: it should be either 'dynrDynamicsFormula' or 'dynrDynamicsMatrix'.")
  }

  if(!class(noise) %in% "dynrNoise"){
    stop("Check to see that noise argument is of the correct class. Hint: it should be 'dynrNoise'.")
  }
  
  # Set ground truth for number/name of states and observations
  nameLatentVars <- measurement$state.names
  numLatentVars <- length(nameLatentVars)
  nameObsVars <- measurement$obs.names
  numObsVars <- length(nameObsVars)
  # numRegimes is defined and checked later in this function
  
  # check the order of the names 
  if (class(dynamics) == "dynrDynamicsFormula"){
    saem <- dynamics$saem

    states.dyn <- lapply(dynamics@formula, function(list){sapply(list, function(fml){as.character(as.list(fml)[[2]])})})
    if(any(sapply(states.dyn, length) != numLatentVars)){
      statePaste <- paste0('(', paste(sapply(states.dyn, length), collapse=', '), ')')
      stop(paste0("Found ", statePaste, " latent states in dynamics formula, but expected (", numLatentVars, ") latent states from measurement model."))
    }
    if(!all(states.dyn[[1]] %in% nameLatentVars)){
      stop(paste0("Latent state names in dynamics (", paste(states.dyn[[1]], collapse=", "), ") do not match those of measurement(", paste(nameLatentVars, collapse=", "), ")."))
    }
    #if (all(sapply(states.dyn, function(x, y){all(x==y)}, y=states.dyn[[1]]))){
    if (!all(sapply(states.dyn, function(x, y){all(x==y)}, y=nameLatentVars))){
      #states.dyn=states.dyn[[1]]#}else{
      stop("The 'state.names' slot of the 'dynrMeasurement' object should match the order \nof the dynamic formulas specified. \nSame order should hold even if you have multiple regimes.")
    }
    
    #Discrepancies in number of regimes in dynamic formula caught below via impliedRegimes function.
  }
  # if class == "dynrDynamicsMatrix" then check that the matrix dynamics is numLatentVars*numLatentVars
  # since prep.matrixDynamics already checks for 1. whether list elements of values.dyn or params.dyn are of the same dimension
  # 2. whether values.dyn and params.dyn are of the same matrix dimension
  # so it should be enough to just check one matrix
  if(class(dynamics) == "dynrDynamicsMatrix"){
    state.dimension <- dim(dynamics@values.dyn[[1]])
    if(state.dimension[1]!=numLatentVars|state.dimension[2]!=numLatentVars){
      stop("The matrix dimensions in prep.matrixDynamics should match the number of latent states in 'dynrMeasurement'.")
    }
  }
  
  # if class == "dynrNoise" then check that 1. the measurement noise is numObsVars by numObsVars
  # 2. the dynamic noise is numLatentVars by numLatentVars
  # Note that prep.noise already checks for 1. whether values.latent, params.latent, values.observed and params.observed are symmetric
  # 2. whether values.latent and params.latent are of the same matrix dimension
  # 3. whether values.observed and params.observed are of the same matrix dimension
  
  if(class(noise) == "dynrNoise"){
    observed.dimension <- dim(noise@values.observed[[1]]) 
    latent.dimension <- dim(noise@values.latent[[1]])
    if(observed.dimension[1]!=numObsVars){
      stop("The dimension of the measurement noise convariance matrix in prep.noise should match the number of observed variables in 'dynrMeasurement'.")
    }
    if(latent.dimension[1]!=numLatentVars){
      stop("The dimension of the dynamic noise covariance matrix in prep.noise should match the number of latent states in 'dynrMeasurement'.")
    }
  }
  
  # if class == "dynrInitial" then 
  ## 1. Check that ini.state is ne by 1
  ## 2. Check that ini.cov is ne by ne
  ## 3. Check that regimep is nr by 1 <- Q: Do separately after numregimes are set?
  if(class(initial) == "dynrInitial"){
    inistate.dim <- dim(initial@values.inistate[[1]])  
    inicov.dim <- dim(initial@values.inicov[[1]])
    if(inistate.dim[1]!=numLatentVars){
      stop("The number of the latent states in 'prep.initial' should match the number of latent states in 'dynrMeasurement'.")
    }
    if(inicov.dim[1]!=numLatentVars){
      msg <- paste0("The dimension of the initial covariance matrix for latent states in prep.noise should correspond to the number of latent states in 'dynrMeasurement'.", "\n", "Ideally: ", numLatentVars, " by ", numLatentVars, "\n", "Current: ", inicov.dim[1], " by ", inicov.dim[2])
      stop(msg, call. = F)
    }
    # iniregime.dim <- dim(initial@values.regimep)
    # if(iniregime.dim[1]!=numRegimes){
    #   stop("The number of the regimes in 'prep.initial' should match the number of regimes in 'dynrRegimes'.")
    # }
    # if(iniregime.dim[2]!= 1){
    #     stop("The initial regime probabilities should be a column vector.")
    # }      
    
  }  
  
  
  if (!all(measurement@obs.names == data$observed.names)){
    stop("The obs.names slot of the 'dynrMeasurement' object should match the 'observed' argument passed to the dynr.data() function.")
  }
  # Check all the covariates
  # Check if any submodel piece has any covariates but data do not
  if(.hasSlot(dynamics, 'covariates')){ dynCovar <- dynamics$covariates} else {dynCovar <- character(0)}
  anyCovariateRecipes <- length(c(measurement$exo.names, initial$covariates, dynCovar)) > 0
  if(anyCovariateRecipes && is.null(data$covariate.names)){
    stop("I found some covariates in your recipes, but not in your data.")
  }
  if (!is.null(data$covariate.names)){
    if(!all(measurement$exo.names %in% data$covariate.names)){
      stop("The 'exo.names' slot of the 'dynrMeasurement' object should match the 'covariates' argument passed to the dynr.data() function.\nA pox on your house if fair Romeo had not found this.")
    }
    if(!all(initial$covariates %in% data$covariate.names)){
      stop("The 'covariates' slot of the 'dynrInitial' object should match the 'covariates' argument passed to the dynr.data() function.\nA pox on your house if fair Romeo had not found this.")
    }
    if(class(dynamics) == "dynrDynamicsMatrix" && !all(dynamics$covariates %in% data$covariate.names)){
      stop("The 'covariates' slot of the 'dynrDynamicsMatrix' object should match the 'covariates' argument passed to the dynr.data() function.\nA pox on your house if fair Romeo had not found this.")
    }
    # Note: formula dynamics covariates are inferred from the formula and swapped in based on the data covariates
  }
  
  #if SAEM is TRUE but noise@values.latent are not zero matrix, warning
  if(saem == TRUE && any(noise@values.latent[[1]] != 0)){
    warning('Currently you can use the SAEM with ordinary differential equations (i.e., null matrix for the "values.latent" in dynrNoise)')
  }
  
  # check and modify the data
  ## For discrete-time models, the time points needs to be equally spaced. 
  if (!dynamics$isContinuousTime){
    time.split = split(data$time, as.factor(data$id))
    time.check = sapply(time.split, function(x) {
      difference = diff(x)
      return(c(spacing = sum(difference%%min(difference)) > 1e-6, #can be a very small positive number
               full = sum(diff(difference)) > 1e-6))
    })
    if(any(time.check["spacing",])){
      stop("Please check the data. The time points are irregularly spaced even with missingness inserted.")
    }else if (any(time.check["full",])){
      if ("covariates" %in% names(data)){
        names(data$covariates) <- data$covariate.names
        names(data$observed) <- data$observed.names
        data.dataframe <- data.frame(id = data$id, time = data$time, data$observed, data$covariates)
        
        data.new.dataframe <- plyr::ddply(data.dataframe, "id", function(df){
          new = data.frame(id = unique(df$id), time = seq(df$time[1], df$time[length(df$time)], by = min(diff(df$time))))
          out = merge(new, df, all.x = TRUE)
        })
        data <- dynr.data(data.new.dataframe, observed = data$observed.names, covariates = data$covariate.names)
      }else{
        names(data$observed) <- data$observed.names
        data.dataframe <- data.frame(id = data$id, time = data$time, data$observed)
        
        data.new.dataframe <- plyr::ddply(data.dataframe, "id", function(df){
          new = data.frame(id = unique(df$id), time = seq(df$time[1], df$time[length(df$time)], by = min(diff(df$time))))
          out = merge(new, df, all.x = TRUE)
        })
        
        data <- dynr.data(data.new.dataframe, observed = data$observed.names)
      }
    }
  }
  
  # gather inputs
  inputs <- list(dynamics=dynamics, measurement=measurement, noise=noise, initial=initial, ...)

  # beginning of new version to actually process the 'options' argument correctly	
  #  # gather inputs	
  #  extraArg <- list(...)	
  #  extraNames <- match.arg(names(extraArg), c('options', 'regimes', 'transform'), several.ok=TRUE)	
  #  inputs <- list(dynamics=dynamics, measurement=measurement, noise=noise, initial=initial, ...)	
  #  if('armadillo' %in% names(inputs)){inputs$armadillo <- NULL}
  
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
  if(saem == FALSE || .hasSlot(inputs$dynamics, 'theta.formula') == FALSE){
	  if(any(sapply(inputs, class) %in% 'dynrTrans')){
		inputs$transform<-createRfun(inputs$transform, param.data, 
									 params.observed=inputs$noise$params.observed, params.latent=inputs$noise$params.latent, params.inicov=inputs$initial$params.inicov,
									 values.observed=inputs$noise$values.observed.inv.ldl, values.latent=inputs$noise$values.latent.inv.ldl, values.inicov=inputs$initial$values.inicov.inv.ldl,
									 values.observed.orig=inputs$noise$values.observed, values.latent.orig=inputs$noise$values.latent, values.inicov.orig=inputs$initial$values.inicov)
		#at this step, the paramnum slot of transform gets populated, which is needed for paramName2Number
	  }else{
		inputs$transform <- createRfun(prep.tfun(), param.data, 
									 params.observed=inputs$noise$params.observed, params.latent=inputs$noise$params.latent, params.inicov=inputs$initial$params.inicov,
									 values.observed=inputs$noise$values.observed.inv.ldl, values.latent=inputs$noise$values.latent.inv.ldl, values.inicov=inputs$initial$values.inicov.inv.ldl,
									 values.observed.orig=inputs$noise$values.observed, values.latent.orig=inputs$noise$values.latent, values.inicov.orig=inputs$initial$values.inicov)
	  }
  }
  # paramName2Number on each recipe (this changes are the params* matrices to contain parameter numbers instead of names
  inputs <- sapply(inputs, paramName2Number, names=param.data$param.name)
  

  #  ------- The following lines obtain the necessary components of SAEM processs ----------------------------
  # Examine variables in formula: if there exist variables that are not in all.params and state.names, issue an error
  #browser()
  legal.variables <- c(all.params, inputs$measurement@state.names, data$covariate.names, inputs$dynamics@random.names, inputs$dynamics@theta.names, 'exp')
  variables <- extractVariablesfromFormula(unlist(dynamics@formula))
  if(any(variables %in% data$observed.names)){
    error.variables <- variables[variables %in% data$observed.names == TRUE]
    stop(paste0("Observed variables, namely, those declared as ``observed'' in dynr.data are not supposed to appear in the dynamic equations. Check to see if they can or should be linked to some state variables: ", paste(error.variables, collapse=", ")))
  }
  
  if(!all(variables %in% legal.variables)){

    error.variables <- variables[variables %in% legal.variables  == FALSE]
    stop(paste0("In formula, there are variables that are not specified in startvals: ", paste(error.variables, collapse=", ")))
  }
  
  #browser()
  # Examine variables in theta.formula: if there exist variables that are not in all.params and state.names, issue an error
  if(.hasSlot(inputs$dynamics, 'theta.formula') && length(inputs$dynamics@theta.formula) > 0){
    variables <- extractVariablesfromFormula(inputs$dynamics@theta.formula)
	if(any(variables %in% data$observed.names)){
      error.variables <- variables[variables %in% data$observed.names == TRUE]
      stop(paste0("Observed variables, namely, those declared as ``observed'' in dynr.data are not supposed to appear in the dynamic equations. Check to see if they can or should be linked to some state variables: ", paste(error.variables, collapse=", ")))
    }
    if(!all(variables %in% legal.variables)){
      error.variables <- variables[variables %in% legal.variables  == FALSE]
      stop(paste0("In theta.formula, there are variables that are not specified in startvals: ", paste(error.variables, collapse=", ")))
    }
  }
  
  if(saem==TRUE){
    # examine whether it is freeIC or fixed IC case
	inistate.names <- unique(as.vector(inputs$initial@params.inistate[[1]]))
	inistate.names <- inistate.names[!inistate.names %in% c(0, 'fixed')]
	inicov.names <- unique(as.vector(inputs$initial@params.inicov[[1]]))
	inicov.names <- inicov.names[!inicov.names %in% c(0, 'fixed')]
	if(length(inistate.names) == 0 && length(inicov.names) > 0){
		warning('The params.inistate in prep.initial are fixed, but params.inicov are free')
	}
	if (length(inistate.names) > 0 || length(inicov.names) > 0){
		freeIC = TRUE
	}
	else{
		freeIC = FALSE
	}
	
	#browser()
	#print("freeIC")
	#print(freeIC)


    #num.x <- length(inputs$dynamics@formula[[1]])
	
	num.x <- length(inputs$measurement$state.names)
	num.theta <- length(inputs$dynamics@theta.formula)
	if(freeIC){
		random.params.inicov = matrix(0L, 
								nrow = num.x + num.theta, 
								ncol = num.x + num.theta)
		random.values.inicov = matrix(0L, 
								nrow = num.x + num.theta, 
								ncol = num.x + num.theta)
	} else {
		random.params.inicov = matrix(0L, 
								nrow = num.theta, 
								ncol = num.theta)
		random.values.inicov = matrix(0L, 
								nrow = num.theta, 
								ncol = num.theta)
	}
	

	# setup InfDS.Sigmab
	# sigmab.names: unique variables in random.params.inicov that needs to be estimated
	sigmab.names <- all.params[unique(as.vector(inputs$initial$params.inicov[[1]]))]
	sigmab.names <- unique(c(as.vector(inputs$dynamics$random.params.inicov), sigmab.names))
	sigmab.names <- sigmab.names[!sigmab.names %in% c("fixed", "0")]
	
	if(length(sigmab.names) > 0){
		#variance/covariance of random.names
		random.params.inicov[1:num.theta,1:num.theta] = inputs$dynamics$random.params.inicov
		random.values.inicov[1:num.theta,1:num.theta] = inputs$dynamics$random.values.inicov
		
		#browser()
		if(freeIC){
			#variance/covariance of states
			random.params.inicov[(num.theta+1):(num.x+num.theta),(num.theta+1):(num.x+num.theta)] = inputs$initial$params.inicov[[1]]
		} 
		
		# LDL transformation
		ret <- symbolicLDLDecomposition(random.params.inicov, random.values.inicov)
		#Solve for numerical values of par0-par9 (the unconstrained parameters)
		# given starting values for the random effect covariance matrix in model@random.params.inicov
		known.vars <- solveStartLDL(ret$ldl, random.values.inicov)
		
		dSigmabdb <- differentiateMatrixOfVariable(ret$ldl, ret$pars)
		dSigmabdb2 <- differentiateMatrixOfVariable(dSigmabdb, ret$pars)
		dSigmabdb2 <- t(dSigmabdb2)
		
		#browser()
		
		#substitute the values in known.vars in
		#comment out [the values should be substituted in dynr.cook]
		#dSigmabdb<-matrix(sapply(dSigmabdb, function(x){eval(x, known.vars)}), nrow=nrow(dSigmabdb), ncol=ncol(dSigmabdb))
		#dSigmabdb2<-matrix(sapply(dSigmabdb2, function(x){eval(x, known.vars)}), nrow=nrow(dSigmabdb2), ncol=ncol(dSigmabdb2))
		#dSigmabdb2 <- t(dSigmabdb2)
		
	}
	
	#browser()
	mu.names <- c()
	if(length(inputs$measurement$params.int) > 0){
	  mu.names <- all.params[unique(as.vector(inputs$measurement$params.int[[1]]))]
	  mu.names <- mu.names[!mu.names %in% c("fixed", "0")]
	  if(length(mu.names) > 0){
	    #dmudparMu <- unList(differentiateMatrixOfVariable(all.params[inputs$measurement$params.int[[1]]]))
	    dmudparMu <- unList(differentiateMatrixOfVariable(matrix(sapply(inputs$measurement$params.int[[1]], function(x, all.params){if(x>0) all.params[x] else x}, all.params), nrow=nrow(inputs$measurement$params.int[[1]])),all.params[inputs$measurement$params.int[[1]]]))
	    dmudparMu2 <- unList(differentiateMatrixOfVariable(dmudparMu, all.params[inputs$measurement$params.int[[1]]]))
        dmudparMu2 <- t(dmudparMu2)
	  }
      else{
	    #when Nmu = 0 (all elements in mu is fixed)
	    dmudparMu <-matrix(0L, nrow=0, ncol=nrow(inputs$measurement$params.int[[1]]))
	    dmudparMu2 <-matrix(0L, nrow=0, ncol=nrow(inputs$measurement$params.int[[1]]))
	  }
	}
	else{
	  # when measurement@prarm.int is absent
      dmudparMu <-matrix(0L, nrow=0, ncol=0)
	  dmudparMu2 <-matrix(0L, nrow=0, ncol=0)
	}
	
	lamdba.names <- c()
	if(length(inputs$measurement$params.load) > 0){
	  lamdba.names <- all.params[unique(as.vector(inputs$measurement$params.load[[1]]))]
	  lamdba.names <- lamdba.names[!lamdba.names %in% c("fixed", "0")]
	  if(length(lamdba.names) > 0){
	    dLambdparLamb <- unList(differentiateMatrixOfVariable(matrix(sapply(inputs$measurement$params.load[[1]], function(x, all.params){if(x>0) all.params[x] else x}, all.params), nrow=nrow(inputs$measurement$params.load[[1]])),all.params[inputs$measurement$params.load[[1]]]))
        dLambdparLamb2 <- unList(differentiateMatrixOfVariable(dLambdparLamb, all.params[inputs$measurement$params.load[[1]]]))
	  }
	  else{
	    # when NLambda = 0 (all elements are fixed)
	    dLambdparLamb <-matrix(0L, nrow=0, ncol=length(inputs$measurement$params.load[[1]]))
	    dLambdparLamb2 <-matrix(0L, nrow=0, ncol=length(inputs$measurement$params.load[[1]]))
	  }
	  dLambdparLamb2 <- t(dLambdparLamb2)
	}
	else{
	  # when measurement@params.load is absent
	  dLambdparLamb <-matrix(0L, nrow=0, ncol=0)
	  dLambdparLamb2 <-matrix(0L, nrow=0, ncol=0)
	}
	
	sigmae.names<- c()
	if(length(inputs$noise$params.observed) > 0){
	  sigmae.names <- all.params[unique(as.vector(inputs$noise$params.observed[[1]]))]
	  sigmae.names <- sigmae.names[!sigmae.names%in% c("fixed", "0")]
	} 
	if(length(sigmae.names) > 0){
	  dSigmaede <-  differentiateMatrixOfVariable(returnExponentialSymbolicTerm(matrix(sapply(inputs$noise$params.observed[[1]], function(x, all.params){if(x>0) all.params[x] else x}, all.params), nrow=nrow(inputs$noise$params.observed[[1]]))), sigmae.names)
	  dSigmaede2<- differentiateMatrixOfVariable(dSigmaede, sigmae.names)
	  dSigmaede2 <- t(dSigmaede2)
	}
	else{
	  dSigmaede <-matrix(0L, nrow=0, ncol=0)
	  dSigmaede2 <-matrix(0L, nrow=0, ncol=0)
	}
	
	
	# [Note] The following part do the formula/theta.formula expansion for freeIC cases to generate the differentiation code, and we are no longer use way. (Instead, we do initial value estimation in dynr.cook). So the following part is commented out on 20/09/09. 
	# For freeIC ones, extend the equations and redo the corresponding differentiation to generate the differentiation code
	# if(freeIC == TRUE){
	  # formula <- inputs$dynamics@formula[[1]]
	  # theta.formula <- inputs$dynamics@theta.formula
      # #zeta_? to be estimated
	  # startval.names <- names(inputs$dynamics@startval)
	  # beta.names <- c()
      # for (i in 1:length(startval.names)){
	    # if(startval.names[[i]] > 0){
          # formula[[length(formula) + 1]] = as.formula(paste0(startval.names[[i]],' ~ 0'))
		  # beta.names = c(beta.names,startval.names[[i]])
		# }
      # }
    
      # # generate the formula for states
      # # note: the formula should be removed in dynr.model if the state in prep.init is fixed
	  # state.names <- inputs$dynamics$state.names
      # for(i in 1:length(state.names)){
        # if(state.names[[i]] > 0){
		  # formula[[length(formula) + 1]] = as.formula(paste0('init_',state.names[[i]],' ~ 0'))
          # beta.names = c(beta.names, paste0('init_',state.names[[i]]))
		# }
      # }
	  # inputs$dynamics@beta.names <- beta.names


      # #process the theta formula
	  # theta.names <- inputs$dynamics$theta.names
      # for (i in 1:length(state.names)){
        # # generate the theta.formula for states
        # # for state x, the corresponding theta.names = x_0
        # #              the corresponding formula is x_0 ~ 0 
        # # (later in dynr.model, the formulas will be modified according to inputs of prep.initial
		# if(state.names[[i]] > 0){
          # theta.formula[[length(theta.formula) + 1]] <- as.formula(paste0(state.names[[i]],'_0 ~ 1 * 0'))
          # theta.names[[length(theta.names) + 1]] <- paste0(state.names[[i]],'_0')
		# }
      # }
	  # inputs$dynamics@theta.names <- theta.names
	  # inputs$dynamics@theta.formula <- theta.formula
	  
	  # inputs$dynamics@jacobianOriginal <- autojacobTry(list(formula))
	  # inputs$dynamics@dfdtheta <- autojacobTry(inputs$dynamics@formulaOriginal, diff.variables=theta.names)
	  # dfdx <- autojacobTry(inputs$dynamics@formulaOriginal, diff.variables=state.names)
	  # inputs$dynamics@dfdx2 <- autojacobTry(dfdx, diff.variables=state.names)
	  # inputs$dynamics@dfdxdtheta <- autojacobTry(dfdx, diff.variables=theta.names)
	  # inputs$dynamics@dfdthetadx <- autojacobTry(inputs$dynamics@dfdtheta, diff.variables=state.names)
	  # inputs$dynamics@dfdtheta2 <- autojacobTry(inputs$dynamics@dfdtheta, diff.variables=theta.names)
	# }
  }
  #  ------- The above lines obtain the necessary components of SAEM processs ----------------------------
  
  if(any(sapply(inputs, class) %in% 'dynrRegimes')){
    numRegimes <- dim(inputs$regimes$values)[1]
    if (!is.null(data$covariate.names) && !all(inputs$regimes$covariates %in% data$covariate.names)){
      stop("The 'covariates' slot of the 'dynrRegimes' object should match the 'covariates' argument passed to the dynr.data() function.\nA pox on your house if fair Romeo had not found this.")
    }
  } else {
    numRegimes <- 1L
  }
  
  #cat('numRegimes = ', numRegimes)
  
  #The following lines basically compare the maximum number of regimes implied by each recipe object
  #against (1 or the maximum regimes across all recipes). Elaborated error messages to be more informative.
  all.regimes <- sapply(inputs[!names(inputs) %in% "data"], impliedRegimes)
  if(!all(all.regimes %in% c(1, max(all.regimes)))){
    likelyNum = as.numeric(names(sort(table(all.regimes[all.regimes > 1]),decreasing=TRUE)[1]))
    deviantR = names(all.regimes[!all.regimes %in% c(likelyNum) & all.regimes > 1])
    stop(paste0("Recipes imply differing numbers of regimes. Here they are:\n", 
                paste(paste0(names(all.regimes), " (", all.regimes, ")"), collapse=", "), 
                "\nNumber of regimes in each recipe must be ",numRegimes," according to prep.regimes, \nor 1 (same configuration automatically extended to all regimes).\nPlease check : ",paste(deviantR,collapse=", ")))
  }
  

  # writeCcode on each recipe
  if(saem == FALSE){
    inputs <- sapply(inputs, writeCcode, data$covariate.names)
  }
  else if (saem == TRUE){
    #browser()
    inputs <- sapply(inputs, writeArmadilloCode, data$covariate.names, as.character(c(param.data$param.name, sigmab.names)), dmudparMu=dmudparMu, dmudparMu2=dmudparMu2, dLambdparLamb=dLambdparLamb, dLambdparLamb2=dLambdparLamb2,dSigmaede=dSigmaede, dSigmaede2=dSigmaede2, dSigmabdb=dSigmabdb, dSigmabdb2=dSigmabdb2, Sigmab=ret$ldl, known.vars=ret$pars)
	#inputs <- sapply(inputs, writeArmadilloCode, covariates=data$covariate.names, param.names=as.character(param.data$param.name), dmudparMu=dmudparMu, dmudparMu2=dmudparMu2)
  } else {stop("Invalid value passed to 'saem' argument. It should be TRUE or FALSE.")}
  all.values <- unlist(sapply(inputs, slot, name='startval'))
  unique.values <- extractValues(all.values, all.params)
  
  if(length(inputs$transform$formula.inv)>0){
    unique.values <- inputs$transform$inv.tfun(unique.values)
  }
  param.data$param.value=unique.values
  
  #initiate a dynrModel object
  if(saem==FALSE){
    # note: random.params.inicov of EstimateRandomAsLV is prepared in dynr.cook. Thus there is no need to attached it here
    obj.dynrModel <- new("dynrModel", c(list(data=data, outfile=outfile, param.names=as.character(param.data$param.name)), inputs))
    obj.dynrModel@dim_latent_var <- dim(inputs$measurement$values.load[[1]])[2] #numbber of columns of the factor loadings
    obj.dynrModel@num_regime <- numRegimes
  }
  else if(saem==TRUE){
  #browser()
    obj.dynrModel <- new("dynrModel", c(list(data=data, outfile=outfile, param.names=as.character(param.data$param.name), random.params.inicov=random.params.inicov, random.values.inicov=random.values.inicov, dmudparMu=dmudparMu, dmudparMu2=dmudparMu2, dLambdparLamb=dLambdparLamb,dLambdparLamb2=dLambdparLamb2, dSigmaede=dSigmaede, dSigmaede2=dSigmaede2, dSigmabdb=dSigmabdb, dSigmabdb2=dSigmabdb2, freeIC=freeIC, Sigmab=ret$ldl, known.vars=known.vars), inputs))
  }
  
  
  
  obj.dynrModel@dim_latent_var <- dim(inputs$measurement$values.load[[1]])[2] #numbber of columns of the factor loadings  
  obj.dynrModel@xstart <- param.data$param.value
  obj.dynrModel@ub <- as.double(rep(NA, length(obj.dynrModel@xstart)))
  obj.dynrModel@lb <- as.double(rep(NA, length(obj.dynrModel@xstart)))
  names(obj.dynrModel@xstart) <- obj.dynrModel@param.names
  names(obj.dynrModel@ub) <- obj.dynrModel@param.names
  names(obj.dynrModel@lb) <- obj.dynrModel@param.names
  
   
  
  
  #write out the C script
  cparts <- unlist(sapply(inputs, slot, name='c.string'))
  if(saem==FALSE){
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
  }
  else if(saem==TRUE){
    glom <- paste(cparts)
    cat(glom, file=obj.dynrModel@outfile)
  } else {stop("Invalid value passed to 'saem' argument. It should be TRUE or FALSE.")}
  
  if (.hasSlot(obj.dynrModel$dynamics, 'theta.formula') && length(obj.dynrModel$dynamics@theta.formula) > 0 && obj.dynrModel$dynamics$saem==FALSE){
    #get the extended model and return it (don't keep the original model)
    extended_model <- ExpandRandomAsLVModel(obj.dynrModel)
	return(extended_model)
  }

  return(obj.dynrModel)
  #modify the object slot, including starting values, etc.
}


addVariableToThetaFormula  <- function(formula, variable.names){
    fml=as.character(formula)
    lhs=fml[[2]]
    rhs=fml[[3]]
    
    rhs = paste0(rhs, "+ 1 *", variable.names)
    
    #print(as.formula(paste0(lhs, ' ~ ', rhs)))
    return(as.formula(paste0(lhs, ' ~ ', rhs)))
}

impliedRegimes <- function(recipe){
	if(inherits(recipe, 'dynrRecipe')){
		sn <- slotNames(recipe)
		if (class(recipe)=="dynrDynamicsFormula"){
		  vs=c("formula","jacobian")
		}
		else{vs <- grep('^values\\.', sn, value=TRUE)}
		if(length(vs) > 0){
			sl <- lapply(vs, FUN=slot, object=recipe)
			nr <- max(sapply(sl, FUN=function(x){if(is.list(x)){return(length(x))} else {return(1)}}))
		} else if(inherits(recipe, 'dynrRegimes')){
			nr <- dim(recipe$values)[1]
		} else {
			nr <- 1
		}
	} else {
		nr <- 1
	}
	return(nr)
}







checkSSMConformable <- function(mat, rows, cols, matname, modname){
  if( nrow(mat) != rows || ncol(mat) != cols ){
    msg <- paste("The ", matname, " matrix is not the correct size",
                 " in the state space expectation of model ", modname,
                 ".  It is ", nrow(mat), " by ", ncol(mat), " and should be ",
                 rows, " by ", cols, ".", sep="")
    stop(msg, call. = FALSE)
  }
}


##' Extend a user-specified model to include random varibles
##' 
##' @param dynrModel a dynrModel object prepared with recipe functions \code{\link{prep.formulaDynamics}}, \code{\link{prep.measurement}}, \code{\link{prep.noise}}, \code{\link{prep.initial}}, \code{\link{dynr.data}}. 

##' 
##' @return an object of dynrModel that is the expanede model.
##' 
##' @details
##' A \code{dynrModel} is a collection of recipes.  The recipes are constructed with the functions unctions \code{\link{prep.formulaDynamics}}, \code{\link{prep.measurement}}, \code{\link{prep.noise}}, \code{\link{prep.initial}}. Additionally, data must be prepared with \code{\link{dynr.data}} and added to the model.
##' 
##' 
##' @examples
##' #model <- dynr.model(dynamics=dynm, measurement=meas, noise=mdcov, initial=initial, data=data, outfile="osc.cpp")
##' #extended_model <- ExpandRandomAsLVModel(model)
##'
##' 
##' #For full demo examples, see:
##' #demo(OscWithRand, package="dynr")
##' #demo(VDPwithRand, package="dynr")
##' @export
EstimateRandomAsLVModel<- function(dynrModel){  
    # Restructure mixed effects structured via theta.formula into an expanded model with 
    # random effects as additional state variables and cook it.
    #browser()
    if(.hasSlot(dynrModel@dynamics,'random.names') && length(dynrModel@dynamics@random.names) > 0){
        user.random.names = setdiff(dynrModel@dynamics@random.names, paste0('b_', dynrModel@measurement@state.names))
        
        # Add the user-specified random effects into states to be estimated
        # state.names2 = state.names + user.random.names
        state.names2 = c(dynrModel@measurement@state.names, user.random.names)
        #print(paste("New state variables to be estimated:", state.names2))
    }
    else{
        stop("There is no random effect variables to be estimated. Initial value estimates are done.")
        #return(list(coefEst=coefEst))
        return(list())
    }
    
	##New cov matrices
    params.latent <- function(){
        icol <- dynrModel@noise@params.latent[[1]]
        ind <- which(icol!=0)
        icol[ind] <- dynrModel$'param.names'[icol[ind]]
        for (i in 1:length(user.random.names)){
        icol <- rbind(icol, 'fixed') 
        icol <- cbind(icol, 'fixed')
        }
        return(icol)
    }
    #browser()
    # If there is random effect to be estimated, set up a new model
    params.latent = diag(c(diag(dynrModel@noise@params.latent[[1]]), rep(0, length(user.random.names))))
    mdcov2 <- prep.noise(
        values.latent=diag(c(diag(dynrModel@noise@values.latent[[1]]), rep(0, length(user.random.names)))),
        params.latent=matrix(mapply(function(x) {if(x > 0){return(dynrModel@param.names[x])} else{return("fixed")}}, params.latent), nrow=nrow(params.latent)),
        #params.latent=diag(rep("fixed",length(state.names2)), length(state.names2)),
        #params.latent=diag(state.names2, length(state.names2)),
        #params.latent=diag(c(diag(dynrModel@noise@params.latent[[1]]), rep('fixed',length(user.random.names)))),
        values.observed=dynrModel@noise@values.observed[[1]],
        params.observed=matrix(mapply(function(x) {if(x > 0){return(dynrModel@param.names[x])} else{return("fixed")}}, dynrModel@noise@params.observed[[1]]), nrow=nrow(dynrModel@noise@params.observed[[1]]))
        #params.observed=dynrModel@noise@params.observed[[1]]
    )
    

    num.y = length(dynrModel@measurement@obs.names)
    #lambda matrix
	if (length(dynrModel@measurement@params.int) > 0){
		meas2 <- prep.measurement(
			values.load = matrix(c(as.vector(dynrModel@measurement@values.load[[1]]), rep(0, num.y * length(user.random.names))), nrow=num.y, ncol= length(state.names2), byrow=FALSE),
			params.load = matrix(c(sapply(dynrModel@measurement@params.load[[1]], function(x) {if(x > 0){return(dynrModel@param.names[x])} else{return("fixed")}}), rep("fixed", num.y * length(user.random.names))), nrow=num.y, ncol= length(state.names2), byrow=FALSE),
			#params.load = matrix(c(dynrModel@measurement@params.load[[1]], rep("fixed", num.y * length(user.random.names))), nrow=num.y, ncol= length(state.names2), byrow=FALSE),
			obs.names = dynrModel@measurement@obs.names,
			state.names = state.names2,
			values.int=dynrModel@measurement@values.int,
			#params.int=dynrModel@measurement@params.int
			params.int=matrix(mapply(function(x) {if(x > 0){return(dynrModel@param.names[x])} else{return("fixed")}}, dynrModel@measurement@params.int[[1]]), nrow=nrow(dynrModel@measurement@params.int[[1]]))
			)
	}
	else{
		meas2 <- prep.measurement(
			values.load = matrix(c(as.vector(dynrModel@measurement@values.load[[1]]), rep(0, num.y * length(user.random.names))), nrow=num.y, ncol= length(state.names2), byrow=FALSE),
			params.load = matrix(c(sapply(dynrModel@measurement@params.load[[1]], function(x) {if(x > 0){return(dynrModel@param.names[x])} else{return("fixed")}}), rep("fixed", num.y * length(user.random.names))), nrow=num.y, ncol= length(state.names2), byrow=FALSE),
			#params.load = matrix(c(dynrModel@measurement@params.load[[1]], rep("fixed", num.y * length(user.random.names))), nrow=num.y, ncol= length(state.names2), byrow=FALSE),
			obs.names = dynrModel@measurement@obs.names,
			state.names = state.names2
			)
	}

    # Generate the new variance/covariance matrix 
    # by adding the user-specified random names into states
    #   - state.names2: c(state.names, random.names)
    #   - values.inicov: initial values of variance/covariance matrix
    #       *state.names: from the first fitted model
    #       *random.names: from user specified in dynrModel@dynamics 
    num.state = length(dynrModel@measurement@state.names)
    num.state2 = length(state.names2) 
    values.inicov = diag(1, length(state.names2))
    values.inicov[1:num.state,1:num.state] = dynrModel@initial@values.inicov[[1]]
    values.inicov[(num.state+1):num.state2,(num.state+1):num.state2] = dynrModel@dynamics@random.values.inicov
    params.inicov = matrix("fixed", nrow=nrow(values.inicov), ncol=ncol(values.inicov))
    params.inicov[1:num.state,1:num.state] = matrix(mapply(function(x) {if(x > 0){return(dynrModel@param.names[x])} else{return("fixed")}}, dynrModel@initial@params.inicov[[1]]), nrow=nrow(dynrModel@initial@params.inicov[[1]]))
    #params.inicov[1:num.state,1:num.state] = matrix(dynrModel@initial@params.inicov[[1]], nrow=nrow(dynrModel@initial@params.inicov[[1]]))
    params.inicov[(num.state+1):num.state2,(num.state+1):num.state2]  = dynrModel@dynamics@random.params.inicov 
    
    initial2 <- prep.initial(
        values.inistate=c(as.vector(dynrModel@initial@values.inistate[[1]]), rep(0, num.state2 - num.state)),
        params.inistate=c(sapply(dynrModel@initial@params.inistate[[1]], function(x) {if(x > 0){return(dynrModel@param.names[x])} else{return("fixed")}}), rep("fixed", num.state2 - num.state)),
        #params.inistate=c(dynrModel@initial@params.inistate[[1]], rep("fixed", num.state2 - num.state)),
        values.inicov=values.inicov, 
        params.inicov=params.inicov)
    
    
    # Formula processing:  
    # 1. If the formula has been already extended to include random.names and mu_x1, mu_x2,
    # only retrieve the formula with state variables as LHS
    formula <- list(dynrModel@dynamics@formulaOriginal[[1]][1:length(dynrModel@measurement@state.names)])
    
    # 2. If theta.formula exists, substitute the content of theta.formula 
    if(length(dynrModel@dynamics@theta.formula) > 0){
        formula <- lapply(formula, function(x){substituteFormula(x, dynrModel@dynamics@theta.formula)})  
    }
    
    #formula <- unlist(dynrModel@dynamics@formula2)[1:length(dynrModel@measurement@state.names)]
    for(i in ((length(dynrModel@measurement@state.names)+1):length(state.names2)) )
        formula[[1]][[i]] <- as.formula(paste0(state.names2[i], '~ 0')) 
    
    dynm2<-prep.formulaDynamics(formula=unlist(formula),
                                startval=dynrModel@dynamics@startval,
                                isContinuousTime=dynrModel@dynamics@isContinuousTime)
    
    model2 <- dynr.model(dynamics=dynm2, measurement=meas2,
                         noise=mdcov2, initial=initial2, data=dynrModel@data)
    #fitted_model2 <- dynr.cook(model2, optimization_flag=optimization_flag, hessian_flag = hessian_flag, verbose=verbose, weight_flag=weight_flag, debug_flag=debug_flag)
    #return(fitted_model2)
	return(model2)
    
}
