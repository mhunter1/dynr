# A Recipe is a helper function that takes user input
#  and produces a C function definition that can be
#  compiled by dynrFuncaddress.

# TODO add check that there is no cross loading on one of the ID variables

# TODO output the starting values for the parameters 
# TODO add documentation


#---------------------------------------------------------------------------------------------
# Create dynrRecipe object class

##' The dynrRecipe Class
##' 
##' @aliases
##' $,dynrRecipe-method
##' 
##' @details
##' This is an internal class structure.  You should not use it directly.
##' The following are all subclasses of this class: \code{\link{dynrMeasurement-class}}, 
##' \code{\link{dynrDynamics-class}}, \code{\link{dynrRegimes-class}}, 
##' \code{\link{dynrInitial-class}}, \code{\link{dynrNoise-class}}, 
##' and \code{\link{dynrTrans-class}}.  Recipes are the things that 
##' go into a \code{\link{dynrModel-class}} using \code{\link{dynr.model}}.
##' Use the recipe prep functions (\code{\link{prep.measurement}}, 
##' \code{\link{prep.formulaDynamics}}, \code{\link{prep.matrixDynamics}}, 
##' \code{\link{prep.regimes}}, \code{\link{prep.initial}}, \code{\link{prep.noise}}, 
##' or \code{\link{prep.tfun}}) to create these classes instead.
setClass(Class =  "dynrRecipe",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnames = "character"
         )
)

##' The dynrMeasurement Class
##' 
##' @details
##' This is an internal class structure.  You should not use it directly.
##' Use \code{\link{prep.measurement}} or \code{\link{prep.loadings}} instead.
setClass(Class = "dynrMeasurement",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnames = "character",
           values.load = "list",
           params.load = "list",
           values.exo = "list",
           params.exo = "list",
           values.int = "list",
           params.int = "list",
           state.names="character",
           obs.names="character",
           exo.names="character"),
         contains = "dynrRecipe"
)

##' The dynrDynamics Class
##' 
##' @aliases
##' dynrDynamicsFormula-class
##' dynrDynamicsMatrix-class
##' 
##' @details
##' This is an internal class structure.  The classes \code{dynrDynamicsFormula-class} 
##' and \code{dynrDynamicsMatrix-class} are subclasses of this.  However, you should 
##' not use it directly.
##' Use \code{\link{prep.matrixDynamics}} or \code{\link{prep.formulaDynamics}} instead.
setClass(Class = "dynrDynamics",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnames = "character",
           isContinuousTime = "logical"
           ),
         contains = "dynrRecipe"
)

setClass(Class = "dynrDynamicsFormula",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnames = "character",
           formula = "list",
           jacobian = "list",
           isContinuousTime = "logical"
           ),
         contains = "dynrDynamics"
)

setClass(Class = "dynrDynamicsMatrix",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnames = "character",
           values.dyn = "list",
           params.dyn = "list",
           values.exo = "list",
           params.exo = "list",
           values.int = "list",
           params.int = "list",
           covariates = "character",
           isContinuousTime = "logical"
           ),
         contains = "dynrDynamics"
)

##' The dynrRegimes Class
##' 
##' @details
##' This is an internal class structure.  You should not use it directly.
##' Use \code{\link{prep.regimes}} instead.
setClass(Class = "dynrRegimes",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnames = "character",
           values = "matrix",
           params = "matrix",
           covariates = "character"),
         contains = "dynrRecipe"
)

##' The dynrInitial Class
##' 
##' @details
##' This is an internal class structure.  You should not use it directly.
##' Use \code{\link{prep.initial}} instead.
setClass(Class = "dynrInitial",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnames = "character",
           values.inistate = "list",
           params.inistate = "list",
           values.inicov = "list",
           values.inicov.inv.ldl = "list",
           params.inicov = "list",
           values.regimep = "matrix",
           params.regimep = "matrix"),
         contains = "dynrRecipe"
)


##' The dynrNoise Class
##' 
##' @details
##' This is an internal class structure.  You should not use it directly.
##' Use \code{\link{prep.noise}} instead.
setClass(Class = "dynrNoise",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnames = "character",
           values.latent = "list",
           params.latent = "list",
           values.observed = "list",
           params.observed = "list",
           values.latent.inv.ldl = "list",
           values.observed.inv.ldl = "list"),
         contains = "dynrRecipe"
)
#TODO we should emphasize that either the full noise covariance structures should be freed or the diagonals because we are to apply the ldl trans  

##' The dynrTrans Class
##' 
##' @details
##' This is an internal class structure.  You should not use it directly.
##' Use \code{\link{prep.tfun}} instead.
setClass(Class = "dynrTrans",
         representation = representation(
           c.string =  "character",
           startval = "numeric",#not sure if needed in dynrTrans
           paramnames = "character",#not sure if needed in dynrTrans
           tfun="function",
           inv.tfun="function",
           inv.tfun.full="function",
           formula.trans="list",
           formula.inv="list",
           transCcode="logical"
           ),
         contains = "dynrRecipe",
         prototype = prototype(
           tfun=function(x){x}
         )
)

setMethod("initialize", "dynrRecipe",
          function(.Object, x){
            for(i in names(x)){
				slot(.Object, name=i, check = TRUE) <- x[[i]]
			}
			return(.Object)
          }
)

setMethod("$", "dynrRecipe",
          function(x, name){slot(x, name)}
)



#------------------------------------------------------------------------------
# printex method definitions

##' The printex Method
##' 
##' @aliases
##' printex,dynrModel-method
##' printex,dynrCook-method
##' printex,dynrMeasurement-method
##' printex,dynrDynamicsFormula-method
##' printex,dynrDynamicsMatrix-method
##' printex,dynrNoise-method
##' printex,dynrInitial-method
##' printex,dynrRegimes-method
##' 
##' @param object The dynr object (recipe, model, cooked, model).
##' @param ParameterAs The parameter values or names to plot. The underscores in parameter names are 
##' saved for use of subscripts.  Greek letters can be specified as corresponding LaTeX symbols without ##' backslashes (e.g., "lambda") and printed as greek letters.
##' @param printDyn logical. Whether or not to print the dynamic model. The default is TRUE.
##' @param printMeas logical. Whether or not to print the measurement model. The default is TRUE.
##' @param printInit logical. Whether or not to print the initial conditions. The default is FALSE.
##' @param printRS logical. Whether or not to print the regime-switching model. The default is FALSE.
##' @param outFile The name of the output tex file.
##' @param show logical indicator of whether or not to show the result in the console. 
##' @param ... Further named arguments, passed to internal method. 
##' \code{AsMatrix} is a logical indicator of whether to put the object in matrix form.
##' 
##' @details
##' This is a general way of getting a LaTeX string for recipes, 
##' models, and cooked models.  It is a great way to check that 
##' you specified the model or recipe you think you did before 
##' estimating its free parameters (cooking).  After the model 
##' is cooked, you can use it to get LaTeX code with the estimated 
##' parameters in it.
##' 
##' Typical inputs to the \code{ParameterAs} argument are (1) the starting values for a model, (2) the final estimated values for a model, and (3) the parameter names.  These are accessible with (1) \code{model$xstart}, (2) \code{coef(cook)}, and (3) \code{model$param.names} or \code{names(coef(cook))}, respectively.
##' 
##' @seealso A way to put this in a plot with \code{\link{plotFormula}}
setGeneric("printex", function(object, ParameterAs, 
	printDyn=TRUE, printMeas=TRUE, printInit=FALSE, printRS=FALSE, outFile, show, ...) { 
	return(standardGeneric("printex")) 
})

setMethod("printex", "dynrMeasurement",
          function(object, ParameterAs, 
			  printDyn=TRUE, printMeas=TRUE, printInit=FALSE, printRS=FALSE, outFile, 
			  show=TRUE, AsMatrix=TRUE){
            if (AsMatrix){
              meas_loadings=lapply(object$values.load, .xtableMatrix, show)
              meas_int=lapply(object$values.int, .xtableMatrix, show)
              meas_exo=lapply(object$values.exo, .xtableMatrix, show)
              meas_exo.names=.xtableMatrix(matrix(object$exo.names,ncol=1), show)
              meas_list = list(meas_loadings=meas_loadings,meas_int=meas_int,
                               meas_exo=meas_exo,meas_exo.names=meas_exo.names)
              return(invisible(meas_list))
            }else{#not to use matrix
              #Matrix to Formula
              meas_loadings <- lapply(object$values.load,matrix2formula,multbyColnames=T)
              meas_int <- lapply((object)$values.int,matrix2formula,multbyColnames=F)
              meas_exo <- lapply((object)$values.exo,matrix2formula,multbyColnames=T)
              meas_fml <- addLLFormulas(list_list_formulae = meas_loadings, 
                                        VecNamesToAdd = c("meas_int","meas_exo"))
              #Formula to LaTex
              meas_tex=lapply(meas_fml, formula2tex, 
                              LHSvarPre = "", LHSvarPost = "", RHSvarPre = "", RHSvarPost = "", 
                              LHSeqPre = "", LHSeqPost = "", RHSeqPre = "", RHSeqPost = "")
              return(meas_tex)
            }
            
          }
)

setMethod("printex", "dynrDynamicsMatrix",
          function(object, ParameterAs, 
			  printDyn=TRUE, printMeas=TRUE, printInit=FALSE, printRS=FALSE, outFile, 
			  show=TRUE, AsMatrix=TRUE){
            if (AsMatrix){
              dyn_tran=lapply((object)$values.dyn,.xtableMatrix, show)
              dyn_int=lapply((object)$values.int,.xtableMatrix, show)
              dyn_exo=lapply((object)$values.exo,.xtableMatrix, show)
              dyn_exo.names=.xtableMatrix(matrix((object)$covariates,ncol=1),show)
              return(invisible(list(dyn_tran = dyn_tran, dyn_int = dyn_int, dyn_exo = dyn_exo, dyn_exo.names = dyn_exo.names) ))
            }else{#not to use matrix
              #Matrix to Formula
              dyn_tran <- lapply(object$values.dyn,matrix2formula,multbyColnames=T)
              dyn_int <- lapply((object)$values.int,matrix2formula,multbyColnames=F)
              dyn_exo <- lapply((object)$values.exo,matrix2formula,multbyColnames=T)
              dyn_fml <- addLLFormulas(list_list_formulae = dyn_tran, 
                                       VecNamesToAdd = c("dyn_int","dyn_exo"))
              #Formula to LaTex
              if (object$isContinuousTime){
                LHSvarPre <- "d("
                LHSvarPost <- "(t))"
                RHSeqPre <- "("
                RHSeqPost <- ")dt"
              }else{
                LHSvarPre <- ""
                LHSvarPost <- "(t+1)"
                RHSeqPre <- ""
                RHSeqPost <- ""
              }
              RHSvarPost="(t)"
              dyn_tex=lapply(dyn_fml, formula2tex, 
                             LHSvarPre = LHSvarPre, LHSvarPost = LHSvarPost, RHSvarPre = "", RHSvarPost = RHSvarPost, 
                             LHSeqPre = "", LHSeqPost = "", RHSeqPre = RHSeqPre, RHSeqPost = RHSeqPost)
              return(dyn_tex)
            }
          }    
)

setMethod("printex", "dynrRegimes",
          function(object, ParameterAs, 
			  printDyn=TRUE, printMeas=TRUE, printInit=FALSE, printRS=FALSE, outFile, 
			  show=TRUE, AsMatrix=TRUE){
            lG <- ifelse(nrow(object$values) != 0, .xtableMatrix(object$values, show), "")
            return(invisible(list(regimes=lG)))
          }
)


setMethod("printex", "dynrInitial",
          function(object, ParameterAs, 
			  printDyn=TRUE, printMeas=TRUE, printInit=FALSE, printRS=FALSE, outFile, 
			  show=TRUE, AsMatrix=TRUE){
            lx0 <- lapply(object$values.inistate, .xtableMatrix, show)
            lP0 <- lapply(object$values.inicov, .xtableMatrix, show)
            lr0 <- .xtableMatrix(object$values.regimep, show)
            return(invisible(list(initial.state=lx0, initial.covariance=lP0, initial.probability=lr0)))
          }
)


setMethod("printex", "dynrNoise",
          function(object, ParameterAs, 
			  printDyn=TRUE, printMeas=TRUE, printInit=FALSE, printRS=FALSE, outFile, 
			  show=TRUE, AsMatrix=TRUE){
            lQ <- lapply(object$values.latent, .xtableMatrix, show=show)
            lR <- lapply(object$values.observed, .xtableMatrix, show=show)
            return(invisible(list(dynamic.noise=lQ, measurement.noise=lR)))
          }
)

setMethod("printex", "dynrDynamicsFormula",
          function(object, ParameterAs, 
			  printDyn=TRUE, printMeas=TRUE, printInit=FALSE, printRS=FALSE, outFile, 
			  show=TRUE, AsMatrix=TRUE){
            if (object$isContinuousTime){
            LHSvarPre <- "d("
            LHSvarPost <- "(t))"
            RHSeqPre <- "("
            RHSeqPost <- ")dt"
            }else{
            LHSvarPre <- ""
            LHSvarPost <- "(t+1)"
            RHSeqPre <- ""
            RHSeqPost <- ""
            }
            RHSvarPost="(t)"
            dyn=lapply(object$formula, formula2tex, 
                       LHSvarPre = LHSvarPre, LHSvarPost = LHSvarPost, RHSvarPre = "", RHSvarPost = RHSvarPost, 
                       LHSeqPre = "", LHSeqPost = "", RHSeqPre = RHSeqPre, RHSeqPost = RHSeqPost)
            return(dyn)
          }
)

#This function transforms a list of formulae to lists of latex code.
#returns a list of left, right, and equation latex code, each a vector.
formula2tex<-function(list_formulae, 
                      LHSvarPre = "", LHSvarPost = "", RHSvarPre = "", RHSvarPost = "", 
                      LHSeqPre = "$", LHSeqPost = "", RHSeqPre = "", RHSeqPost = "$"){
  eq.char=lapply(list_formulae, as.character)
  str.left=sapply(eq.char,"[",2)
  str.right=sapply(eq.char,"[",3)
  neq=length(list_formulae)
  mulpatn<-"([[:print:]]*)"
  sigpatn<-"([0-9A-Za-z_. ^*]*)"
  for (j in 1:neq){
    #The order of gsub matters.
    str.right[j]=gsub(paste0("\\(",mulpatn,"\\)/\\(",mulpatn,"\\)"),"\\\\frac{\\1}{\\2}",str.right[j])
    str.right[j]=gsub(paste0("\\(",mulpatn,"\\)/",sigpatn),"\\\\frac{\\1}{\\2}",str.right[j])
    str.right[j]=gsub(paste0(sigpatn,"/\\(",mulpatn,"\\)"),"\\\\frac{\\1}{\\2}",str.right[j])
    str.right[j]=gsub(paste0(sigpatn,"/",sigpatn),"\\\\frac{\\1}{\\2}",str.right[j])
    str.right[j]=gsub("\\bexp\\(", "\\\\exp\\(", str.right[j])
    str.right[j]=gsub("\\blog\\(", "\\\\log\\(", str.right[j])
    str.right[j]=gsub("\\bsin\\(", "\\\\sin\\(", str.right[j])
    str.right[j]=gsub("\\bcos\\(", "\\\\cos\\(", str.right[j])
    str.right[j]=gsub("\\*","\\\\times",str.right[j])
  }
  
  #Modify the presentation of the variables on RHS
  for (j in 1:neq){
    for (i in 1:neq){
      str.right[j]=gsub(paste0("\\<",str.left[i],"\\>"),paste0(RHSvarPre,str.left[i],RHSvarPost),str.right[j])
    }
  }
  #Modify the presentation of the variables on LHS
  for (j in 1:neq){
    str.left[j]=gsub(paste0(sigpatn),paste0(LHSvarPre,"\\1",LHSvarPost),str.left[j])
  }
  
  return(invisible(#list(left=str.left,right=str.right,equation=
    paste0(LHSeqPre, str.left, LHSeqPost," = ", RHSeqPre, str.right, RHSeqPost))
    )#)
}

.xtableMatrix <- function(m, show){
	x <- xtable::xtable(m, align=rep("", ncol(m)+1))
	out <- print(x, floating=FALSE, tabular.environment="bmatrix", 
		hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE,
		print.results=show)
	return(out)
}


#------------------------------------------------------------------------------
# paramName2Number method definitions

setGeneric("paramName2Number", function(object, names) { 
	return(standardGeneric("paramName2Number")) 
})

setGeneric("paramName2NumericNumber", function(object, paramList) { 
  return(standardGeneric("paramName2NumericNumber")) 
})

.exchangeNamesAndNumbers <- function(params, names){
	matrix(match(params, names, nomatch=0), nrow(params), ncol(params))
}

.exchangeformulaNamesAndNumbers <- function(formula, paramnames, names){
	string<-paste0(deparse(formula,width.cutoff = 500L),collapse="")
	pattern <- paramnames
	#pattern=gsub("\\{","\\\\\\{",paramnames)
	#pattern=gsub("\\}","\\\\\\}",pattern)
	pattern=gsub("\\\\","\\\\\\\\",pattern)
	
	for(i in 1:length(paramnames)){
		string<-gsub(paste0("\\<",pattern[i],"\\>"),paste0("param[",match(paramnames[i], names, nomatch=0)-1,"]"), string)
	}
	eval(parse(text=string))
}

.replaceNamesbyNumbers <- function(formula, paramtoPlot){
  string<-paste0(deparse(formula,width.cutoff = 500L),collapse="")
  names <-   paste0("param\\[",(1:length(paramtoPlot))-1,"\\]")
  #names(paramtoPlot)
  #pattern=gsub("\\{","\\\\\\{",names)
  #pattern=gsub("\\}","\\\\\\}",pattern)
  #pattern=gsub("\\\\","\\\\\\\\",pattern)
  for (i in 1:length(names)){
  string<-gsub(paste0("\\<",names[i],"\\>"),paramtoPlot[i], string)
  }
  eval(parse(text=string))
}


setMethod("paramName2Number", "dynrMeasurement",
	function(object, names){
		object@params.load <- lapply(object$params.load, .exchangeNamesAndNumbers, names=names)
		object@params.exo <- lapply(object$params.exo, .exchangeNamesAndNumbers, names=names)
		object@params.int <- lapply(object$params.int, .exchangeNamesAndNumbers, names=names)
		return(object)
	}
)

.exchangeNamesInFormula <- function(formula, names){
	sall <- formula2string(formula)
	sr <- sall$rhs
	for(i in 1:length(sr)){
		for(j in 1:length(names)){
			pat <- paste0('\\<', names[j], '\\>')
			sr[i] <- gsub(pat, paste0('param[', j, ']'), sr[i])
		}
		formula[[i]] <- as.formula(paste(sall$lhs[i], '~', sr[i]))
	}
	return(formula)
}

# Example
## dynamics
#formula=list(
#  list(x1~a*x1,
#       x2~b*x2),
#  list(x1~a*x1+e*(exp(abs(x2))/(1+exp(abs(x2))))*x2,
#       x2~b*x2+f*(exp(abs(x1))/(1+exp(abs(x1))))*x1) 
#)
#paramnames <- c('a', 'b', 'e', 'f')

#.exchangeNamesInFormula(formula[[1]], paramnames)
#lapply(formula, .exchangeNamesInFormula, names=paramnames)


setMethod("paramName2Number", "dynrDynamicsFormula",
	function(object, names){
	  object@formula = .exchangeformulaNamesAndNumbers(object@formula, object@paramnames, names)
	  object@jacobian = .exchangeformulaNamesAndNumbers(object@jacobian, object@paramnames, names)
	  return(object)
	}
)

setMethod("paramName2NumericNumber", "dynrDynamicsFormula",
          function(object, paramList){
            object@formula = .replaceNamesbyNumbers(object@formula, paramList)
            object@jacobian =.replaceNamesbyNumbers(object@jacobian, paramList)
            return(object)
          }
)


setMethod("paramName2Number", "dynrDynamicsMatrix",
	function(object, names){
		object@params.dyn <- lapply(object$params.dyn, .exchangeNamesAndNumbers, names=names)
		object@params.exo <- lapply(object$params.exo, .exchangeNamesAndNumbers, names=names)
		object@params.int <- lapply(object$params.int, .exchangeNamesAndNumbers, names=names)
		return(object)
	}
)


setMethod("paramName2Number", "dynrRegimes",
	function(object, names){
		object@params <- .exchangeNamesAndNumbers(object$params, names)
		return(object)
	}
)


setMethod("paramName2Number", "dynrInitial",
	function(object, names){
		object@params.inistate <- lapply(object$params.inistate, .exchangeNamesAndNumbers, names=names)
		object@params.inicov <- lapply(object$params.inicov, .exchangeNamesAndNumbers, names=names)
		object@params.regimep <- .exchangeNamesAndNumbers(object$params.regimep, names=names)
		return(object)
	}
)


setMethod("paramName2Number", "dynrNoise",
	function(object, names){
		object@params.latent <- lapply(object$params.latent, .exchangeNamesAndNumbers, names=names)
		object@params.observed <- lapply(object$params.observed, .exchangeNamesAndNumbers, names=names)
		return(object)
	}
)

setMethod("paramName2Number", "dynrTrans",
  function(object, names){
    if (length(object@formula.trans)>0){
      object@formula.trans = .exchangeformulaNamesAndNumbers(object@formula.trans, object@paramnames, names)
    }
    return(object)
  }
)


#------------------------------------------------------------------------------
# writeCcode method definitions

setGeneric("writeCcode", 
           function(object, covariates, show=TRUE) { 
             return(standardGeneric("writeCcode")) 
           })



setMethod("writeCcode", "dynrMeasurement",
	function(object, covariates){
		values.load <- object$values.load
		params.load <- object$params.load
		values.exo <- object$values.exo
		params.exo <- object$params.exo
		values.int <- object$values.int
		params.int <- object$params.int
		hasCovariates <- length(values.exo) > 0
		hasIntercepts <- length(values.int) > 0
		nregime <- length(values.load)
		ret <- "void function_measurement(size_t t, size_t regime, double *param, const gsl_vector *eta, const gsl_vector *co_variate, gsl_matrix *Ht, gsl_vector *y){\n\n"
		if(hasCovariates){
		  ret <- paste0(ret, gslcovariate.front(object@exo.names, covariates))
		  ret <- paste(ret, createGslMatrix(nrow(params.exo[[1]]), ncol(params.exo[[1]]), "Bmatrix"), sep="\n\t")
		  }#create covariates matrix
		if(hasIntercepts){ret <- paste(ret, createGslVector(nrow(params.int[[1]]), "intVector"), sep="\n\t")}#create intercepts vector
		if(nregime > 1) {
			ret <- paste(ret, "\tswitch (regime) {\n", sep="\n")
			for(reg in 1:nregime){
				ret <- paste0(ret, paste0("\t\tcase ", reg-1, ":"), "\n")
				ret <- paste0(ret, "\t\t\t", setGslMatrixElements(values.load[[reg]], params.load[[reg]], "Ht"))
				if(hasCovariates){ #set covariates matrix
					ret <- paste0(ret, "\t\t\t", setGslMatrixElements(values.exo[[reg]], params.exo[[reg]], "Bmatrix"))
				}
				if(hasIntercepts){ret <- paste0(ret, "\t\t\t", setGslVectorElements(values.int[[reg]], params.int[[reg]], "intVector"))}#set intercepts vector
				ret <- paste0(ret, "\t\t", "break;", "\n") 
			}
			ret <- paste0(ret, "\t", "}\n\n")
		} else {
			ret <- paste(ret, setGslMatrixElements(values.load[[1]], params.load[[1]], "Ht"), sep="\n")
			if(hasCovariates){ret <- paste(ret, setGslMatrixElements(values.exo[[1]], params.exo[[1]], "Bmatrix"), sep="\n")}#set covariates matrix
			if(hasIntercepts){ret <- paste(ret, setGslVectorElements(values.int[[1]], params.int[[1]], "intVector"), sep="\n")}#set intercepts vector
		}
		ret <- paste(ret, "\n\tgsl_blas_dgemv(CblasNoTrans, 1.0, Ht, eta, 0.0, y);\n")
		if(hasCovariates){
		  ret <- paste(ret, "\n\tgsl_blas_dgemv(CblasNoTrans, 1.0, Bmatrix, covariate_local, 1.0, y);\n", destroyGslMatrix("Bmatrix"))
		  ret <- paste0(ret, destroyGslVector("covariate_local"))
		  }#multiply, add, and destroy covariates matrix
		if(hasIntercepts){ret <- paste(ret, "\n\tgsl_vector_add(y, intVector);\n", destroyGslVector("intVector"))}#add and destroy intercepts vector
		ret <- paste(ret, "\n}\n\n")
		object@c.string <- ret
		return(object)
	}
)


# not sure what to do here yet
setMethod("writeCcode", "dynrDynamicsFormula",
	function(object, covariates){
	  formula <- object$formula
	  jacob <- object$jacobian
	  nregime=length(formula)
	  n=sapply(formula,length)
	  
	  fml=lapply(formula,processFormula)
	  lhs=lapply(fml,function(x){lapply(x,"[[",1)})
	  rhs=lapply(fml,function(x){lapply(x,"[[",2)})
	  
	  fmlj=lapply(jacob,processFormula)
	  row=lapply(fmlj,function(x){lapply(x,"[[",1)})
	  col=lapply(fmlj,function(x){lapply(x,"[[",2)})
	  rhsj=lapply(fmlj,function(x){lapply(x,"[[",3)})
	  
	  for (i in 1:length(covariates)){
	    selected <- covariates[i]
	    get <- paste0("gsl_vector_get(co_variate, ", which(covariates == selected)-1,")")
	    rhs <- lapply(gsub(paste0("\\<",selected,"\\>"), get, rhs), function(x){eval(parse(text=x))})
	    rhsj <- lapply(gsub(paste0("\\<",selected,"\\>"), get, rhsj), function(x){eval(parse(text=x))})
	  }
	  
	  if (object@isContinuousTime){
	    #function_dx_dt
	    ret="void function_dx_dt(double t, size_t regime, const gsl_vector *x, double *param, size_t n_param, const gsl_vector *co_variate, gsl_vector *F_dx_dt){"
	    
	    if (nregime>1){
	      ret=paste(ret,"switch (regime) {",sep="\n\t")
	      for (r in 1:nregime){
	        ret=paste(ret,paste0("\tcase ",r-1,":"),sep="\n\t")
	        for (i in 1:n[r]){
	          for (j in 1:length(lhs[[r]])){
	            rhs[[r]][[i]]=gsub(paste0("\\<",lhs[[r]][[j]],"\\>"),paste0("gsl_vector_get(x,",j-1,")"),rhs[[r]][[i]])
	          }
	          ret=paste(ret,paste0("\tgsl_vector_set(F_dx_dt,",i-1,",",rhs[[r]][[i]],");"),sep="\n\t")    
	        }
	        ret=paste(ret,paste0("break;\n"),sep="\n\t")
	        
	      }
	      ret=paste(ret,paste0("\t}"),sep="\n\t")
	      
	    }else{
	      for (i in 1:n){
	        for (j in 1:length(lhs[[1]])){
	          rhs[[1]][[i]]=gsub(paste0("\\<",lhs[[1]][[j]],"\\>"),paste0("gsl_vector_get(x,",j-1,")"),rhs[[1]][[i]])
	        }
	        ret=paste(ret,paste0("\tgsl_vector_set(F_dx_dt,",i-1,",",rhs[[1]][[i]],");"),sep="\n\t")    
	      }
	    }
	    
	    ret=paste0(ret,"\n\t}")
	    
	    #function_dF_dx
	    ret=paste0(ret,"\n\n/**\n* The dF/dx function\n* The partial derivative of the jacobian of the DE function with respect to the variable x\n* @param param includes at the end the current state estimates in the same order as the states following the model parameters\n*/void function_dF_dx(double t, size_t regime, double *param, const gsl_vector *co_variate, gsl_matrix *F_dx_dt_dx){")
	    if (nregime>1){
	      ret=paste(ret,"switch (regime) {",sep="\n\t")
	      for (r in 1:nregime){
	        ret=paste(ret,paste0("case ",r-1,":"),sep="\n\t")
	        for (i in 1:length(jacob[[r]])){
	          for (j in 1:length(lhs[[r]])){
	            rhsj[[r]][[i]]=gsub(paste0("\\<",lhs[[r]][[j]],"\\>"),paste0("param[NUM_PARAM+",j-1,"]"),rhsj[[r]][[i]])
	          }
	          
	          ret=paste(ret,paste0("\tgsl_matrix_set(F_dx_dt_dx,",which(lhs[[r]]==row[[r]][[i]])-1,",",which(lhs[[r]]==col[[r]][[i]])-1,",",rhsj[[r]][[i]],");"),sep="\n\t")    
	        }
	        ret=paste(ret,paste0("break;\n"),sep="\n\t")
	        
	      }
	      ret=paste(ret,paste0("\t}"),sep="\n\t")
	      
	    }else{
	      for (i in 1:length(jacob[[1]])){
	        for (j in 1:length(lhs[[1]])){
	          rhsj[[1]][[i]]=gsub(paste0("\\<",lhs[[1]][[j]],"\\>"),paste0("param[NUM_PARAM+",j-1,"]"),rhsj[[1]][[i]])
	        }
	        
	        ret=paste(ret,paste0("\tgsl_matrix_set(F_dx_dt_dx,",which(unlist(lhs[[1]])==row[[1]][[i]])-1,",",which(unlist(lhs[[1]])==col[[1]][[i]])-1,",",rhsj[[1]][[i]],");"),sep="\n\t")    
	      }
	    }
	    
	    ret=paste0(ret,"\n\t}")
	    
	  }else{
	    #function_dynam
	    ret="void function_dynam(const double tstart, const double tend, size_t regime, const gsl_vector *xstart,\n\tdouble *param, size_t n_gparam,const gsl_vector *co_variate,\n\tvoid (*g)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),\n\tgsl_vector *x_tend){"
	    
	    if (nregime>1){
	      ret=paste(ret,"switch (regime) {",sep="\n\t")
	      for (r in 1:nregime){
	        ret=paste(ret,paste0("\tcase ",r-1,":"),sep="\n\t")
	        for (i in 1:n[r]){
	          for (j in 1:length(lhs[[r]])){
	            rhs[[r]][[i]]=gsub(paste0("\\<",lhs[[r]][[j]],"\\>"),paste0("gsl_vector_get(xstart,",j-1,")"),rhs[[r]][[i]])
	          }
	          ret=paste(ret,paste0("\tgsl_vector_set(x_tend,",i-1,",",rhs[[r]][[i]],");"),sep="\n\t")    
	        }
	        ret=paste(ret,paste0("break;\n"),sep="\n\t")
	        
	      }
	      ret=paste(ret,paste0("\t}"),sep="\n\t")
	      
	    }else{
	      for (i in 1:n){
	        for (j in 1:length(lhs[[1]])){
	          rhs[[1]][[i]]=gsub(paste0("\\<",lhs[[1]][[j]],"\\>"),paste0("gsl_vector_get(xstart,",j-1,")"),rhs[[1]][[i]])
	        }
	        ret=paste(ret,paste0("\tgsl_vector_set(x_tend,",i-1,",",rhs[[1]][[i]],");"),sep="\n\t")    
	      }
	    }
	    
	    ret=paste0(ret,"\n\t}")
	    
	    #function_jacob_dynam
	    ret=paste0(ret,"\n\nvoid function_jacob_dynam(const double tstart, const double tend, size_t regime, const gsl_vector *xstart,\n\tdouble *param, size_t num_func_param, const gsl_vector *co_variate,\n\tvoid (*g)(double, size_t, double *, const gsl_vector *, gsl_matrix *),\n\tgsl_matrix *Jx){")
	    if (nregime>1){
	      ret=paste(ret,"switch (regime) {",sep="\n\t")
	      for (r in 1:nregime){
	        ret=paste(ret,paste0("case ",r-1,":"),sep="\n\t")
	        for (i in 1:length(jacob[[r]])){
	          for (j in 1:length(lhs[[r]])){
	            rhsj[[r]][[i]]=gsub(paste0("\\<",lhs[[r]][[j]],"\\>"),paste0("gsl_vector_get(xstart,",j-1,")"),rhsj[[r]][[i]])
	          }
	          
	          ret=paste(ret,paste0("\tgsl_matrix_set(Jx,",which(lhs[[r]]==row[[r]][[i]])-1,",",which(lhs[[r]]==col[[r]][[i]])-1,",",rhsj[[r]][[i]],");"),sep="\n\t")    
	        }
	        ret=paste(ret,paste0("break;\n"),sep="\n\t")
	        
	      }
	      ret=paste(ret,paste0("\t}"),sep="\n\t")
	      
	    }else{
	      for (i in 1:length(jacob[[1]])){
	        for (j in 1:length(lhs[[1]])){
	          rhsj[[1]][[i]]=gsub(paste0("\\<",lhs[[1]][[j]],"\\>"),paste0("gsl_vector_get(xstart,",j-1,")"),rhsj[[1]][[i]])
	        }
	        
	        ret=paste(ret,paste0("\tgsl_matrix_set(Jx,",which(unlist(lhs[[1]])==row[[1]][[i]])-1,",",which(unlist(lhs[[1]])==col[[1]][[i]])-1,",",rhsj[[1]][[i]],");"),sep="\n\t")    
	      }
	    }
	    
	    ret=paste0(ret,"\n\t}")
	  }
	  object@c.string <- ret
	  return(object)
	}
)
#Example
# formula=list(
#   list(x1~a1*x1,
#        x2~a2*x2),
#   list(x1~a1*x1+c12*(exp(x2)/(1+exp(x2)))*x2,
#        x2~a2*x2+c21*(exp(x1)/(1+exp(x1)))*x1) 
# )
# dynam<-prep.formulaDynamics(formula=formula,startval=c(a1=.3,a2=.4,c12=-.5,c21=-.5),isContinuousTime=FALSE)
# cat(writeCcode(dynam)$c.string)

setMethod("writeCcode", "dynrDynamicsMatrix",
	function(object, covariates){
		isContinuousTime <- object$isContinuousTime
		params.dyn <- object$params.dyn
		values.dyn <- object$values.dyn
		params.exo <- object$params.exo
		values.exo <- object$values.exo
		params.int <- object$params.int
		values.int <- object$values.int
		
		nregime <- length(values.dyn)
		hasCovariates <- length(values.exo) > 0
		hasIntercepts <- length(values.int) > 0
		
		time <- ifelse(isContinuousTime, 'continuous', 'discrete')
		if(time == 'continuous'){
			# Construct matrices for A and B with A ~ dyn, B ~ exo
			# dx/dt ~ A %*% x + B %*% u
			# x is the latent state vector, u is the covariate vector
			# F_dx_dt = A %*% x + B %*% u
			# F_dx_dt_dx = A
			dynHead <- "void function_dx_dt(double t, size_t regime, const gsl_vector *x, double *param, size_t n_param, const gsl_vector *co_variate, gsl_vector *F_dx_dt){"
			jacHead <- "void function_dF_dx(double t, size_t regime, double *param, const gsl_vector *co_variate, gsl_matrix *F_dx_dt_dx){"
			inName <- "x"
			outName <- "F_dx_dt"
			jacName <- "F_dx_dt_dx"
		} else if(time == 'discrete'){
			# Construct matrices for A and B with A ~ dyn, B ~ exo
			# x[t] ~ A %*% x[t-1] + B %*% u
			# x is the latent state vector, u is the covariate vector
			dynHead <- "void function_dynam(const double tstart, const double tend, size_t regime, const gsl_vector *xstart, double *param, size_t n_gparam, const gsl_vector *co_variate, void (*g)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *), gsl_vector *x_tend){"
			jacHead <- "void function_jacob_dynam(const double tstart, const double tend, size_t regime, const gsl_vector *xstart, double *param, size_t num_func_param, const gsl_vector *co_variate, void (*g)(double, size_t, double *, const gsl_vector *, gsl_matrix *), gsl_matrix *Jx){"
			inName <- "xstart"
			outName <- "x_tend"
			jacName <- "Jx"
		}
		
		#If numRegimes > 1
		if(nregime > 1){
			# Create dynamics (state-transition or drift matrix) with covariate effects
			ret <- paste(dynHead,
				createGslMatrix(nrow(params.dyn[[1]]), ncol(params.dyn[[1]]), "Amatrix"),
				ifelse(hasCovariates, 
				       paste0(gslcovariate.front(object@covariates, covariates), createGslMatrix(nrow(params.exo[[1]]), ncol(params.exo[[1]]), "Bmatrix")),
				       ""),
				ifelse(hasIntercepts, createGslVector(nrow(params.int[[1]]), "intVector"), ""),
				sep="\n")
			
			#Loop through regimes for dynamics
			ret <- paste(ret, "\tswitch (regime) {\n", sep="\n")
			for(reg in 1:nregime){
				ret <- paste0(ret, paste0("\t\tcase ", reg-1, ":"), "\n")
				ret <- paste(ret,
					setGslMatrixElements(values=values.dyn[[reg]], params=params.dyn[[reg]], name="Amatrix"),
					ifelse(hasCovariates, paste0("\t\t\t", setGslMatrixElements(values=values.exo[[reg]], params=params.exo[[reg]], name="Bmatrix")), ""),
					ifelse(hasIntercepts, paste0("\t\t\t", setGslVectorElements(values=values.int[[reg]], params=params.int[[reg]], name="intVector")), ""),
					sep="\n")
				ret <- paste0(ret, "\t\t", "break;", "\n") 
			}
			ret <- paste0(ret, "\t", "}\n\n")
			
			ret <- paste(ret,
				blasMV(FALSE, "1.0", "Amatrix", inName, "0.0", outName),
				ifelse(hasCovariates, blasMV(FALSE, "1.0", "Bmatrix", "covariate_local", "1.0", outName), ""),
				ifelse(hasIntercepts, paste0("\tgsl_vector_add(", outName, ", intVector);\n"), ""),
				destroyGslMatrix("Amatrix"),
				ifelse(hasCovariates, paste0(destroyGslMatrix("Bmatrix"), destroyGslVector("covariate_local")), ""),
				ifelse(hasIntercepts, destroyGslMatrix("intVector"), ""),
				"}\n\n", sep="\n")
			
			
			# Create jacobian function
			ret <- paste(ret, jacHead)
			#Loop through regimes for jacobian
			ret <- paste(ret, "\tswitch (regime) {\n", sep="\n")
			for(reg in 1:nregime){
				ret <- paste0(ret, paste0("\t\tcase ", reg-1, ":"), "\n")
				ret <- paste(ret,
					setGslMatrixElements(values=values.dyn[[reg]], params=params.dyn[[reg]], name=jacName),
					sep="\n")
				ret <- paste0(ret, "\t\t", "break;", "\n") 
			}
			ret <- paste0(ret, "\t", "}\n\n")
			ret <- paste(ret, "}\n\n")
		} else { #Else If numRegimes == 1
			# Create dynamics (state-transition or drift matrix) with covariate effects
			ret <- paste(dynHead,
				createGslMatrix(nrow(params.dyn[[1]]), ncol(params.dyn[[1]]), "Amatrix"),
				setGslMatrixElements(values=values.dyn[[1]], params=params.dyn[[1]], name="Amatrix"),
				blasMV(FALSE, "1.0", "Amatrix", inName, "0.0", outName),
				destroyGslMatrix("Amatrix"),
				sep="\n")
			
			if(hasCovariates){
				ret <- paste(ret,
				  gslcovariate.front(object@covariates, covariates),
					createGslMatrix(nrow(params.exo[[1]]), ncol(params.exo[[1]]), "Bmatrix"),
					setGslMatrixElements(values=values.exo[[1]], params=params.exo[[1]], name="Bmatrix"),
					blasMV(FALSE, "1.0", "Bmatrix", "covariate_local", "1.0", outName),
					destroyGslMatrix("Bmatrix"),
					destroyGslVector("covariate_local"),
					sep="\n")
			}
			
			if(hasIntercepts){
				ret <- paste0(ret, "\n",
					createGslVector(nrow(params.int[[1]]), "intVector"), "\n",
					setGslVectorElements(values=values.int[[1]], params=params.int[[1]], name="intVector"), "\n",
					"\tgsl_vector_add(", outName, ", intVector);", "\n",
					destroyGslMatrix("intVector"))
			}
			
			ret <- paste(ret, "}\n\n", sep="\n")
			
			
			# Create jacobian function
			ret <- paste(ret,
				jacHead,
				setGslMatrixElements(values=values.dyn[[1]], params=params.dyn[[1]], name=jacName),
				"}\n\n",
				sep="\n")
		}
		object@c.string <- ret
		return(object)
	}
)


setMethod("writeCcode", "dynrRegimes",
	function(object, covariates){
		values <- object$values
		params <- object$params
		covariates_local <- object$covariates
		numCovariates <- length(covariates_local)
		numRegimes <- nrow(values)
		
		#Restructure values matrix for row-wise
		if(nrow(values)!=0 && nrow(params)!=0){
			values <- matrix(t(values), nrow=numRegimes*numRegimes, ncol=numCovariates+1, byrow=TRUE)
			params <- matrix(t(params), nrow=numRegimes*numRegimes, ncol=numCovariates+1, byrow=TRUE)
			rowBeginSeq <- seq(1, nrow(values), by=numRegimes)
			rowEndSeq <- seq(numRegimes, nrow(values), by=numRegimes)
		}
		
		ret <- "void function_regime_switch(size_t t, size_t type, double *param, const gsl_vector *co_variate, gsl_matrix *regime_switch_mat){"
		if(prod(nrow(values), nrow(params), numCovariates) != 0){
			ret <- paste(ret,
			             gslcovariate.front(covariates_local, covariates),
			             createGslMatrix(numRegimes, numCovariates, "Gmatrix"),
			             createGslVector(numRegimes, "Pvector"),
			             createGslVector(numRegimes, "Presult"),
			             sep="\n")
			for(reg in 1L:numRegimes){
				selRows <- rowBeginSeq[reg]:rowEndSeq[reg]
				ret <- paste(ret,
					setGslVectorElements(values=values[selRows, 1, drop=FALSE], params=params[selRows, 1, drop=FALSE], name="Pvector"),
					setGslMatrixElements(values=values[selRows, -1, drop=FALSE], params=params[selRows, -1, drop=FALSE], name="Gmatrix"),
					blasMV(FALSE, "1.0", "Gmatrix", "covariate_local", "1.0", "Pvector"),
					"\tmathfunction_softmax(Pvector, Presult);",
					gslVector2Column("regime_switch_mat", reg-1, "Presult", 'row'),
					"\tgsl_vector_set_zero(Pvector);",
					"\tgsl_matrix_set_zero(Gmatrix);",
					sep="\n")
			}
			ret <- paste(ret,
				destroyGslMatrix("Gmatrix"),
				destroyGslVector("Pvector"),
				destroyGslVector("Presult"),
				destroyGslVector("covariate_local"),
				sep="\n")
		} else if(nrow(values) != 0 && nrow(params) != 0 && numCovariates == 0){
			# same as above processing, but without the Gmatrix for the covariates
			ret <- paste(ret,
				createGslVector(numRegimes, "Pvector"),
				createGslVector(numRegimes, "Presult"),
				sep="\n")
			for(reg in 1L:numRegimes){
				selRows <- rowBeginSeq[reg]:rowEndSeq[reg]
				ret <- paste(ret,
					setGslVectorElements(values=values[selRows, 1, drop=FALSE], params=params[selRows, 1, drop=FALSE], name="Pvector"),
					"\tmathfunction_softmax(Pvector, Presult);",
					gslVector2Column("regime_switch_mat", reg-1, "Presult", 'row'),
					"\tgsl_vector_set_zero(Pvector);",
					sep="\n")
			}
			ret <- paste(ret,
				destroyGslVector("Pvector"),
				destroyGslVector("Presult"),
				sep="\n")
		} else {
			ret <- paste(ret, "\tgsl_matrix_set_identity(regime_switch_mat);", sep="\n")
		}
		ret <- paste(ret, "}\n\n", sep="\n")
		object@c.string <- ret
		return(object)
	}
)


setMethod("writeCcode", "dynrInitial",
	function(object, covariates){
		values.inistate <- object$values.inistate
		params.inistate <- object$params.inistate
		values.inicov <- object@values.inicov.inv.ldl
		params.inicov <- object$params.inicov
		values.regimep <- object$values.regimep
		params.regimep <- object$params.regimep
		
		nregime <- max(length(values.inistate), length(values.inicov))
		someStatesNotZero <- sapply(values.inistate, function(x){!all(x == 0)})
		someStatesNotZero2 <- sapply(params.inistate, function(x){!all(x == 0)})
		someStatesNotZero <- someStatesNotZero | someStatesNotZero2
		
		ret <- "void function_initial_condition(double *param, gsl_vector **co_variate, gsl_vector *pr_0, gsl_vector **eta_0, gsl_matrix **error_cov_0){\n"
		ret <- paste(ret, setGslVectorElements(values.regimep,params.regimep, "pr_0"), sep="\n")
		ret <- paste0(ret,"\tsize_t num_regime=pr_0->size;\n\tsize_t dim_latent_var=error_cov_0[0]->size1;")
		if (any(someStatesNotZero)){
			ret <- paste0(ret,"\n\tsize_t num_sbj=(eta_0[0]->size)/(dim_latent_var);\n\tsize_t i;")
		}
		ret <- paste0(ret,"\n\tsize_t j;\n")
		
		if(nregime > 1){
			ret <- paste(ret, "\tswitch (regime) {\n", sep="\n")
			for(reg in 1:nregime){
				ret <- paste0(ret, paste0("\t\tcase ", reg-1, ":"), "\n")
				ret <- paste0(ret, "\t\t\tfor(j=0;j<num_regime;j++){")
				if (any(someStatesNotZero[reg])){
					ret <- paste0(ret,"\n\t\t\t\tfor(i=0;i<num_sbj;i++){\n")
					ret <- paste0(ret, setGslVectorElements(values=values.inistate[[reg]], params=params.inistate[[reg]], name='(eta_0)[j]', fill="i*dim_latent_var+", depth=5))
					ret <- paste0(ret,"\t\t\t\t}") # close i loop
				}
				ret <- paste(ret, setGslMatrixElements(values.inicov[[reg]], params.inicov[[reg]], "(error_cov_0)[j]", depth=4), sep="\n")
				ret <- paste0(ret, "\t\t\t}\n") # close j loop
				ret <- paste0(ret, "\t\t", "break;", "\n")
			}
			ret <- paste0(ret, "\t", "}\n\n") # close case switch
		} else{
			ret <- paste0(ret, "\tfor(j=0;j<num_regime;j++){")
			if (any(someStatesNotZero)){
				ret <- paste0(ret,"\n\t\tfor(i=0;i<num_sbj;i++){\n")
				ret <- paste0(ret, setGslVectorElements(values=values.inistate[[1]], params=params.inistate[[1]], name='(eta_0)[j]', fill="i*dim_latent_var+", depth=3))
				ret <- paste0(ret,"\t\t}") # close i loop
			}
			ret <- paste(ret, setGslMatrixElements(values.inicov[[1]], params.inicov[[1]], "(error_cov_0)[j]"), sep="\n")
			ret <- paste0(ret, "\t}") # close j loop
		}
		
		ret <- paste0(ret, "\n}\n") #Close function definition
		object@c.string <- ret
		
		return(object)
	}
)

setMethod("writeCcode", "dynrNoise",
	function(object, covariates){
		params.latent <- object$params.latent
		params.observed <- object$params.observed
		#Note: use mutated values, instead of raw values
		values.latent <- object@values.latent.inv.ldl
		values.observed <- object@values.observed.inv.ldl
		
		nregime <- length(values.latent)
		
		ret <- "void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov){\n\n"
		if(nregime > 1){
			ret <- paste(ret, "\tswitch (regime) {\n", sep="\n")
			for(reg in 1:nregime){
				ret <- paste0(ret, paste0("\t\tcase ", reg-1, ":"), "\n")
				ret <- paste(ret, setGslMatrixElements(values.latent[[reg]], params.latent[[reg]], "eta_noise_cov", depth=3), sep="\n")
				ret <- paste(ret, setGslMatrixElements(values.observed[[reg]], params.observed[[reg]], "y_noise_cov", depth=3), sep="\n")
				ret <- paste0(ret, "\t\t", "break;", "\n")
			} # close for loop
			ret <- paste0(ret, "\t", "}\n\n") # close case switch
		} else {
			ret <- paste(ret, setGslMatrixElements(values.latent[[1]], params.latent[[1]], "eta_noise_cov"), sep="\n")
			ret <- paste(ret, setGslMatrixElements(values.observed[[1]], params.observed[[1]], "y_noise_cov"), sep="\n")
		}
		ret <- paste(ret, "\n}\n\n") # close C function
		
		object@c.string <- ret
		
		return(object)
	}
)

setMethod("writeCcode", "dynrTrans",
          function(object, covariates){
            #function_transform
            ret="/**\n * This function modifies some of the parameters so that it satisfies the model constraint.\n * Do not include parameters in noise_cov matrices \n */\nvoid function_transform(double *param){"
            n=length(object$formula.trans)
            if(n == 0){
              object@c.string <- paste0(ret, "\n}\n\n")
              return(object)
            }else if (object@transCcode){
              formula.trans=object$formula.trans
              fml=processFormula(formula.trans)
              lhs=lapply(fml,function(x){x[1]})
              rhs=lapply(fml,function(x){x[2]})
              
              for (i in 1:n){
                ret=paste(ret,paste0(lhs[i],"=",rhs[i],";"),sep="\n\t") 
              }
            }
            
            ret=paste0(ret,"\n\t}\n")
            object@c.string <- ret
            return(object)
          }
)

setGeneric("createRfun", function(object, param.data, params.observed, params.latent, params.inicov,
                                  values.observed, values.latent, values.inicov,
                                  values.observed.orig, values.latent.orig, values.inicov.orig, show=TRUE) { 
  return(standardGeneric("createRfun")) 
})

setMethod("createRfun", "dynrTrans",
          function(object, param.data, params.observed, params.latent, params.inicov,
                   values.observed, values.latent, values.inicov,
                   values.observed.orig, values.latent.orig, values.inicov.orig){
            #inv.tfun
            inv.tf <- NULL
            tf <- NULL
            f.string<-"inv.tf<-function(vec){"
            if (length(object@formula.inv) > 0){
              #TODO If formula.inv is missing (i.e. length == 0), point-wise inverse functions that are based on the uniroot function will be used to calculate starting values.
              #inverse = function (f, lower = -100, upper = 100) {function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)$root}
              #eval(parse(text="inv.tf<-function(namedvec){return(c(a=inverse(exp)(namedvec[\"a\"]),b=inverse(function(x)x^2,0.0001,100)(namedvec[\"b\"])))}"))
              fml.str=formula2string(object@formula.inv)
              lhs=fml.str$lhs
              rhs=fml.str$rhs
              sub=sapply(lhs,function(x){paste0("vec[",param.data$param.number[param.data$param.name==x],"]")})
              for (j in 1:length(rhs)){
                for (i in 1:length(lhs)){
                  rhs[j]=gsub(paste0("\\<",lhs[i],"\\>"),sub[i],rhs[j])
                }
                eq=paste0(sub[j],"=",rhs[j])
                f.string<-paste(f.string,eq,sep="\n\t")
              }
            }
            f.string2 <- f.string
            #observed
            if (sum(param.data$ldl.observed)>0){
            f.string2 <- paste(f.string2, makeldlchar(param.data, lapply(values.observed.orig, replaceDiagZero), params.observed, reverse=TRUE),sep="\n\t")
            }
            #latent
            if (sum(param.data$ldl.latent)>0){
            f.string2 <- paste(f.string2, makeldlchar(param.data, lapply(values.latent.orig, replaceDiagZero), params.latent, reverse=TRUE),sep="\n\t")
            }
            #inicov
            if (sum(param.data$ldl.inicov)>0){
            f.string2 <- paste(f.string2, makeldlchar(param.data, lapply(values.inicov.orig, replaceDiagZero), params.inicov, reverse=TRUE),sep="\n\t")
            }
            f.string2<-paste0(f.string2,"\n\treturn(vec)}")
            f.string<-paste0(f.string,"\n\treturn(vec)}")
            
            eval(parse(text=f.string2))
            object@inv.tfun.full <- inv.tf
            
            eval(parse(text=f.string)) 
            object@inv.tfun <- inv.tf
            
            f.string<-"tf<-function(vec){"
              
              if (length(object@formula.trans)>0){
              fml.str=formula2string(object@formula.trans)
              lhs=fml.str$lhs
              rhs=fml.str$rhs
              sub=sapply(lhs,function(x){paste0("vec[",param.data$param.number[param.data$param.name==x],"]")})
              
              for (j in 1:length(rhs)){
                for (i in 1:length(lhs)){
                  rhs[j]=gsub(paste0("\\<",lhs[i],"\\>"),sub[i],rhs[j])
                }
                eq=paste0(sub[j],"=",rhs[j])
                f.string<-paste(f.string,eq,sep="\n\t")
              }
              object@paramnames<-lhs
              }
            
            #observed
            if (sum(param.data$ldl.observed)>0){
              f.string<-paste(f.string, makeldlchar(param.data, values.observed, params.observed),sep="\n\t")
            }
            #latent
            if (sum(param.data$ldl.latent)>0){
              f.string<-paste(f.string, makeldlchar(param.data, values.latent, params.latent),sep="\n\t")
            }
            #inicov
            if (sum(param.data$ldl.inicov)>0){
              f.string<-paste(f.string, makeldlchar(param.data, values.inicov, params.inicov),sep="\n\t")
            }
              
            f.string<-paste0(f.string,"\n\treturn(vec)}")
            eval(parse(text=f.string))   
            object@tfun <- tf
            
            return(object)
          }
)

makeldlchar<-function(param.data, values, params, reverse=FALSE){
  #Note: Because we need to gurantee the uniqueness of the parameter values in estimation after applying reverseldl transformations, 
  #there should be no intersections between the parameter names used in different matrices that are subject to the ldl transformations
  #See below for an example where the (2,2) parameters are the same in two original matrices but the reverseldl-ed parameters differ.
  #This is problematic in estimation. We should treat the parameters in matrix(c(2,1,2,5) and matrix(c(3,1,2,5) as two different sets of parameters.
  #Also, the noise structure can only be either fully freed or with diagonals freed.
  #   > dynr:::reverseldl(matrix(c(2,1,2,5),ncol=2))
  #   [,1]     [,2]
  #   [1,] 0.6931472 0.000000
  #   [2,] 1.0000000 1.098612
  #   > dynr:::reverseldl(matrix(c(3,1,2,5),ncol=2))
  #   [,1]     [,2]
  #   [1,] 1.0986123 0.000000
  #   [2,] 0.6666667 1.299283
  params.numbers <- lapply(params, .exchangeNamesAndNumbers, param.data$param.name)
  values <- PopBackMatrix(values, params.numbers, paste0("vec[",param.data$param.number,"]"))
  char<-character(0)
  if (length(values)>0){
    for (i in 1:length(values)){
      param.number=unique(params.numbers[[i]][which(params.numbers[[i]]!=0,arr.ind = TRUE)])
      if (length(param.number)!=0){
        vec.noise=paste0("vec[",paste0("c(",paste(param.number,collapse=","),")"),"]")
        vec.sub <- as.vector(values[[i]])
        mat.index=sapply(param.number,function(x){min(which(params.numbers[[i]]==x))})
        char.i=paste0(vec.noise,"=as.vector(", ifelse(reverse, 'reverseldl', 'transldl'), "(matrix(", paste0("c(",paste(vec.sub,collapse=","),")"),",ncol=",ncol(values[[i]]),")))[",
                      paste0("c(",paste(mat.index,collapse=","),")"),"]", ifelse(i!=length(values),"\n\t",""))
        char=paste0(char, char.i)
      }
    }
  }
  return(char)
}

vec2mat<-function(vectr,dimension){
  n.vec=length(vectr)
  if (n.vec==dimension){
    #diagonal
    mat=diag(vectr)
  }else if (dimension==sqrt(2*n.vec+.25) - .5){
    #symmetric
    mat=matrix(0,dimension,dimension)
    diag(mat)<-vectr[1:dimension]
    for (i in 1:(dimension-1)){
      for (j in (i+1):dimension){
        mat[i,j]<-vectr[i+j+dimension-2]
        mat[j,i]<-vectr[i+j+dimension-2]
      }
    }
  }else{
    cat('Length of the vector does not match the matrix dimension!\n')
  }	
  return(mat)	
}

mat2vec<-function(mat,n.vec){
  dimension=dim(mat)[1]
  vec=numeric(n.vec)
  if (n.vec==dimension){
    #diagonal
    vec=diag(mat)
  }else if (dimension==sqrt(2*n.vec+.25) - .5){
    #symmetric
    vec[1:dimension]<-diag(mat)
    for (i in 1:(dimension-1)){
      for (j in (i+1):dimension){
        vec[i+j+dimension-2]<-mat[i,j]
      }
    }
  }else{
    cat('Length of the vector does not match the matrix dimension!\n')
  }	
  return(vec)	
}
#transldl function for caluaclating the LDL values
transldl <- function(mat){
  L <- mat
  diag(L) <- 1
  L[upper.tri(L)] <- 0
  D <- diag(exp(diag(mat)), nrow=nrow(mat))
  # final caluclation
  outldl <- L %*% D %*% t(L)
  return(outldl)
}


##' LDL Decomposition for Matrices
##' @param x a numeric matrix
##'
##' This is a wrapper function around the \code{\link{chol}} function.
##' The goal is to factor a square, symmetric, positive (semi-)definite matrix into the product of a lower triangular matrix, a diagonal matrix, and the transpose of the lower triangular matrix.
##' The value returned is a lower triangular matrix with the elements of D on the diagonal.
dynr.ldl <- function(x){
	ret <- t(chol(x))
	d <- diag(ret)
	ret <- ret %*% diag(1/d, nrow=length(d))
	diag(ret) <- d^2
	return(ret)
}

reverseldl <- function(values){
	if(dim(values)[1]==1){
		return(log(values))
	} else if(any(is.na(values))){
		warning("Avast ye swarthy dog! NA was passed to LDL. Unset bounds might be fine. Values might be wrong.")
		return(values)
	} else{
		mat <- dynr.ldl(values)
		diag(mat) <- log(diag(mat))
		return(mat)
	}
}

#------------------------------------------------------------------------------
# Some usefull helper functions
#

##' Create a diagonal matrix from a character vector
##' @aliases diag.character diag
##' 
##' @param x Character vector used to create the matrix
##' @param nrow Numeric. Number of rows for the resulting matrix.
##' @param ncol Numeric. Number of columns for the resulting matrix.
##' 
##' @details
##' The default behavior for missing \code{nrow} and/or \code{ncol} arguments is the same
##' as for the \code{\link{diag}} function in the base package.  Off-diagonal entries
##' are filled with "0".
##' 
##' @examples
##' diag(letters[1:3])
# a new method for diag with character input
setMethod("diag", "character", #diag.character <-
	function(x=1, nrow, ncol){
		n <- length(x)
		if(!missing(nrow)) n <- nrow
		if(missing(ncol)) ncol <- n
		r <- matrix("0", n, ncol)
		diag(r) <- x
		return(r)
	}
)


# allow free/fixed to be used in values/params arguments
preProcessValues <- function(x){
	x[is.na(x)] <- 'freed'
	if(is.null(dim(x))){
		numRow <- length(x)
		numCol <- 1
		rowNam <- NULL
		colNam <- NULL
	} else {
		numRow <- nrow(x)
		numCol <- ncol(x)
		rowNam <- rownames(x)
		colNam <- colnames(x)
	}
	x <- c(x)
	xl <- tolower(x)
	sel <- pmatch(xl, "freed", duplicates.ok=TRUE)
	notAnumber <- is.na(sapply(x, function(x){suppressWarnings(try(as.numeric(x), silent=TRUE))}))
	if(any( pchars <- (is.na(sel) & notAnumber))){
		msg <- paste0("The following are not numbers or some version of 'freed': ", paste(x[pchars], collapse=', '), "\nEnter numerical starting values for free parameters, and 0 or other meaningful values for fixed parameters.")
		stop(msg)
	}
	# Set starting values to zero for 'free'
	x[sel %in% 1] <- 0
	x <- matrix(as.numeric(x), numRow, numCol, dimnames=list(rowNam, colNam))
	return(x)
}

preProcessParams <- function(x){
	x[is.na(x)] <- 'fixed'
	if(is.null(dim(x))){
		numRow <- length(x)
		numCol <- 1
		rowNam <- NULL
		colNam <- NULL
	} else {
		numRow <- nrow(x)
		numCol <- ncol(x)
		rowNam <- rownames(x)
		colNam <- colnames(x)
	}
	x <- c(x)
	xl <- tolower(x)
	sel <- xl %in% c("0", "fix", "fixed")
	x[sel] <- "fixed"
	validName <- sapply(x, function(x){x == make.names(x)})
	if(any(!validName)){
		msg <- paste0("Whoa nelly! The following are not valid free parameter names: ", paste(x[!validName], collapse=", "), ".\n  Valid names are things that could be variable names in R.  See ?make.names")
		stop(msg)
	}
	x <- matrix(as.character(x), numRow, numCol, dimnames=list(rowNam, colNam))
	return(x)
}

coProcessValuesParams <- function(values=NULL, params=NULL, missingOK=FALSE){
	if(is.null(values) || missing(values)){
		if(missingOK){
			values <- list()
			params <- list()
		} else {
			stop("values missing and it's not okay.")
		}
	}
	if(!is.list(values)){
		values <- list(values)
	}
	if(is.null(params) || missing(params)){
		params <- rep(list(matrix(0, nrow(values[[1]]), ncol(values[[1]]))), length(values))
	}
	if(!is.list(params)){
		params <- list(params)
	}
	if(length(values) != length(params)){
		stop(paste0("Mismatch between values and params.  'values' argument indicates ", length(values), " regimes but 'params' argument indicates ", length(params), " regimes.  Get your mind right."))
	}
	return(list(values=values, params=params))
}


preProcessNames <- function(x,rnames=character(0),cnames=character(0)){
  rownames(x) <- rnames
  if (!length(cnames)==0) colnames(x) <- cnames
  return(x)
}

extractWhichParams <- function(p){
	p!="fixed" & !duplicated(p, MARGIN=0) & (p!=0)
}

extractParams <- function(p){
	if(is.list(p) && length(p) > 0){
		ret <- c()
		for(i in 1:length(p)){
			ret <- c(ret, extractParams(p[[i]]))
		}
		return(ret)
	} else if(length(p) == 0){
		return(character(0))
	}else {
		return(p[extractWhichParams(p)])
	}
}

extractValues <- function(v, p){
	if(is.list(v) && is.list(p) && length(v) == length(p) && length(v) > 0){
		ret <- c()
		for(i in 1:length(v)){
			ret <- c(ret, extractValues(v[[i]], p[[i]]))
		}
		return(ret)
	} else if(length(v) == 0){
		return(numeric(0))
	}else {
		return(v[extractWhichParams(p)])
	}
}

#------------------------------------------------------------------------------
# Create recipe for function measurement

# TODO add ability to use covariates in these functions

#--------------------------------------
# brief input version

##' Recipe function to quickly create factor loadings
##'
##' @param map list giving how the latent variables map onto the observed variables
##' @param params parameter numbers
##' @param idvar Names of the variables used to identify the factors
##' @param exo.names names of the exogenous covariates
##' 
##' @details
##' The default pattern for 'idvar' is to fix the first factor loading 
##' for each factor to one.  The variable names listed in 'idvar' have 
##' their factor loadings fixed to one.  However, if the names of the 
##' latent variables are used for 'idvar', then all the factor loadings
##' will be freely estimated and you should fix the factor variances 
##' in the noise part of the model (e.g. \code{\link{prep.noise}}).
##'
##' This function does not have the full set of features possible in 
##' the dynr package. In particular, it does not have regime-swtiching 
##' factor loadings, any intercepts, or any covariates.  For complete 
##' functionality use \code{\link{prep.measurement}}.
##' 
##' @examples
##' #Single factor model with one latent variable fixing first loading
##' prep.loadings(list(eta1=paste0('y', 1:4)), paste0("lambda_",2:4))
##'
##' #Single factor model with one latent variable fixing the fourth loading
##' prep.loadings(list(eta1=paste0('y', 1:4)), paste0("lambda_",1:3), idvar='y4')
##' 
##' #Single factor model with one latent variable freeing all loadings
##' prep.loadings(list(eta1=paste0('y', 1:4)), paste0("lambda_", 1:4), idvar='eta1')
##' 
##' # Two factor model with simple structure
##' prep.loadings(list(eta1=paste0('y', 1:4), eta2=paste0('y', 5:7)), 
##' paste0("lambda_", c(2:4, 6:7)))
##' 
##' #Two factor model with repeated use of a free parameter
##' prep.loadings(list(eta1=paste0('y', 1:4), eta2=paste0('y', 5:8)), 
##' paste0("lambda_", c(2:4, 6:7, 4)))
##' 
##' #Two factor model with a cross loading
##' prep.loadings(list(eta1=paste0('y', 1:4), eta2=c('y5', 'y2', 'y6')), 
##' paste0("lambda_", c("21", "31", "41", "22", "62")))
prep.loadings <- function(map, params, idvar, exo.names=character(0)){
	if(missing(idvar)){
		idvar <- sapply(map, '[', 1)
	}
	
	allVars <- unique(unlist(map))
	
	nx <- length(allVars)
	ne <- length(map)
	
	if(!all(idvar %in% c(names(map), unlist(map)))){
		stop("The 'idvar' must all be either in the names of the 'map' argument or parts of the 'map' argument.")
	}
	
	paramsNeeded <- length(unlist(map)) - sum(idvar %in% allVars)
	if(length(params) != paramsNeeded){
		stop(paste0("Number of free parameters provided (", length(params), ") does not match number needed (", paramsNeeded, ")."))
	}
	
	
	k <- 1
	valuesMat <- matrix(0, nx, ne)
	paramsMat <- matrix(0, nx, ne)
	for(j in 1:ne){
		for(i in 1:nx){
			if(allVars[i] %in% map[[j]] & !(allVars[i] %in% idvar)){
				paramsMat[i, j] <- params[k]
				valuesMat[i, j] <- .8
				k <- k+1
			} else if(allVars[i] %in% map[[j]] & allVars[i] %in% idvar){
				valuesMat[i, j] <- 1
			}
		}
	}
	rownames(valuesMat) <- allVars
	rownames(paramsMat) <- allVars
	colnames(valuesMat) <- names(map)
	colnames(paramsMat) <- names(map)
	x <- prep.measurement(values.load=valuesMat, params.load=paramsMat,state.names=names(map),obs.names=allVars,exo.names=exo.names)
	return(x)
}



#--------------------------------------
# matrix input version

# values, and params are all MxN matrices
# a zero param is taken to be fixed.

##' Prepare the measurement recipe
##'
##' @param values.load matrix of starting or fixed values for factor loadings. 
##' For models with regime-specific factor loadings provide a list of matrices of factor loadings.
##' @param params.load matrix or list of matrices. Contains parameter names of the factor loadings. 
##' @param values.exo matrix or list of matrices. Contains starting/fixed values of the covariate regression slopes.
##' @param params.exo matrix or list of matrices. Parameter names of the covariate regression slopes.
##' @param values.int vector of intercept values specified as matrix or list of matrices. Contains starting/fixed values of the intercepts.
##' @param params.int vector of names for intercept parameters specified as a matrix or list of matrices. 
##' @param obs.names  vector of names for the observed variables in the order they appear in the measurement model.
##' @param state.names  vector of names for the latent variables in the order they appear in the measurement model.
##' @param exo.names  (optional) vector of names for the exogenous variables in the order they appear in the measurement model.
##'
##' The values.* arguments give the starting and fixed values for their respective matrices.
##' The params.* arguments give the free parameter labels for their respective matrices.
##' Numbers can be used as labels.
##' The number 0 and the character 'fixed' are reserved for fixed parameters.
##'
##' When a single matrix is given to values.*, that matrix is not regime-switching.
##' Correspondingly, when a list of length r is given, that matrix is regime-switching with values and params for the r regimes in the elements of the list.
##' 
##' @examples
##' prep.measurement(diag(1, 5), diag("lambda", 5))
##' prep.measurement(matrix(1, 5, 5), diag(paste0("lambda_", 1:5)))
##' prep.measurement(diag(1, 5), diag(0, 5)) #identity measurement model
##' 
##' #Regime-switching measurement model where the first latent variable is
##' # active for regime 1, and the second latent variable is active for regime 2
##' # No free parameters are present.
##' prep.measurement(values.load=list(matrix(c(1,0), 1, 2), matrix(c(0, 1), 1, 2)))
prep.measurement <- function(values.load, params.load=NULL, values.exo=NULL, params.exo=NULL, values.int=NULL, params.int=NULL,
                             obs.names, state.names, exo.names){
	# Handle load
	r <- coProcessValuesParams(values.load, params.load)
	values.load <- r$values
	params.load <- r$params
	
	# Handle exo
	r <- coProcessValuesParams(values.exo, params.exo, missingOK=TRUE)
	values.exo <- r$values
	params.exo <- r$params
	
	# Handle int
	r <- coProcessValuesParams(values.int, params.int, missingOK=TRUE)
	values.int <- r$values
	params.int <- r$params
	
	
	if(missing(obs.names)){
		obs.names = paste0('y',1:nrow(values.load[[1]]))
	}
	
	if(missing(state.names)){
		state.names = paste0('state',1:ncol(values.load[[1]]))
	}
	
	if(missing(exo.names)){
		exo.names = character(0)
	}
	
	values.load <- lapply(values.load, preProcessValues)
	params.load <- lapply(params.load, preProcessParams)
	values.exo <- lapply(values.exo, preProcessValues)
	params.exo <- lapply(params.exo, preProcessParams)
	values.int <- lapply(values.int, preProcessValues)
	params.int <- lapply(params.int, preProcessParams)
	
	values.load <- lapply(values.load, preProcessNames,obs.names,state.names)
	#params.load <- lapply(params.load, preProcessNames,obs.names,state.names)
	values.exo <- lapply(values.exo, preProcessNames,obs.names,exo.names)
	#params.exo <- lapply(params.exo, preProcessNames,obs.names,exo.names)
	values.int <- lapply(values.int, preProcessNames,obs.names)
	#params.int <- lapply(params.int, preProcessNames,obs.names)
	
	sv <- c(extractValues(values.load, params.load), extractValues(values.exo, params.exo), extractValues(values.int, params.int))
	pn <- c(extractParams(params.load), extractParams(params.exo), extractParams(params.int))
	sv <- extractValues(sv, pn)
	pn <- extractParams(pn)
	x <- list(startval=sv, paramnames=pn, values.load=values.load, params.load=params.load,
		values.exo=values.exo, params.exo=params.exo, values.int=values.int, params.int=params.int,
		obs.names=obs.names, state.names=state.names, exo.names=exo.names)
	return(new("dynrMeasurement", x))
}


#------------------------------------------------------------------------------
# Error covariance matrix
# N.B. This function produces BOTH the latent and observed error covariance matrices.


##' Recipe function for specifying the measurement error and process noise covariance structures
##' 
##' @param values.latent a positive definite matrix or a list of positive definite matrices of the starting or fixed values of the process noise covariance structure(s) in one or more regimes. If only one matrix is specified for a regime-switching dynamic model, the process noise covariance structure stays the same across regimes. To ensure the matrix is positive definite in estimation, we apply LDL transformation to the matrix. Values are hence automatically adjusted for this purpose.
##' @param params.latent a matrix or list of matrices of the parameter names that appear in the process noise covariance(s) in one or more regimes. If an element is 0 or "fixed", the corresponding element is fixed at the value specified in the values matrix; Otherwise, the corresponding element is to be estimated with the starting value specified in the values matrix. If only one matrix is specified for a regime-switching dynamic model, the process noise structure stays the same across regimes. If a list is specified, any two sets of the parameter names as in two matrices should be either the same or totally different to ensure proper parameter estimation.
##' @param values.observed a positive definite matrix or a list of positive definite matrices of the starting or fixed values of the measurement error covariance structure(s) in one or more regimes. If only one matrix is specified for a regime-switching measurement model, the measurement noise covariance structure stays the same across regimes. To ensure the matrix is positive definite in estimation, we apply LDL transformation to the matrix. Values are hence automatically adjusted for this purpose. 
##' @param params.observed a matrix or list of matrices of the parameter names that appear in the measurement error covariance(s) in one or more regimes. If an element is 0 or "fixed", the corresponding element is fixed at the value specified in the values matrix; Otherwise, the corresponding element is to be estimated with the starting value specified in the values matrix. If only one matrix is specified for a regime-switching dynamic model, the process noise structure stays the same across regimes. If a list is specified, any two sets of the parameter names as in two matrices should be either the same or totally different to ensure proper parameter estimation.
##' 
##' @examples 
##' prep.noise(values.latent=diag(c(0.8, 1)), params.latent=diag(c('fixed', "e_x")), 
##' values.observed=diag(1.5,1), params.observed=diag("e_y", 1))
prep.noise <- function(values.latent, params.latent, values.observed, params.observed){
	# Handle latent covariance
	r <- coProcessValuesParams(values.latent, params.latent)
	values.latent <- r$values
	params.latent <- r$params
	
	# Handle observed covariance
	r <- coProcessValuesParams(values.observed, params.observed)
	values.observed <- r$values
	params.observed <- r$params
	
	if(length(values.observed) != length(values.latent)){
		msg <- paste0("Stuff and nonsense: number of regimes for observed covariance (", length(values.observed),
			") is not equal to that for the latent covariance(", length(values.latent), ").")
		stop(msg)
	}
	
	values.latent <- lapply(values.latent, preProcessValues)
	params.latent <- lapply(params.latent, preProcessParams)
	values.observed <- lapply(values.observed, preProcessValues)
	params.observed <- lapply(params.observed, preProcessParams)
	
	values.latent.inv.ldl <- lapply(values.latent, replaceDiagZero)
	values.latent.inv.ldl <- lapply(values.latent.inv.ldl, reverseldl)
	values.observed.inv.ldl <- lapply(values.observed, replaceDiagZero)
	values.observed.inv.ldl <- lapply(values.observed.inv.ldl, reverseldl)
	
	sv <- c(extractValues(values.latent.inv.ldl, params.latent), extractValues(values.observed.inv.ldl, params.observed))
	pn <- c(extractParams(params.latent), extractParams(params.observed))
	sv <- extractValues(sv, pn)
	pn <- extractParams(pn)
	x <- list(startval=sv, paramnames=pn, values.latent=values.latent, values.observed=values.observed, params.latent=params.latent, params.observed=params.observed, values.latent.inv.ldl=values.latent.inv.ldl, values.observed.inv.ldl=values.observed.inv.ldl)
	return(new("dynrNoise", x))
}


replaceDiagZero <- function(x){
	diag(x)[diag(x) == 0] <- 1e-6
	return(x)
}


#------------------------------------------------------------------------------
# Regime switching matrix/function

# For regime switching function, each element of the transition probability
# matrix should be the result of a multinomial logistic regression.
# That is, for a 2x2 matrix we have
# Behind the scenes it does this
#   p11 ~ softmax(1 + x1 + x2 + ... + xn)
#   p21 ~ softmax(1 + x1 + x2 + ... + xn)
#   p12 ~ softmax(1 + x1 + x2 + ... + xn)
#   p22 ~ softmax(1 + x1 + x2 + ... + xn)
#
# Users specify this
#   p11 ~ 1 + x1 + x2
#   p21 ~ 1 + x2 + x3

# Perhaps have two regime functions
# One is analogous to the linear dynamics function, and has more limited interface
# The other may be more general.
# Even the limited one should be able to handle covariates as above via multinomial logistic regression

# add identification constraints
# 

# Each element of the transition probability matrix (TPM) has a linear predictor (lp).
# LPM = 
# lp(p11) ~ 1 + x1 + x2 + ... + xn,   lp(p12) ~ 1 + x1 + x2 + ... + xn
# lp(p21) ~ 1 + x1 + x2 + ... + xn,   lp(p22) ~ 1 + x1 + x2 + ... + xn
# 
# The softmax function then operates on the columns of the linear predictor matrix (LPM).
# TPM = 
#       cbind( softmax(LPM[,1]), softmax(LPM[,2]), ..., softmax(LPM[,k]) )
# for k regimes

##' Recipe function for creating regime switching (Markov transition) functions
##' 
##' @param values matrix giving the values. Should have (number of Regimes) rows 
##' and (number of regimes x number of covariates) columns
##' @param params matrix of the same size as "values" consisting of the names of the parameters
##' @param covariates a vector of the names of the covariates to be used in the regime-switching functions
##' 
##' @details
##' Note that each row of the transition probability matrix must sum to one. To accomplish this
##' fix at least one transition log odds parameter in each row of "values" (including its intercept 
##' and the regression slopes of all covariates) to 0.
##' 
##' @examples
##' #Regime-switching with no covariates (self-transition ID)
##' b <- prep.regimes(values=matrix(0, 3, 3), 
##' params=matrix(c(0, 'p1', 'p2', 'p3', 0, 'p4', 'p5', 'p6', 0), 3, 3))
##' 
##' #Regime switching with no covariates (second regime ID)
##' b <- prep.regimes(values=matrix(0, 3, 3), 
##' params=matrix(c('p1', 'p2', 'p3', 0, 0, 0, 'p4', 'p5', 'p6'), 3, 3))
##' 
##' #2 regimes with three covariates
##' b <- prep.regimes(values=matrix(c(0), 2, 8), 
##' params=matrix(c(paste0('p', 8:15), rep(0, 8)), 2, 8), 
##' covariates=c('x1', 'x2', 'x3'))
prep.regimes <- function(values, params, covariates){
	if(!missing(values)){
		values <- preProcessValues(values)
	}
	if(!missing(params)){
		params <- preProcessParams(params)
	}
	if(missing(values)){
		values <- preProcessValues(matrix(0, 0, 0))
		params <- preProcessParams(matrix(0, 0, 0))
	}
	if(missing(covariates)){
		covariates <- character(0)
	}
	
	if(!all(dim(values))==all(dim(params))) {
	  stop('The dimensions of values and params matrices must match.')
	}else if(!ncol(values )== nrow(values)*(length(covariates)+1)){
	  stop(paste0("The matrix values should have ",nrow(values)*(length(covariates)+1), " columns."))
	}
	
	sv <- extractValues(values, params)
	pn <- extractParams(params)
	
	if(length(sv) > nrow(values)*(nrow(values)-1)*(length(covariates)+1)){
		stop("Regime transition probability parameters not uniquely identified.\nFix all parameters in at least one cell of each row of the \ntransition probability matrix to be zero.")
	}
	x <- list(startval=sv, paramnames=pn, values=values, params=params, covariates=covariates)
	return(new("dynrRegimes", x))
}



autojacob<-function(formula,n){
  #formula=list(x1~a*x1,x2~b*x2)
  tuple=lapply(formula,as.list)
  lhs=sapply(tuple,function(x){deparse(x[[2]])})
  rhs=sapply(tuple,function(x){x[[3]]})
  rhsj=vector("list", n*n)
  jacob=vector("list", n*n)
  for (i in 1:n){
    for (j in 1:n){
      rhsj[[(i-1)*n+j]]=paste0(deparse(D(rhs[[i]],lhs[[j]]),width.cutoff = 500L),collapse="")
      jacob[[(i-1)*n+j]]=as.formula(paste0(lhs[[i]],"~",lhs[[j]],"~",rhsj[[(i-1)*n+j]],collapse=""))
    }
  }
  return(list(row=as.list(rep(lhs,each=2)),col=as.list(rep(lhs,2)),rhsj=rhsj,jacob=jacob))
}

#------------------------------------------------------------------------------
# "Dynamics" functions

##' Recipe function for specifying dynamic functions using formulas
##' 
##' @param formula a list of formulas specifying the drift or state-transition 
##' equations for the latent variables in continuous or discrete time, respectively.
##' @param startval a named vector of starting values of the parameters in the 
##' formulas for estimation with parameter names as its name.
##' @param isContinuousTime if True, the left hand side of the formulas represent 
##' the first-order derivatives of the specified variables; if False, the left hand 
##' side of the formulas represent the current state of the specified variable while 
##' the same variable on the righ hand side is its previous state.  
##' @param jacobian (optional) a list of formulas specifying the analytic jacobian matrices 
##' containing the analytic differentiation function of the dynamic functions with respect to
##' the latent variables. If this is not provided, dynr will invoke an automatic differentiation
##' procedure to compute the jacobian functions.
##' 
##' @examples
##' #Not run: 
##' #For a full demo example that uses user-supplied analytic jacobian functions see:
##' #demo(RSNonlinearDiscrete, package="dynr")
##' formula=list(list(x1~a1*x1,x2~a2*x2),list(x1~a1*x1+c12*(exp(abs(x2)))/(1+exp(abs(x2)))*x2,
##'     x2~a2*x2+c21*(exp(abs(x1)))/(1+exp(abs(x1)))*x1))
##' jacob=list(
##'   list(x1~x1~a1,
##'       x2~x2~a2),
##'   list(x1~x1~a1,
##'       x1~x2~c12*(exp(abs(x2))/(exp(abs(x2))+1)+x2*sign(x2)*exp(abs(x2))/(1+exp(abs(x2))^2)),
##'       x2~x2~a2,
##'       x2~x1~c21*(exp(abs(x1))/(exp(abs(x1))+1)+x1*sign(x1)*exp(abs(x1))/(1+exp(abs(x1))^2))))
##' dynm<-prep.formulaDynamics(formula=formula,startval=c(a1=.3,a2=.4,c12=-.5,c21=-.5),
##' isContinuousTime=FALSE,jacobian=jacob)
##' 
##' #For a full demo example that uses automatic jacobian functions see:
##' #demo(RSNonlinearODE , package="dynr")
##' formula=list(prey~ a*prey - b*prey*predator, predator~ -c*predator + d*prey*predator)
##' dynm<-prep.formulaDynamics(formula=formula,
##'                           startval=c(a = 2.1, c = 0.8, b = 1.9, d = 1.1),
##'                           isContinuousTime=TRUE)
prep.formulaDynamics <- function(formula, startval, isContinuousTime=FALSE, jacobian){
  if(is.null(names(startval))){
    stop('startval must be a named vector')
  }
  # e.g. for the one-regime case, if we get a list of formula, make a list of lists of formula
  if(is.list(formula) && plyr::is.formula(formula[[1]])){
    formula <- list(formula)
  }
  x <- list(formula=formula, startval=startval, paramnames=c(preProcessParams(names(startval))), isContinuousTime=isContinuousTime)
  if (missing(jacobian)){
    autojcb=try(lapply(formula,autojacob,length(formula[[1]])))
    if (class(autojcb) == "try-error") {
      stop("Automatic differentiation is not supported by part of the dynamic functions.\n 
           Please provide the analytic jacobian functions.")
    }else{
      jacobian=lapply(autojcb,"[[","jacob")
    }
  }
  x$jacobian <- jacobian
  x$paramnames<-names(x$startval)
  return(new("dynrDynamicsFormula", x))
}

##' Recipe function for creating Linear Dynamcis using matrices
##' 
##' @param params.dyn the matrix of parameter names for the transition matrix in the 
##' specified linear dynamic model
##' @param values.dyn the matrix of starting/fixed values for the transition matrix in the 
##' specified linear dynamic model
##' @param params.exo the matrix of parameter names for the regression slopes of covariates on the latent variables (see details)
##' @param values.exo matrix of starting/fixed values for the regression slopes of covariates on the latent variables (see details)
##' @param params.int vector of names for intercept parameters in the dynamic model specified as a matrix or list of matrices. 
##' @param values.int vector of intercept values in the dynamic model specified as matrix or list of matrices. Contains starting/fixed values of the intercepts.
##' @param covariates the names or the index numbers of the covariates used in the dynamic model
##' @param isContinuousTime logical. When TRUE, use a continuous time model.  When FALSE use a discrete time model.
##' 
##' @details
##' The right-hand-side of the dynamic model consists of a vector of latent variables for the next time point in the discrete time case,
##' and the vector of derivatives for the latent variables at the current time point in the continuous time case.
prep.matrixDynamics <- function(params.dyn=NULL, values.dyn, params.exo=NULL, values.exo=NULL, params.int=NULL, values.int=NULL, 
                                covariates, isContinuousTime){
	# Handle numerous cases of missing or non-list arguments
	# General idea
	# If they give us a non-list argument, make it a one-element list
	# If they don't give us params, assume it's all fixed params
	# If they don't give us values, assume they don't want that part of the model
	
	# Handle dyn
	r <- coProcessValuesParams(values.dyn, params.dyn)
	values.dyn <- r$values
	params.dyn <- r$params
	
	# Handle exo
	r <- coProcessValuesParams(values.exo, params.exo, missingOK=TRUE)
	values.exo <- r$values
	params.exo <- r$params
	
	# Handle int
	r <- coProcessValuesParams(values.int, params.int, missingOK=TRUE)
	values.int <- r$values
	params.int <- r$params
	
	
	if(missing(covariates)){
		covariates <- character(0)
	}
	values.dyn <- lapply(values.dyn, preProcessValues)
	params.dyn <- lapply(params.dyn, preProcessParams)
	values.exo <- lapply(values.exo, preProcessValues)
	params.exo <- lapply(params.exo, preProcessParams)
	values.int <- lapply(values.int, preProcessValues)
	params.int <- lapply(params.int, preProcessParams)

	
	sv <- c(extractValues(values.dyn, params.dyn), extractValues(values.exo, params.exo), extractValues(values.int, params.int))
	pn <- c(extractParams(params.dyn), extractParams(params.exo), extractParams(params.int))
	sv <- extractValues(sv, pn)
	pn <- extractParams(pn)
	
	x <- list(startval=sv, paramnames=pn, params.dyn=params.dyn, values.dyn=values.dyn, params.exo=params.exo, values.exo=values.exo, params.int=params.int, values.int=values.int, isContinuousTime=isContinuousTime,covariates=covariates)
	return(new("dynrDynamicsMatrix", x))
}


# nonlinear functions for lookup
# logit, logistic, softmax

#https://en.wikipedia.org/wiki/C_mathematical_functions

#carl ~ param[5]*carl + param[7]*logistic(abs(bob))
#bob ~ param[4]*bob + param[6]*logistic(abs(carl))

#cform <- carl ~ param[5] * carl + param[7] * logistic(abs(bob)) + dan**2
#rhs.vars(cform)

isSymbolNumberFunction <- function(x){is.symbol(x) || is.numeric(x) || is.function(x)}

# Recursive function for formula parsing
parseFormula <- function(formula, debug=FALSE){
  tuple <- as.list(formula)
  op<-tuple[[1]]
  left <- tuple[[2]]
  right <- tuple[[3]]
	if(debug){
		print("LEFT")
		print(left)
		print("RIGHT")
		print(right)
	}
	if(!isSymbolNumberFunction(left)){
	  leftTuple <- as.list(left)
		if(length(leftTuple)==3){
			left <- parseFormula(left,debug)
		} else {
			left <- parseNested(left,debug)
		}
	}
	if(!isSymbolNumberFunction(right)){
	  rightTuple <- as.list(right)
		if(length(rightTuple)==3){
			right <- parseFormula(right,debug)
		} else{
			right <- parseNested(right,debug)
		}
	}

  op <- trans2CFunction(op)
  
  tuple.new <- list(op,left,right)
  outFormula <- as.call(tuple.new)
	return(outFormula)
}

parseNested <- function(formula,debug=FALSE){
	tuple <- as.list(formula)
	outer <- tuple[[1]]
	inner <- tuple[[2]]
	if(debug){
	  print(tuple)
	  print("OUTER")
	  print(outer)
	  print("INNER")
	  print(inner)
	}
	if(!isSymbolNumberFunction(inner)){
	  innerTuple=as.list(inner)
	  if(length(innerTuple)==3){
	    inner <- parseFormula(inner,debug)
	  } else{
	    inner <- parseNested(inner,debug)
	  }
	}
	
	outer <- trans2CFunction(outer)
	
	tuple.new <-as.list(c(outer,inner))
	outFormula <- as.call(tuple.new)
	return(outFormula)
}

#TODO check type casting int->double
#TODO check the parameter indices
trans2CFunction<-function(op.symbol){
  op.char=deparse(op.symbol)
  if (op.char %in% c("d","abs","^","**","mod","max","min","sign")) {
    op.symbol<-switch(op.char,
              d = NULL,       
            abs = as.name("fabs"),
            "^" = as.name("pow"),
            "**"= as.name("pow"),
            mod = as.name("fmod"),
            max = as.name("fmax"),
            min = as.name("fmin"),
            sign = c(as.name("copysign"),1))
    }
  return(op.symbol)
}

processFormula<-function(formula.list){
  formula.char=sapply(formula.list,FUN=function(f){paste0(deparse(parseFormula(f),width.cutoff = 500L),collapse="")})
  out=strsplit(formula.char," ~ ")
  return(out)
}
# Example usage
#require(dynr)
#logistic <- function(x){x}
#
#cform <- carl ~ param[5] * carl + param[7] * logistic(abs(bob)) + dan**2
#cform2 <- parseFormula(cform)
#cform2
#str(cform2)
#
# Note also some helpful functions in the pryr package
#pryr::ast(carl ~ param[5] * carl + param[7] * logistic(abs(bob)) + dan**2)
#pryr::call_tree(cform)
#
# The source code for 
#pryr:::tree
# is useful in this regard.
# We are essentially making our own version of pryr:::tree
# that does some substitutions before collapsing back to a 
# string that is composed of C code.


##' Recipe function for preparing the initial conditions for the model. 
##' 
##' @param values.inistate a vector or list of vectors of the starting or fixed values of the initial state vector in one or more regimes.
##' @param params.inistate a vector or list of vectors of the parameter names that appear in the initial state vector in one or more regimes. If an element is 0 or "fixed", the corresponding element is fixed at the value specified in the values vector; Otherwise, the corresponding element is to be estimated with the starting value specified in the values vector.
##' @param values.inicov a positive definite matrix or a list of positive definite matrices of the starting or fixed values of the initial error covariance structure(s) in one or more regimes. If only one matrix is specified for a regime-switching dynamic model, the initial error covariance structure stays the same across regimes. To ensure the matrix is positive definite in estimation, we apply LDL transformation to the matrix. Values are hence automatically adjusted for this purpose.
##' @param params.inicov a matrix or list of matrices of the parameter names that appear in the initial error covariance(s) in one or more regimes. If an element is 0 or "fixed", the corresponding element is fixed at the value specified in the values matrix; Otherwise, the corresponding element is to be estimated with the starting value specified in the values matrix. If only one matrix is specified for a regime-switching dynamic model, the process noise structure stays the same across regimes. If a list is specified, any two sets of the parameter names as in two matrices should be either the same or totally different to ensure proper parameter estimation.
##' @param values.regimep a vector of the starting or fixed values of the initial probabilities of being in each regime. By default, the initial probability of being in the first regime is fixed at 1.
##' @param params.regimep a vector of the parameter indices of the initial probabilities of 
##' being in each regime. If an element is 0 or "fixed", the corresponding element is fixed at the value 
##' specified in the "values" vector; Otherwise, the corresponding element is to be estimated 
##' with the starting value specified in the values vector.
##' 
##' @details
##' The initial condition model includes specifications for the intial state vector, 
##' initial error covariance matrix, initial probabilities of 
##' being in each regime and all associated parameter specifications.
##' 
prep.initial <- function(values.inistate, params.inistate, values.inicov, params.inicov, values.regimep=1, params.regimep=0){
	# Handle initial state
	r <- coProcessValuesParams(values.inistate, params.inistate)
	values.inistate <- r$values
	params.inistate <- r$params
	
	# Handle initial covariance
	r <- coProcessValuesParams(values.inicov, params.inicov)
	values.inicov <- r$values
	params.inicov <- r$params
	
	if(length(values.inistate) != length(values.inicov)){
		stop('Initial state and covariance imply a different number of regimes.')
	}
	
	values.inistate <- lapply(values.inistate, preProcessValues)
	params.inistate <- lapply(params.inistate, preProcessParams)
	
	values.inicov <- lapply(values.inicov, preProcessValues)
	params.inicov <- lapply(params.inicov, preProcessParams)
	values.inicov.inv.ldl <- lapply(values.inicov, replaceDiagZero)
	values.inicov.inv.ldl <- lapply(values.inicov.inv.ldl, reverseldl)
	
	values.regimep <- preProcessValues(values.regimep)
	params.regimep <- preProcessParams(params.regimep)
	
	sv <- c(extractValues(values.inistate, params.inistate), extractValues(values.inicov.inv.ldl, params.inicov), extractValues(values.regimep, params.regimep))
	pn <- c(extractParams(params.inistate), extractParams(params.inicov), extractParams(params.regimep))
	sv <- extractValues(sv, pn)
	pn <- extractParams(pn)
	
	x <- list(startval=sv, paramnames=pn, values.inistate=values.inistate, params.inistate=params.inistate, values.inicov=values.inicov, values.inicov.inv.ldl=values.inicov.inv.ldl, params.inicov=params.inicov, values.regimep=values.regimep, params.regimep=params.regimep)
	return(new("dynrInitial", x))
}

##' Create a dynrTrans object to handle the transformations and inverse 
##' transformations of model paramters
##' 
##' @param formula.trans a list of formulae for transforming freed parameters 
##' other than variance-covariance parameters during the optimization process. 
##' These transformation functions may be helpful for transforming parameters 
##' that would normally appear on a constrained scale to an unconstrained 
##' scale (e.g., parameters that can only take on positive values can be 
##' subjected to exponential transformation to ensure positivity.)
##' @param formula.inv a list of formulae that inverse the transformation 
##' on the free parameters and will be used to calculate the starting values 
##' of the parameters.
##' @param transCcode a logical value indicating whether the functions in 
##' formula.trans need to be transformed to functions in C. The default 
##' for transCcode is TRUE, which means that the formulae will be translated 
##' to C functions and utilized during the optimization process. 
##' If transCcode = FALSE, the transformations are only performed at the end 
##' of the optimization process for standard error calculations but not 
##' during the optimization process. 
prep.tfun<-function(formula.trans, formula.inv, transCcode = TRUE){
  #input: formula.trans=list(a~exp(a),b~b^2)
  #input: formula.inv=list(a~log(a),b~sqrt(b))
  if (missing(formula.trans)&missing(formula.inv)){
    #TODO modify this part if needed
    x<-list(transCcode = FALSE)
  }else{
    if (missing(formula.inv)){
      x <- list(formula.trans=formula.trans, transCcode = transCcode)
    }else if (missing(formula.trans)){
      x <- list(formula.inv=formula.inv, transCcode = FALSE)
    }else{
      x <- list(formula.trans=formula.trans, formula.inv=formula.inv, transCcode = transCcode)
    }
  }
  return(new("dynrTrans", x))
}

#------------------------------------------------------------------------------

formula2string<-function(formula.list){
  tuple=lapply(formula.list,as.list)
  lhs=sapply(tuple,function(x){paste0(deparse(x[[2]],width.cutoff = 500L),collapse="")})
  rhs=sapply(tuple,function(x){paste0(deparse(x[[3]],width.cutoff = 500L),collapse="")})
  return(list(lhs=lhs,rhs=rhs))
}


matrix2formula <- function(x,multbyColnames=TRUE){
	if(!is.matrix(x)){
		stop("Dude! You have to give me a matrix. If you do, I'll give you a formula. Seriously.")
	}
	if(is.null(rownames(x))){
		rownames(x) <- paste0('y', 1:nrow(x))
	}
	if(is.null(colnames(x))){
		colnames(x) <- paste0('x', 1:ncol(x))
	}
	preds <- character(nrow(x))
	for(i in 1:nrow(x)){
	  if (multbyColnames==FALSE){
	    preds[i] <- paste(rownames(x)[i], paste(x[i,]), sep=' ~ ')
	  }else{
		  preds[i] <- paste(rownames(x)[i], paste(x[i,], 
		                       colnames(x), sep='*', collapse=' + '), sep=' ~ ')
		  } 
	}
	if(is.numeric(x)){
		preds <- gsub(' 1\\*', ' ', preds)
		preds <- gsub(paste(" 0\\*(", paste(colnames(x), collapse = "|"), ") \\+", sep=""), '', preds)
		preds <- gsub(paste(" \\+ 0\\*(", paste(colnames(x), collapse = "|"), ")", sep=""), '', preds)
	} else {
		preds <- gsub(' 1\\*', ' ', preds)
		preds <- gsub(paste(" 0\\*(", paste(colnames(x), collapse = "|"), ") \\+", sep=""), '', preds)
		preds <- gsub(paste(" \\+ 0\\*(", paste(colnames(x), collapse = "|"), ")", sep=""), '', preds)
	}
	form <- lapply(preds, formula, env=.GlobalEnv)
	#eq.char=lapply(form, as.character)
	#str.left=sapply(eq.char,"[",2)
	#str.right=sapply(eq.char,"[",3)
	#if (bothSides) {return(form)}else{return(str.right)}
	return(form)
}

# # Examples
# meas <- prep.loadings(
#  map=list(
#    eta1=paste0('y', 1:3),
#    eta2=paste0('y', 4:6)),
#    params = paste0("lambda",c(1:4)))
# matrix2formula(meas$values.load[[1]])
# matrix2formula(meas$params.load[[1]])
# lapply(meas$values.load,matrix2formula)
# #TODO check if fixed in formula/tex is substituted by fixed values.

addFormulas <- function(f1, f2){
	charf1 <- as.character(f1)
	charf2 <- as.character(f2)
	if((length(charf1) != 3) || (length(charf2) != 3)){
		msg <- paste("Formula 1 and/or 2 look(s) strange, and not in a good way.", "A formula used here should only have one tilde (~) in it", "I don't know what to do with them")
		stop(msg)
	}
	if(charf1[2] != charf2[2]){
		stop('The two formulas do not have the same outcome variable. I simply refuse to add them.')
	}
	formula(paste(charf1[2], "~", charf1[3], "+", charf2[3]), env=.GlobalEnv)
}

## Examples
#f1 <- y~9*x1+x2
#f2 <- y~1.2
#f3 <- y~5*u1
#addFormulas(f1, f2)
#addFormulas(f2, f1)

addLLFormulas <- function(list_list_formulae, VecNamesToAdd){
  nRegime <- length(list_list_formulae)
  for (j in 1:nRegime){
    neq <- length(list_list_formulae[[j]])
    for (k in 1:neq){
      AddedFml <- list_list_formulae[[j]][[k]]
      for (l in 1:length(VecNamesToAdd)){
        list_list_formulae_add=get(VecNamesToAdd[l], parent.frame())#eval(parse(text=VecNamesToAdd[l]),environment())
        if (length(list_list_formulae_add)>0){
          if (length(list_list_formulae_add[[j]])>0){
            AddedFml <- addFormulas(list_list_formulae_add[[j]][[k]],AddedFml)
          }
        }
      }
      list_list_formulae[[j]][[k]] <- AddedFml
    }
  }
  return(list_list_formulae)
}

#------------------------------------------------------------------------------
dP_dt <- "/**\n * The dP/dt function: depend on function_dF_dx, needs to be compiled on the user end\n * but user does not need to modify it or care about it.\n */\nvoid mathfunction_mat_to_vec(const gsl_matrix *mat, gsl_vector *vec){\n\tsize_t i,j;\n\tsize_t nx=mat->size1;\n\t/*convert matrix to vector*/\n\tfor(i=0; i<nx; i++){\n\t\tgsl_vector_set(vec,i,gsl_matrix_get(mat,i,i));\n\t\tfor (j=i+1;j<nx;j++){\n\t\t\tgsl_vector_set(vec,i+j+nx-1,gsl_matrix_get(mat,i,j));\n\t\t\t/*printf(\"%lu\",i+j+nx-1);}*/\n\t\t}\n\t}\n}\nvoid mathfunction_vec_to_mat(const gsl_vector *vec, gsl_matrix *mat){\n\tsize_t i,j;\n\tsize_t nx=mat->size1;\n\t/*convert vector to matrix*/\n\tfor(i=0; i<nx; i++){\n\t\tgsl_matrix_set(mat,i,i,gsl_vector_get(vec,i));\n\t\tfor (j=i+1;j<nx;j++){\n\t\t\tgsl_matrix_set(mat,i,j,gsl_vector_get(vec,i+j+nx-1));\n\t\t\tgsl_matrix_set(mat,j,i,gsl_vector_get(vec,i+j+nx-1));\n\t\t}\n\t}\n}\nvoid function_dP_dt(double t, size_t regime, const gsl_vector *p, double *param, size_t n_param, const gsl_vector *co_variate, gsl_vector *F_dP_dt){\n\t\n\tsize_t nx;\n\tnx = (size_t) floor(sqrt(2*(double) p->size));\n\tgsl_matrix *P_mat=gsl_matrix_calloc(nx,nx);\n\tmathfunction_vec_to_mat(p,P_mat);\n\tgsl_matrix *F_dx_dt_dx=gsl_matrix_calloc(nx,nx);\n\tfunction_dF_dx(t, regime, param, co_variate, F_dx_dt_dx);\n\tgsl_matrix *dFP=gsl_matrix_calloc(nx,nx);\n\tgsl_matrix *dP_dt=gsl_matrix_calloc(nx,nx);\n\tgsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, F_dx_dt_dx, P_mat, 0.0, dFP);\n\tgsl_matrix_transpose_memcpy(dP_dt, dFP);\n\tgsl_matrix_add(dP_dt, dFP);\n\tsize_t n_Q_vec=(1+nx)*nx/2;\n\tgsl_vector *Q_vec=gsl_vector_calloc(n_Q_vec);\n\tsize_t i;\n\tfor(i=1;i<=n_Q_vec;i++){\n\t\t\tgsl_vector_set(Q_vec,n_Q_vec-i,param[n_param-i]);\n\t}\n\tgsl_matrix *Q_mat=gsl_matrix_calloc(nx,nx);\n\tmathfunction_vec_to_mat(Q_vec,Q_mat);\n\tgsl_matrix_add(dP_dt, Q_mat);\n\tmathfunction_mat_to_vec(dP_dt, F_dP_dt);\n\tgsl_matrix_free(P_mat);\n\tgsl_matrix_free(F_dx_dt_dx);\n\tgsl_matrix_free(dFP);\n\tgsl_matrix_free(dP_dt);\n\tgsl_vector_free(Q_vec);\n\tgsl_matrix_free(Q_mat);\n}\n"

#N.B. The formula 
#	nx = (size_t) floor(sqrt(2*(double) p->size));
# is used instead of the analytic formula
#	nx = (size_t) sqrt(2*double p->size + .25) - .5
# due to issues with rounding.



#------------------------------------------------------------------------------
# Utility functions for writing GSL code

setGslMatrixElements <- function(values, params, name, depth=1){
	ret <- ""
	numRow <- nrow(values)
	numCol <- ncol(values)
	for(j in 1:numCol){
		for(i in 1:numRow){
			if(params[i, j] > 0){
				ret <- paste(ret,
					paste(rep('\t', depth), collapse=""), 'gsl_matrix_set(', name, ', ', i-1, ', ', j-1,
					', param[', params[i, j] - 1, ']);\n', sep='')
			} else if(values[i, j] != 0){
				ret <- paste(ret,
					paste(rep('\t', depth), collapse=""), 'gsl_matrix_set(', name, ', ', i-1, ', ', j-1,
					', ', values[i, j], ');\n', sep='')
			}
		}
	}
	return(ret)
}

setGslVectorElements <- function(values, params, name, fill="", depth=1){
  ret <- ""
  numLength <- length(values)
  for(i in 1:numLength){
      if(params[i] > 0){
        ret <- paste(ret,
                     paste(rep('\t', depth), collapse=""), 'gsl_vector_set(', name, ', ', fill, i-1,
                     ', param[', params[i] - 1, ']);\n', sep='')
      } else if(values[i] != 0){
        ret <- paste(ret,
                     paste(rep('\t', depth), collapse=""), 'gsl_vector_set(', name, ', ', fill, i-1,
                     ', ', values[i], ');\n', sep='')
      }
    }
  return(ret)
}


createGslMatrix <- function(nrow, ncol, name){
	paste0("\tgsl_matrix *", name, " = gsl_matrix_calloc(", nrow, ", ", ncol, ");\n")
}

destroyGslMatrix <- function(name){
	paste0("\tgsl_matrix_free(", name, ");\n")
}

createGslVector <- function(size, name){
	paste0("\tgsl_vector *", name, " = gsl_vector_calloc(", size, ");\n")
}

destroyGslVector <- function(name){
	paste0("\tgsl_vector_free(", name, ");\n")
}

#y <- alpha * transA(A) %*% x + beta * y
blasMV <- function(transA, alpha, A, x, beta, y){
	transA <- ifelse(transA, "CblasTrans", "CblasNoTrans")
	paste0("\tgsl_blas_dgemv(", paste(transA, alpha, A, x, beta, y, sep=", "), ");\n")
}

# C <- alpha * transA(A) %*% transB(B) + beta * C
blasMM <- function(transA, transB, alpha, A, B, beta, C){
	transA <- ifelse(transA, "CblasTrans", "CblasNoTrans")
	transB <- ifelse(transB, "CblasTrans", "CblasNoTrans")
	paste0("\tgsl_blas_dgemm(", paste(transA, transB, alpha, A, B, beta, C, sep=", "), ");\n")
}

gslVector2Column <- function(matrix, index, vector, which){
	if(which=='column'){
		paste0("\tgsl_matrix_set_col(", paste(matrix, index, vector, sep=", "), ");\n")
	} else if(which=='row'){
		paste0("\tgsl_matrix_set_row(", paste(matrix, index, vector, sep=", "), ");\n")
	}
}

gslcovariate.front<-function(selected, covariates){
  ret <- createGslVector(length(selected), "covariate_local")
  for (i in 1:length(selected)){
    get <- paste0("gsl_vector_get(co_variate, ", which(covariates == selected[i])-1,")")
    set <- paste0("\tgsl_vector_set(covariate_local, ", i-1,", ", get,");\n")
    ret <- paste0(ret, set)
  }
  return(ret)
}
