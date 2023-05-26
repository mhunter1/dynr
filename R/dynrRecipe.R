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
##' print,dynrRecipe-method
##' show,dynrRecipe-method
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
           state.names = "character",
           obs.names = "character",
           exo.names="character",
           dmudparMu="matrix",
           dmudparMu2="matrix"
           ),
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
           formulaOriginal = "list",
           jacobianOriginal = "list",
           state.names = "character",
           theta.names = "character",
           beta.names = "character",
           random.names = "character",
           random.lb = "numeric",
           random.ub = "numeric",
           intercept.names = "character",
           theta.formula = "list",
           isContinuousTime = "logical",
           dfdtheta= "list",
           dfdx2= "list",
           dfdxdtheta= "list", 
           dfdthetadx= "list", 
           dfdtheta2= "list",
           formula2 = "list",
           saem = "logical",
           random.params.inicov="matrix",
           random.values.inicov="matrix",
		   covariate.formula = "list",
		   covariate.names = "character"
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
           isContinuousTime = "logical",
           random.names = "character",
           random.lb = "numeric",
           random.ub = "numeric",
           random.params.inicov="matrix",
           random.values.inicov="matrix"
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
           covariates = "character",
           deviation = "logical",
           refRow = "numeric"),
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
           params.regimep = "matrix",
           covariates = "character",
           deviation = "logical",
           refRow = "numeric",
           y0 = "list"),
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
           values.observed.inv.ldl = "list",
		   covariates = "character",
		   latent.formula = "list",
		   state.names = "character",
		   latent.startval = "numeric",
		   ldl.transformed = "matrix",
		   is_eta_cov_formula = "logical",
		   known.vars = "list"),
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

# initial class for variables of random effects (currently not used)
setClass(Class =  "dynrRandom",
         representation = representation(
           random.names = "character",
           random.lb = "numeric",
           random.ub = "numeric",
           num.subj = "numeric",
           values.inicov = "matrix",
           params.inicov = "matrix",
           b="matrix"
         ),
         contains = "dynrRecipe"
)

setClass(Class =  "dynrSaemParameter",
         representation = representation(
           MAXGIB = "numeric",
           MAXITER = "numeric",
           maxIterStage1 = "numeric",
           gainpara = "numeric",
           gainparb = "numeric",
           gainpara1 = "numeric",
           gainparb1 = "numeric",
           bAdaptParams = "vector", 
           KKO = "numeric",
		   scaleb = "numeric",
		   setScaleb = "numeric",
		   setAccept = "numeric",
		   seed = "numeric",
		   trueb = "matrix",
		   errtrol1 = "numeric",
		   errtrol = "numeric"
         ),
         contains = "dynrRecipe"
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

printRecipeOrModel <- function(x, ...){
	for(aslot in slotNames(x)){
		if( !(aslot %in% c("c.string", "startval", "paramnames", "formulaOriginal", "jacobianOriginal")) ){
			cat(" $", aslot, "\n", sep="")
			print(slot(x, aslot))
			cat("\n")
		}
	}
	return(invisible(x))
}

setMethod("print", "dynrRecipe", printRecipeOrModel)
setMethod("show", "dynrRecipe", function(object){printRecipeOrModel(object)})


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
##' @param object The dynr object (recipe, model, or cooked model).
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
##' 
##' @return character text suitable for use fiel LaTeX
##' 
##' @seealso A way to put this in a plot with \code{\link{plotFormula}}
setGeneric("printex", function(object, ParameterAs, 
	printDyn=TRUE, printMeas=TRUE, printInit=FALSE, printRS=FALSE, outFile, show, ...) { 
	return(standardGeneric("printex")) 
})

setMethod("printex", "dynrMeasurement",
	function(object, ParameterAs, printDyn=TRUE, printMeas=TRUE, printInit=FALSE, printRS=FALSE, outFile, show=TRUE, AsMatrix=TRUE){
		if (AsMatrix){
			meas_loadings <- lapply(object$values.load, .xtableMatrix, show)
			meas_int <- lapply(object$values.int, .xtableMatrix, show)
			meas_exo <- lapply(object$values.exo, .xtableMatrix, show)
			meas_exo.names <- .xtableMatrix(matrix(object$exo.names, ncol=1), show)
			meas_list <- list(meas_loadings=meas_loadings, meas_int=meas_int,
				meas_exo=meas_exo, meas_exo.names=meas_exo.names)
			return(invisible(meas_list))
		} else { #not to use matrix
			#Matrix to Formula
			meas_loadings <- lapply(object$values.load, matrix2formula, multbyColnames=TRUE)
			meas_int <- lapply(object$values.int, matrix2formula, multbyColnames=FALSE)
			meas_exo <- lapply(object$values.exo, matrix2formula, multbyColnames=TRUE)
			meas_fml <- addLLFormulas(list_list_formulae = meas_loadings,
				VecNamesToAdd = c("meas_int", "meas_exo"))
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


mvpaste <- function(m, v, a){
	nr <- nrow(m)
	res <- matrix("", nrow=nr, ncol=1)
	for(r in 1:nr){
		mat1 <- m[r,]
		mat2 <- v
		mat2 <- mat2[mat1 !=0 ]
		mat1 <- mat1[mat1 !=0 ]
		mat3 <- paste(mat1, mat2, sep="*")
		mat3 <- gsub("*1", "", mat3, fixed=TRUE)
		a <- a[ a != 0]
		b <- implode(c(mat3, a), sep=" + ")
		b[b %in% ""] <- "0"
		res[r,] <- b
	}
	return(res)
}

setMethod("printex", "dynrInitial",
	function(object, ParameterAs, printDyn=TRUE, printMeas=TRUE,
	printInit=FALSE, printRS=FALSE, outFile, show=TRUE, AsMatrix=TRUE){
			values.regimep <- object$values.regimep
			params.regimep <- object$params.regimep
			numRegimes <- nrow(values.regimep)
			covariates <- object$covariates
			numCovariates <- length(covariates)
			deviation <- object$deviation
			refRow <- object$refRow
			
			nx <- nrow(object$values.inistate[[1]])
			covar <- c("1", object$covariates)
			x0val <- lapply(object$values.inistate, mvpaste, v=covar, a=rep("0", nx))
			lx0 <- lapply(x0val, .xtableMatrix, show)
			lP0 <- lapply(object$values.inicov, .xtableMatrix, show)
			
			if(deviation){
				if(nrow(values.regimep)!=0 && nrow(params.regimep)!=0){
					values.regIntercept <- matrix(values.regimep[refRow, 1], nrow=numRegimes, 1)
					params.regIntercept <- matrix(params.regimep[refRow, 1], nrow=numRegimes, 1)
					values.regimep[refRow, 1] <- 0
					params.regimep[refRow, 1] <- 0
				}
			} else {
				values.regIntercept <- matrix(0, numRegimes, 1)
				params.regIntercept <- matrix(0, numRegimes, 1)
			}
			
			p0val <- mvpaste(values.regimep, covar, values.regIntercept)
			lr0 <- .xtableMatrix(p0val, show)
			
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

setGeneric("paramName2Number", function(object, names, invert=FALSE) { 
    return(standardGeneric("paramName2Number")) 
	return(standardGeneric("paramName2Number")) 
})

setGeneric("paramName2NumericNumber", function(object, paramList) { 
  return(standardGeneric("paramName2NumericNumber")) 
})

# params is a matrix with parameter names in it
# names is a vector of the parameter names
# returns matrix with the parameters numbers instead of the names
.exchangeNamesAndNumbers <- function(params, names){
	matrix(match(params, names, nomatch=0), nrow(params), ncol(params), dimnames=dimnames(params))
}

# params is a matrix with parameter numbers in it
# names is a vector of the parameter names
# returns matrix with the parameters names instead of the numbers
.exchangeNumbersAndNames <- function(params, names){
	matrix(c(0, names)[params + 1], nrow(params), ncol(params), dimnames=dimnames(params))
}

# formula is the list of dynamics formula
# paramnames is the character vector of names of parameters used in the formula
# names is the character vector of the names of all the free parameters
# returns list of formula with e.g. aParamName swapped to param[2]
.exchangeformulaNamesAndNumbers <- function(formula, paramnames, names){
	string <- paste0(deparse(formula, width.cutoff = 500L), collapse="")
	pattern <- paramnames
	pattern <- gsub("\\\\", "\\\\\\\\", pattern)
	
	for(i in 1:length(paramnames)){
		string <- gsub(paste0("\\<", pattern[i], "\\>"), paste0("param[", match(paramnames[i], names, nomatch=0)-1, "]"), string)
	}
	eval(parse(text=string))
}

# formula is the list of dynamics formula using C-style param[2] in place of character names
# paramnames is the character vector of names of parameters used in the formula
# names is the character vector of the names of all the free parameters
# returns list of formula with e.g. param[2] swapped to aParamName
.exchangeformulaNumbersAndNames <- function(formula, paramnames, names){
	string <- paste0(deparse(formula, width.cutoff = 500L), collapse="")
	pattern <- paramnames
	pattern <- gsub("\\\\", "\\\\\\\\", pattern)
	
	for(i in 1:length(paramnames)){
		string <- gsub(paste0("param\\[", match(paramnames[i], names, nomatch=0)-1, "\\]"), paste0("", pattern[i], ""), string)
	}
	eval(parse(text=string))
}

# Don't know why this exists
.replaceNamesbyNumbers <- function(formula, paramtoPlot){
  string <- paste0(deparse(formula,width.cutoff = 500L),collapse="")
  names <- paste0("param\\[",(1:length(paramtoPlot))-1,"\\]")
  #names(paramtoPlot)
  #pattern=gsub("\\{","\\\\\\{",names)
  #pattern=gsub("\\}","\\\\\\}",pattern)
  #pattern=gsub("\\\\","\\\\\\\\",pattern)
  for (i in 1:length(names)){
    string <- gsub(paste0("\\<", names[i], "\\>"), paramtoPlot[i], string)
  }
  eval(parse(text=string))
}


setMethod("paramName2Number", "dynrMeasurement",
	function(object, names, invert=FALSE){
		if(!invert){
			exchangeFUN <- .exchangeNamesAndNumbers
		} else {
			exchangeFUN <- .exchangeNumbersAndNames
		}
		object@params.load <- lapply(object$params.load, exchangeFUN, names=names)
		object@params.exo <- lapply(object$params.exo, exchangeFUN, names=names)
		object@params.int <- lapply(object$params.int, exchangeFUN, names=names)
		return(object)
	}
)

# Don't know why this exists
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
	function(object, names, invert=FALSE){
		if(!invert){
			exchangeFUN <- .exchangeformulaNamesAndNumbers
		} else {
			exchangeFUN <- .exchangeformulaNumbersAndNames
		}
		object@formula <- exchangeFUN(object@formula, object@paramnames, names)
		object@jacobian <- exchangeFUN(object@jacobian, object@paramnames, names)
		return(object)
	}
)

# Don't know why this exists
setMethod("paramName2NumericNumber", "dynrDynamicsFormula",
          function(object, paramList){
            object@formulaOriginal <- object@formula
            object@jacobianOriginal <- object@jacobian
            object@formula = .replaceNamesbyNumbers(object@formula, paramList)
            object@jacobian =.replaceNamesbyNumbers(object@jacobian, paramList)
            return(object)
          }
)


setMethod("paramName2Number", "dynrDynamicsMatrix",
	function(object, names, invert=FALSE){
		if(!invert){
			exchangeFUN <- .exchangeNamesAndNumbers
		} else {
			exchangeFUN <- .exchangeNumbersAndNames
		}
		object@params.dyn <- lapply(object$params.dyn, exchangeFUN, names=names)
		object@params.exo <- lapply(object$params.exo, exchangeFUN, names=names)
		object@params.int <- lapply(object$params.int, exchangeFUN, names=names)
		return(object)
	}
)


setMethod("paramName2Number", "dynrRegimes",
	function(object, names, invert=FALSE){
		if(!invert){
			exchangeFUN <- .exchangeNamesAndNumbers
		} else {
			exchangeFUN <- .exchangeNumbersAndNames
		}
		object@params <- exchangeFUN(object$params, names)
		return(object)
	}
)


setMethod("paramName2Number", "dynrInitial",
	function(object, names, invert=FALSE){
		if(!invert){
			exchangeFUN <- .exchangeNamesAndNumbers
		} else {
			exchangeFUN <- .exchangeNumbersAndNames
		}
		object@params.inistate <- lapply(object$params.inistate, exchangeFUN, names=names)
		object@params.inicov <- lapply(object$params.inicov, exchangeFUN, names=names)
		object@params.regimep <- exchangeFUN(object$params.regimep, names=names)
		return(object)
	}
)


setMethod("paramName2Number", "dynrNoise",
	function(object, names, invert=FALSE){
		if(!invert){
			exchangeFUN <- .exchangeNamesAndNumbers
		} else {
			exchangeFUN <- .exchangeNumbersAndNames
		}
		object@params.latent <- lapply(object$params.latent, exchangeFUN, names=names)
		object@params.observed <- lapply(object$params.observed, exchangeFUN, names=names)
		return(object)
	}
)

setMethod("paramName2Number", "dynrTrans",
	function(object, names, invert=FALSE){
		if(!invert){
			exchangeFUN <- .exchangeformulaNamesAndNumbers
		} else {
			exchangeFUN <- .exchangeformulaNumbersAndNames
		}
		if (length(object@formula.trans) > 0){
			# TODO is this invert behavior correct?
			object@formula.trans <- exchangeFUN(object@formula.trans, object@paramnames, names)
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
		nregime <- length(formula)
		n <- sapply(formula, length)
		nj <- sapply(jacob, length)
		
		fml <- lapply(formula, processFormula)
		lhs <- lapply(fml, function(x){lapply(x, "[[", 1)})
		rhs <- lapply(fml, function(x){lapply(x, "[[", 2)})
		
		fmlj <- lapply(jacob, processFormula)
		row <- lapply(fmlj, function(x){lapply(x, "[[", 1)})
		col <- lapply(fmlj, function(x){lapply(x, "[[", 2)})
		rhsj <- lapply(fmlj, function(x){lapply(x, "[[", 3)})
		
		for (i in 1:length(covariates)){
			selected <- covariates[i]
			get <- paste0("gsl_vector_get(co_variate, ", which(covariates == selected)-1, ")")
			rhs <- lapply(gsub(paste0("\\<", selected, "\\>"), get, rhs), function(x){eval(parse(text=x))})
			rhsj <- lapply(gsub(paste0("\\<", selected, "\\>"), get, rhsj), function(x){eval(parse(text=x))})
		}
		
		if (object@isContinuousTime){
			#function_dx_dt
			ret <- "void function_dx_dt(double t, size_t regime, const gsl_vector *x, double *param, size_t n_param, const gsl_vector *co_variate, gsl_vector *F_dx_dt){"
			
			ret <- paste0(ret, cswapDynamicsFormulaLoop(nregime=nregime, n=n, lhs=lhs, rhs=rhs, target1='gsl_vector_get(x, ', close=')', target2='gsl_vector_set(F_dx_dt, ', vector=TRUE))
			
			#function_dF_dx
			ret <- paste0(ret, "\n\n/**\n* The dF/dx function\n* The partial derivative of the jacobian of the DE function with respect to the variable x\n* @param param includes at the end the current state estimates in the same order as the states following the model parameters\n*/void function_dF_dx(double t, size_t regime, double *param, const gsl_vector *co_variate, gsl_matrix *F_dx_dt_dx){")
			
			ret <- paste0(ret, cswapDynamicsFormulaLoop(nregime=nregime, n=nj, lhs=lhs, rhs=rhsj, target1='param[NUM_PARAM+', close=']', target2='gsl_matrix_set(F_dx_dt_dx, ', row=row, col=col))
			
		} else{ # is Discrete Time
			#function_dynam
			ret <- "void function_dynam(const double tstart, const double tend, size_t regime, const gsl_vector *xstart,\n\tdouble *param, size_t n_gparam, const gsl_vector *co_variate,\n\tvoid (*g)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),\n\tgsl_vector *x_tend){"
			
			ret <- paste0(ret, cswapDynamicsFormulaLoop(nregime=nregime, n=n, lhs=lhs, rhs=rhs, target1='gsl_vector_get(xstart, ', close=')', target2='gsl_vector_set(x_tend, ', vector=TRUE))
			
			#function_jacob_dynam
			ret <- paste0(ret, "\n\nvoid function_jacob_dynam(const double tstart, const double tend, size_t regime, const gsl_vector *xstart,\n\tdouble *param, size_t num_func_param, const gsl_vector *co_variate,\n\tvoid (*g)(double, size_t, double *, const gsl_vector *, gsl_matrix *),\n\tgsl_matrix *Jx){")
			
			ret <- paste0(ret, cswapDynamicsFormulaLoop(nregime=nregime, n=nj, lhs=lhs, rhs=rhsj, target1='gsl_vector_get(xstart, ', close=')', target2='gsl_matrix_set(Jx, ', row=row, col=col))
			
		} # end discrete time ifelse
	  #browser()
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

cswapDynamicsFormula <- function(lhs, index, rhs, vtarget, close){
	return(gsub(paste0("\\<", lhs, "\\>"), paste0(vtarget, index, close), rhs))
}

cswapDynamicsFormulaLoop <- function(nregime, n, lhs, rhs, target1, close, target2, vector=FALSE, row=NULL, col=NULL){
	ans <- paste0("\n\t", "switch (regime) {")
	for (r in 1:nregime){
		ans <- paste(ans, paste0("\tcase ", r-1, ":"), sep="\n\t")
		for (i in 1:n[r]){
			for (j in 1:length(lhs[[r]])){
				rhs[[r]][[i]] <- cswapDynamicsFormula(lhs=lhs[[r]][[j]], index=j-1, rhs=rhs[[r]][[i]], vtarget=target1, close=close)
			}
			if(vector){
				ans <- paste(ans, paste0("\t\t", target2, i-1, ", ", rhs[[r]][[i]], ");"), sep="\n\t")
			} else {
				ans <- paste(ans, paste0("\t\t", target2, which(lhs[[r]] == row[[r]][[i]])-1, ", ", which(lhs[[r]] == col[[r]][[i]])-1, ", ", rhs[[r]][[i]], ");"), sep="\n\t")
			}
		}
		ans <- paste(ans, paste0("\t\tbreak;"), sep="\n\t")
	}
	ans <- paste(ans, paste0("}"), sep="\n\t")
	
	ans <- paste0(ans, "\n}")
	return(ans)
}

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
				ifelse(hasIntercepts, destroyGslVector("intVector"), ""),
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
					destroyGslVector("intVector"))
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
		deviation <- object$deviation
		refRow <- object$refRow
		
		if(deviation){
			if(nrow(values)!=0 && nrow(params)!=0){
				# I need to look into this one
				interceptSel <- seq(1, ncol(values), by=numCovariates+1)
				intercept.values <- matrix(values[refRow, interceptSel], numRegimes, 1)
				intercept.params <- matrix(params[refRow, interceptSel], numRegimes, 1)
				values[refRow, interceptSel] <- 0
				params[refRow, interceptSel] <- 0
			}
		} else {
			intercept.values <- matrix(0, numRegimes, 1)
			intercept.params <- matrix(0, numRegimes, 1)
		}
		
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
			             createGslVector(numRegimes, "Pintercept"),
			             setGslVectorElements(values=intercept.values, params=intercept.params, name="Pintercept"),
			             sep="\n")
			for(reg in 1L:numRegimes){
				selRows <- rowBeginSeq[reg]:rowEndSeq[reg]
				ret <- paste(ret,
					setGslVectorElements(values=values[selRows, 1, drop=FALSE], params=params[selRows, 1, drop=FALSE], name="Pvector"),
					setGslMatrixElements(values=values[selRows, -1, drop=FALSE], params=params[selRows, -1, drop=FALSE], name="Gmatrix"),
					blasMV(FALSE, "1.0", "Gmatrix", "covariate_local", "1.0", "Pvector"),
					"\tgsl_vector_add(Pvector, Pintercept);",
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
				destroyGslVector("Pintercept"),
				destroyGslVector("covariate_local"),
				sep="\n")
		} else if(nrow(values) != 0 && nrow(params) != 0 && numCovariates == 0){
			# same as above processing, but without the Gmatrix for the covariates
			ret <- paste(ret,
				createGslVector(numRegimes, "Pvector"),
				createGslVector(numRegimes, "Presult"),
				createGslVector(numRegimes, "Pintercept"),
				setGslVectorElements(values=intercept.values, params=intercept.params, name="Pintercept"),
				sep="\n")
			for(reg in 1L:numRegimes){
				selRows <- rowBeginSeq[reg]:rowEndSeq[reg]
				ret <- paste(ret,
					setGslVectorElements(values=values[selRows, 1, drop=FALSE], params=params[selRows, 1, drop=FALSE], name="Pvector"),
					"\tgsl_vector_add(Pvector, Pintercept);",
					"\tmathfunction_softmax(Pvector, Presult);",
					gslVector2Column("regime_switch_mat", reg-1, "Presult", 'row'),
					"\tgsl_vector_set_zero(Pvector);",
					sep="\n")
			}
			ret <- paste(ret,
				destroyGslVector("Pvector"),
				destroyGslVector("Presult"),
				destroyGslVector("Pintercept"),
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
		covariates_local <- object$covariates
		numCovariates <- length(covariates_local)
		hasCovariates <- numCovariates > 0
		deviation <- object$deviation
		refRow <- object$refRow
		
		
		nregime <- max(length(values.inistate), length(values.inicov), nrow(values.regimep))
		if(nregime != 1 && length(values.inistate) == 1){
			values.inistate <- rep(values.inistate, nregime)
			params.inistate <- rep(params.inistate, nregime)
		}
		if(nregime != 1 && length(values.inicov) == 1){
			values.inicov <- rep(values.inicov, nregime)
			params.inicov <- rep(params.inicov, nregime)
		}
		someStatesNotZero <- sapply(values.inistate, function(x){!all(x == 0)})
		someStatesNotZero2 <- sapply(params.inistate, function(x){!all(x == 0)})
		someStatesNotZero <- someStatesNotZero | someStatesNotZero2
		someStatesNotZero <- someStatesNotZero | rep(TRUE, length(someStatesNotZero))
		
		numLatent <- nrow(values.inistate[[1]])
		
		values.etaIntercept <- lapply(values.inistate, function(x){x[,1, drop=FALSE]})
		params.etaIntercept <- lapply(params.inistate, function(x){x[,1, drop=FALSE]})
		values.covEffects <- lapply(values.inistate, function(x){x[,-1, drop=FALSE]})
		params.covEffects <- lapply(params.inistate, function(x){x[,-1, drop=FALSE]})
		
		if(deviation){
			if(nrow(values.regimep)!=0 && nrow(params.regimep)!=0){
				values.regIntercept <- matrix(values.regimep[refRow, 1], nrow=nregime, 1)
				params.regIntercept <- matrix(params.regimep[refRow, 1], nrow=nregime, 1)
				values.regimep[refRow, 1] <- 0
				params.regimep[refRow, 1] <- 0
			}
		} else {
			values.regIntercept <- matrix(0, nregime, 1)
			params.regIntercept <- matrix(0, nregime, 1)
		}
		
		ret <- "void function_initial_condition(double *param, gsl_vector **co_variate, gsl_vector **pr_0, gsl_vector **eta_0, gsl_matrix **error_cov_0, size_t *index_sbj){\n\t"
		ret <- paste0(ret, "\n", createGslVector(nregime, "Pvector"))
		ret <- paste0(ret, createGslVector(nregime, "Pintercept"))
		ret <- paste0(ret, createGslVector(nregime, "Padd"))
		ret <- paste0(ret, createGslVector(nregime, "Preset"))
		ret <- paste0(ret, setGslVectorElements(values.regimep[ , 1, drop=FALSE], params.regimep[ , 1, drop=FALSE], "Pvector"))
		ret <- paste0(ret, setGslVectorElements(values.regIntercept, params.regIntercept, "Pintercept"))
		ret <- paste0(ret, "\tgsl_vector_add(Padd, Pvector);\n")
		ret <- paste0(ret, "\tgsl_vector_add(Padd, Pintercept);\n")
		ret <- paste0(ret, "\tgsl_vector_add(Preset, Pvector);\n")
		ret <- paste0(ret, "\tgsl_vector_add(Preset, Pintercept);\n")
		
		ret <- paste0(ret, createGslVector(numLatent, "eta_local"))
		if(hasCovariates){
			ret <- paste0(ret, createGslVector(numCovariates, "covariate_local"))
			ret <- paste0(ret, createGslMatrix(numLatent, numCovariates, "Cmatrix"))
			ret <- paste0(ret, createGslMatrix(nregime, numCovariates, "Pmatrix"))
			ret <- paste0(ret, setGslMatrixElements(values=values.regimep[ , -1, drop=FALSE], params=params.regimep[ , -1, drop=FALSE], name="Pmatrix"))
			ret <- paste0(ret, "\t\n")
		}
		
		ret <- paste0(ret,"\tsize_t num_regime=pr_0[0]->size;\n\tsize_t dim_latent_var=error_cov_0[0]->size1;")
		if (any(someStatesNotZero)){
			ret <- paste0(ret,"\n\tsize_t num_sbj=(eta_0[0]->size)/(dim_latent_var);\n\tsize_t i;")
		}
		ret <- paste0(ret,"\n\tsize_t regime;\n")
		
		if(nregime > 1){
			ret <- paste0(ret, "\tfor(regime=0; regime < num_regime; regime++){")
			ret <- paste(ret, "\t\tswitch (regime) {\n", sep="\n")
				for(reg in 1:nregime){
					ret <- paste0(ret, paste0("\t\t\tcase ", reg-1, ":"), "\n")
					if (any(someStatesNotZero[reg])){
						ret <- paste0(ret,"\t\t\t\tfor(i=0; i < num_sbj; i++){\n")
						ret <- paste0(ret, setGslVectorElements(values=values.etaIntercept[[reg]], params=params.etaIntercept[[reg]], name='eta_local', depth=5))
						if(hasCovariates){
							ret <- paste0(ret, gslVectorCopy("co_variate[index_sbj[i]]", "covariate_local", fromLoc=match(covariates_local, covariates), toLoc=1:numCovariates, depth=5))
							ret <- paste0(ret, setGslMatrixElements(values=values.covEffects[[reg]], params=params.covEffects[[reg]], name="Cmatrix", depth=5))
							ret <- paste0(ret, "\t\t\t\t", blasMV(FALSE, "1.0", "Cmatrix", "covariate_local", "1.0", "eta_local"))
						}
						ret <- paste0(ret, gslVectorCopy("eta_local", "eta_0[regime]", 1:numLatent, 1:numLatent, toFill="i*dim_latent_var+", depth=5))
						ret <- paste0(ret, "\t\t\t\t\tgsl_vector_set_zero(eta_local);\n")
						if(hasCovariates){
							ret <- paste0(ret, "\t\t\t\t\tgsl_matrix_set_zero(Cmatrix);\n")
						}
						ret <- paste0(ret,"\t\t\t\t}") # close i loop
					}
					ret <- paste(ret, setGslMatrixElements(values.inicov[[reg]], params.inicov[[reg]], "(error_cov_0)[regime]", depth=4), sep="\n")
					ret <- paste0(ret, "\t\t\t", "break;", "\n")
				}
			ret <- paste0(ret, "\t\t", "}\n") # close case switch
			ret <- paste0(ret, "\t}") # close regime loop
		} else{
			ret <- paste0(ret, "\tfor(regime=0; regime < num_regime; regime++){")
			if (any(someStatesNotZero)){
				ret <- paste0(ret,"\n\t\tfor(i=0; i < num_sbj; i++){\n")
				ret <- paste0(ret, setGslVectorElements(values=values.etaIntercept[[1]], params=params.etaIntercept[[1]], name='eta_local', depth=3))
				if(hasCovariates){
					ret <- paste0(ret, gslVectorCopy("co_variate[index_sbj[i]]", "covariate_local", fromLoc=match(covariates_local, covariates), toLoc=1:numCovariates, depth=3))
					ret <- paste0(ret, setGslMatrixElements(values=values.covEffects[[1]], params=params.covEffects[[1]], name="Cmatrix", depth=3))
					ret <- paste0(ret, "\t\t", blasMV(FALSE, "1.0", "Cmatrix", "covariate_local", "1.0", "eta_local"))
				}
				ret <- paste0(ret, gslVectorCopy("eta_local", "eta_0[regime]", 1:numLatent, 1:numLatent, toFill="i*dim_latent_var+", depth=3))
				ret <- paste0(ret, "\t\t\tgsl_vector_set_zero(eta_local);\n")
				if(hasCovariates){
					ret <- paste0(ret, "\t\t\tgsl_matrix_set_zero(Cmatrix);\n")
				}
				ret <- paste0(ret,"\t\t}") # close i loop
			}
			ret <- paste(ret, setGslMatrixElements(values.inicov[[1]], params.inicov[[1]], "(error_cov_0)[regime]"), sep="\n")
			ret <- paste0(ret, "\t}") # close regime loop
		}
		
		ret <- paste0(ret,"\n\tfor(i=0; i < num_sbj; i++){\n")
		if(hasCovariates){
			ret <- paste0(ret, gslVectorCopy("co_variate[index_sbj[i]]", "covariate_local", fromLoc=match(covariates_local, covariates), toLoc=1:numCovariates, depth=2))
			ret <- paste0(ret, "\t", blasMV(FALSE, "1.0", "Pmatrix", "covariate_local", "1.0", "Padd"))
		}
		ret <- paste0(ret, "\t\tmathfunction_softmax(Padd, pr_0[i]);\n")
		if(hasCovariates){
			ret <- paste0(ret, "\t\tgsl_vector_memcpy(Padd, Preset);\n")
		}
		ret <- paste0(ret,"\t}") # close i loop
		
		ret <- paste0(ret, "\n", destroyGslVector("Pvector"), destroyGslVector("Pintercept"), destroyGslVector("Padd"), destroyGslVector("Preset"), destroyGslVector("eta_local"))
		if(hasCovariates){
			ret <- paste0(ret, destroyGslVector("covariate_local"))
			ret <- paste0(ret, destroyGslMatrix("Cmatrix"))
			ret <- paste0(ret, destroyGslMatrix("Pmatrix"))
		}
		ret <- paste0(ret, "}\n") #Close function definition
		object@c.string <- ret
		
		return(object)
	}
)

setMethod("writeCcode", "dynrNoise",
	function(object, covariates){
	    #browser()
		params.latent <- object$params.latent
		params.observed <- object$params.observed
		#Note: use mutated values, instead of raw values
		values.latent <- object@values.latent.inv.ldl
		values.observed <- object@values.observed.inv.ldl
		
		nregime <- length(values.latent)
		
		#browser()
		# not formula --------------
		if (length(object@latent.formula) == 0){
			ret <- "void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov, const gsl_vector *co_variate){\n\n"
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
		}
		# formula --------------
		else{	
		    # this chunk may not work for multi-regime
		    ldl_list <- vectorizeMatrix(object$ldl.transformed, byrow= FALSE)
			fml=processFormula(ldl_list)
			lhs=lapply(fml,function(x){x[1]})
            rhs=lapply(fml,function(x){x[2]})
			
			fkv = processFormula(object$known.vars)
			rhskv = lapply(fkv,function(x){x[2]})
			lhskv = names(object$known.vars)
			
			
			for (i in 1:length(object$paramnames)){
			  pattern <- object$paramnames[i]
			  rhs  <- lapply(rhs, function(x){gsub(pattern, paste0('param[', i-1, ']'),x, fixed = TRUE)})
			  rhskv <- lapply(rhskv, function(x){gsub(pattern, paste0('param[', i-1, ']'),x, fixed = TRUE)})
			  lhskv <- gsub(pattern, paste0('param[', i-1, ']'), lhskv , fixed = TRUE)
			}
			
			for (i in 1:length(covariates)){
			  pattern <- covariates[i]
			  rhs  <- lapply(rhs, function(x){gsub(pattern, paste0('gsl_vector_get(co_variate,', i-1, ')'),x, fixed = TRUE)})
			  rhskv  <- lapply(rhskv, function(x){gsub(pattern, paste0('gsl_vector_get(co_variate,', i-1, ')'),x, fixed = TRUE)})
			}
			
			#browser()
			for(i in 1:nrow(params.latent[[1]])){
			  for(j in 1:ncol(params.latent[[1]])){
			    pattern <- lhs[(i-1)*nrow(params.latent[[1]])+j][[1]]
			    rhs  <- lapply(rhs, function(x){gsub(pattern, paste0('noise_formula[', i-1,'][', j-1, ']'),x, fixed = TRUE)})
			    rhskv  <- lapply(rhskv, function(x){gsub(pattern, paste0('noise_formula[', i-1,'][', j-1, ']'),x, fixed = TRUE)})
			  }
			}

			# writeC code here
			ret <- "void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov, const gsl_vector *co_variate){\n\n"
			ret <- paste0(ret, "\tdouble noise_formula [", nrow(object$ldl.transformed), "][", ncol(object$ldl.transformed) ,"];\n\n")
			for(i in 1:nrow(object$ldl.transformed)){
			  for(j in 1:ncol(object$ldl.transformed)){
			    ret <- paste0(ret, "\tnoise_formula[", i-1,'][', j-1, ']=', rhs[(j-1)*2+i] ,";\n")
			  }
			}
			ret <- paste0(ret, '\n')
			ret <- paste0(ret, 'printf("noise_formula\\n%.4lf %.4lf\\n%.4lf %.4lf\\n", noise_formula[0][0], noise_formula[0][1], noise_formula[1][0], noise_formula[1][1]);\nprintf("param3-5 %lf %lf %lf\\n", param[3], param[4], param[5]);')
			#for(i in 1:length(lhskv)){
			#  ret <- paste0(ret, "\t", lhskv[i], '=', rhskv[i] ,";\n")
			#}
			
			
			
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
			    ret<- paste(ret, "\n")
				for(i in 1:nrow(object$ldl.transformed)){
				  for(j in 1:ncol(object$ldl.transformed)){
					#gsl_matrix_set(eta_noise_cov, 0, 0, param[0]);
					ret <- paste0(ret, "\tgsl_matrix_set(eta_noise_cov,", i-1, ",",  j-1, ", noise_formula[", i-1, "][",  j-1,"]);\n")
				  }
				}
				ret <- paste(ret, setGslMatrixElements(values.observed[[1]], params.observed[[1]], "y_noise_cov"), sep="\n")
			}
			ret <- paste(ret, "\n}\n\n") # close C function
		
		}
		
		object@c.string <- ret
		
		#browser()
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
			#browser()
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
			#browser()
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
  #browser()
  params.numbers <- lapply(params, .exchangeNamesAndNumbers, param.data$param.name)
  values <- PopBackMatrix(values, params.numbers, paste0("vec[",param.data$param.number,"]"))
  char <- character(0)
  char.vec <- character(0)
  if (length(values) > 0){
    for (i in 1:length(values)){
      param.number <- unique(params.numbers[[i]][which(params.numbers[[i]]!=0,arr.ind = TRUE)])
      if (length(param.number)!=0){
        vec.noise <- paste0("vec[",paste0("c(",paste(param.number,collapse=","),")"),"]")
        vec.sub <- as.vector(values[[i]])
        mat.index <- sapply(param.number,function(x){min(which(params.numbers[[i]]==x))})
        char.i <- paste0(vec.noise,"=as.vector(", ifelse(reverse, 'reverseldl', 'transldl'), "(matrix(", paste0("c(",paste(vec.sub,collapse=","),")"),",ncol=",ncol(values[[i]]),")))[",
                      paste0("c(",paste(mat.index,collapse=","),")"),"]", ifelse(i!=length(values),"\n\t",""))
        char.vec <- c(char.vec, char.i)
      }
    }
    # Check for and remove any duplicates
    dups <- duplicated(gsub('\n\t$', '', char.vec))
    if(length(char.vec) > 1 && any(dups)){
      char.vec <- char.vec[!dups]
    }
    char <- paste0(char.vec, collapse='')
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
#transldl function for calculating the LDL values
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
##' 
##' @return A matrix
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
		if(all(is.na(values))){
			# if everything is NA, then we're giving up.
			return(values)
		} else if(is.matrix(values) && all(is.na(values[lower.tri(values, diag=FALSE)]))){
			# if it's a matrix and all the lower triangular parts are NA
			# then we're just setting bounds on the variances
			values[lower.tri(values, diag=FALSE)] <- 0
			values[upper.tri(values, diag=FALSE)] <- 0
			mat <- dynr.ldl(values)
			diag(mat) <- log(diag(mat))
			mat[lower.tri(mat, diag=FALSE)] <- NA
			return(mat)
		} else if(is.matrix(values) && all(is.na(diag(values))) && all( values[lower.tri(values, diag=FALSE)] == 0)){
			# if the matrix has NA diagonal and is otherwise 0
			# then leave it alone
			return(values)
		} else{
			# only warn when values is
			#  1 not all missing
			#  2 not missing everywhere except the diagonal
			#  3 not missing on the diagonal with zero everywhere else
			warning("Avast ye lowly dog! NA was passed to LDL in confusing way. Not doing LDL.")
			return(values)
		}
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
##' We create a new method for \code{diag} with character input.  The default behavior for missing \code{nrow} and/or \code{ncol} arguments is the same
##' as for the \code{\link{diag}} function in the base package.  Off-diagonal entries
##' are filled with "0".
##' 
##' @return A matrix
##' 
##' @examples
##' diag(letters[1:3])
setMethod("diag", "character",
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
preProcessValues <- function(x, rowNam, colNam){
	x[is.na(x)] <- 'freed'
	if(is.null(dim(x))){
		numRow <- length(x)
		numCol <- 1
		if(missing(rowNam)) rowNam <- NULL
		if(missing(colNam)) colNam <- NULL
	} else {
		numRow <- nrow(x)
		numCol <- ncol(x)
		if(missing(rowNam)) rowNam <- rownames(x)
		if(missing(colNam)) colNam <- colnames(x)
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

preProcessParams <- function(x, rowNam, colNam){
	if (!is.null(x)) {x[is.na(x)] <- 'fixed'}
	if(is.null(dim(x))){
		numRow <- length(x)
		numCol <- 1
		if(missing(rowNam)) rowNam <- NULL
		if(missing(colNam)) colNam <- NULL
	} else {
		numRow <- nrow(x)
		numCol <- ncol(x)
		if(missing(rowNam)) rowNam <- rownames(x)
		if(missing(colNam)) colNam <- colnames(x)
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
	vdim <- sapply(lapply(values, as.matrix), dim)
	pdim <- sapply(lapply(params, as.matrix), dim)
	if(length(vdim) > 0 && any(apply(vdim, 1, function(x) length(unique(x)) > 1))){
		stop("Some of the 'values' list elements are not the same size as each other\nNot cool, Donny.")
	}
	if(length(pdim) > 0 && any(apply(pdim, 1, function(x) length(unique(x)) > 1))){
		stop("Some of the 'params' list elements are not the same size as each other\nNo-go for launch.")
	}
	if(any(vdim != pdim)){
		stop("'values' and 'params' are not all the same size.\nWalter Sobchak says you can't do that.")
	}
	return(list(values=values, params=params))
}

symmExtract <- function(x, s){if(s){x[lower.tri(x, diag=TRUE)]}else{x}}

checkMultipleStart <- function(values, params, symmetric=FALSE){
	if(is.list(values) & is.list(params)){
		values <- unlist(lapply(values, symmExtract, s=symmetric))
		params <- unlist(lapply(params, symmExtract, s=symmetric))
	} else if(is.matrix(values) && is.matrix(params)){
		values <- c(symmExtract(values, symmetric))
		params <- c(symmExtract(params, symmetric))
	} else if(xor(is.list(values), is.list(params)) || xor(is.matrix(values), is.matrix(params))){
		stop("Invalid input to checkMultipleStart() function.\nMust be a pair of vectors, a pair of lists, or a pair of matrices.\nFound a mix.")
	}
	names(values) <- params
	chparam <- setdiff(unique(params), c("0", "fixed")) # Don't check parameters labeled "fixed" or "0"
	for(i in chparam){
		ch <- values[ names(values) %in% i]
		if(sum(!duplicated(ch)) > 1){
			stop(paste0("Found multiple (transformed) start values for parameter '", i, "': ", paste(ch, sep="", collapse=", ")), call.=FALSE)
		}
	}
}


preProcessNames <- function(x, rnames=character(0), cnames=character(0), xname=character(0), rarg=character(0), carg=character(0)){
	d <- dim(x)
	if(d[1] != length(rnames)){
		stop(paste0("Matrix ", xname, " has ", d[1], " rows and ", length(rnames), " ", rarg, " (i.e., row names).\nHint: they should match.\nNot even the King of Pop could set these row names."), call.=FALSE)
	}
	rownames(x) <- rnames
	if(!(length(cnames)==0 & length(carg)==0) && d[2] != length(cnames)){
	stop(paste0("Matrix ", xname, " has ", d[2], " columns and ", length(cnames), " ", carg, " (i.e., column names).\nHint: they should match.\nNot even the King of Soul could set these column names."), call.=FALSE)
	}
	if(!length(cnames)==0 & length(carg)==0) colnames(x) <- cnames
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

extractValues <- function(v, p, symmetric=FALSE){
	if(is.list(v) && is.list(p) && length(v) == length(p) && length(v) > 0){
		ret <- c()
		for(i in 1:length(v)){
			ret <- c(ret, extractValues(v[[i]], p[[i]], symmetric))
		}
		return(ret)
	} else if(length(v) == 0){
		return(numeric(0))
	}else {
		checkMultipleStart(v, p, symmetric)
		return(v[extractWhichParams(p)])
	}
}

#------------------------------------------------------------------------------
# Create recipe for function measurement


#--------------------------------------
# brief input version

##' Recipe function to quickly create factor loadings
##'
##' @param map list giving how the latent variables map onto the observed variables
##' @param params parameter numbers
##' @param idvar names of the variables used to identify the factors
##' @param exo.names names of the exogenous covariates
##' @param intercept logical. Whether to include freely esimated intercepts
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
##' the dynr package. In particular, it does not have any regime-swtiching.
##' Covariates can be included with the \code{exo.names} argument, but
##' all covariate effects are freely estimated and the starting values
##' are all zero.  Likewise, intercepts can be included with the \code{intercept}
##' logical argument, but all intercept terms are freely estimated with 
##' zero as the starting value.
##' For complete functionality use \code{\link{prep.measurement}}.
##' 
##' @return Object of class 'dynrMeasurement'
##' 
##' @examples
##' #Single factor model with one latent variable fixing first loading
##' prep.loadings(list(eta1=paste0('y', 1:4)), paste0("lambda_", 2:4))
##'
##' #Single factor model with one latent variable fixing the fourth loading
##' prep.loadings(list(eta1=paste0('y', 1:4)), paste0("lambda_", 1:3), idvar='y4')
##' 
##' #Single factor model with one latent variable freeing all loadings
##' prep.loadings(list(eta1=paste0('y', 1:4)), paste0("lambda_", 1:4), idvar='eta1')
##' 
##' #Single factor model with one latent variable fixing first loading
##' # and freely estimated intercept
##' prep.loadings(list(eta1=paste0('y', 1:4)), paste0("lambda_", 2:4),
##'  intercept=TRUE)
##' 
##' #Single factor model with one latent variable fixing first loading
##' # and freely estimated covariate effects for u1 and u2
##' prep.loadings(list(eta1=paste0('y', 1:4)), paste0("lambda_", 2:4),
##'  exo.names=paste0('u', 1:2))
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
prep.loadings <- function(map, params=NULL, idvar, exo.names=character(0), intercept=FALSE){
	if(missing(idvar)){
		idvar <- sapply(map, '[', 1)
	}
	
	allVars <- unique(unlist(map))
	
	nx <- length(allVars)
	ne <- length(map)
	nu <- length(exo.names)
	
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
	
	if(nu > 0){
		exoVal <- matrix(0, nx, nu)
		exoPar <- outer(paste0("b_", allVars), exo.names, paste0)
	} else {
		exoVal <- NULL
		exoPar <- NULL
	}
	
	if(intercept){
		intVal <- matrix(0, nx, 1)
		intPar <- matrix(paste0("int_", allVars), nx, 1)
	}else{
		intVal <- NULL
		intPar <- NULL
	}
	
	rownames(valuesMat) <- allVars
	rownames(paramsMat) <- allVars
	colnames(valuesMat) <- names(map)
	colnames(paramsMat) <- names(map)
	x <- prep.measurement(values.load=valuesMat, params.load=paramsMat, values.exo=exoVal, params.exo=exoPar, values.int=intVal, params.int=intPar, state.names=names(map), obs.names=allVars, exo.names=exo.names)
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
##' @details
##' The values.* arguments give the starting and fixed values for their respective matrices.
##' The params.* arguments give the free parameter labels for their respective matrices.
##' Numbers can be used as labels.
##' The number 0 and the character 'fixed' are reserved for fixed parameters.
##'
##' When a single matrix is given to values.*, that matrix is not regime-switching.
##' Correspondingly, when a list of length r is given, that matrix is regime-switching with values and params for the r regimes in the elements of the list.
##' 
##' @return Object of class 'dynrMeasurement'
##' 
##' @seealso 
##' Methods that can be used include: \code{\link{print}}, \code{\link{printex}}, \code{\link{show}} 
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
	if(missing(values.int) && !missing(params.int)){
		warning('Intercept parameters are specified in params.int without any initial values. Using 0s as the default starting values. To change, please use values.int.')
		values.int = matrix(rep(0, times = nrow(params.int)), ncol=1)
	}
    r <- coProcessValuesParams(values.int, params.int, missingOK=TRUE)
    values.int <- r$values
    params.int <- r$params
    

	if(missing(obs.names)){
		obs.names <- paste0('y',1:nrow(values.load[[1]]))
	}
	
	if(missing(state.names)){
		state.names <- paste0('state',1:ncol(values.load[[1]]))
	}
	
	if(missing(exo.names)){
		exo.names = character(0)
	}
	if(any(duplicated(obs.names))){
		stop("'obs.names' uses some of the same names more than once.\n  I repeat myself when under stress.  I repeat myself when under stress.  I repeat myself when under stress.\nI repeat myself when under stress.  I repeat myself when under stress.")
	}
	if(any(duplicated(state.names))){
		stop("'state.names' uses some of the same names more than once.\n  Sensei says this lacks discipline.")
	}
	if(any(duplicated(exo.names))){
		stop("'exo.names' uses some of the same names more than once.\n  You cannot be the King of Crimson with this indiscipline.")
	}
	
	values.load <- lapply(values.load, preProcessNames, obs.names, state.names, 'values.load', 'obs.names', 'state.names')
	#params.load <- lapply(params.load, preProcessNames, obs.names, state.names, 'values.load', 'obs.names', 'state.names')
	values.exo <- lapply(values.exo, preProcessNames, obs.names, exo.names, 'values.exo', 'obs.names', 'exo.names')
	#params.exo <- lapply(params.exo, preProcessNames, obs.names, exo.names, 'values.exo', 'obs.names', 'exo.names')
	values.int <- lapply(values.int, preProcessNames, obs.names, character(0), 'values.int', 'obs.names')
	#params.int <- lapply(params.int, preProcessNames, obs.names, 'values.int', 'obs.names')
	
	values.load <- lapply(values.load, preProcessValues, rowNam=obs.names, colNam=state.names)
	params.load <- lapply(params.load, preProcessParams, rowNam=obs.names, colNam=state.names)
	values.exo <- lapply(values.exo, preProcessValues, rowNam=obs.names, colNam=exo.names)
	params.exo <- lapply(params.exo, preProcessParams, rowNam=obs.names, colNam=exo.names)
	values.int <- lapply(values.int, preProcessValues, rowNam=obs.names, colNam='one')
	params.int <- lapply(params.int, preProcessParams, rowNam=obs.names, colNam='one')
	
	# Check that all 'values' imply 0, 1, or the same number of regimes.
	# Note that the 'values' and 'params' have already been checked to imply this.
	nregs <- sapply(list(values.load=values.load, values.exo=values.exo, values.int=values.int), length)
	msg <- paste0("Y'all iz trippin! Different numbers of regimes implied:\n'load' has ", nregs[1], ", 'exo' has ",  nregs[2], ", and 'int' has ", nregs[3], " regimes.")
	mr <- max(nregs)
	if(!all(nregs %in% c(0, 1, mr))){
		stop(msg)
	} else if (!all(nregs %in% c(0, mr))){
		if(nregs[1] == 1){
			ae <- autoExtendSubRecipe(values.load, params.load, 'values.load', 'loadings', mr)
			values.load <- ae[[1]]
			params.load <- ae[[2]]
		}
		if(nregs[2] == 1){
			ae <- autoExtendSubRecipe(values.exo, params.exo, 'values.exo', 'covariate regressions', mr)
			values.exo <- ae[[1]]
			params.exo <- ae[[2]]
		}
		if(nregs[3] == 1){
			ae <- autoExtendSubRecipe(values.int, params.int, 'values.int', 'intercepts', mr)
			values.int <- ae[[1]]
			params.int <- ae[[2]]
		}
	}
	
	sv <- c(extractValues(values.load, params.load), extractValues(values.exo, params.exo), extractValues(values.int, params.int))
	pn <- c(extractParams(params.load), extractParams(params.exo), extractParams(params.int))
	sv <- extractValues(sv, pn)
	pn <- extractParams(pn)
	
	x <- list(startval=sv, paramnames=pn, values.load=values.load, params.load=params.load,
		values.exo=values.exo, params.exo=params.exo, values.int=values.int, params.int=params.int,
		obs.names=obs.names, state.names=state.names, exo.names=exo.names)
	return(new("dynrMeasurement", x))
}

regime1msg <- "\nEven non-regime-switching parts of a recipe must match in their numbers of regimes.\nE.g., use rep(list(blah), 3) to make 'blah' repeat 3 times in a list."

autoExtendSubRecipe <- function(values, params, formalName, informalName, maxReg){
	v <- rep(values, maxReg)
	p <- rep(params, maxReg)
	message(paste0("Oi, Chap! I found 1 regime for '",  formalName, "'  but ", maxReg, " regimes elsewhere, so I extended the ", informalName, " to match.\nIf this is what you wanted, all is sunshine and puppy dogs."))
	return(list(v, p))
}

#------------------------------------------------------------------------------
# Error covariance matrix
# N.B. This function produces BOTH the latent and observed error covariance matrices.


##' Recipe function for specifying the measurement error and process noise covariance structures
##' 
##' @param values.latent a positive definite matrix or a list of positive definite matrices of the starting or fixed values of the process noise covariance structure(s) in one or more regimes. If only one matrix is specified for a regime-switching dynamic model, the process noise covariance structure stays the same across regimes. To ensure the matrix is positive definite in estimation, we apply LDL transformation to the matrix. Values are hence automatically adjusted for this purpose.
##' @param params.latent a matrix or list of matrices of the parameter names that appear in the process noise covariance(s) in one or more regimes. If an element is 0 or "fixed", the corresponding element is fixed at the value specified in the values matrix; Otherwise, the corresponding element is to be estimated with the starting value specified in the values matrix. If only one matrix is specified for a regime-switching dynamic model, the process noise structure stays the same across regimes. If a list is specified, any two sets of the parameter names as in two matrices should be either the same or totally different to ensure proper parameter estimation.  See Details.
##' @param values.observed a positive definite matrix or a list of positive definite matrices of the starting or fixed values of the measurement error covariance structure(s) in one or more regimes. If only one matrix is specified for a regime-switching measurement model, the measurement noise covariance structure stays the same across regimes. To ensure the matrix is positive definite in estimation, we apply LDL transformation to the matrix. Values are hence automatically adjusted for this purpose. 
##' @param params.observed a matrix or list of matrices of the parameter names that appear in the measurement error covariance(s) in one or more regimes. If an element is 0 or "fixed", the corresponding element is fixed at the value specified in the values matrix; Otherwise, the corresponding element is to be estimated with the starting value specified in the values matrix. If only one matrix is specified for a regime-switching dynamic model, the process noise structure stays the same across regimes. If a list is specified, any two sets of the parameter names as in two matrices should be either the same or totally different to ensure proper parameter estimation.  See Details.
##' @param ... Further named arguments.  Currently we only accept 'covariates' and 'latent.formula'.
##' 
##' @details
##' The arguments of this function should generally be either matrices or lists of matrices.  Lists of matrices are used for regime-switching models with each list element corresponding to a regime.  Thus, a list of three matrices implies a three-regime model.  Single matrices are for non-regime-switching models.  Some checking is done to ensure that the number of regimes implied by one part of the model matches that implied by the others.  For example, the noise model (\code{prep.noise}) cannot suggest three regimes when the measurement model (\code{\link{prep.measurement}}) suggests two regimes. An exception to this rule is single-regime (i.e. non-regime-switching) components.  For instance, the noise model can have three regimes even though the measurement model implies one regime.  The single-regime components are simply assumed to be invariant across regimes.
##' 
##' Care should be taken that the parameters names for the latent covariances do not overlap with the parameters in the observed covariances.  Likewise, the parameter names for the latent covariances in each regime should either be identical or completely distinct. Because the LDL' transformation is applied to the covariances, sharing a parameter across regimes may cause problems with the parameter estimation.
##' 
##' Use $ to show specific arguments from a dynrNoise object (see examples).
##' 
##' @return Object of class 'dynrNoise'
##' 
##' @seealso 
##' \code{\link{printex}} to show the covariance matrices in latex.
##'  
##' @examples 
##' # Two latent variables and one observed variable in a one-regime model
##' Noise <- prep.noise(values.latent=diag(c(0.8, 1)),
##'     params.latent=diag(c('fixed', "e_x")), 
##'     values.observed=diag(1.5,1), params.observed=diag("e_y", 1))
##' # For matrices that can be import to latex:
##' printex(Noise, show=TRUE)
##' # If you want to check specific arguments you've specified, for example,
##' # values for variance structure of the latent variables
##' Noise$values.latent
##' 
##' # Two latent variables and one observed variable in a two-regime model
##' Noise <- prep.noise(values.latent=list(diag(c(0.8, 1)), diag(c(0.8, 1))),
##'     params.latent=list(diag(c('fixed', "e_x1")), diag(c('fixed', "e_x2"))),
##'     values.observed=list(diag(1.5,1), diag(0.5,1)),
##'     params.observed=list(diag("e_y1", 1), diag("e_y2",1)))
##' # If the error and noise structures are assumed to be the same across regimes,
##' #  it is okay to use matrices instead of lists.
prep.noise <- function(values.latent, params.latent, values.observed, params.observed, ...){
  #browser()
  # ---- To incorporate covariates and formulas into covariance functions ----
  if (!missing(...)){
     dots <- list(...)
  }
  if ((missing(values.latent) || missing(params.latent) || missing(values.observed) || missing(params.observed)) && 
      (length(dots==0)))
    stop("You have to provide the noise structure as lists of matrices or formulas. Neither is available now.")
  #if(length(dots) > 0 && is.list(dots$latent.formula) && length(dots$latent.formula) > 0){
  if (!missing(...)){
    if(!all(names(dots) %in% c('covariates', 'latent.formula', 'state.names', 'latent.startval'))){
      stop("You passed some invalid names to the ... argument. Check the ?prep.noise help page for more information.")
    }
    
    # if the argument formula is not a list 
    if(!is.list(dots$latent.formula)){
      msg <- paste0(ifelse(plyr::is.formula(dots$latent.formula), "'latent.formula' argument is a formula but ", ""),
                    "'formula'",
                    ifelse(plyr::is.formula(dots$latent.formula), " ", " argument "),
                    "should be a list of formulas.\nCan't nobody tell me nothin'")
      stop(msg)
    }
	#browser()
	x_constant <- processCovConstant(values.latent, params.latent, values.observed, params.observed, exclulde_latent = TRUE)
	x_formula <- processCovFormula(dots, values.latent, params.latent)
    x <- c(x_constant, x_formula)
	x$paramnames <- c(x_constant$paramnames, x_formula$paramnames)
	x$startval <- c(x_constant$startval, x_formula$startval)
	x$is_eta_cov_formula <- TRUE

    # ---- End of Cov formula modifications ----  
  } else{
    # Else process cov-related things the usual way
    x <- processCovConstant(values.latent, params.latent, values.observed, params.observed)
	x$is_eta_cov_formula <- FALSE
  }
  # Handle latent covariance
  return(new("dynrNoise", x))
}


#SMC 9/28/22
#Generally for a symmetric cov matrix, A (e.g., process noise covariance matrix), A = LDL'
#we want to do optimization of parameters in the L and D. Capitalizing on
#https://handwiki.org/wiki/Cholesky_decomposition#:~:text=LDL%20decomposition%20A%20closely%20related%20variant%20of%20the,an%20additional%20diagonal%20matrix%20D%20in%20the%20decomposition.
#L = CS^{-1}, where C = the cholesky of A, and S S is a diagonal matrix that contains the main diagonal of C
#D = S^2
#This operation is implemented in dynr.ldl so it returns the values of L and D.
#In classical dynr, here is the flow of actions:
#Before cooking the model:
#(1a) checkSymmetric
#(1b) autoExtendSubRecipe to all regimes
#(1c) replaceDiagZero - replace diagonal values that are = 0 with a small constant
#(2) reverseldl - calls dynr.ldl, which provides values for elements in L and D given A,
#                 and then apply log transformation to elements in D so the optimization
#                 algorithm can optimize these elements on the unconstrained scale.
#(3) writeCcode - makeldlchar does something to ensure uniqueness of parameter names (YET TO GO THROUGH)
#After cooking: take unconstrained parameters returned by optimization algorithm and:
#(1)transformation - calls transldl, which retrieves A = LDL'
#
#In SAEM: we generate starting values for L and D element-wise, following the symmetric indefinite factorization in:
#https://en.wikipedia.org/wiki/Cholesky_decomposition
#, which does not require taking the square root of A
#See symbolicLDLDecomposition and solveStartLDL
#These two are equivalent, as verified in the first chunk of code in Brekfis/CovarianceTesting/CompareLDL.R
#The advantages of these functions are that symbolicLDLDecomposition outputs the symbolic expressions
#of elements of L, D, and LDL', which are now being extended to allow more general covariance function
#modeling involving covariates (needed e.g., for exact discrete-time models)
processCovConstant <- function(values.latent, params.latent, values.observed, params.observed, exclulde_latent = FALSE){
  r <- coProcessValuesParams(values.latent, params.latent)
  values.latent <- r$values
  params.latent <- r$params
  
  # Handle observed covariance
  r <- coProcessValuesParams(values.observed, params.observed)
  values.observed <- r$values
  params.observed <- r$params
  
  values.latent <- lapply(values.latent, preProcessValues)
  params.latent <- lapply(params.latent, preProcessParams)
  values.observed <- lapply(values.observed, preProcessValues)
  params.observed <- lapply(params.observed, preProcessParams)
  
  lapply(values.latent, checkSymmetric, name="values.latent")
  lapply(params.latent, checkSymmetric, name="params.latent")
  lapply(values.observed, checkSymmetric, name="values.observed")
  lapply(params.observed, checkSymmetric, name="params.observed")
  
  # Check that all 'values' imply 0, 1, or the same number of regimes.
  # Note that the 'values' and 'params' have already been checked to imply this.
  nregs <- sapply(list(values.latent=values.latent, values.observed=values.observed), length)
  msg <- paste0("Different numbers of regimes implied:\n'latent' has ", nregs[1], " and 'observed' has ",  nregs[2], " regimes.\nCardi B don't like it like that!")
  mr <- max(nregs)
  if(!all(nregs %in% c(1, mr))){
    stop(msg)
  } else if(!all(nregs %in% mr)) {
    if(nregs[1] == 1){
      ae <- autoExtendSubRecipe(values.latent, params.latent, 'values.latent', 'latent covariances', mr)
      values.latent <- ae[[1]]
      params.latent <- ae[[2]]
    }
    if(nregs[2] == 1){
      ae <- autoExtendSubRecipe(values.observed, params.observed, 'values.observed', 'observed covariances', mr)
      values.observed <- ae[[1]]
      params.observed <- ae[[2]]
    }
  }
  
  values.latent.inv.ldl <- lapply(values.latent, replaceDiagZero)
  values.latent.inv.ldl <- lapply(values.latent.inv.ldl, reverseldl)
  values.observed.inv.ldl <- lapply(values.observed, replaceDiagZero)
  values.observed.inv.ldl <- lapply(values.observed.inv.ldl, reverseldl)
  
  
  if(exclulde_latent == FALSE){
    sv <- c(extractValues(values.latent.inv.ldl, params.latent, symmetric=TRUE), extractValues(values.observed.inv.ldl, params.observed, symmetric=TRUE))
	pn <- c(extractParams(params.latent), extractParams(params.observed))
  }
  else{
    sv <- c(extractValues(values.observed.inv.ldl, params.observed, symmetric=TRUE))
	pn <- c(extractParams(params.observed))
  }
  sv <- extractValues(sv, pn)
  pn <- extractParams(pn)
  x <- list(startval=sv, paramnames=pn, values.latent=values.latent, values.observed=values.observed, params.latent=params.latent, params.observed=params.observed, values.latent.inv.ldl=values.latent.inv.ldl, values.observed.inv.ldl=values.observed.inv.ldl)
  return(x)
}

# Who wrote this?
# It's filled with undefined global variables.
# You can't define a variable inside a conditional, not define it if the condition is not met, and then use it regardless throughout the function.

processCovFormula <- function(dots, values.latent, params.latent){
  latent.formula <- NULL
  state.regimes <- NULL
  latent.startval <- NULL
  if('latent.formula' %in% names(dots))
    latent.formula <- dots$latent.formula
  if('covariates' %in% names(dots))
    covariates <- dots$covariates
  if('state.names' %in% names(dots)){
    state.names <- dots$state.names
  } else {
    state.names <- c()
  }
  
  # e.g. for the one-regime case, if we get a list of formula, make a list of lists of formula
  if(is.list(latent.formula) && plyr::is.formula(latent.formula[[1]])){
    latent.formula <- list(latent.formula)
  }
  latent.names <- lapply(latent.formula, function(fml){sapply(fml, function(x){as.character(x[[2]])})})
  latent.regimes <- paste(paste0('Regime ', 1:length(latent.names), ': ', lapply(latent.names, paste, collapse=', ')), collapse='\n')
  if(any(sapply(lapply(latent.names, duplicated), any))){
    stop(paste0("Found duplicated var-cov names:\n", latent.regimes))
  }
  if(length(unique(sapply(latent.names, length))) != 1){
    stop(paste0("Found different number of var-cov names for different regimes:\n", var.regimes))
  }
  if(!all(sapply(latent.names, function(x, y){all(x==y)}, y=latent.names[[1]]))){
    stop(paste0("Found different var-cov names or different ordering of var-cov names across regimes:\n", var.regimes))
  }
  x <- list(latent.formula=latent.formula, latent.startval=dots$latent.startval)

  #browser() 	  
  sv <- c(extractValues(x$latent.startval, names(x$latent.startval)))
  pn <- c(extractParams(names(x$latent.startval)))
  sv <- extractValues(sv, pn)
  pn <- extractParams(pn)
  
  # Check for var-cov names with same name as free parameter
  if(any(sapply(lapply(latent.names, "%in%", x$paramnames), any))){
    stop(paste0("See no evil, but var-cov parameters had the same names as other free parameters.",
                "\nParameters that are both var-cov and other free parameter names: ",
                paste(unique(c(sapply(latent.names, function(nam){x$paramnames[x$paramnames %in% nam]}))), collapse=', ')))
  }
  
  fml=lapply(latent.formula[[1]],as.character)
  lhs=lapply(fml,function(x){x[[2]]})
  rhs=lapply(fml,function(x){x[[3]]})

  #replace covariates with cov_i
  # if(length(covariates)>0){
    # ind_order <- order(nchar(covariates), covariates, decreasing= TRUE)
    # var_order <- covariates[ind_order]
    # for (i in 1:length(covariates)){
      # pattern <- var_order[i]
	  # ind <- ind_order[i]
      # for (e in 1:length(unlist(var.formula))){
         # rhs[[e]]<- gsub(pattern, paste0('cov_',ind), rhs[[e]], fixed = TRUE)
	  # }
    # }
  # }
  
  #replace state.names with state_i
  # if(length(state.names)>0){
    # ind_order <- order(nchar(state.names), state.names, decreasing= TRUE)
    # var_order <- state.names[ind_order]
    # for (i in 1:length(state.names)){
      # pattern <- var_order[i]
	  # ind <- ind_order[i]
      # for (e in 1:length(unlist(var.formula))){
         # rhs[[e]]<- gsub(pattern, paste0('state_',ind), rhs[[e]], fixed = TRUE)
	  # }
    # }
  # }

  # for (e in 1:length(unlist(var.formula))){
    # x$var.formula[[1]][[e]] <- as.formula(paste(lhs[[e]], '~', rhs[[e]]))
  # }
  #TODO: functions I need to call later:
  #processFormula -> parseFormula -> # reverseldl (which calls dynr.ldl) -> trans2ArmadilloFunction (as opposed to trans2CFunction):
  
  #browser()
  formula.latent <- params.latent
  for (i in 1:nrow(params.latent)){
	for (j in 1:ncol(params.latent)){
	  for (e in 1:length(unlist(latent.formula))){
	    formula.latent[i,j] <- gsub(lhs[[e]], rhs[[e]], formula.latent[i,j], fixed = TRUE)
	  }
	}
  }
  
  #browser()
  out <- symbolicLDLDecomposition2(formula.latent, values.latent)
  #known.vars <- solveStartLDLSymbolically(out$ldl, params.latent)  
  ldl.transformed <- out$ldl
  for (i in 1:nrow(params.latent)){
	for (j in 1:ncol(params.latent)){
	  ldl.transformed[i, j][[1]] <- as.formula(paste0( params.latent[i, j], '~', out$ldl[i, j]))
	  ldl.transformed[i, j][[1]] <- as.formula(paste0( params.latent[i, j], '~', formula.latent[i,j]))
	}
  }
  known.vars = list(x~0)
  
  
  #sampleCovformula=
  #  list(Sigma11 ~ par1*delta_t^3,
  #       Sigma12 ~ par2*deltat, 
  #       Sigma21 ~ par2*deltat,
  #       Sigma22 ~ par3*delta_t^3)
  
  # SymbolicLDLDecomposition (dynrModel.R)-> L, D -> unique parameters in L & D -> InfDS.par
  # a.params = sampleCovformula, which is a list of formulas
  # a.values =  c(par1 = 3, par2 = .1, par3 = 1)
  # covariate.names = c(delta_t)
  # Call SymbolicLDLDecomposition and have it return L and D as matrices of expressions
  # solveStartLDL = evaluate the whole formula and find starting values for the parameters in L and D
  # Write out Armadillo code
  # D = exp(par11*delta_t)
  
  x <- list(startval=sv, paramnames=pn, ldl.transformed=ldl.transformed, latent.formula=latent.formula, latent.startval=dots$latent.startval, known.vars=known.vars)
  return(x)
}  


replaceDiagZero <- function(x){
	diag(x)[diag(x) == 0] <- 1e-6
	return(x)
}

checkSymmetric <- function(m, name="matrix"){
	if(any(m != t(m))){
		msg <- paste0("Covariance matrix called '", name, "' is not symmetric.\n",
			"Covariance matrices should be equal to their transposes.")
		stop(msg)
	}
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
##' @param deviation logical. Whether to use the deviation form or not.  See Details.
##' @param refRow numeric. Which row is treated at the reference.  See Details.
##' 
##' @details
##' Note that each row of the transition probability matrix must sum to one. To accomplish this
##' fix at least one transition log odds parameter in each row of "values" (including its intercept 
##' and the regression slopes of all covariates) to 0.
##' 
##' When \code{deviation=FALSE}, the non-deviation form of the multinomial logistic regression is used. This form has a separate intercept term for each entry of the transition probability matrix (TPM). When \code{deviation=TRUE}, the deviation form of the multinomial logistic regression is used. This form has an intercept term that is common to each column of the TPM. The rows are then distinguished by their own individual deviations from the common intercept. The deviation form requires the same reference column constraint as the non-deviation form; however, the deviation form also requires one row to be indicated as the reference row (described below). By default the reference row is taken to be the same as the reference column.
##' 
##' The \code{refRow} argument determines which row is used as the intercept row. It is only
##' used in the deviation form (i.e. \code{deviation=TRUE}). In the deviation form, one row of \code{values} and \code{params} contains the intercepts, other rows contain deviations from these intercepts. The \code{refRow} argument says which row contains the intercept terms. The default behavior for \code{refRow} is to be the same as the reference column.  The reference column is automatically detected. If we have problems detecting which is the reference column, then we provide error messages that are as helpful as we can make them.
##' 
##' @return Object of class 'dynrRegimes'
##' 
##' @seealso 
##' Methods that can be used include: \code{\link{print}}, \code{\link{printex}}, \code{\link{show}} 
##'
##' @examples
##' #Two-regime example with a covariate, x; log odds (LO) parameters represented in default form,
##' #2nd regime set to be the reference regime (i.e., have LO parameters all set to 0).
##' #The values and params matrices are of size 2 (numRegimes=2) x 4 (numRegimes*(numCovariates+1)).
##' #    The LO of staying within the 1st regime (corresponding to the (1,1) entry in the
##' #              2 x 2 transition probability matrix for the 2 regimes) = a_11 + d_11*x
##' #    The log odds of switching from the 1st to the 2nd regime (the (1,2) entry in the
##' #              transition probability matrix) = 0
##' #    The log odds of moving from regime 2 to regime 1 (the (2,1) entry) = a_21 + d_21*x
##' #    The log odds of staying within the 2nd regime (the (2,2) entry) = 0
##' b <- prep.regimes(
##' values=matrix(c(8,-1,rep(0,2),
##'                -4,.1,rep(0,2)),
##'              nrow=2, ncol=4, byrow=TRUE), 
##' params=matrix(c("a_11","d_11x",rep("fixed",2),
##'                "a_21","d_21x",rep("fixed",2)), 
##'              nrow=2, ncol=4, byrow=TRUE), covariates=c("x"))
##'  
##' # Same example as above, but expressed in deviation form by specifying 'deviation = TRUE'
##' #    The LO of staying within the 1st regime (corresponding to the (1,1) entry in the
##' #              2 x 2 transition probability matrix for the 2 regimes) = a_21 + a_11 + d_11*x
##' #    The log odds of switching from the 1st to the 2nd regime (the (1,2) entry in the
##' #              transition probability matrix) = 0
##' #    The log odds of moving from regime 2 to regime 1 (the (2,1) entry) = a_21 + d_21*x
##' #    The log odds of staying within the 2nd regime (the (2,2) entry) = 0            
##' b <- prep.regimes(
##' values=matrix(c(8,-1,rep(0,2),
##'                -4,.1,rep(0,2)),
##'              nrow=2, ncol=4, byrow=TRUE), 
##' params=matrix(c("a_11","d_11x",rep("fixed",2),
##'                "a_21","d_21x",rep("fixed",2)), 
##'              nrow=2, ncol=4, byrow=TRUE), covariates=c("x"), deviation = TRUE)
##'              
##' #An example of regime-switching with no covariates. The diagonal entries are fixed
##' #at zero for identification purposes
##' b <- prep.regimes(values=matrix(0, 3, 3), 
##' params=matrix(c('fixed', 'p12', 'p13', 
##'                 'p21', 'fixed', 'p23', 
##'                 'p31', 'p32', 'fixed'), 3, 3, byrow=TRUE))
##' 
##' #An example of regime-switching with no covariates. The parameters for the second regime are 
##' #  fixed at zero for identification purposes, making the second regime the reference regime.
##' b <- prep.regimes(values=matrix(0, 3, 3), 
##' params=matrix(c('p11', 'fixed', 'p13',
##'                 'p21', 'fixed', 'p23', 
##'                 'p31', 'fixed', 'p33'), 3, 3, byrow=TRUE))
##' 
##' #2 regimes with three covariates
##' b <- prep.regimes(values=matrix(c(0), 2, 8), 
##' params=matrix(c(paste0('p', 8:15), rep(0, 8)), 2, 8), 
##' covariates=c('x1', 'x2', 'x3'))
##' 
prep.regimes <- function(values, params, covariates, deviation=FALSE, refRow){
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
	if(missing(refRow)){ # or other default values???
		refRow <- numeric(0)
	}
	
	if(!all(dim(values))==all(dim(params))) {
	  stop('The dimensions of values and params matrices must match.')
	}else if(!ncol(values ) == nrow(values)*(length(covariates)+1)){
	  stop(paste0("The matrix values should have ", nrow(values)*(length(covariates)+1), " columns."))
	}
	
	sv <- extractValues(values, params)
	pn <- extractParams(params)
	
	if(length(sv) > nrow(values)*(nrow(values)-1)*(length(covariates)+1)){
		stop("Regime transition probability parameters not uniquely identified.\nFix all parameters in at least one cell of each row of the \ntransition probability matrix to be zero.")
	}
	# Do lots or processing and checking for the reference row and column
	if(deviation == TRUE){
		# if needed set refRow
		if(length(refRow) == 0){
			# check valid reference column (i.e. there is a single column as reference, no diagonal identification)
			# detect reference column
			colIsZero <- apply(params, 2, function(x){all(x %in% c('fixed'))})
			blockSize <- length(covariates)+1
			colIsZeroM <- matrix(colIsZero, nrow=blockSize, ncol=nrow(values))
			refCol <- apply(colIsZeroM, 2, function(x){all(x==TRUE)})
			if(sum(refCol) < 1){ #no single reference column found
				stop(paste("No single reference column found. Identification along the diagonal is not allowed without you specifiy the reference row. Set e.g. refRow=1."))
			}
			if(sum(refCol) > 1){#found more than one reference column
				stop(paste0("Found multiple possible reference columns: ", paste(which(refCol==1), collapse=', '), ".  Reconsider your model specification or just set e.g. refRow=1."))
			}
			# set reference row to reference column
			refRow <- which(refCol)
		}
		if(length(refRow) > 1){
			# error reference row must have length 0 or 1
			stop(paste("refRow must being a single number. That is, have length 0 or 1, we found length", length(refRow)))
		}
		if(refRow > nrow(values)){
			# error reference row must be between 1 and nrow(values)
			stop(paste("refRow must be between 1 and", nrow(values), "but we found", refRow))
		}
	} else {
		if(length(refRow) > 0){
			#warning refRow is ignored when for the non-deviation case (deviation=FALSE)
			warning("refRow argument is ignored in the non-deviation case (i.e. when deviation=FALSE)")
		}
	}
	# Create the list for the object to return
	x <- list(startval=sv, paramnames=pn, values=values, params=params, covariates=covariates, deviation=deviation, refRow=refRow)
	return(new("dynrRegimes", x))
}


#autojacob: obtain differentiation of formula by diff.variables); if diff.variables are not given, do differentiation by left-hand-side of variables by formula
autojacob <- function(formula, n, diff.variables){
  tuple=lapply(formula,as.list)
  lhs=sapply(tuple,function(x){deparse(x[[2]])})
  rhs=sapply(tuple,function(x){x[[3]]})
  
  if(missing(diff.variables)){
    # default setting, differentiate formula by the varibles in lhs
    n.diff = n
    diff.variables = lhs
  }
  else{
    # differentiate by the varialbes by diff.variables
    n.diff = length(diff.variables)
  }
  
  rhsj=vector("list", n*n.diff)
  jacob=vector("list", n*n.diff)
  for (i in 1:n){
    for (j in 1:n.diff){
      #jacob[[(i-1)*n.par+j]] = df[i]/dpar[j], column-major order   
      rhsj[[(i-1)*n.diff+j]]=paste0(deparse(D(rhs[[i]],diff.variables[[j]]),width.cutoff = 500L),collapse="")
      jacob[[(i-1)*n.diff+j]]=as.formula(paste0(lhs[[i]],"~",diff.variables[[j]],"~",rhsj[[(i-1)*n.diff+j]],collapse=""))
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
##' formulas for estimation with parameter names as its name. If there are no free parameters in 
##' the dynamic functions, leave startval as the default \code{numeric(0)}.
##' @param isContinuousTime if True, the left hand side of the formulas represent 
##' the first-order derivatives of the specified variables; if False, the left hand 
##' side of the formulas represent the current state of the specified variable while 
##' the same variable on the righ hand side is its previous state.  
##' @param jacobian (optional) a list of formulas specifying the analytic jacobian matrices 
##' containing the analytic differentiation function of the dynamic functions with respect to
##' the latent variables. If this is not provided, dynr will invoke an automatic differentiation
##' procedure to compute the jacobian functions.
##' @param ... further named arguments. Some of these arguments may include:
##' 
##' \code{theta.formula} specifies a list consisting of formula(s) of the form 
##' \code{list (par ~ 1 * b_0  + covariate_1 * b_1 + ... + covariate_p * b_p 
##'  + 1 * rand_par)}, where \code{par} is a parameter is a unit- (e.g., person-) 
##'  specific that appears in a dynamic formula and is assumed to follow
##'  a linear mixed effects structure. Here, \code{b_p} are fixed effects 
##'  parameters; \code{covariate_1}, ..., \code{covariate_p} are known covariates as predeclared in
##'  \code{dynr.data}, and \code{rand_par} is a random effect component representing unit i's random deviation
##'  in \code{par} value from that predicted by \code{b_0 + covariate_1*b_1 + ... + covariate_p*b_p}. 
##'
##' \code{random.names} specifies names of random effect components in the \code{theta.formula}
##'
##' \code{random.params.inicov} specifies names of elements in the covariance matrix of the random effect components
##'
##' \code{random.values.inicov} specifies starting values of elements in the covariance matrix of the random effect components
##' 
##' \code{random.lb}  and \code{random.ub} specify the lower and upper bound of the random effect estimates
##' 
##' 
##' @details
##' This function defines the dynamic functions of the model either in discrete time or in continuous time.
##' The function can be either linear or nonlinear, with free or fixed parameters, numerical constants, 
##' covariates, and other mathematical functions that define the dynamics of the latent variables.
##' Every latent variable in the model needs to be defined by a differential (for continuous time model), or
##' difference (for discrete time model) equation.  The names of the latent variables should match 
##' the specification in \code{prep.measurement()}.
##' For nonlinear models, the estimation algorithm generally needs a Jacobian matrix that contains
##' elements of first differentiations of the dynamic functions with respect to the latent variables
##' in the model. For most nonlinear models, such differentiations can be handled automatically by
##' dynr. However, in some cases, such as when the absolute function (\code{abs}) is used, the automatic
##' differentiation would fail and the user may need to provide his/her own Jacobian functions.
##' When \code{theta.formula} and other accompanying elements in "\code{...}" are provided, the program
##' automatically inserts the random effect components specified in random.names as additional
##' latent (state) variables in the model, and estimate (cook) this expanded model. Do check
##' that the expanded model satisfies conditions such as observability for the estimation to work.
##' 
##' @return Object of class 'dynrDynamicsFormula'
##'
##' @examples
##' # In this example, we present how to define the dynamics of a bivariate dual change score model
##' # (McArdle, 2009). This is a linear model and the user does not need to worry about 
##' # providing any jacobian function (the default). 
##'  
##' # We start by creating a list of formula that describes the model. In this model, we have four 
##' # latent variables, which are "readLevel", "readSlope", "mathLevel", and "math Slope".  The right-
##' # hand side of each formula gives a function that defines the dynamics.   
##'  
##'  formula <- list(
##'           list(readLevel~ (1+beta.read)*readLevel + readSlope + gamma.read*mathLevel,
##'           readSlope~ readSlope,
##'           mathLevel~ (1+beta.math)*mathLevel + mathSlope + gamma.math*readLevel, 
##'           mathSlope~ mathSlope
##'           ))
##'
##' # Then we use prep.formulaDynamics() to define the formula, starting value of the parameters in
##' # the model, and state the model is in discrete time by setting isContinuousTime=FALSE.
##'  
##' dynm  <- prep.formulaDynamics(formula=formula,
##'                              startval=c(beta.read = -.5, beta.math = -.5, 
##'                                         gamma.read = .3, gamma.math = .03
##'                              ), isContinuousTime=FALSE)
##' 
##' 
##' # For a full demo example of regime switching nonlinear discrete time model, you
##' # may refer to a tutorial on 
##' # \url{https://quantdev.ssri.psu.edu/tutorials/dynr-rsnonlineardiscreteexample}
##' 
##' #Not run: 
##' #For a full demo example that uses user-supplied analytic jacobian functions see:
##' #demo(RSNonlinearDiscrete, package="dynr")
##' formula <- list(
##'     list(
##'       x1 ~ a1*x1,
##'       x2 ~ a2*x2),
##'     list(
##'       x1 ~ a1*x1 + c12*(exp(abs(x2)))/(1+exp(abs(x2)))*x2,
##'       x2 ~ a2*x2 + c21*(exp(abs(x1)))/(1+exp(abs(x1)))*x1)
##'   )
##' jacob <- list(
##'   list(x1~x1~a1,
##'       x2~x2~a2),
##'   list(x1~x1~a1,
##'       x1~x2~c12*(exp(abs(x2))/(exp(abs(x2))+1)+x2*sign(x2)*exp(abs(x2))/(1+exp(abs(x2))^2)),
##'       x2~x2~a2,
##'       x2~x1~c21*(exp(abs(x1))/(exp(abs(x1))+1)+x1*sign(x1)*exp(abs(x1))/(1+exp(abs(x1))^2))))
##' dynm <- prep.formulaDynamics(formula=formula, startval=c( a1=.3, a2=.4, c12=-.5, c21=-.5),
##'                              isContinuousTime=FALSE, jacobian=jacob)
##' 
##' #For a full demo example that uses automatic jacobian functions (the default) see:
##' #demo(RSNonlinearODE , package="dynr")
##' formula=list(prey ~ a*prey - b*prey*predator, predator ~ -c*predator + d*prey*predator)
##' dynm <- prep.formulaDynamics(formula=formula,
##'                           startval=c(a = 2.1, c = 0.8, b = 1.9, d = 1.1),
##'                           isContinuousTime=TRUE)
##'
##' #For a full demo example that includes unit-specific random effects in theta.formula see:
##' #demo(OscWithRand, package="dynr")
##' formula <- list(x ~ dx,
##'                dx ~ eta_i * x + zeta*dx)
##' theta.formula  = list (eta_i ~ 1 * eta0  + u1 * eta1 + u2 * eta2 + 1 * b_eta)
##' dynm <- prep.formulaDynamics(formula=formula,
##'                            startval=c(eta0=-1, eta1=.1, eta2=-.1,zeta=-.02),
##'                            isContinuousTime=TRUE,
##'                            theta.formula=theta.formula,
##'                            random.names=c('b_eta'),
##'                            random.params.inicov=matrix(c('sigma2_b_eta'), ncol=1,byrow=TRUE),
##'                            random.values.inicov=matrix(c(0.1), ncol=1,byrow=TRUE))
prep.formulaDynamics <- function(formula, startval = numeric(0), isContinuousTime=FALSE, jacobian, saem=FALSE, ...){
	#browser()
	dots <- list(...)
	
	# if the argument formula is not a list 
	if(!is.list(formula)){
		msg <- paste0(ifelse(plyr::is.formula(formula), "'formula' argument is a formula but ", ""),
			"'formula'",
			ifelse(plyr::is.formula(formula), " ", " argument "),
			"should be a list of formulas.\nCan't nobody tell me nothin'")
		stop(msg)
	}
	
	
	if(length(startval) == 0){
		warning("You provided no start values: length(startval)==0. If you have no free parameters, keep calm and carry on.")
	}
	
	# parsing the parameters in ...
	if(length(dots) > 0){
		if(!all(names(dots) %in% c('theta.formula', 'random.names',  'random.params.inicov', 'random.values.inicov', 'random.ub', 'random.lb'))){
			stop("You passed some invalid names to the ... argument. Check with US Customs or the ?prep.formulaDynamics help page.")
		}
		if('theta.formula' %in% names(dots)){
			theta.formula <- dots$theta.formula
			# retrieve theta.names from theta.formula instead of asking users to give inputs
			fml=lapply(theta.formula, as.character)
			theta.names=unlist(lapply(fml,function(x){x[[2]]}))
		}
	}
	if(length(startval) > 0){
		beta.names = names(startval)
	}

	if(length(startval) > 0 & is.null(names(startval))){
		stop('startval must be a named vector.')
	}



  
	# e.g. for the one-regime case, if we get a list of formula, make a list of lists of formula
	if(is.list(formula) && plyr::is.formula(formula[[1]])){
		formula <- list(formula)
	}
	state.names <- lapply(formula, function(fml){sapply(fml, function(x){as.character(x[[2]])})})
	state.regimes <- paste(paste0('Regime ', 1:length(state.names), ': ', lapply(state.names, paste, collapse=', ')), collapse='\n')
	if(any(sapply(lapply(state.names, duplicated), any))){
		stop(paste0("Found duplicated latent state names:\n", state.regimes))
	}
	if(length(unique(sapply(state.names, length))) != 1){
		stop(paste0("Found different number of latent states for different regimes:\n", state.regimes))
	}
	if(!all(sapply(state.names, function(x, y){all(x==y)}, y=state.names[[1]]))){
		stop(paste0("Found different latent states or different ordering of latent states across regimes:\n", state.regimes))
	}
	x <- list(formula=formula, startval=startval, paramnames=c(preProcessParams(names(startval))), isContinuousTime=isContinuousTime)
	
	#browser()
	# if (missing(jacobian)){
	  # if('theta.formula' %in% names(dots)){
	    # # This line has bug.
        # #jacobian <- autojacobTry(lapply(formula, function(x){parseFormulaTheta(x, theta.formula)}))
		# jacobian <- autojacobTry(formula)
	    # #jacobianOriginal <- autojacobTry(formula)
	  # }
	  # else{   
      # jacobian <- autojacobTry(formula)
      # autojcb=try(lapply(formula,autojacob,length(formula[[1]])))
		# if (class(autojcb) == "try-error") {
			# stop("Automatic differentiation is not supported by part of the dynamic functions.\n 
					 # Please provide the analytic jacobian functions.")
		# }else{
			# jacobian=lapply(autojcb,"[[","jacob") 
		# }
	  # }
	# }
	#jacobian <- autojacob(formula)
	if (missing(jacobian)){
		autojcb=try(lapply(formula,autojacob,length(formula[[1]])))
		if ("try-error" %in% class(autojcb)) {
			stop("Automatic differentiation is not supported by part of the dynamic functions.\n 
					 Please provide the analytic jacobian functions.")
		}else{
			jacobian=lapply(autojcb,"[[","jacob") 
		}
	}
  
  


	# Check that all 'values' imply 0, 1, or the same number of regimes.
	# Note that the 'values' and 'params' have already been checked to imply this.
	nregs <- sapply(list(formula=formula, jacobian=jacobian), length)
	if(nregs[1] != nregs[2]){
		stop(paste0("Don't bring that trash up in my house!\nDifferent numbers of regimes implied:\n'formula' has ", nregs[1], " but 'jacobian' has ", nregs[2], " regimes."))
	}
	x$jacobian <- jacobian
	x$formulaOriginal <- x$formula
	x$jacobianOriginal <- jacobian
	x$paramnames <- names(x$startval)
	
	# Check for latent states with same name as free parameter
	if(any(sapply(lapply(state.names, "%in%", x$paramnames), any))){
		stop(paste0("See no evil, but latent states had the same names as free parameters.",
			"\nParameters that are both latent states and free parameters: ",
			paste(unique(c(sapply(state.names, function(nam){x$paramnames[x$paramnames %in% nam]}))), collapse=', ')))
	}
	if('theta.formula' %in% names(dots)){
		x$theta.formula <- theta.formula
		x$theta.names <- theta.names
	}

	
	# parsing the parameters in ...
	if(length(dots) > 0){
		if('random.names' %in% names(dots))
			x$random.names <- dots$random.names
		if('random.params.inicov' %in% names(dots))
			x$random.params.inicov <- dots$random.params.inicov
		if('random.values.inicov' %in% names(dots))
			x$random.values.inicov <- dots$random.values.inicov
		
		if('random.ub' %in% names(dots))
			x$random.ub <- dots$random.ub
		#else if ('random.values.inicov' %in% names(dots))
		#	x$random.ub <- diag(dots$random.values.inicov)[1] * 10
		else
			x$random.ub <- 10
		
		 
		if('random.lb' %in% names(dots))
			x$random.lb <- dots$random.lb
		#else if ('random.values.inicov' %in% names(dots))
		#	x$random.lb <- diag(dots$random.values.inicov)[1] * -10
		else
			x$random.lb <- -10#dots$random.lb


	}
	

	startval.names <- names(startval)
	x$saem <- saem
	
	if(saem == FALSE){
		return(new("dynrDynamicsFormula", x))
	}

	
	# The following is only for saem	

		

	
	if(is.list(state.names)){
		state.names = unlist(state.names)
	}
	
	
	#dfdtheta <- autojacobTry(formula_onlystate, diff.variables=theta.names)
	dfdtheta <- autojacob(formula[[1]], n = length(formula[[1]]), diff.variables=theta.names)
	x$dfdtheta <- dfdtheta$jacob
	x$theta.names<-theta.names
	
	dfdx <- autojacob(formula[[1]], n = length(formula[[1]]), diff.variables=state.names)
	dfdx2 <- autojacob(dfdx$jacob, n = length(dfdx$jacob),  diff.variables=state.names)
	x$dfdx2 <- dfdx2$jacob
	
	dfdxdtheta <- autojacob(dfdx$jacob, n = length(dfdx$jacob), diff.variables=theta.names)
	x$dfdxdtheta <- dfdxdtheta$jacob
	
	dfdthetadx <- autojacob(dfdtheta$jacob, n = length(dfdtheta$jacob), diff.variables=state.names)
	x$dfdthetadx <- dfdthetadx$jacob
	
	dfdtheta2 <- autojacob(dfdtheta$jacob, n = length(dfdtheta$jacob), diff.variables=theta.names)
	x$dfdtheta2 <- dfdtheta2$jacob
	x$beta.names <- startval.names
	
	
	#transfer differentiation into matrix
	#d(i,j): the differention of i-th matrix to j-th variable
	#browser()
	x$jacobian <- list(matrix(unlist(x$jacobian), ncol =length(state.names), byrow= FALSE))
	x$jacobianOriginal <- list(matrix(unlist(x$jacobianOriginal), ncol = length(state.names), byrow= FALSE))
	x$dfdtheta <- list(matrix(unlist(x$dfdtheta), nrow =length(theta.names), byrow = FALSE))
	x$dfdx2 <- list(matrix(unlist(x$dfdx2), ncol = length(state.names),byrow= TRUE))
	x$dfdxdtheta <- list(matrix(unlist(x$dfdxdtheta), ncol =length(theta.names),byrow= TRUE))
	x$dfdthetadx <- list(matrix(unlist(x$dfdthetadx), ncol = length(state.names),byrow= TRUE))
	x$dfdtheta2 <- list(matrix(unlist(x$dfdtheta2), ncol =length(theta.names),byrow= TRUE))
	
	
	
	#x$theta.formula <- theta.formula
	x$theta.names <- theta.names
	x$state.names <- state.names
	
	#x$covariate.names <- dots$covariate.names
	#x$covariate.formula <- dots$covariate.formula

	
	return(new("dynrDynamicsFormula", x))
}


autojacobTry <- function(formula, formula2, ...){
    if(missing(formula2)) formula2 <- formula
    autojcb <- try( lapply(formula2, autojacob, length(formula[[1]]), ...) )
    if (class(autojcb) == "try-error") {
        # TODO Add more specific feedback about which formula(s) failed automatic differentiation
        stop("Automatic differentiation is not supported for these specific dynamic functions.\n 
            You may want to provide the analytic jacobian functions instead.")
    }else{
        return(lapply(autojcb,"[[", "jacob"))
    }
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
##' A recipe function for specifying the deterministic portion of a set of linear dynamic functions as:
##'
##' Discrete-time model: eta(t+1) = int + dyn*eta(t) + exo*x(t), 
##' where eta(t) is a vector of latent variables, x(t) is a vector of covariates,
##' int, dyn, and exo are vectors and matrices specified via the arguments *.int, *.dyn, and *.exo. 
##'
##' Continuous-time model: d/dt eta(t) = int + dyn*eta(t) + exo*x(t), 
##' where eta(t) is a vector of latent variables, x(t) is a vector of covariates,
##' int, dyn, and exo are vectors and matrices specified via the arguments *.int, *.dyn, and *.exo.
##' 
##' The left-hand side of the dynamic model consists of a vector of latent variables for the next time point in the discrete-time case,
##' and the vector of derivatives for the latent variables at the current time point in the continuous-time case.
##'
##' For models with regime-switching dynamic functions, the user will need to provide a list of the *.int, *.dyn, and *.exo arguments. 
##' (when they are specified to take on values other than the default of zero vectors and matrices), or if a single set of vectors/matrices are provided, the same 
##' vectors/matrices are assumed to hold across regimes.
##' 
##' \code{prep.matrixDynamics} serves as an alternative to \code{\link{prep.formulaDynamics}}.
##' 
##' @return Object of class 'dynrDynamicsMatrix'
##' 
##' @seealso 
##' Methods that can be used include: \code{\link{print}}, \code{\link{show}} 
##' 
##' @examples 
##' #Single-regime, continuous-time model. For further details run: 
##' #demo(RSNonlinearDiscrete, package="dynr"))
##' dynamics <- prep.matrixDynamics(
##'     values.dyn=matrix(c(0, -0.1, 1, -0.2), 2, 2),
##'     params.dyn=matrix(c('fixed', 'spring', 'fixed', 'friction'), 2, 2),
##'     isContinuousTime=TRUE)
##' 
##' #Two-regime, continuous-time model. For further details run: 
##' #demo(RSNonlinearDiscrete, package="dynr"))
##' dynamics <- prep.matrixDynamics(
##'     values.dyn=list(matrix(c(0, -0.1, 1, -0.2), 2, 2),
##'                     matrix(c(0, -0.1, 1, 0), 2, 2)),
##'     params.dyn=list(matrix(c('fixed', 'spring', 'fixed', 'friction'), 2, 2),
##'                     matrix(c('fixed', 'spring', 'fixed', 'fixed'), 2, 2)),
##'     isContinuousTime=TRUE) 
prep.matrixDynamics <- function(params.dyn=NULL, values.dyn, params.exo=NULL, values.exo=NULL, params.int=NULL, values.int=NULL, 
                                covariates, isContinuousTime, ...){
	# Handle numerous cases of missing or non-list arguments
	# General idea
	# If they give us a non-list argument, make it a one-element list
	# If they don't give us params, assume it's all fixed params
	# If they don't give us values, assume they don't want that part of the model
	
	# Handle dyn
  dots <- list(...)
  
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
	
	# Check that the number of covariates implied by the 'covariates' arg is the same as that
	#  implied by the number of columns in the 'values.exo' arg.
	matCovariates <- lapply(lapply(values.exo, dim), "[[", 2)
	if(length(matCovariates) == 0){matCovariates <- 0}
	argCovariates <- length(covariates)
	if(!all(matCovariates == argCovariates)){
		msg <- paste0("Mind your teaspoons and tablespoons.  The 'exo.values' argument says there are\n (", paste(matCovariates[[1]], collapse=", "), ") covariates, but the 'covariates' arg says there are (", argCovariates, ").")
		stop(msg)
	}
	
	# Check that all 'values' imply 0, 1, or the same number of regimes.
	# Note that the 'values' and 'params' have already been checked to imply this.
	nregs <- sapply(list(values.dyn=values.dyn, values.exo=values.exo, values.int=values.int), length)
	msg <- paste0("Different numbers of regimes implied:\n'dyn' has ", nregs[1], ", 'exo' has ",  nregs[2], ", and 'int' has ", nregs[3], " regimes.\nWhat do you want from me? I'm not America's Sweetheart!")
	mr <- max(nregs)
	if(!all(nregs %in% c(0, 1, mr))){
		stop(msg)
	} else if(!all(nregs %in% c(0, mr))) {
		if(nregs[1] == 1){
			ae <- autoExtendSubRecipe(values.dyn, params.dyn, 'values.dyn', 'dynamics', mr)
			values.dyn <- ae[[1]]
			params.dyn <- ae[[2]]
		}
		if(nregs[2] == 1){
			ae <- autoExtendSubRecipe(values.exo, params.exo, 'values.exo', 'covariate regessions', mr)
			values.exo <- ae[[1]]
			params.exo <- ae[[2]]
		}
		if(nregs[3] == 1){
			ae <- autoExtendSubRecipe(values.int, params.int, 'values.int', 'intercepts', mr)
			values.int <- ae[[1]]
			params.int <- ae[[2]]
		}
	}
	
	sv <- c(extractValues(values.dyn, params.dyn), extractValues(values.exo, params.exo), extractValues(values.int, params.int))
	pn <- c(extractParams(params.dyn), extractParams(params.exo), extractParams(params.int))
	sv <- extractValues(sv, pn)
	pn <- extractParams(pn)
	
	x <- list(startval=sv, paramnames=pn, params.dyn=params.dyn, values.dyn=values.dyn, params.exo=params.exo, values.exo=values.exo, params.int=params.int, values.int=values.int, isContinuousTime=isContinuousTime,covariates=covariates)
	# TEMP: Add in dots and additional arguments needed for saem
	if(length(dots) > 0){
	  if('random.names' %in% names(dots))
	    x$random.names <- dots$random.names
	  if('random.params.inicov' %in% names(dots))
	    x$random.params.inicov <- dots$random.params.inicov
	  if('random.values.inicov' %in% names(dots))
	    x$random.values.inicov <- dots$random.values.inicov
	  
	  if('random.ub' %in% names(dots))
	    x$random.ub <- dots$random.ub
	  #else if ('random.values.inicov' %in% names(dots))
	  #	x$random.ub <- diag(dots$random.values.inicov)[1] * 10
	  else
	    x$random.ub <- 10
	  
	  
	  if('random.lb' %in% names(dots))
	    x$random.lb <- dots$random.lb
	  #else if ('random.values.inicov' %in% names(dots))
	  #	x$random.lb <- diag(dots$random.values.inicov)[1] * -10
	  else
	    x$random.lb <- -10#dots$random.lb
	}
	# END OF TEMP
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
##' @param values.inistate a vector or list of vectors of the starting or fixed values of the initial state vector in one or more regimes.  May also be a matrix or list of matrices.
##' @param params.inistate a vector or list of vectors of the parameter names that appear in the initial state vector in one or more regimes. If an element is 0 or "fixed", the corresponding element is fixed at the value specified in the values vector; Otherwise, the corresponding element is to be estimated with the starting value specified in the values vector.  May also be a matrix or list of matrices.
##' @param values.inicov a positive definite matrix or a list of positive definite matrices of the starting or fixed values of the initial error covariance structure(s) in one or more regimes. If only one matrix is specified for a regime-switching dynamic model, the initial error covariance structure stays the same across regimes. To ensure the matrix is positive definite in estimation, we apply LDL transformation to the matrix. Values are hence automatically adjusted for this purpose.
##' @param params.inicov a matrix or list of matrices of the parameter names that appear in the initial error covariance(s) in one or more regimes. If an element is 0 or "fixed", the corresponding element is fixed at the value specified in the values matrix; Otherwise, the corresponding element is to be estimated with the starting value specified in the values matrix. If only one matrix is specified for a regime-switching dynamic model, the process noise structure stays the same across regimes. If a list is specified, any two sets of the parameter names as in two matrices should be either the same or totally different to ensure proper parameter estimation.
##' @param values.regimep a vector/matrix of the starting or fixed values of the initial probabilities of being in each regime. By default, the initial probability of being in the first regime is fixed at 1.
##' @param params.regimep a vector/matrix of the parameter indices of the initial probabilities of 
##' being in each regime. If an element is 0 or "fixed", the corresponding element is fixed at the value 
##' specified in the "values" vector/matrix; Otherwise, the corresponding element is to be estimated 
##' with the starting value specified in the values vector/matrix.
##' @param covariates character vector of the names of the (person-level) covariates
##' @param deviation logical. Whether to use the deviation form or not.  See Details.
##' @param refRow numeric. Which row is treated at the reference.  See Details.
##' 
##' @details
##' The initial condition model includes specifications for the intial state vector, 
##' initial error covariance matrix, initial probabilities of 
##' being in each regime and all associated parameter specifications.
##' The initial probabilities are specified in multinomial logistic regression form.  When there are no covariates, this implies multinomial logistic regression with intercepts only.  In particular, the initial probabilities not not specified on a 0 to 1 probability scale, but rather a negative infinity to positive infinity log odds scale.  Fixing an initial regime probability to zero does not mean zero probability.  It translates to a comparison log odds scale against which other regimes will be judged.
##' 
##' The structure of the initial state vector and the initial probability vector depends on the presence of covariates.  When there are no covariates these should be vectors, or equivalently single-column matrices.  When there are covariates they should have \eqn{c+1} columns for \eqn{c} covariates.  For \code{values.regimep} and \code{params.regimep} the number of rows should be the number of regimes.  For \code{inistate} and \code{inicov} the number of rows should be the number of latent states.  Of course, \code{inicov} is a square and symmetric so its number of rows should be the same as its number of columns.
##' 
##' When \code{deviation=FALSE}, the non-deviation form of the multinomial logistic regression is used. This form has a separate intercept term for each entry of the initial probability vector. When \code{deviation=TRUE}, the deviation form of the multinomial logistic regression is used. This form has an intercept term that is common to all rows of the initial probability vector. The rows are then distinguished by their own individual deviations from the common intercept. The deviation form requires the same reference row constraint as the non-deviation form (described below). By default the reference row is taken to be the row with all zero covariate effects.  Of course, if there are no covariates and the deviation form is desired, then the user must provide the reference row.
##' 
##' The \code{refRow} argument determines which row is used as the intercept row. It is only
##' used in the deviation form (i.e. \code{deviation=TRUE}). In the deviation form, one row of \code{values.regimep} and \code{params.regimep} contains the intercepts, other rows contain deviations from these intercepts. The \code{refRow} argument says which row contains the intercept terms. The default behavior for \code{refRow} is to detect the reference row automatically based on which parameters are \code{fixed}.  If we have problems detecting which is the reference row, then we provide error messages that are as helpful as we can make them.
##' 
##' @seealso 
##' Methods that can be used include: \code{\link{print}}, \code{\link{printex}}, \code{\link{show}}
##' 
##' @return Object of class 'dynrInitial'
##' 
##' @examples
##' #### No-covariates
##' # Single regime, no covariates
##' # latent states are position and velocity
##' # initial position is free and called 'inipos'
##' # initial slope is fixed at 1
##' # initial covariance is fixed to a diagonal matrix of 1s
##' initialNoC <- prep.initial(
##' 	values.inistate=c(0, 1),
##' 	params.inistate=c('inipos', 'fixed'),
##' 	values.inicov=diag(1, 2),
##' 	params.inicov=diag('fixed', 2))
##' 
##' #### One covariate
##' # Single regime, one covariate on the inital mean
##' # latent states are position and velocity
##' # initial covariance is fixed to a diagonal matrix of 1s
##' # initial latent means have
##' #   nrow = numLatentState, ncol = numCovariates + 1
##' # initial position has free intercept and free u1 effect
##' # initial slope is fixed at 1
##' initialOneC <- prep.initial(
##' 	values.inistate=matrix(
##' 		c(0, .5,
##' 		  1,  0), byrow=TRUE,
##' 		nrow=2, ncol=2),
##' 	params.inistate=matrix(
##' 		c('iniPosInt', 'iniPosSlopeU1',
##' 		'fixed', 'fixed'), byrow=TRUE,
##' 		nrow=2, ncol=2),
##' 	values.inicov=diag(1, 2),
##' 	params.inicov=diag('fixed', 2),
##' 	covariates='u1')
##' 
##' #### Regime-switching, one covariate
##' # latent states are position and velocity
##' # initial covariance is fixed to a diagonal matrix of 1s
##' # initial latent means have
##' #   nrow = numLatentState, ncol = numCovariates + 1
##' # initial position has free intercept and free u1 effect
##' # initial slope is fixed at 1
##' # There are 3 regimes but the mean and covariance
##' #   are not regime-switching.
##' initialRSOneC <- prep.initial(
##' 	values.regimep=matrix(
##' 		c(1, 1,
##' 		  0, 1,
##' 		  0, 0), byrow=TRUE,
##' 		nrow=3, ncol=2),
##' 	params.regimep=matrix(
##' 		c('r1int', 'r1slopeU1',
##' 		  'r2int', 'r2slopeU2',
##' 		  'fixed', 'fixed'), byrow=TRUE,
##' 		nrow=3, ncol=2),
##' 	values.inistate=matrix(
##' 		c(0, .5,
##' 		  1,  0), byrow=TRUE,
##' 		nrow=2, ncol=2),
##' 	params.inistate=matrix(
##' 		c('iniPosInt', 'iniPosSlopeU1',
##' 		'fixed', 'fixed'), byrow=TRUE,
##' 		nrow=2, ncol=2),
##' 	values.inicov=diag(1, 2),
##' 	params.inicov=diag('fixed', 2),
##' 	covariates='u1')
##' 
prep.initial <- function(values.inistate, params.inistate, values.inicov, params.inicov, values.regimep=1, params.regimep=0, covariates, deviation=FALSE, refRow){
    if(missing(covariates)){
        covariates <- character(0)
    }
    if(missing(refRow)){ # or other default values???
        refRow <- numeric(0)
    }
    
    # Handle initial state
    r <- coProcessValuesParams(values.inistate, params.inistate)
    values.inistate <- r$values
    params.inistate <- r$params
    
    # Handle initial covariance
    r <- coProcessValuesParams(values.inicov, params.inicov)
    values.inicov <- r$values
    params.inicov <- r$params
    
    
    values.inistate <- lapply(values.inistate, preProcessValues)
    params.inistate <- lapply(params.inistate, preProcessParams)
    
    
    values.inicov <- lapply(values.inicov, preProcessValues)
    params.inicov <- lapply(params.inicov, preProcessParams)
    lapply(values.inicov, checkSymmetric, name="values.inicov")
    lapply(params.inicov, checkSymmetric, name="params.inicov")
    
    if(nrow(values.inistate[[1]]) != nrow(values.inicov[[1]])){
        stop(paste0('Number of latent variables implied by initial state and initial covariance differ:\n',
            'initial state (', nrow(values.inistate[[1]]), '), initial cov (', nrow(values.inicov[[1]]), ')'))
    }
    if(ncol(values.inistate[[1]]) != (length(covariates) + 1)){
		stop(paste0('Incorrect dimensions for initial state\nFound ', paste0(dim(values.inistate[[1]]), collapse=' by '),
			' but should be ', nrow(values.inistate[[1]]), ' by ', length(covariates) + 1, '\n',
			'k by (c+1) for k=number of latent variables, c=number of covariates\n',
			"Maybe you forgot to add your covariates to the 'covariates' argument of prep.initial()\n",
			"or forgot to add columns for the covariates to inistate"))
	}
    
    if(identical(values.regimep, 1) && identical(params.regimep, 0)){
        values.regimep <- matrix(c(1, rep(0, length(covariates))), 1, length(covariates)+1)
        params.regimep <- matrix(c(0), nrow(values.regimep), length(covariates)+1)
    }
    values.regimep <- preProcessValues(values.regimep)
    params.regimep <- preProcessParams(params.regimep)
    
    if(ncol(values.regimep) != (length(covariates) + 1)){
        stop(paste0('Incorrect dimensions for initial probabilities\nFound ', paste0(dim(values.regimep), collapse=' by '),
            ' but should be ', nrow(values.regimep), ' by ', length(covariates) + 1, '\n',
            'r by (c+1) for r=number of regimes, c=number of covariates'))
    }
    
    # Do lots or processing and checking for the reference row and column
    if(deviation == TRUE){
        # if needed set refRow
        if(length(refRow) == 0){
            # check valid reference row (i.e. there is a single row as reference)
            # detect reference row
            # reference row has all zero covariate effects
            covAreZero <- apply(params.regimep, 1, function(x){all(x[-1] %in% c('fixed'))})
            blockSize <- max(length(covariates), 1)
            covAreZeroM <- matrix(covAreZero, nrow=blockSize, ncol=nrow(values.regimep))
            refRow <- apply(covAreZeroM, 2, function(x){all(x==TRUE)})
            if(sum(refRow) < 1){ #no single reference row found
                stop(paste("No single reference row found. Initial probabilities might not be identified. Set e.g. refRow=1."))
            }
            if(sum(refRow) > 1){#found more than one reference row
                stop(paste0("Found multiple possible reference rows: ", paste(which(refRow==1), collapse=', '), ".  Reconsider your model specification or just set e.g. refRow=1."))
            }
            # set reference row to reference column
            refRow <- which(refRow)
        }
        if(length(refRow) > 1){
            # error reference row must have length 0 or 1
            stop(paste("'refRow' must being a single number. That is, have length 0 or 1, we found length", length(refRow)))
        }
        if(refRow > nrow(values.regimep)){
            # error reference row must be between 1 and nrow(values)
            stop(paste("'refRow' must be between 1 and", nrow(values.regimep), "but we found", refRow))
        }
    } else {
        if(length(refRow) > 0){
            #warning refRow is ignored when for the non-deviation case (deviation=FALSE)
            warning("'refRow' argument is ignored in the non-deviation case (i.e. when deviation=FALSE)")
        }
    }
    
    # Check that all 'values' imply 0, 1, or the same number of regimes.
    # Note that the 'values' and 'params' have already been checked to imply this.
    nregs <- sapply(list(values.inistate=values.inistate, values.inicov=values.inicov), length)
    nregs <- c(nregs, nrow(values.regimep))
    mr <- max(nregs)
    msg <- paste0("Initial state means, initial state covariance matrix, and initial regime probabilities imply different numbers of regimes:\n'inistate' has ",
            nregs[1], ", 'inicov' has ",  nregs[2], ", and 'regimep' has ", nregs[3], " regimes.\nEven Black Eyed Peas know that's not how you get it started.")
    if(!all(nregs %in% c(1, mr)) || (nregs[3]==1 & any(nregs[-3] > 1))){
        stop(msg)
    } else if(!all(nregs %in% mr)) {
        if(nregs[1] == 1){
            ae <- autoExtendSubRecipe(values.inistate, params.inistate, 'values.inistate', 'initial states', mr)
            values.inistate <- ae[[1]]
            params.inistate <- ae[[2]]
        }
        if(nregs[2] == 1){
            ae <- autoExtendSubRecipe(values.inicov, params.inicov, 'values.inicov', 'initial covariances', mr)
            values.inicov <- ae[[1]]
            params.inicov <- ae[[2]]
        }
    }
    
    # check if values.inicov is positive definite
    if (min(diag(MASS::ginv(values.inicov[[1]])))<0){
    stop("Starting values for entries in random.values.inicov
        and values.inicov generate a non-positive 
         definite covariance matrix. Please check
         the starting values for these entries.")
    }
  
    values.inicov.inv.ldl <- lapply(values.inicov, replaceDiagZero)
    values.inicov.inv.ldl <- lapply(values.inicov.inv.ldl, reverseldl)
    
    sv <- c(extractValues(values.inistate, params.inistate), extractValues(values.inicov.inv.ldl, params.inicov, symmetric=TRUE), extractValues(values.regimep, params.regimep))
    pn <- c(extractParams(params.inistate), extractParams(params.inicov), extractParams(params.regimep))
    sv <- extractValues(sv, pn)
    pn <- extractParams(pn)
    
    x <- list(startval=sv, paramnames=pn, values.inistate=values.inistate, params.inistate=params.inistate, values.inicov=values.inicov, values.inicov.inv.ldl=values.inicov.inv.ldl, params.inicov=params.inicov, values.regimep=values.regimep, params.regimep=params.regimep, covariates=covariates, deviation=deviation, refRow=refRow)
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
##' ##' 
##' @details
##' Prepares a dynr recipe that specifies the names of the parameters that are 
##' to be subjected to user-supplied transformation functions and the 
##' corresponding transformation and reverse-transformation functions. 
##' This can be very handy in fitting dynamic models in which certain parameters can 
##' only take on permissible values in particular ranges (e.g., a parameter may 
##' have to positive). Note that all variance-covariance parameters in the model
##' are automatically subjected to transformation functions to ensure that
##' the resultant covariance matrices are positive-definite. Thus, no additional
##' transformation functions are needed for variance-covariance parameters.
##' 
##' @return Object of class 'dynrTrans'
##' 
##' @examples
##' #Specifies a transformation recipe, r20, that subjects the parameters
##' #'r10' and 'r20' to exponential transformation to ensure that they are positive.
##' trans <-prep.tfun(formula.trans=list(r10~exp(r10), r20~exp(r20)),
##'                   formula.inv=list(r10~log(r10),r20~log(r20)))
##'
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

formula2string <- function(formula.list){
	tuple <- lapply(formula.list, as.list)
	lhs <- sapply(tuple, function(x){paste0(deparse(x[[2]], width.cutoff = 500L), collapse="")})
	rhs <- sapply(tuple, function(x){paste0(deparse(x[[3]], width.cutoff = 500L), collapse="")})
	return(list(lhs=lhs, rhs=rhs))
}


matrix2formula <- function(x, multbyColnames=TRUE, convertWholeMatrix = FALSE){
	if(!is.matrix(x)){
		stop("Dude! You have to give me a matrix. If you do, I'll give you a formula. Seriously.")
	}
	if(is.null(rownames(x))){
		rownames(x) <- paste0('y', 1:nrow(x))
	}
	if(is.null(colnames(x))){
		colnames(x) <- paste0('x', 1:ncol(x))
	}
	
	if (convertWholeMatrix == FALSE){
	  preds <- character(nrow(x))
	for(i in 1:nrow(x)){
		if (multbyColnames==FALSE){
			preds[i] <- paste(rownames(x)[i], paste(x[i,]), sep=' ~ ')
		}else{
			preds[i] <- paste(rownames(x)[i], paste(x[i,],
		                      colnames(x), sep='*', collapse=' + '), sep=' ~ ')
			}
	}}else{ #If wanting to convert the whole matrix into formula
	  preds <- character(nrow(x)*ncol(x))
	  for(i in 1:nrow(x)){
	    for (j in 1:ncol(x)){
	      preds[j+(i-1)*ncol(x)] <- paste(rownames(x)[i], paste(x[i,j]), sep=' ~ ')
	    }}  
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
				list_list_formulae_add <- get(VecNamesToAdd[l], parent.frame()) #eval(parse(text=VecNamesToAdd[l]),environment())
				if (length(list_list_formulae_add) > 0){
					if (length(list_list_formulae_add[[j]]) > 0){
						AddedFml <- addFormulas(list_list_formulae_add[[j]][[k]], AddedFml)
					}
				}
			}
			list_list_formulae[[j]][[k]] <- AddedFml
		}
	}
	return(list_list_formulae)
}



#------------------------------------------------------------------------------
# Formula to matrix functions

# TODO interpret strict as character with match.arg.  strict=c('allowOne' let's people fly with x without writing 1*x, 'convent' let's no one get away with anything, 'hippie' is very loose
formula2matrix <- function(formula, variables, strict=TRUE, col.match=TRUE, process=TRUE){
    if(process){
        pf <- dynr:::processFormula(list(formula))[[1]]
        lhs <- pf[1]
        rhs <- pf[2]
        lhs2 <- sapply(lhs, strsplit, split=' + ', fixed=TRUE)[[1]]
        rhs2 <- sapply(rhs, strsplit, split=' + ', fixed=TRUE)[[1]]
        rhs3 <- strsplit(rhs2, split=" * ", fixed=TRUE)
        if(strict == FALSE){
            rhs4 <- lapply(rhs3, function(x){if(length(x) == 1) c("1", x) else x})
        } else {rhs4 <- rhs3}
    } else {
        lhs2 <- formula[[1]]
        rhs4 <- formula[[2]]
    }
    lens <- sapply(rhs4, length)
    if(any(lens > 2)){
        stop(paste0("I spy with my little eye terms in the right hand side of your formula with more than two parts.\n",
            "This function can only handle linear formulas.\nOffending term has: ",
            paste(rhs4[lens > 2], collapse=", ")))
    }
    if(any(lens < 2)){
        stop(paste0("I spy with my little eye terms in the right hand side of your formula with less than two parts.\n",
            "Perhaps you forgot to multiply by 1.\nOffending term has: ",
            paste(rhs4[lens < 2], collapse=", ")))
    }
    
    eleNames <- sapply(rhs4, function(x) {x[!(x %in% variables)]} )
    ncol <- ifelse(col.match, length(variables), length(eleNames))
    cnam <- if(col.match) variables else eleNames
    
    rmat <- matrix("0", nrow=length(lhs2), ncol=ncol, dimnames=list(lhs2, cnam))
    for(aterm in 1:length(rhs4)){
        x <- rhs4[[aterm]]
        if(sum(x %in% colnames(rmat)) == 2){stop("Both parts of term are in the 'variables'")}
        # TODO warning/error if strict = TRUE and trying to overwrite nonzero element
        # if strict = FALSE try to write next nonzero element that matches
        if(x[2] %in% colnames(rmat)){
            if(strict==TRUE & !all(rmat[, x[2]] %in% "0")){
                warning(paste0("Overwriting element in column ", x[2], ". It was ", rmat[1, x[2]], " and now will be ", x[1], "."), call.=FALSE)
            }
            rmat[, colnames(rmat) %in% x[2]] <- x[1]
        } else if(x[1] %in% colnames(rmat)){
            if(strict==TRUE & !all(rmat[, x[1]] %in% "0")){
                warning(paste0("Overwriting element in column ", x[1], ". It was ", rmat[1, x[1]], " and now will be ", x[2], "."), call.=FALSE)
            }
            rmat[, colnames(rmat) %in% x[1]] <- x[2]
        }
    }
    return(rmat)
}


# Example useage
#formula2matrix(theta1 ~ x1 + x2 + mylabel*x3, variables=c('x1', 'x2', 'x3'), strict=FALSE)

#formula2matrix(theta1 ~ x1 + x2 + mylabel*x3, variables=c('1', 'mylabel'), strict=FALSE)

#formula2matrix(theta1 + theta2 ~ x1 + x2 + mylabel*x3, variables=c('x1', 'x2', 'x3'), strict=FALSE)

## Factor model
#formula2matrix(x1 + x2 + x3 + x4 + x5 ~ F, variables="F", strict=FALSE)
#t(formula2matrix(F ~ x1 + x2 + x3 + x4 + x5, variables=paste0('x', 1:5), strict=FALSE))
#t(formula2matrix(F ~ x1 + load2*x2 + load3*x3 + load4*x4 + load5*x5, variables=paste0('x', 1:5), strict=FALSE))

## Regression among latent variables
#formula2matrix(F1 ~ F2 + F3, variables=c("F2", "F3"), strict=FALSE)

# original version
# formula2design <- function(dots, covariates, random.names, beta.names){
    # #dots <- list(...)  
    
    # # Extending the theta formula in dots to add variables in beta. names
	# # beta.names include variables without random effects but needs to be estimated.
	# if(length(beta.names) > 0){
		# for (i in 1:length(beta.names)){
			# dots = c(dots, as.formula(paste(beta.names[i], '~ 1 *', beta.names[i])))
		# }
	# }
    
    
    # pf <- list(dynr:::processFormula(dots))
    # lhs <- unlist(lapply(pf,function(x){lapply(x,"[[",1)})[[1]])
    # rhs <- lapply(pf,function(x){lapply(x,"[[",2)})[[1]]
    # rhs2 <- sapply(rhs, strsplit, split=' + ', fixed=TRUE)
    # rhs3 <- lapply(rhs2, strsplit, split=" * ", fixed=TRUE)
    # rhs4 <- lapply(rhs3, function(x){if(length(x) == 1) c("1", x) else x})
    # eleNames <- unique(unlist(sapply(rhs3, function(rlist){ sapply(rlist, function(x) {x[!(x %in% covariates) & !(x %in% random.names)]} ) } )))
	
    
    # dots <- cbind(lhs, rhs3)
    
    # fmat <- t(apply(dots, 1, formula2matrix, variables=eleNames, strict=TRUE, process=FALSE))
    # rmat <- t(apply(dots, 1, formula2matrix, variables=random.names, strict=TRUE, process=FALSE))
    # dimnames(fmat) <- list(lhs, eleNames)
    # dimnames(rmat) <- list(lhs, random.names)
    
    # #print(fmat) #H
    # #print(rmat) #Z 
    # return(list(fixed=as.matrix(fmat), random=as.matrix(rmat)))
# }

# formula2design <- function(..., covariates, random.names, beta.names){
  # #browser()
  # dots <- unlist(list(...))  
  # len_dots = length(dots)
  
  # # Extending the theta formula in dots to add variables in beta. names
  # # beta.names include variables without random effects but needs to be estimated.
  # if(length(beta.names) > 0){
    # for (i in 1:length(beta.names)){
      # dots = c(dots, as.formula(paste(beta.names[i], '~ 1 *', beta.names[i])))
    # }
  # }
  
  
  # pf <- list(dynr:::processFormula(dots))
  # lhs <- unlist(lapply(pf,function(x){lapply(x,"[[",1)})[[1]])
  # rhs <- lapply(pf,function(x){lapply(x,"[[",2)})[[1]]
  # rhs2 <- sapply(rhs, strsplit, split=' + ', fixed=TRUE)
  # rhs3 <- lapply(rhs2, strsplit, split=" * ", fixed=TRUE)
  # rhs4 <- lapply(rhs3, function(x){if(length(x) == 1) c("1", x) else x})
  # eleNames <- unique(unlist(sapply(rhs3, function(rlist){ sapply(rlist, function(x) {x[!(x %in% covariates) & !(x %in% random.names)]} ) } )))
  
  
  # dots <- cbind(lhs, rhs3)
  # fmat <- t(apply(dots, 1, formula2matrix, variables=eleNames, strict=TRUE, process=FALSE))
  # rmat <- t(apply(dots, 1, formula2matrix, variables=random.names, strict=TRUE, process=FALSE))
  # fmat = fmat[seq_len(len_dots), , drop=FALSE]
  # rmat = rmat[seq_len(len_dots), , drop=FALSE]
  # dimnames(fmat) <- list(lhs[1:len_dots], eleNames)
  # dimnames(rmat) <- list(lhs[1:len_dots], random.names)
  
  # #browser()
  # #fmat = fmat[1:len_dots, ]
  # #rmat = rmat[1:len_dots, ]
  # return(list(fixed=as.matrix(fmat), random=as.matrix(rmat)))
# }

#formula2design_symiin
formula2design <- function(..., covariates, random.names, startval.names){
  #browser()
  dots <- unlist(list(...))   
  len_dots = length(dots) # = length of theta_formula
  
  pf <- list(dynr:::processFormula(dots))
  lhs <- unlist(lapply(pf,function(x){lapply(x,"[[",1)})[[1]])
  rhs <- lapply(pf,function(x){lapply(x,"[[",2)})[[1]]
  rhs2 <- sapply(rhs, strsplit, split=' + ', fixed=TRUE)
  rhs3 <- lapply(rhs2, strsplit, split=" * ", fixed=TRUE)
  rhs4 <- lapply(rhs3, function(x){if(length(x) == 1) c("1", x) else x})
  thetaFormula_beta.names <- unique(unlist(sapply(rhs4, function(rlist){ 
    sapply(rlist, function(x) {x[!(x %in% covariates) & !(x %in% random.names)]} ) } )))
  
  
  # Extending the theta formula in dots to add variables in startvals that are not in theta formula
  # beta.names include variables without random effects but needs to be estimated.
  variablesToAdd = startval.names[!startval.names %in% thetaFormula_beta.names]
  if(length(variablesToAdd) > 0){
    for (i in 1:length(variablesToAdd)){
      toAddRandom = NULL
      if (length(random.names)>0){
      toAddRandom = paste0(" + 0 * ",random.names[1])}
      dots = c(dots, as.formula(paste0(variablesToAdd[i], '~ 1 *', variablesToAdd[i],toAddRandom)))
    }}
  
  pf <- list(dynr:::processFormula(dots))
  lhs <- unlist(lapply(pf,function(x){lapply(x,"[[",1)})[[1]])
  rhs <- lapply(pf,function(x){lapply(x,"[[",2)})[[1]]
  rhs2 <- sapply(rhs, strsplit, split=' + ', fixed=TRUE)
  rhs3 <- lapply(rhs2, strsplit, split=" * ", fixed=TRUE)
  #rhs4 <- lapply(rhs3, function(x){if(length(x) == 1) c("1", x) else x})
  eleNames <- unique(unlist(sapply(rhs3, function(rlist){ 
    sapply(rlist, function(x) {x[!(x %in% covariates) & !(x %in% random.names)&
           (x!="0") & (x!= "1")]})})))
  
  dots <- cbind(lhs, rhs3)
  fmat <- t(apply(dots, 1, formula2matrix, variables=eleNames, strict=TRUE, process=FALSE))
  rmat <- t(apply(dots, 1, formula2matrix, variables=random.names, strict=FALSE, process=FALSE))
  #fmat = fmat[seq_len(len_dots), , drop=FALSE]
  rmat = matrix(rmat,ncol=length(random.names),byrow=TRUE)
  dimnames(fmat) <- list(lhs, eleNames)
  dimnames(rmat) <- list(lhs, random.names)
  
  #browser()
  #fmat = fmat[1:len_dots, ]
  #rmat = rmat[1:len_dots, ]
  return(list(fixed=as.matrix(fmat), random=as.matrix(rmat)))
}


# Example usage
# formula2design(
    # theta1 ~ 1*zeta_0 + u1*zeta_1 + u2*zeta_2 + 1*b_zeta,
    # theta2 ~ mu_x1*1 + 1*b_x1,
    # theta3 ~ mu_x2*1 + 1*b_x2,
    # covariates=c('1', 'u1', 'u2'),
    # random.names=c('b_zeta', 'b_x1', 'b_x2'))
# lme4::lmer style interface
#theta1 ~ u1 + u2 + (1 + u2 | Subject)
#theta2 ~ 1 + (1 | Subject)
#theta3 ~ 1 + (1 | Subject)




multiformula2matrix <- function(..., variables){
    dots <- list(...)
    
    pf <- list(dynr:::processFormula(dots))
    lhs <- unlist(lapply(pf,function(x){lapply(x,"[[",1)})[[1]])
    rhs <- lapply(pf,function(x){lapply(x,"[[",2)})[[1]]
    rhs2 <- sapply(rhs, strsplit, split=' + ', fixed=TRUE)
    rhs3 <- lapply(rhs2, strsplit, split=" * ", fixed=TRUE)
    #rhs4 <- lapply(rhs3, function(x){if(length(x) == 1) c("1", x) else x})
    
    dots <- cbind(lhs, rhs3)
    
    fmat <- t(apply(dots, 1, formula2matrix, variables=variables, strict=TRUE, process=FALSE))
    dimnames(fmat) <- list(lhs, variables)
    return(fmat)
}

# Example useage
# Regression among latent variables
#multiformula2matrix(F1 ~ a*F2 + b*F3, F2 ~ 0*F2, F3 ~ c*F2, variables=c("F1", "F2", "F3"))


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

# fromName character name of the from vector
# toName character name of the to vector
# fromLoc numeric integer vector of from locations
# toLoc numeric integer vector of to locations
# depth number of tabs to indent
# create logical whether to create the toName before copying.  The vector created will be the length of toLoc
gslVectorCopy <- function(fromName, toName, fromLoc, toLoc, fromFill="", toFill="", depth=1, create=FALSE){
	# check lengths match
	if(length(fromLoc) != length(toLoc)){stop("'fromLoc' and 'toLoc' lengths must match")}
	tabs <- paste(rep('\t', depth), collapse="")
	tabsm1 <- paste(rep('\t', depth-1), collapse="")
	ret <- character(0)
	if(create){
		ret <- paste0(ret, tabsm1, createGslVector(length(toLoc), toName))
	}
	for(i in 1:length(toLoc)){
		get <- paste0("gsl_vector_get(", fromName, ", ", fromFill, fromLoc[i]-1, ")")
		set <- paste0(tabs, "gsl_vector_set(", toName, ", ", toFill, toLoc[i]-1, ", ", get, ");\n")
		ret <- paste0(ret, set)
	}
	return(ret)
}

gslcovariate.front <- function(selected, covariates){
	gslVectorCopy("co_variate", "covariate_local", match(selected, covariates), 1:length(selected), create=TRUE)
}

#------------------------------------------------------------------------------
# functions that are used in only SAEM

#retriveStateFormula: remains only formulas whose LHS is in state.names 
retriveStateFormula <- function(formula, state.names){
    #fun
    fml=lapply(formula, as.character)
    lhs=lapply(fml,function(x){x[[2]]})
    rhs=lapply(fml,function(x){x[[3]]})
    
    isStateVariables <- 1:length(lhs)
    for(i in 1:length(formula[[1]])){
        isStateVariables[[i]] = FALSE;
        formula[[i]]=as.character(formula[[i]])
        for (j in 1:length(state.names)){
            if(lhs[[i]] == state.names[[j]]){
                isStateVariables[[i]] = TRUE;
                break;
            }
        }
        #formula[[i]]=as.formula(paste0(lhs[[i]], ' ~ ', rhs[[i]]))
    }
    
        
    j = 1;
    formula2= list();
    for(i in 1:length(formula)){
        if(isStateVariables[[i]] == TRUE){
            formula2[[j]]=as.formula(paste0(lhs[[i]], ' ~ ', rhs[[i]]))   
            j = j +1
        }
    }
    return(list(formula2))
}

#parseFormulaTheta: retrieve the LHS and RHS of the theta.formula
parseFormulaTheta <- function(formula, theta.formula){
    #fun
    fml=lapply(formula, as.character)
    lhs=lapply(fml,function(x){x[[2]]})
    rhs=lapply(fml,function(x){x[[3]]})
    
    #thetai
    fmlt=lapply(theta.formula, as.character)
    lhst=lapply(fmlt,function(x){x[[2]]})
    rhst=lapply(fmlt,function(x){x[[3]]})

    for(i in 1:length(formula)){
        formula[[i]]=as.character(formula[[i]])
        for (j in 1:length(lhst)){
            # gsub (a, b, c) : in c replace a with b
            rhs[[i]]=gsub(paste0(lhst[[j]]),paste0("(",rhst[[j]],")"),rhs[[i]], fixed = TRUE)
        }
        formula[[i]]=as.formula(paste0(lhs[[i]], ' ~ ', rhs[[i]]))
    }
    

    return(formula)
}


#function for parsing theta.formula (remove intercept terms and random names) 
prep.thetaFormula <- function(formula, intercept.names, random.names){
    if(length(formula) == 0)
      return(list())
 
    fml=lapply(formula, as.character)
    lhs=lapply(fml,function(x){x[[2]]})
    rhs=lapply(fml,function(x){x[[3]]})
    
    for(i in 1:length(formula)){
        formula[[i]]=as.character(formula[[i]])
        if(length(intercept.names) > 0){
            for (j in 1:length(intercept.names)){
                rhs[[i]]=gsub(paste0(intercept.names[j]),paste0("0"),rhs[[i]], fixed = TRUE)
            }
        }
        
        if(length(random.names) > 0){
            for (j in 1:length(random.names)){
                rhs[[i]]=gsub(paste0(random.names[j]),paste0("0"),rhs[[i]], fixed = TRUE)
            }
        }
        
        formula[[i]]=as.formula(paste0(lhs[[i]], ' ~ ', rhs[[i]]))
    }
    return(formula)
}


#---
#cook for random effect b in SAEM (currently not used)
# prep.random<- function(random.names, random.lb=rep(-5, length(random.names)), random.ub=rep(5, length(random.names)), num.subj, params.inicov, values.inicov){
    # if(length(random.names) != length(random.lb) || length(random.names) != length(random.ub)){stop("The number of variables is different from the number of elements in lower/upper bound")}
    
    # b <- matrix(rnorm(num.subj*length(random.names)), num.subj, length(random.names))
    # for(i in 1:num.subj){
        # b[i, b[i, ] < random.lb | b[i, ] > random.ub] <-  0
    # }

    # if(nrow(params.inicov) != nrow(values.inicov)){
        # stop('Unmatched dimension in random effect variance-covariance matrix')
    # }
    

    # x <- list(random.names= random.names, random.lb = random.lb, random.ub = random.ub, b=b, params.inicov=params.inicov, values.inicov=values.inicov)
    # return(new("dynrRandom", x))
# }



##' Recipe function for specifying the control parameters used in the SAEM estimation
##' 
##' @param MAXGIB A positive integer indicates the number of iterations to be run in the Gibbs sampler (Default: 200). 
##' @param MAXITER A positive integer indicates the maximum number of iterations of SAEM (Default: 200).
##' @param maxIterStage1 A positive integer indicates the maximum number of iterations of SAEM Stage 1 (Default: 100).
##' @param gainpara A floating number; SAEM control parameter (Default: 0.6).
##' @param gainparb A floating number; SAEM control parameter (Default: 3).
##' @param gainpara1 A floating number; SAEM control parameter (Default: 0.9).
##' @param gainparb1 A floating number; SAEM control parameter (Default: 1).
##' @param bAdaptParams A vector of two floating number \code{[min, max]} to specify the control parameters in the MCMC sampler of b (Default: c(0.5, 2.5)).
##' @param KKO A integer. The SAEM process only starts to evaluate whether to transit to stage 2 after \code{KKO} iterations (Default: 20).
##' @param scaleb A floating number indicating the initial value of \code{scaleb}, which is a weighting parameter to control the b acceptance rate (Default: 1). 
##' @param setScaleb A booling value. If \code{setScaleb} is set to 0, \code{scaleb} will remain in its initial value. If \code{setScaleb} is set to 1, \code{scaleb} will be updated in the MCMC sampler of b (Default: 1).
##' @param setAccept A floating number, which is a control parameter of the b acceptance rate (Default: 0.8).
##' @param seed An integer, the random seed. If \code{seed} is not specified, a random number will be given as seed (Default: NA).
##' 
##' @details
##' Use @ to show specific arguments from a dynrSaemParameter object, see the example.
##' 
##' @return Object of class 'dynrSaemParameter'
##' 
##'  
##' @examples 
##' saemp <- prep.saemParameter(MAXGIB = 200, MAXITER = 200, maxIterStage1 = 100)
##' print(saemp@MAXGIB)
##' ## 200
##' print(saemp@setAccept)
##' ## 0.8
##'
prep.saemParameter<- function(MAXGIB = 50, MAXITER = 200, 
                              maxIterStage1 = 100, 
                              gainpara = 0.600000, 
                              gainparb = 3.000000, 
                              gainpara1 = 0.900000, 
                              gainparb1 = 1.000000, 
                              bAdaptParams = c(1, 2), KKO = 10, scaleb = 1, 
                              setScaleb = 0, setAccept = 0.8, seed = NA, trueb = matrix(),
                              errtrol1 = 1, errtrol = 1){ #errtrol1 and errtrol are stages 1 and 2 error tolerance levels in saem
    if(is.numeric(seed) == TRUE)
		seed <- as.integer(seed)
	else
		seed <- as.integer((100000-1) * runif(1) + 1)
    x <- list(MAXGIB = MAXGIB, MAXITER = MAXITER, maxIterStage1 = maxIterStage1, gainpara = gainpara, 
              gainparb = gainparb, gainpara1 = gainpara1, gainparb1 = gainparb1,  
              bAdaptParams = bAdaptParams, KKO = KKO, scaleb=scaleb, 
              setScaleb = setScaleb, setAccept=setAccept, seed = seed, 
              trueb=trueb, errtrol1 = errtrol1, errtrol = errtrol)
    return(new("dynrSaemParameter", x))
}

#------------------------------------------------------------------------------
# functions that are used in only SAEM




# A internal-use only function for substituting formula
# @param formula a list of original formulas
# @param term.formula a list of term formulas
# 
# @details
# If the RHS of \code{formula} has terms in the LHS of \code{term.formula}, this function replaces any appearance with the RHS of \code{term.formula}.
#
# @return a list of formulas after the replacement
# @examples
# #substitutedformula <- substituteFormula(formula, term.formula)

substituteFormula <- function(formula, term.formula){
	#parseFormulaTheta <- function(formula, theta.formula){
    #fun
	fml=lapply(formula, as.character)
    lhs=lapply(fml,function(x){x[[2]]})
    rhs=lapply(fml,function(x){x[[3]]})
    
    #thetai
    fmlt=lapply(term.formula, as.character)
    lhst=lapply(fmlt,function(x){x[[2]]})
    rhst=lapply(fmlt,function(x){x[[3]]})
	
	#deciding the substitute order: variables with longer name first (e.g., substitute zeta_i then eta_i
	sub_order = order(1:length(lhst), key=sapply(1:length(lhst), function(i)nchar(lhst[[i]])), decreasing= TRUE)
	
    for(i in 1:length(formula)){
        formula[[i]]=as.character(formula[[i]])
		for (j in sub_order){
			# gsub (a, b, c) : in c replace a with b
			rhs[[i]]=gsub(paste0(lhst[[j]]),paste0("(",rhst[[j]],")"),rhs[[i]], fixed = TRUE)
		}
		formula[[i]]=as.formula(paste0(lhs[[i]], ' ~ ', rhs[[i]]))
	}
	
	
	return(formula)
}


##' A internal-use function for removing intercept terms and random names in the formula 
##' @param formula the formula that are going to be processed
##' @param intercept.names variables that are going to be removed if it appears in formula
##' @param random.names variables that are going to be removed if it appears in formula
#prep.thetaFormula <- function(formula, intercept.names, random.names){
#    if(length(formula) == 0)
#     return(list())
# 
#    fml=lapply(formula, as.character)
#    lhs=lapply(fml,function(x){x[[2]]})
#    rhs=lapply(fml,function(x){x[[3]]})
#    
#    for(i in 1:length(formula)){
#        formula[[i]]=as.character(formula[[i]])
#       if(length(intercept.names) > 0){
#           for (j in 1:length(intercept.names)){
#               rhs[[i]]=gsub(paste0(intercept.names[j]),paste0("0"),rhs[[i]], fixed = TRUE)
#           }
#       }
#        
#       if(length(random.names) > 0){
#           for (j in 1:length(random.names)){
#               rhs[[i]]=gsub(paste0(random.names[j]),paste0("0"),rhs[[i]], fixed = TRUE)
#           }
#       }
#        
#        formula[[i]]=as.formula(paste0(lhs[[i]], ' ~ ', rhs[[i]]))
#    }
#    return(formula)
#}

# returnExponentialSymbolicTerm <- function(inputs){
	# #browser()
	# if(is.character(inputs[1,1])){
		# ret <- sapply(inputs, function(term){
								# as.list(as.formula(paste0('y ~ exp(', term, ')')))[[3]]})
		
	# } else {
		# #ret <- sapply(inputs, function(term){
		# #						return(as.list(as.formula(paste0('y ~ exp(', deparse(term[[1]], width.cutoff = 500), ')')))[[3]])})
		# ret = rep(list(), ncol(inputs) *nrow(inputs))
		# for(i in 1:nrow(inputs)){
			# for(j in 1:ncol(inputs)){
				# ret[nrow(inputs)*(j-1)+i] <- list(as.list(as.formula(paste0('y ~ exp(', deparse(inputs[i,j][[1]], width.cutoff = 500), ')')))[[3]])
		# }}
	# }
	# ret <- matrix(ret, ncol = ncol(inputs), nrow=nrow(inputs))
	# return(ret)
# }

# only do exp to diag ones
returnExponentialSymbolicTerm <- function(inputs){
	#browser()
	
	ret = rep(list(), ncol(inputs) *nrow(inputs))
	for(i in 1:nrow(inputs)){
		for(j in 1:ncol(inputs)){
			if(i == j){
			if(is.character(inputs[i,j])){
				ret[nrow(inputs)*(j-1)+i] <- list(as.list(as.formula(paste0('y ~ exp(', inputs[i,j], ')')))[[3]])
			} else {
				ret[nrow(inputs)*(j-1)+i] <- list(as.list(as.formula(paste0('y ~ exp(', deparse(inputs[i,j][[1]], width.cutoff = 500), ')')))[[3]])
			}
			} else {
			#if(is.character(inputs[i,j])){
				ret[nrow(inputs)*(j-1)+i] <- list(as.list(as.formula(paste0('y ~ (', inputs[i,j], ')')))[[3]])
			#} else {
				#ret[nrow(inputs)*(j-1)+i] <- list(as.list(as.formula(paste0('y ~ (', deparse(inputs[i,j][[1]], width.cutoff = 500), ')')))[[3]])
			#}
			}
		}
	}
	
	ret <- matrix(ret, ncol = ncol(inputs), nrow=nrow(inputs))
	return(ret)
}


# differentiateMatrixOfVariable <- function(inputs, variable.names=character(0)){
	# #browser()
	# if(is.vector(inputs)){
		# inputs <- as.matrix(inputs, ncol=1)
		# #rownames(inputs) <- inputs
	# }
	
	# if(length(variable.names) > 0){
		# variable.names <- unique(as.vector(variable.names))
		# variable.names <- variable.names[!variable.names%in% c("fixed", "0")]
	# } else {
		# variable.names <- unique(as.vector(inputs))
		# variable.names <- variable.names[!variable.names%in% c("fixed", "0")]
	# }

	# ret <- matrix(0, nrow= length(inputs), ncol= length(variable.names))
	# #rownames(ret) <- rep(rownames(inputs), length(inputs)/nrow(inputs))
	# rownames(ret) <- as.vector(inputs)
	# colnames(ret) <- variable.names
	# for(i in 1:length(inputs)){
		# for(j in 1:length(variable.names)){
			# ret[i,j] <- D(as.symbol(inputs[i]),variable.names[j])
			# #print(paste(matrix[i],variable.names[j], D(as.symbol(matrix[i]),variable.names[j])))
		# }
	# }
	

	# return(t(ret))
# }


# merage differentiateMatrixOfVariable and differentiateMatrixOfVariable2
differentiateMatrixOfVariable <- function(inputs, variable.names=character(0)){
	#inputs in call
	#browser()
	if(is.vector(inputs)){
		inputs <- as.matrix(inputs, ncol=1)
		#rownames(inputs) <- inputs
	}
	else{
		inputs <- as.matrix(as.vector(inputs), ncol=1)
	}
	
	if(length(variable.names) > 0){
		variable.names <- unique(as.vector(variable.names))
		variable.names <- variable.names[!variable.names%in% c("fixed", "0", 0)]
	} else {
		variable.names <- unique(as.vector(inputs))
		variable.names <- variable.names[!variable.names%in% c("fixed", "0", 0)]
	}
	
	if(length(variable.names) == 0){
		return(matrix(list(), nrow= 0, ncol= 0))
	}
	

	ret <- matrix(list(), nrow= length(inputs), ncol= length(variable.names))
	#rownames(ret) <- rep(rownames(inputs), length(inputs)/nrow(inputs))
	rownames(ret) <- as.vector(inputs)
	colnames(ret) <- variable.names
	for(i in 1:length(inputs)){
		for(j in 1:length(variable.names)){
			if(is.character(inputs[[i]])){
				ret[i,j] <- D(as.symbol(inputs[i]),variable.names[j])
			} else{
				ret[i,j] <- list(D(inputs[[i]],variable.names[j]))
				#ret[i,j] <- D(inputs[[i]],variable.names[j])
			}
			#print(paste(matrix[i],variable.names[j], D(as.symbol(matrix[i]),variable.names[j])))
		}
	}

	#ret = matrix(sapply(ret, function(term){term[[1]]}), nrow(ret), ncol(ret))
	return(t(ret))
}



##' Obtain symbolic expressions for \code{L\%*\%D\%*\%t(L)}, and unique parameter names for parameters that are on the unconstrained scale (e.g., par0-par9)
##' 
##' @param a the random effect covariance matrix in \code{dynrModel@random.params.inicov}
##' @return a list of \code{ldl} and \code{pars}, where \code{ldl} is a matrix of the symbolic expressions of \code{L\%*\%D\%*\%t(L)}, and \code{pars} is a list of unique parameter names in \code{ldl}
##' @examples 
##' #dynrModel@random.params.inicov
##' a.params = matrix(c(
##'     "a11","a12",
##'     "a12","a22"
##'     ),byrow=TRUE,ncol=2)
##'
##' #dynrModel@random.values.inicov
##' a.values = matrix(c(
##'            .6, .3,
##'            .3, .6
##'           ),byrow=TRUE,ncol=2)
##'
##' # Obtain symbolic expressions for L%*%D%*%t(L), and unique parameter names for parameters
##' # that are on the unconstrained scale (e.g., par0-par2). 
##' r = symbolicLDLDecomposition(a.params, a.values)
##'
##' # Solve for numerical values of par0-par2 (the unconstrained parameters)
##' # given starting values for the random effect covariance matrix in model@random.params.inicov
##' par.values <- solveStartLDL(r$ldl, a.values)
##'
##' # Evaluate the LDL' expressions in r at the values of par.values and return
##' # the "reassembled" L%*%D%*%t(L) numerical matrices based on the unconstrained par.values
##' a.values2 = matrix(sapply(r$ldl, function(x){eval(x, par.values)}), nrow=2)
##' 
##' # compare the difference between the original LDL and the reassembled LDL
##' a.values2 - a.values
##'
##' @details Exponential constraints are imposed on the diagonal elements of \code{D} to ensure the positive definiteness of the product LDL'.
##' @seealso \code{\link{solveStartLDL}}
symbolicLDLDecomposition <- function(a.params, a.values){
	if(!is.matrix(a.params) || nrow(a.params) != ncol(a.params))
		stop("The variable 'a.params' must be a square matrix")
		
	if(!isSymmetric(a.params))
		stop("The variable 'a.params' must be a symmetric matrix")
	
	
	
	# repalce elements equal to NA or 'NA' as 0
	a.params[ is.na(a.params) ] = 0
	a.params[ a.params == 'NA'] = 0
	a.params[ a.params == 'fixed'] = 0
	a.params[ a.params == '0'] = 0
	# a.params stores the input (with severl inputs processing as above)
	a = a.params
	
	
	# relace elements 'fixed' w/ the corresponding values in a.values
	a = matrix(mapply(function(a, v){if(a == 0) a = v else a = a}, as.list(a), as.list(a.values)), nrow = nrow(a))
	
	n <- nrow(a)
	L <- rep(list(0), n*n)
	D <- rep(list(1), n)
	ret <- rep(list(0), n*n)
	temp <- rep(list(0), n*n)
	
	if(is.character(a[1,1])){
		# the input variable, in the type of "character" 
		a = matrix(sapply(a, function(x){list(x)}), nrow=nrow(a), ncol=ncol(a))
	} else {
		# the input is expression, in the type of "call"
		a = matrix(sapply(a, function(x){return(deparse(unlist(x)))}), nrow=nrow(a), ncol=ncol(a))
	}
	# the element in a should be a list with a element, which contains the expression in the form of character
	
	for(i in 1:n){
		#D[i,i]
		#if(!is.character(a[i,i][[1]]))
		term <- a[i,i][[1]]
		if(i > 1){
			for(k in 1:(i-1)){
				#term <- paste(term, '-', as.character(eval(L[i,k][[1]]*L[i,k][[1]]*D[k,k][[1]])))
				e_l <- evaluateExpression(L[i+(k-1)*n][[1]])
				e_d <- evaluateExpression(D[k][[1]])
				if(!is.na(e_l)){
					if(e_l == 0)
						term <- term
					else
						term <- as.character(paste0(term, '-(', e_l, ')^2*(', D[k][[1]], ')'))
				}
				else if(!is.na(e_d)){
					if(e_d == 0)
						term <- term
					else
						term <- as.character(paste0(term, '-(', L[i+(k-1)*n][[1]], ')^2*(', e_d, ')'))
				}
				else{
					term <- as.character(paste0(term, '-(', L[i+(k-1)*n][[1]], ')^2*(', D[k][[1]], ')'))
				}
			}
		}
		#D[i][[1]] <- as.list(as.formula(paste0('x ~ ',term)))[[3]]
		D[i][[1]] <- term
		if(i < n){
			for(j in (i+1):n){
				term <- a[j,i][[1]]
				
				if(i > 1){
					for(k in 1:(i-1)){
						e_l <- evaluateExpression(L[i+(k-1)*n][[1]])
						e_d <- evaluateExpression(D[k][[1]])
						if(!is.na(e_l)){
							if(e_l == 0)
								term <- term
							else
								term <- as.character(paste0(term, '-(', L[j+(k-1)*n][[1]], ')*(', e_l, ')*(', D[k][[1]], ')'))
						}
						else if(!is.na(e_d) && e_d == 0){
							if(e_d == 0)
								term <- term
							else
								term <- as.character(paste0(term, '-(', L[j+(k-1)*n][[1]], ')*(', L[i+(k-1)*n][[1]], ')*(', e_d, ')'))
						}
						else{
							term <- as.character(paste0(term, '-(', L[j+(k-1)*n][[1]], ')*(', L[i+(k-1)*n][[1]], ')*(', D[k][[1]], ')'))
						}
					}
				}
				#L[j+(i-1)*n][[1]] <- as.list(as.formula(paste0('x ~ ',term)))[[3]]
				e_t <- evaluateExpression(term)
				if(is.na(e_t) || e_t != 0){
					term <- paste0('(', term, ')/(', D[i][[1]],')')
				}
				L[j+(i-1)*n][[1]] <- term
			}
		}
		L[i+(i-1)*n][[1]] <- "1"
		
	}
	
	
	# examine whether we need to have a pari to estimate i
	par_list <- list()
	solved_a <- list() # for each solved element, we assign an arbitrary value (here is its value) for it temporarily
	for(i in 1:n){
	
		e_t <- evaluateExpression(D[i][[1]], solved_a) 
	
		if(is.na(e_t) && a.params[i,i] != 0){
			new_par <- paste0('par', length(par_list))
			D[i][[1]] <- new_par
			#par_list <- c(par_list, new_par)
			par_list[new_par] <- a.params[i,i]
			solved_a[a.params[i,i]] = a.values[i,i]
		}
		else if(!is.na(e_t)){
			D[i][[1]] <- as.character(e_t)
		}
		else{ #this paramter needs to be rewritten in the form of pars
			for(par_ in names(par_list)){
				D[i][[1]] <- gsub(par_list[par_], par_, D[i][[1]])
			}
		}
		
		if(i > 1){
			for(j in 1:(i-1)){
				e_t <- evaluateExpression(L[i+(j-1)*n][[1]], solved_a) 
				
				if(is.na(e_t) && a.params[i,j] != 0){
					new_par <- paste0('par', length(par_list))
					L[i+(j-1)*n][[1]] <- new_par
					#par_list <- c(par_list, new_par)
					par_list[new_par] <- a.params[i,j]
					solved_a[a.params[i,j]] = a.values[i,j]
				}
				else if(!is.na(e_t)){
					L[i+(j-1)*n][[1]] <- as.character(e_t)
				}
				else{ #this paramter needs to be rewritten in the form of pars
					for(par_ in names(par_list)){
						L[i+(j-1)*n][[1]] <- gsub(par_list[par_], par_, L[i+(j-1)*n][[1]])
					}
				}
			}
		}
	}
	
	#browser()
	#browser()
	# exponential tranformation of D
	D <- sapply(D, function(term){ if(is.na(evaluateExpression(term[[1]]))){term[[1]]<- paste0('exp(', term[[1]], ')')} else{ term[[1]]<- paste0('(', term[[1]], ')')}})
	
	#L*D
	for(i in 1:n){
		for(j in 1:n){
			e_l <- evaluateExpression(L[i+(j-1)*n][[1]])
			e_d <- evaluateExpression(D[j][[1]])
			if(!is.na(e_l) && e_l == 0){
				term <- term
			}
			else if(!is.na(e_d) && e_d == 0){
				term <- term
			}
			else{
				temp[i+(j-1)*n][[1]] <- paste0('(', L[i+(j-1)*n][[1]], ')*(', D[j][[1]],')')
			}
		}
	}
	
	#D*t(L)
	for(i in 1:n){
		for(j in 1:n){
			term <- paste0('(', temp[i][[1]],')*(', L[j][[1]],')')
			if(n>1){
				for(k in 2:n){
					e_l <- evaluateExpression(temp[i+(k-1)*n][[1]])
					e_d <- evaluateExpression(L[j+(k-1)*n][[1]])
					if(!is.na(e_l) && e_l == 0){
						term <- term
					}
					else if(!is.na(e_d) && e_d == 0){
						term <- term
					}
					else{
						term <- paste0('(', term, ')+(', temp[i+(k-1)*n][[1]],')*(', L[j+(k-1)*n][[1]],')')
					}
				}
			}
			ret[i+(j-1)*n][[1]]<- term
		}
	}
	
	
	L = matrix(sapply(L, function(x){as.list(as.formula(paste0('x ~ ',x)))[[3]]}), nrow=nrow(a), ncol=ncol(a))
	D = matrix(sapply(D, function(x){as.list(as.formula(paste0('x ~ ',x)))[[3]]}), nrow=nrow(a), ncol=ncol(a))
	for(i in 1:nrow(D)){
	  for(j in 1:ncol(D)){
	    if(i != j)
	      D[i,j] = expression(0)
	  }
	}
	
	# examine if there is any discrepancy between a.values and expressions in the multiplication of LDL' 
	#if(any(mapply(function(x, value){
	#					e_l <- evaluateExpression(x)
	#					return(!is.na(e_l) && e_l != value)}, ret, as.list(a.values)))){
	#	stop('unreasonable setting of dynrModel@random.values.inicov')
	#}
	
	#
	ret = matrix(mapply(function(term, ret, value){
			if(term == 0)
				as.list(as.formula(paste0('x ~ ',value)))[[3]]
			else
				as.list(as.formula(paste0('x ~ ',ret)))[[3]]}, as.list(a.params), ret, as.list(a.values)), nrow=nrow(a), ncol=ncol(a))
				
	#ret = matrix(sapply(ret, function(x){as.list(as.formula(paste0('x ~ ',x)))[[3]]}), nrow=nrow(a), ncol=ncol(a))
	
	
	
	#return(list(L=L, D=D, r=ret))
	
	return(list(ldl=ret, pars=names(par_list), L=L, D=D))

}

unList <- function(ret){
	ret = matrix(sapply(ret, function(term){term[[1]]}), nrow(ret), ncol(ret))
	return(ret)
}


# evaluate whether the expression can be solved after substitutes the values in par.values
# a: the expression in the form of character
# par.values: the values to be substituted in (e.g. list(par0 = 1, par0 = 2)
# return NA if the expression can not be solved, return the value if it is solved
evaluateExpression <- function(a, par.values = list()){
	if(is.numeric(a) == TRUE)
		return(a)
	if(is.character(a) == FALSE)
		stop('evaluateExpression: term should be in the type of character')
	
	value=NA
	tryCatch(
		{value=eval(as.list(as.formula(paste0('x ~ ',a)))[[3]], par.values)},
		error=function(e){value=NA})
	
	return(value)
}

##' A sub-function called by solveStartLDL
##' 
##' The goal is to solve equation call = b 
##' There may be multiple variables in call, but there must be at most one unknown variables in call
##'
##' @param  call the LHS expression given in the type of call
##' @param  b the RHS value
##' @param known.vars the list of values of already known variables (e.g., list(par1= 1, par2 = 2))
##' @return the concatenation list of known.vars and values of new-solved variable (e.g., list(par1= 1, par2 = 2, par2=3))
solveOneElementEquation <- function(call, b, known.vars){
  #substitute the known values in known.vars in
  expression <- do.call('substitute',list(call, known.vars))
  
  #find the unknown var is in exp
  mode <- NA
  str <- deparse(expression, ,width.cutoff = 500L)
  pos <- regexpr('par[0-9]+', str)
  if (pos != -1){
    len <- attr(regexpr('par[0-9]+', str), "match.length")
    unknown.variable <- substr(str, pos, pos + len -1)
    mode <- 'normal'
	# constant in LHS in character type
	const.str <- gsub(unknown.variable, '0', str)
	const <- evaluateExpression(const.str)
  }
  else{
    return(known.vars)}
  
  # further examine whehter it is a exp(par), if yes, set mode to be 'exp'
  pos <- regexpr(paste0('exp\\(', unknown.variable, '\\)'), str)
  if (pos != -1){
    len <- attr(regexpr(paste0('exp\\(', unknown.variable, '\\)'), str), "match.length")
	#unknown.variable <- substr(str, pos, pos + len -1)
	mode <- 'exp'
	const.str <- gsub(paste0('exp\\(', unknown.variable, '\\)'), '0', str)
	const <- evaluateExpression(const.str)
  }
  
  # now we get the form a * par + const = b (or a * exp(par) + const = b)
  # the following lines obtain a
  if(mode == 'exp'){
    a.str <- gsub(paste0('exp\\(', unknown.variable, '\\)'), '1', paste0(str , '-', const))
  } else if(mode == 'normal'){
    a.str <- gsub(unknown.variable, '1', paste0(str , '-', const))
  }
  a <- evaluateExpression(a.str)
  par <- (b - const)/a 
  
  
  #log transformation if the unknown var is in the form of exp(par)
  if (mode == 'exp'){
    par <- log(par)
  }

  # add it to the list of known variables
  known.vars[unknown.variable] = par
  
  # return a new list of known vars
  return (known.vars)
}

# Given the expression and values of L%*%D%*%t(L) in the contrained scale, obtain the numerical values of the parameeters in the unconstrained scale.
##' Solve for numerical values of the unconstrained parameters (e.g.,par0-par2),
##' given starting values for the random effect covariance matrix (i.e. model@random.values.inicov)  in the contrained scale
##' 
##' @param ldl the expression of the random effect covariance matrix (which can be obtained by symbolicLDLDecomposition)
##' @param values.ldl the starting values for the random effect covariance matrix in model@random.values.inicov
##' @return the list of numerical values of the unconstrained parameters (e.g., list(par0= 1, par1 =2, par2 =3))
##' @seealso \code{\link{symbolicLDLDecomposition}}
solveStartLDL <- function(ldl, values.ldl){
  known.vars <- list()
  #browser()
  
  if (min(diag(MASS::ginv(values.ldl)))<0){
    stop("Starting values for entries in random.values.inicov
        and values.inicov generate a non-positive 
         definite covariance matrix. Please check
         the starting values for these entries.")
  }
  
  for(i in 1:nrow(ldl)){
    for(j in 1:nrow(ldl)){
       known.vars <- solveOneElementEquation(ldl[i,j][[1]], values.ldl[i,j], known.vars) 
    }
  }
  return(known.vars)
}


# retrieve all variables from RHS of formulas
# input example:
# formula = list(x ~ dx,
#               dx ~ eta_i * x + zeta*dx)
# output example:
# [1] "dx"    "eta_i" "x"     "zeta" 
extractVariablesfromFormula<- function(formula){
    variables <-c()
    for (i in 1:length(formula))
        variables <- c(variables, extractVariablesfromSingleFormula(formula[[i]]))
	
	# exclude the functions cos, sin, log, exp
	variables <- setdiff(unique(variables), c('exp', 'log', 'sin', 'cos'))
	
    return(variables)
}

# retrieve all variables from RHS of a single formula
# variables is the combination of A-Z, a-z, 0-9, and '_' (at least one A-Z/a-z)
# \\W: all characters excluding A-z, a-z, and '_'
# \\w: A-z, a-z, and '_'
extractVariablesfromSingleFormula <- function(formula){
    #browser()
    rhs <- as.character(formula)[[3]]
    #a <- unlist(stringr::str_split(rhs, '[^A-z0-9_.]+'))
	a <- unlist(stringr::str_split(rhs, '[^[:alnum:]_.]+'))
    variables <- a[grep('[\\w.]*[a-zA-Z]+[\\w.]*', a)]
    return(unique(variables))
}

# The helper function that vectorize a matrix
# Input: a matrix in column-major order (byrow=FALSE) or in row-major order (byrow=TRUe)
# output: The vectorization result
# Example: x = matrix(c(1:9), nrow=3) 
# x =   [ 1    4    7
#         2    5    8
#         3    6    9]
# vectorizeMatrix(x, byrow=FALSE) Result: [1] 1 2 3 4 5 6 7 8 9
# vectorizeMatrix(x, byrow=TRUE)  Result: [1] 1 4 7 2 5 8 3 6 9
vectorizeMatrix <- function(matrix, byrow= FALSE){
    if(byrow == FALSE)
        return(c(matrix))
    else
        return(c(t(matrix)))
}


symbolicLDLDecomposition2 <- function(a.params, a.values){
  #browser()
  if(!is.matrix(a.params) || nrow(a.params) != ncol(a.params))
    stop("The variable 'a.params' must be a square matrix")
  
  if(!isSymmetric(a.params))
    stop("The variable 'a.params' must be a symmetric matrix")
  
  
  
  # repalce elements equal to NA or 'NA' as 0
  a.params[ is.na(a.params) ] = 0
  a.params[ a.params == 'NA'] = 0
  a.params[ a.params == 'fixed'] = 0
  a.params[ a.params == '0'] = 0
  # a.params stores the input (with severl inputs processing as above)
  a = a.params
  
  
  # repalce elements 'fixed' w/ the corresponding values in a.values
  a = matrix(mapply(function(a, v){if(a == 0) a = v else a = a}, as.list(a), as.list(a.values)), nrow = nrow(a))
  
  n <- nrow(a)
  L <- rep(list(0), n*n)
  D <- rep(list(1), n)
  ret <- rep(list(0), n*n)
  temp <- rep(list(0), n*n)
  
  if(is.character(a[1,1])){
    # the input variable, in the type of "character" 
    a = matrix(sapply(a, function(x){list(x)}), nrow=nrow(a), ncol=ncol(a))
  } else {
    # the input is expression, in the type of "call"
    a = matrix(sapply(a, function(x){return(deparse(unlist(x)))}), nrow=nrow(a), ncol=ncol(a))
  }
  # the element in a should be a list with a element, which contains the expression in the form of character
  
  for(i in 1:n){
    #D[i,i]
    #if(!is.character(a[i,i][[1]]))
    term <- a[i,i][[1]]
    if(i > 1){
      for(k in 1:(i-1)){
        #term <- paste(term, '-', as.character(eval(L[i,k][[1]]*L[i,k][[1]]*D[k,k][[1]])))
        e_l <- evaluateExpression(L[i+(k-1)*n][[1]])
        e_d <- evaluateExpression(D[k][[1]])
        if(!is.na(e_l)){
          if(e_l == 0)
            term <- term
          else
            term <- as.character(paste0(term, '-(', e_l, ')^2*(', D[k][[1]], ')'))
        }
        else if(!is.na(e_d)){
          if(e_d == 0)
            term <- term
          else
            term <- as.character(paste0(term, '-(', L[i+(k-1)*n][[1]], ')^2*(', e_d, ')'))
        }
        else{
          term <- as.character(paste0(term, '-(', L[i+(k-1)*n][[1]], ')^2*(', D[k][[1]], ')'))
        }
      }
    }
    #D[i][[1]] <- as.list(as.formula(paste0('x ~ ',term)))[[3]]
    D[i][[1]] <- term
    if(i < n){
      for(j in (i+1):n){
        term <- a[j,i][[1]]
        
        if(i > 1){
          for(k in 1:(i-1)){
            e_l <- evaluateExpression(L[i+(k-1)*n][[1]])
            e_d <- evaluateExpression(D[k][[1]])
            if(!is.na(e_l)){
              if(e_l == 0)
                term <- term
              else
                term <- as.character(paste0(term, '-(', L[j+(k-1)*n][[1]], ')*(', e_l, ')*(', D[k][[1]], ')'))
            }
            else if(!is.na(e_d) && e_d == 0){
              if(e_d == 0)
                term <- term
              else
                term <- as.character(paste0(term, '-(', L[j+(k-1)*n][[1]], ')*(', L[i+(k-1)*n][[1]], ')*(', e_d, ')'))
            }
            else{
              term <- as.character(paste0(term, '-(', L[j+(k-1)*n][[1]], ')*(', L[i+(k-1)*n][[1]], ')*(', D[k][[1]], ')'))
            }
          }
        }
        #L[j+(i-1)*n][[1]] <- as.list(as.formula(paste0('x ~ ',term)))[[3]]
        e_t <- evaluateExpression(term)
        if(is.na(e_t) || e_t != 0){
          term <- paste0('(', term, ')/(', D[i][[1]],')')
        }
        L[j+(i-1)*n][[1]] <- term
      }
    }
    L[i+(i-1)*n][[1]] <- "1"
    
  }
  
  
  # examine whether we need to have a pari to estimate i
  # par_list <- list()
  # solved_a <- list() # for each solved element, we assign an arbitrary value (here is its value) for it temporarily
  # for(i in 1:n){
  #   
  #   e_t <- evaluateExpression(D[i][[1]], solved_a) 
  #   
  #   if(is.na(e_t) && a.params[i,i] != 0){
  #     new_par <- paste0('par', length(par_list), '_star')
  #     D[i][[1]] <- new_par
  #     #par_list <- c(par_list, new_par)
  #     par_list[new_par] <- a.params[i,i]
  #     solved_a[a.params[i,i]] = a.values[i,i]
  #   }
  #   else if(!is.na(e_t)){
  #     D[i][[1]] <- as.character(e_t)
  #   }
  #   else{ #this paramter needs to be rewritten in the form of pars
  #     for(par_ in names(par_list)){
  #       D[i][[1]] <- gsub(par_list[par_], par_, D[i][[1]])
  #     }
  #   }
  #   
  #   if(i > 1){
  #     for(j in 1:(i-1)){
  #       e_t <- evaluateExpression(L[i+(j-1)*n][[1]], solved_a) 
  #       
  #       if(is.na(e_t) && a.params[i,j] != 0){
  #         new_par <- paste0('par', length(par_list), '_star')
  #         L[i+(j-1)*n][[1]] <- new_par
  #         #par_list <- c(par_list, new_par)
  #         par_list[new_par] <- a.params[i,j]
  #         solved_a[a.params[i,j]] = a.values[i,j]
  #       }
  #       else if(!is.na(e_t)){
  #         L[i+(j-1)*n][[1]] <- as.character(e_t)
  #       }
  #       else{ #this paramter needs to be rewritten in the form of pars
  #         for(par_ in names(par_list)){
  #           L[i+(j-1)*n][[1]] <- gsub(par_list[par_], par_, L[i+(j-1)*n][[1]])
  #         }
  #       }
  #     }
  #   }
  # }
  
  #browser()
  #browser()
  # exponential tranformation of D
  D <- sapply(D, function(term){ if(is.na(evaluateExpression(term[[1]]))){term[[1]]<- paste0('exp(', term[[1]], ')')} else{ term[[1]]<- paste0('(', term[[1]], ')')}})
  
  #L*D
  for(i in 1:n){
    for(j in 1:n){
      e_l <- evaluateExpression(L[i+(j-1)*n][[1]])
      e_d <- evaluateExpression(D[j][[1]])
      if(!is.na(e_l) && e_l == 0){
        term <- term
      }
      else if(!is.na(e_d) && e_d == 0){
        term <- term
      }
      else{
		if(!is.na(e_l) && !is.na(e_d) && e_l ==  1 && e_d == 1){
		  temp[i+(j-1)*n][[1]] <- '(1)'
		}
		else if (!is.na(e_l) && e_l ==  1){
		  temp[i+(j-1)*n][[1]] <- paste0('(', D[j][[1]],')')
		}
		else if (!is.na(e_d) && e_d ==  1){
		  temp[i+(j-1)*n][[1]] <- paste0('(', L[i+(j-1)*n][[1]], ')')
		}
		else{
		  temp[i+(j-1)*n][[1]] <- paste0('(', L[i+(j-1)*n][[1]], ')*(', D[j][[1]],')')
		}
      }
    }
  }
  
  #D*t(L)
  for(i in 1:n){
    for(j in 1:n){
	  e_t <- evaluateExpression(temp[i][[1]])
	  e_l <- evaluateExpression(L[j][[1]])
	  if(!is.na(e_l) && !is.na(e_t) && e_l ==  1 && e_t == 1){
	    term <- '(1)'
	  }
	  else if (!is.na(e_l) && e_l ==  1){
	    term <- paste0('(', temp[i][[1]], ')')
	  }
	  else if (!is.na(e_t) && e_t ==  1){
	    term <- paste0('(', L[j][[1]],')')
	  }
	  else{
        term <- paste0('(', temp[i][[1]],')*(', L[j][[1]],')')
	  }
      if(n>1){
        for(k in 2:n){
          e_t <- evaluateExpression(temp[i+(k-1)*n][[1]])
          e_l <- evaluateExpression(L[j+(k-1)*n][[1]])
          if(!is.na(e_l) && e_l == 0){
            term <- term
          }
          else if(!is.na(e_t) && e_t == 0){
            term <- term
          }
          else{
            #term <- paste0('(', term, ')+(', temp[i+(k-1)*n][[1]],')*(', L[j+(k-1)*n][[1]],')')
			if(!is.na(e_l) && !is.na(e_t) && e_l ==  1 && e_t == 1){
			  term <- paste0('(',term, ')+ (1)')
			}
			else if (!is.na(e_l) && e_l ==  1){
			  term <- paste0('(', term, ')+(', temp[i+(k-1)*n][[1]],')')
			}
			else if (!is.na(e_t) && e_t ==  1){
			  term <- paste0('(', term, ')+(', L[j+(k-1)*n][[1]],')')
			}
			else{
			  term <- paste0('(', term, ')+(', temp[i+(k-1)*n][[1]],')*(', L[j+(k-1)*n][[1]],')')
			}
          }
        }
      }
      ret[i+(j-1)*n][[1]]<- term
    }
  }
  
  
  L = matrix(sapply(L, function(x){as.list(as.formula(paste0('x ~ ',x)))[[3]]}), nrow=nrow(a), ncol=ncol(a))
  D = matrix(sapply(D, function(x){as.list(as.formula(paste0('x ~ ',x)))[[3]]}), nrow=nrow(a), ncol=ncol(a))
  for(i in 1:nrow(D)){
    for(j in 1:ncol(D)){
      if(i != j)
        D[i,j] = expression(0)
    }
  }
  
  # examine if there is any discrepancy between a.values and expressions in the multiplication of LDL' 
  #if(any(mapply(function(x, value){
  #					e_l <- evaluateExpression(x)
  #					return(!is.na(e_l) && e_l != value)}, ret, as.list(a.values)))){
  #	stop('unreasonable setting of dynrModel@random.values.inicov')
  #}
  
  #
  
  ret = matrix(mapply(function(term, ret, value){
    if(term == 0)
      as.list(as.formula(paste0('x ~ ',value)))[[3]]
    else
      as.list(as.formula(paste0('x ~ ',ret)))[[3]]}, as.list(a.params), ret, as.list(a.values)), nrow=nrow(a), ncol=ncol(a))
  
  #ret = matrix(sapply(ret, function(x){as.list(as.formula(paste0('x ~ ',x)))[[3]]}), nrow=nrow(a), ncol=ncol(a))
  
  
  
  return(list(L=L, D=D, ldl=ret))
  
  #return(list(ldl=ret, pars=names(par_list), L=L, D=D))
  
}

solveStartLDLSymbolically <- function(ldl, params.ldl){
  known.vars <- list()
  #browser()
  
  # if (min(diag(MASS::ginv(values.ldl)))<0){
    # stop("Starting values for entries in random.values.inicov
        # and values.inicov generate a non-positive 
         # definite covariance matrix. Please check
         # the starting values for these entries.")
  # }
  
  for(i in 1:nrow(ldl)){
    for(j in 1:nrow(ldl)){
	   #equ <- paste(deparse(ldl[i,j][[1]], width.cutoff = 256L), '-', values.ldl[i,j], '== 0')
	   equ <- paste(deparse(ldl[i,j][[1]], width.cutoff = 256L), '-', params.ldl[i,j], '== 0')
	   equ <- gsub('exp', 'Exp', equ)
	   equ <- gsub('log', 'Ln', equ)
	   equ <- gsub(' _', '_', equ)
	   if(i == 1 && j == 1)
         cmd <- equ %>% y_fn("Solve", "sigma2beta")
	   else if(i == 1 && j == 2)
         cmd <- equ %>% y_fn("Solve", "sigma2beta")
	   else if(i == 2 && j == 2)
         cmd <- equ %>% y_fn("Solve", "sigma2beta")
	   else
	     next
	   sol <- yac_str(cmd)
	   sol<- substr(sol, 2, nchar(sol)-1)
	   sol <- gsub('Exp', 'exp', sol)
	   sol <- gsub('Ln', 'log', sol)
	   sol <- gsub(' _', '_', sol)
	   lhs_sol <- strsplit(sol, "==")[[1]][1]
	   rhs_sol <- strsplit(sol, "==")[[1]][2]
	   #browser()
	   known.vars[lhs_sol] <- list(as.formula(paste0('x ~ ',rhs_sol))) 
    }
  }
  return(known.vars)
}