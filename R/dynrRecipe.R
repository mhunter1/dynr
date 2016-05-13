# A Recipe is a helper function that takes user input
#  and produces a C function definition that can be
#  compiled by dynrFuncaddress.

# TODO add check that there is no cross loading on one of the ID variables

# TODO add assignment of first loading to 1 in C
# TODO rename function to more like fastLoadings
# TODO create new mid-level function (for which this is a wrapper) that takes three
#  matrices (parameter staring values, free-fixed, and parameter numbers)

# TODO output the starting values for the parameters 
# TODO add documentation


#------------------------------------------------------------------------------
# Create dynrRecipe object class

##DONE##
setClass(Class =  "dynrRecipe",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnames = "character"
         )
)

##DONE##
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

##DONE##
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

##DONE##
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

##DONE##
setClass(Class = "dynrInitial",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnames = "character",
           values.inistate = "matrix",
           params.inistate = "matrix",
           values.inicov = "matrix",
           values.inicov.inv.ldl = "matrix",
           params.inicov = "matrix",
           values.regimep = "matrix",
           params.regimep = "matrix"),
         contains = "dynrRecipe"
)


##DONE##
setClass(Class = "dynrNoise",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnames = "character",
           values.latent = "matrix",
           params.latent = "matrix",
           values.observed = "matrix",
           params.observed = "matrix",
           values.latent.inv.ldl = "matrix",
           values.observed.inv.ldl = "matrix"),#TODO we should emphasize that either the full noise covariance structures should be freed or the diagonals because we are to apply the ldl trans  
         contains = "dynrRecipe"
)

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

#setGeneric("printex", function(object, observed, latent, covariates, show=TRUE) { 
#  return(standardGeneric("printex")) 
#})
setGeneric("printex", function(object, show=TRUE) { 
	return(standardGeneric("printex")) 
})



setMethod("printex", "dynrMeasurement",
          function(object, show=TRUE){
            meas_loadings=lapply((object)$values.load,.xtableMatrix,show)
            meas_int=lapply((object)$values.int,.xtableMatrix,show)
            meas_exo=lapply((object)$values.exo,.xtableMatrix,show)
            meas_exo.names=.xtableMatrix(matrix((object)$exo.names,ncol=1),show)
            meas_list = list(meas_loadings=meas_loadings,meas_int=meas_int,
                             meas_exo=meas_exo,meas_exo.names=meas_exo.names)
            return(invisible(meas_list))
          }
)


setMethod("printex", "dynrDynamicsMatrix",
          function(object, show=TRUE){
            dyn_tran=lapply((object)$values.dyn,.xtableMatrix,show)   
            dyn_int=lapply((object)$values.int,.xtableMatrix,show)   
            dyn_exo=lapply((object)$values.exo,.xtableMatrix,show)   
            dyn_exo.names=.xtableMatrix(matrix((object)$covariates,ncol=1),show)
            return(invisible(list(dyn_tran = dyn_tran, dyn_int = dyn_int, dyn_exo = dyn_exo, dyn_exo.names = dyn_exo.names) ))
          }
)

setMethod("printex", "dynrRegimes",
          function(object, show=TRUE){
            lG <- ifelse(nrow(object$values) != 0, .xtableMatrix(object$values, show), "")
            return(invisible(list(regimes=lG)))
          }
)


setMethod("printex", "dynrInitial",
          function(object, show=TRUE){
            lx0 <- .xtableMatrix(object$values.inistate, show)
            lP0 <- .xtableMatrix(object$values.inicov, show)
            lr0 <- .xtableMatrix(object$values.regimep, show)
            return(invisible(list(initial.state=lx0, initial.covariance=lP0, initial.probability=lr0)))
          }
)


setMethod("printex", "dynrNoise",
          function(object, show=TRUE){
            lQ <- .xtableMatrix(object$values.latent, show)
            lR <- .xtableMatrix(object$values.observed, show)
            return(invisible(list(dynamic.noise=lQ, measurement.noise=lR)))
          }
)

setMethod("printex", "dynrDynamicsFormula",
          function(object, show=TRUE){
            dyn=lapply(object$formula,dynfmltex,object$isContinuousTime)
            return(invisible(dyn))
          }
)



#setGeneric("printmath", 
#           function(object, observed, latent, covariates,show=TRUE) { 
#  return(standardGeneric("printmath")) 
#})



dynfmltex<-function(eqregime,isContinuousTime){
  eq.char=lapply(eqregime, as.character)
  str.left=sapply(eq.char,"[",2)
  str.right=sapply(eq.char,"[",3)
  neq=length(eqregime)
  mulpatn<-"([[:print:]]*)"
  sigpatn<-"([0-9A-Za-z ^*]*)"
  for (j in 1:neq){
    if (!isContinuousTime){
      for (i in 1:neq){
        str.right[j]=gsub(str.left[i],paste0(str.left[i],"_{t-1}"),str.right[j])
      }
    }
    str.right[j]=gsub(paste0("\\(",mulpatn,"\\)/\\(",mulpatn,"\\)"),"\\\\frac{\\1}{\\2}",str.right[j])
    str.right[j]=gsub(paste0("\\(",mulpatn,"\\)/",sigpatn),"\\\\frac{\\1}{\\2}",str.right[j])
    str.right[j]=gsub(paste0(sigpatn,"/\\(",mulpatn,"\\)"),"\\\\frac{\\1}{\\2}",str.right[j])
    str.right[j]=gsub(paste0(sigpatn,"/",sigpatn),"\\\\frac{\\1}{\\2}",str.right[j])
    str.right[j]=gsub("exp","\\\\exp",str.right[j])
    str.right[j]=gsub("log","\\\\log",str.right[j])
    str.right[j]=gsub("\\*","\\\\times",str.right[j])
  }
  
  for (j in 1:neq){
    str.left[j]=gsub(paste0(sigpatn),ifelse(isContinuousTime,"d(\\1)","\\1_t"),str.left[j])
  }
  
  return(list(left=str.left,right=str.right))
}



dynfm_math<-function(eqregime,isContinuousTime){
  eq.char=lapply(eqregime, as.character)
  str.left=sapply(eq.char,"[",2)
  str.right=sapply(eq.char,"[",3)
  neq=length(eqregime)
  pretty.list = vector("list", neq) 
  #mulpatn<-"([[:print:]]*)"
  #sigpatn<-"([0-9A-Za-z ^*]*)"
  for (j in 1:neq){
    if (isContinuousTime){
      str.left[j]=paste0("frac(italic(d)",str.left[j],", italic(dt))")
    }else{
      str.left[j]=paste0(str.left[j],"_{t}")      
    }
    str.right[j]=gsub(str.left[j],paste0(str.left[j],"_{t-1}"),str.right[j])
    pretty.list[[j]] = c(paste0(str.left[j]," = ", str.right[j]))
    }
  return(pretty.list)
}


#setMethod("printmath", "dynrDynamicsFormula",
#          function(object, observed, latent, covariates, show=TRUE){
#            dyn=lapply(object$formula,dynfm_math,object$isContinuousTime)
#            return(invisible(dyn))
#          }
#)

#TODO: Incorporate process noise into dynamic formula
.formulatoTex <- function(object,model){
  nregime = length(object)
  exp1 = NULL
  for (r in 1:nregime){
    ne = length(object[[r]])
    for (j in 1:ne){
      if ((model$dynamics)$isContinuousTime){
        exp1 = c(exp1,paste0("$\\frac{",paste(object[[r]]$left[j]),"}{dt} = $", 
                             paste0("$",paste(object[[r]]$right[j]),"$")))
      }else{
        exp1 = c(exp1,paste0("$",paste(object[[r]]$left[j])," = $", 
                             paste0("$",paste(object[[r]]$right[j]),"$")))
      }
    }
  }
  return(invisible(exp1))
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
    pattern=gsub("\\{","\\\\\\{",paramnames)
    pattern=gsub("\\}","\\\\\\}",pattern)
    pattern=gsub("\\\\","\\\\\\\\",pattern)
    
    for (i in 1:length(paramnames)){
      
      string<-gsub(pattern[i],paste0("param[",match(paramnames[i], names, nomatch=0)-1,"]"), string, perl = TRUE)
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
  string<-gsub(names[i],paramtoPlot[i], string)
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
		object@params.inistate <- .exchangeNamesAndNumbers(object$params.inistate, names)
		object@params.inicov <- .exchangeNamesAndNumbers(object$params.inicov, names)
		object@params.regimep <- .exchangeNamesAndNumbers(object$params.regimep, names)
		return(object)
	}
)


setMethod("paramName2Number", "dynrNoise",
	function(object, names){
		object@params.latent <- .exchangeNamesAndNumbers(object$params.latent, names)
		object@params.observed <- .exchangeNamesAndNumbers(object$params.observed, names)
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

setGeneric("writeCcode", function(object, show=TRUE) { 
	return(standardGeneric("writeCcode")) 
})



setMethod("writeCcode", "dynrMeasurement",
	function(object){
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
		if(hasCovariates){ret <- paste(ret, createGslMatrix(nrow(params.exo[[1]]), ncol(params.exo[[1]]), "Bmatrix"), sep="\n\t")}#create covariates matrix
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
			if(hasCovariates){paste(ret, setGslMatrixElements(values.exo[[1]], params.exo[[1]], "Bmatrix"), sep="\n")}#set covariates matrix
			if(hasIntercepts){paste(ret, setGslVectorElements(values.int[[1]], params.int[[1]], "intVector"), sep="\n")}#set intercepts vector
		}
		ret <- paste(ret, "\n\tgsl_blas_dgemv(CblasNoTrans, 1.0, Ht, eta, 0.0, y);\n")
		if(hasCovariates){ret <- paste(ret, "\n\tgsl_blas_dgemv(CblasNoTrans, 1.0, Bmatrix, co_variate, 1.0, y);\n", destroyGslMatrix("Bmatrix"))}#multiply, add, and destroy covariates matrix
		if(hasIntercepts){ret <- paste(ret, "\n\tgsl_vector_add(y, intVector);\n", destroyGslVector("intVector"))}#add and destroy intercepts vector
		ret <- paste(ret, "\n}\n\n")
		object@c.string <- ret
		return(object)
	}
)


# not sure what to do here yet
setMethod("writeCcode", "dynrDynamicsFormula",
	function(object){
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
	  
	  #TODO add covariate
	  
	  if (object@isContinuousTime){
	    #function_dx_dt
	    ret="void function_dx_dt(double t, size_t regime, const gsl_vector *x, double *param, size_t n_param, const gsl_vector *co_variate, gsl_vector *F_dx_dt){"
	    
	    if (nregime>1){
	      ret=paste(ret,"switch (regime) {",sep="\n\t")
	      for (r in 1:nregime){
	        ret=paste(ret,paste0("\tcase ",r-1,":"),sep="\n\t")
	        for (i in 1:n[r]){
	          for (j in 1:length(lhs[[r]])){
	            rhs[[r]][[i]]=gsub(lhs[[r]][[j]],paste0("gsl_vector_get(x,",j-1,")"),rhs[[r]][[i]])
	          }
	          ret=paste(ret,paste0("\tgsl_vector_set(F_dx_dt,",i-1,",",rhs[[r]][[i]],");"),sep="\n\t")    
	        }
	        ret=paste(ret,paste0("break;\n"),sep="\n\t")
	        
	      }
	      ret=paste(ret,paste0("\t}"),sep="\n\t")
	      
	    }else{
	      for (i in 1:n){
	        for (j in 1:length(lhs[[1]])){
	          rhs[[1]][[i]]=gsub(lhs[[1]][[j]],paste0("gsl_vector_get(x,",j-1,")"),rhs[[1]][[i]])
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
	            rhsj[[r]][[i]]=gsub(lhs[[r]][[j]],paste0("param[NUM_PARAM+",j-1,"]"),rhsj[[r]][[i]])
	          }
	          
	          ret=paste(ret,paste0("\tgsl_matrix_set(F_dx_dt_dx,",which(lhs[[r]]==row[[r]][[i]])-1,",",which(lhs[[r]]==col[[r]][[i]])-1,",",rhsj[[r]][[i]],");"),sep="\n\t")    
	        }
	        ret=paste(ret,paste0("break;\n"),sep="\n\t")
	        
	      }
	      ret=paste(ret,paste0("\t}"),sep="\n\t")
	      
	    }else{
	      for (i in 1:length(jacob[[1]])){
	        for (j in 1:length(lhs[[1]])){
	          rhsj[[1]][[i]]=gsub(lhs[[1]][[j]],paste0("param[NUM_PARAM+",j-1,"]"),rhsj[[1]][[i]])
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
	            rhs[[r]][[i]]=gsub(lhs[[r]][[j]],paste0("gsl_vector_get(xstart,",j-1,")"),rhs[[r]][[i]])
	          }
	          ret=paste(ret,paste0("\tgsl_vector_set(x_tend,",i-1,",",rhs[[r]][[i]],");"),sep="\n\t")    
	        }
	        ret=paste(ret,paste0("break;\n"),sep="\n\t")
	        
	      }
	      ret=paste(ret,paste0("\t}"),sep="\n\t")
	      
	    }else{
	      for (i in 1:n){
	        for (j in 1:length(lhs[[1]])){
	          rhs[[1]][[i]]=gsub(lhs[[1]][[j]],paste0("gsl_vector_get(xstart,",j-1,")"),rhs[[1]][[i]])
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
	            rhsj[[r]][[i]]=gsub(lhs[[r]][[j]],paste0("gsl_vector_get(xstart,",j-1,")"),rhsj[[r]][[i]])
	          }
	          
	          ret=paste(ret,paste0("\tgsl_matrix_set(Jx,",which(lhs[[r]]==row[[r]][[i]])-1,",",which(lhs[[r]]==col[[r]][[i]])-1,",",rhsj[[r]][[i]],");"),sep="\n\t")    
	        }
	        ret=paste(ret,paste0("break;\n"),sep="\n\t")
	        
	      }
	      ret=paste(ret,paste0("\t}"),sep="\n\t")
	      
	    }else{
	      for (i in 1:length(jacob[[1]])){
	        for (j in 1:length(lhs[[1]])){
	          rhsj[[1]][[i]]=gsub(lhs[[1]][[j]],paste0("gsl_vector_get(xstart,",j-1,")"),rhsj[[1]][[i]])
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
	function(object){
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
				ifelse(hasCovariates, createGslMatrix(nrow(params.exo[[1]]), ncol(params.exo[[1]]), "Bmatrix"), ""),
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
				ifelse(hasCovariates, blasMV(FALSE, "1.0", "Bmatrix", "co_variate", "1.0", outName), ""),
				ifelse(hasIntercepts, paste0("\tgsl_vector_add(", outName, ", intVector);\n"), ""),
				destroyGslMatrix("Amatrix"),
				ifelse(hasCovariates, destroyGslMatrix("Bmatrix"), ""),
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
					createGslMatrix(nrow(params.exo[[1]]), ncol(params.exo[[1]]), "Bmatrix"),
					setGslMatrixElements(values=values.exo[[1]], params=params.exo[[1]], name="Bmatrix"),
					blasMV(FALSE, "1.0", "Bmatrix", "co_variate", "1.0", outName),
					destroyGslMatrix("Bmatrix"),
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
	function(object){
		values <- object$values
		params <- object$params
		covariates <- object$covariates
		numCovariates <- length(covariates)
		numRegimes <- nrow(values)
		
		#Restructure values matrix for row-wise
		if(nrow(values)!=0 && nrow(params)!=0){
			values <- matrix(t(values), nrow=numRegimes*numRegimes, ncol=numCovariates+1, byrow=TRUE)
			params <- matrix(t(params), nrow=numRegimes*numRegimes, ncol=numCovariates+1, byrow=TRUE)
			rowBeginSeq <- seq(1, nrow(values), by=numRegimes)
			rowEndSeq <- seq(numRegimes, nrow(values), by=numRegimes)
		}
		
		ret <- "void function_regime_switch(size_t t, size_t type, double *param, const gsl_vector *co_variate, gsl_matrix *regime_switch_mat){"
		if(prod(nrow(values), nrow(params), length(covariates)) != 0){
			ret <- paste(ret,
				createGslMatrix(numRegimes, numCovariates, "Gmatrix"),
				createGslVector(numRegimes, "Pvector"),
				createGslVector(numRegimes, "Presult"),
				sep="\n")
			for(reg in 1L:numRegimes){
				selRows <- rowBeginSeq[reg]:rowEndSeq[reg]
				ret <- paste(ret,
					setGslVectorElements(values=values[selRows, 1, drop=FALSE], params=params[selRows, 1, drop=FALSE], name="Pvector"),
					setGslMatrixElements(values=values[selRows, -1, drop=FALSE], params=params[selRows, -1, drop=FALSE], name="Gmatrix"),
					blasMV(FALSE, "1.0", "Gmatrix", "co_variate", "1.0", "Pvector"),
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
				sep="\n")
		} else if(nrow(values) != 0 && nrow(params) != 0 && length(covariates) == 0){
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
	function(object){
		values.inistate <- object$values.inistate
		params.inistate <- object$params.inistate
		values.inicov <- object@values.inicov.inv.ldl
		params.inicov <- object$params.inicov
		values.regimep <- object$values.regimep
		params.regimep <- object$params.regimep
		ret <- "void function_initial_condition(double *param, gsl_vector **co_variate, gsl_vector *pr_0, gsl_vector **eta_0, gsl_matrix **error_cov_0){\n"
		ret <- paste(ret, setGslVectorElements(values.regimep,params.regimep, "pr_0"), sep="\n")
		ret <- paste0(ret,"\tsize_t num_regime=pr_0->size;\n\tsize_t dim_latent_var=error_cov_0[0]->size1;")
		if (sum(values.inistate!=0)!=0){
		  ret <- paste0(ret,"\n\tsize_t num_sbj=(eta_0[0]->size)/(dim_latent_var);\n\tsize_t i;")
		}
		ret <- paste0(ret,"\n\tsize_t j;\n\tfor(j=0;j<num_regime;j++){")
		if (sum(values.inistate!=0)!=0){
		  ret <- paste0(ret,"\n\t\tfor(i=0;i<num_sbj;i++){\n")
		  for(i in 1:length(values.inistate)){
		    if(params.inistate[i] > 0){
		      ret <- paste(ret,
		                   '\t\t\tgsl_vector_set((eta_0)[j],i*dim_latent_var+', i-1,
		                   ', param[', params.inistate[i] - 1, ']);\n', sep='')
		    } else if(values.inistate[i] != 0){
		      ret <- paste(ret,
		                   '\t\t\tgsl_vector_set((eta_0)[j],i*dim_latent_var+', i-1,
		                   ', ', values.inistate[i], ');\n', sep='')
		    }
		  }
		  ret <- paste0(ret,"\t\t}")
		}

		ret <- paste(ret, setGslMatrixElements(values.inicov,params.inicov, "(error_cov_0)[j]"), sep="\n")
		ret <- paste(ret, "\t}\n}\n")
		
		object@c.string <- ret
		
		return(object)
	}
)

setMethod("writeCcode", "dynrNoise",
	function(object){
		params.latent <- object$params.latent
		params.observed <- object$params.observed
		#Note: should we mutate these other slots of the object or leave them alone?
		values.latent <- object@values.latent.inv.ldl
		values.observed <- object@values.observed.inv.ldl
		
		ret <- "void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov){\n\n"
		ret <- paste(ret, setGslMatrixElements(values.latent, params.latent, "eta_noise_cov"), sep="\n")
		ret <- paste(ret, setGslMatrixElements(values.observed, params.observed, "y_noise_cov"), sep="\n")
		ret <- paste(ret, "\n}\n\n")

		object@c.string <- ret
		

		return(object)
	}
)

setMethod("writeCcode", "dynrTrans",
          function(object){
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
                                  values.observed, values.latent, values.inicov, show=TRUE) { 
  return(standardGeneric("createRfun")) 
})

setMethod("createRfun", "dynrTrans",
          function(object, param.data, params.observed, params.latent, params.inicov,
                   values.observed, values.latent, values.inicov){
            #inv.tfun
            if (length(object@formula.inv)==0){
              #TODO If formula.inv is missing, point-wise inverse functions that are based on the uniroot function will be used to calculate starting values.
              #inverse = function (f, lower = -100, upper = 100) {function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)$root} 
              #eval(parse(text="inv.tf<-function(namedvec){return(c(a=inverse(exp)(namedvec[\"a\"]),b=inverse(function(x)x^2,0.0001,100)(namedvec[\"b\"])))}")) 
            }else{
              fml.str=formula2string(object@formula.inv)
              lhs=fml.str$lhs
              rhs=fml.str$rhs
              sub=sapply(lhs,function(x){paste0("vec[",param.data$param.number[param.data$param.name==x],"]")})
              f.string<-"inv.tf<-function(vec){"
              for (j in 1:length(rhs)){
                #TODO modify the sub pattern, and make sure "a" in "abs" will not be substituted
                for (i in 1:length(lhs)){
                  rhs[j]=gsub(lhs[i],sub[i],rhs[j])
                }
                eq=paste0(sub[j],"=",rhs[j])
                f.string<-paste(f.string,eq,sep="\t\n")
              }
              f.string2 <- f.string
              #observed
              if (sum(param.data$ldl.observed)>0){
                f.string2 <- paste(f.string2, makeldlchar(param.data, "ldl.observed", values.observed, params.observed, reverse=TRUE),sep="\t\n")
              }
              #latent
              if (sum(param.data$ldl.latent)>0){
                f.string2 <- paste(f.string2, makeldlchar(param.data, "ldl.latent", values.latent, params.latent, reverse=TRUE),sep="\t\n")
              }
              #inicov
              if (sum(param.data$ldl.inicov)>0){
                f.string2 <- paste(f.string2, makeldlchar(param.data, "ldl.inicov", values.inicov, params.inicov, reverse=TRUE),sep="\t\n")
              }
              f.string2<-paste0(f.string2,"\t\nreturn(vec)}")
              f.string<-paste0(f.string,"\t\nreturn(vec)}")
              
              eval(parse(text=f.string2))
              object@inv.tfun.full <- inv.tf
              
              eval(parse(text=f.string)) 
              object@inv.tfun <- inv.tf
            }
            
            f.string<-"tf<-function(vec){"
              
              if (length(object@formula.trans)>0){
              fml.str=formula2string(object@formula.trans)
              lhs=fml.str$lhs
              rhs=fml.str$rhs
              sub=sapply(lhs,function(x){paste0("vec[",param.data$param.number[param.data$param.name==x],"]")})
              
              for (j in 1:length(rhs)){
                #TODO modify the sub pattern, and make sure "a" in "abs" will not be substituted
                for (i in 1:length(lhs)){
                  rhs[j]=gsub(lhs[i],sub[i],rhs[j])
                }
                eq=paste0(sub[j],"=",rhs[j])
                f.string<-paste(f.string,eq,sep="\t\n")
              }
              object@paramnames<-lhs
              }
            
            #observed
            if (sum(param.data$ldl.observed)>0){
              f.string<-paste(f.string, makeldlchar(param.data, "ldl.observed", values.observed, params.observed),sep="\t\n")
            }
            #latent
            if (sum(param.data$ldl.latent)>0){
              f.string<-paste(f.string, makeldlchar(param.data, "ldl.latent", values.latent, params.latent),sep="\t\n")
            }
            #inicov
            if (sum(param.data$ldl.inicov)>0){
              f.string<-paste(f.string, makeldlchar(param.data, "ldl.inicov", values.inicov, params.inicov),sep="\t\n")
            }
              
            f.string<-paste0(f.string,"\t\nreturn(vec)}")
            eval(parse(text=f.string))   
            object@tfun <- tf
              
			      
            
            return(object)
          }
)
makeldlchar<-function(param.data, ldl.char, values, params, reverse=FALSE){
  vec.noise=paste0("vec[",paste0("c(",paste(param.data$param.number[param.data[,ldl.char]],collapse=","),")"),"]")
  param.name=param.data$param.name[param.data[,ldl.char]]

  vec.sub=as.vector(params)
  mat.index=sapply(param.name,function(x){min(which(vec.sub==x))})
  for (i in 1:length(param.name)){
    vec.sub=gsub(param.name[i], paste0("vec[",param.data$param.number[param.data$param.name==param.name[i]],"]"), vec.sub)
  }
  vec.sub[which(vec.sub=="fixed")]<-as.vector(values)[which(vec.sub=="fixed")]
  
  char=paste0(vec.noise,"=as.vector(", ifelse(reverse, 'reverseldl', 'transldl'), "(matrix(",paste0("c(",paste(vec.sub,collapse=","),")"),",ncol=",ncol(params),")))[",paste0("c(",paste(mat.index,collapse=","),")"),"]")
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
##' The goal is for factor a square, symmetrix, positive (semi-)definite matrix into the product of a lower triangular matrix, a diagonal matrix, and the transpose of the lower triangular matrix.
##' The value returned is a lower triangular matrix with the elements of D on the diagonal.
dynr.ldl <- function(x){
	ret <- t(chol(x))
	d <- diag(ret)
	ret <- ret %*% diag(1/d, nrow=length(d))
	diag(ret) <- d^2
	return(ret)
}

reverseldl<-function(values){
	if(dim(values)[1]==1){
		return(log(values))
	}else{
		mat <- dynr.ldl(values)
		diag(mat) <- log(diag(mat))
		return(mat)
	}
}

#------------------------------------------------------------------------------
# Some usefull helper functions
#

# a new method for diag with character input
setMethod("diag", "character", #diag.character <-
	function(x=1, nrow, ncol){
		n <- length(x)
		if(!missing(nrow)) n <- nrow
		if(missing(ncol)) ncol <- n
		r <- matrix("0", nrow, ncol)
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
	x <- tolower(c(x))
	sel <- pmatch(x, "freed", duplicates.ok=TRUE)
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
	x <- tolower(c(x))
	sel <- x %in% c("0", "fix", "fixed")
	x[sel] <- "fixed"
	x <- matrix(as.character(x), numRow, numCol, dimnames=list(rowNam, colNam))
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
##' 
##' @details
##' The default pattern for 'idvar' is to fix the first factor loading for each factor to one.
##' The variable names listed in 'idvar' have their factor loadings fixed to one.
##' 
##' @examples
##' #Single factor model with one latent variable
##' prep.loadings( list(eta1=paste0('y', 1:4)), 4:6)
##' 
##' # Two factor model with simple structure
##' prep.loadings( list(eta1=paste0('y', 1:4), eta2=paste0('y', 5:7)), c(4:6, 1:2))
##' 
##' #Two factor model with repeated use of a free parameter
##' prep.loadings( list(eta1=paste0('y', 1:4), eta2=paste0('y', 5:8)), c(4:6, 1:2, 4))
##' 
##' #Two factor model with a cross loading
##' prep.loadings( list(eta1=paste0('y', 1:4), eta2=c('y5', 'y2', 'y6')), c(4:6, 1:2))
prep.loadings <- function(map, params, idvar,exo.names=NULL){
	if(missing(idvar)){
		idvar <- sapply(map, '[', 1)
	}
	
	allVars <- unique(unlist(map))
	
	nx <- length(allVars)
	ne <- length(map)
	
	if(!all(idvar %in% c(names(map), unlist(map)))){
		stop("The 'idvar' must all be either in the names of the 'map' argument or parts of the part 'map' argument.")
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
##' @param values.load matrix or list of matrices. Values of the factor loadings.
##' @param params.load matrix or list of matrices. Params of the factor loadings.
##' @param values.exo matrix or list of matrices. Values of the covariate effects.
##' @param params.exo matrix or list of matrices. Params of the covariate effects.
##' @param values.int matrix or list of matrices. Values of the intercepts.
##' @param params.int matrix or list of matrices. Params of the intercerpts.
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
##' prep.measurement(diag(1, 5), diag(1:5))
##' prep.measurement(matrix(1, 5, 5), diag(1:5))
##' prep.measurement(diag(1, 5), diag(0, 5)) #identity measurement model
##' 
##' #Regime-switching measurement model where the first latent variable is
##' # active for regime 1, and the second latent variable is active for regime 2
##' # No free parameters are present.
##' prep.measurement(values.load=list(matrix(c(1,0), 1, 2), matrix(c(0,1), 1, 2)))
prep.measurement <- function(values.load, params.load, values.exo, params.exo, values.int, params.int,
                             obs.names,state.names,exo.names){

  
  if(!is.list(values.load)){
		values.load <- list(values.load)
	}
	if(missing(params.load)){
		params.load <- rep(list(matrix(0, nrow(values.load[[1]]), ncol(values.load[[1]]))), length(values.load))
	}
	if(!is.list(params.load)){
		params.load <- list(params.load)
	}
	if(missing(values.exo)){
		values.exo <- list()
		params.exo <- list()
	}
	if(missing(params.exo)){
		params.exo <- list()
	}
  
  if(missing(obs.names)){
    obs.names = paste0('y',1:nrow(values.load[[1]]))
  }
  
  if(missing(state.names)){
    state.names = paste0('state',1:ncol(values.load[[1]]))
  }
  
  if(missing(exo.names)){
    exo.names = character(0)
  }
  
	if(missing(values.int)){
		values.int <- list()
		params.int <- list()
	}
	if(missing(params.int)){
		params.int <- rep(list(matrix(0, nrow(values.int[[1]]), ncol(values.int[[1]]))), length(values.int))
	}
	values.load <- lapply(values.load, preProcessValues)
	params.load <- lapply(params.load, preProcessParams)
	values.exo <- lapply(values.exo, preProcessValues)
	params.exo <- lapply(params.exo, preProcessParams)
	values.int <- lapply(values.int, preProcessValues)
	params.int <- lapply(params.int, preProcessParams)
	sv <- c(extractValues(values.load, params.load), extractValues(values.exo, params.exo), extractValues(values.int, params.int))
	pn <- c(extractParams(params.load), extractParams(params.exo), extractParams(params.int))
	sv <- extractValues(sv, pn)
	pn <- extractParams(pn)
	x <- list(startval=sv, paramnames=pn, values.load=values.load, params.load=params.load,
		values.exo=values.exo, params.exo=params.exo, values.int=values.int, params.int=params.int,
		obs.names=obs.names, state.names=state.names,exo.names=exo.names)
	return(new("dynrMeasurement", x))
}

# Examples
# a <- prep.loadings( list(eta1=paste0('y', 1:4), eta2=c('y5', 'y2', 'y6')), c(4:6, 1:2))
# prep.measurement(a@values, a@params)
#
# prep.measurement(diag(1, 5), diag(1:5))
# prep.measurement(matrix(1, 5, 5), diag(1:5))
# prep.measurement(diag(1, 5), diag(0, 5)) #identity measurement model


#------------------------------------------------------------------------------
# Error covariance matrix
# N.B. This function produces BOTH the latent and observed error covariance matrices.

# TODO break below into two functions that populate into the same C function
#  at cook time.-- maybe we do not need this

##' The translation function for measurement error and process noise covariances
##' Output a C function to set up measurement error and process noise covariances and the starting values of the related parameters.
##' 
##' @param values.latent a positive definite matrix of the starting or fixed values of the process noise covariance matrix. To ensure the matrix is positive definite in estimation, we apply LDL transformation to the matrix. Values are hence automatically adjusted for this purpose. If theorectically an element is of value 0, please adjust it to some small number (e.g., 0.000001).
##' @param params.latent a matrix of the parameter indices of the measurement error covariance. If an element is 0, the corresponding element is fixed at the value specified in the values matrix; Otherwise, the corresponding element is to be estimated with the starting value specified in the values matrix.
##' @param values.observed a positive definite matrix of the starting or fixed values of the measurement error covariance matrix. To ensure the matrix is positive definite in estimation, we apply LDL transformation to the matrix. Values are hence automatically adjusted for this purpose. If theorectically an element is of value 0, please adjust it to some small number (e.g., 0.000001).
##' @param params.observed a matrix of the parameter indices of the process noise covariance. If an element is 0, the corresponding element is fixed at the value specified in the values matrix; Otherwise, the corresponding element is to be estimated with the starting value specified in the values matrix.
prep.noise <- function(values.latent, params.latent, values.observed, params.observed){
	values.latent <- preProcessValues(values.latent)
	params.latent <- preProcessParams(params.latent)
	values.observed <- preProcessValues(values.observed)
	params.observed <- preProcessParams(params.observed)
	
	values.latent.inv.ldl <- replaceDiagZero(values.latent)
	values.latent.inv.ldl <- reverseldl(values.latent.inv.ldl)
	values.observed.inv.ldl <- replaceDiagZero(values.observed)
	values.observed.inv.ldl <- reverseldl(values.observed.inv.ldl)
	
	sv <- c(extractValues(values.latent.inv.ldl, params.latent), extractValues(values.observed.inv.ldl, params.observed))
	pn <- c(extractParams(params.latent), extractParams(params.observed))
	sv <- extractValues(sv, pn)
	pn <- extractParams(pn)
	x <- list(startval=sv, paramnames=pn, values.latent=values.latent, values.observed=values.observed, params.latent=params.latent, params.observed=params.observed, values.latent.inv.ldl=values.latent.inv.ldl, values.observed.inv.ldl=values.observed.inv.ldl)
	return(new("dynrNoise", x))
}

# Examples
#prep.noise(values.latent=diag(c('Free', 1)), params.latent=diag(c('fixed', 3)),
#	values.observed=diag(1.5,1), params.observed=diag(4, 1))


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

##' Recipe function for creating Regime Switching
##' 
##' @param values matrix giving the values. Should have Number of Regimes rows and Number of Regimes time Number of Covariates columns
##' @param param matrix of the same size as values giving the free parameters
##' 
##' @details
##' Note that the ROW sums for the transition probability matrix must be one.
##' 
##' @examples
##' #Regime-switching with no covariates (self-transition ID)
##' b <- prep.regimes(values=matrix(0, 3, 3), params=matrix(c(0, 1, 2, 3, 0, 4, 5, 6, 0), 3, 3))
##' 
##' #Regime switching with no covariates (second regime ID)
##' b <- prep.regimes(values=matrix(0, 3, 3), params=matrix(c(1, 2, 3, 0, 0, 0, 4, 5, 6), 3, 3))
##' 
##' #2 regimes with three covariates
##' b <- prep.regimes(values=matrix(c(0), 2, 8), params=matrix(c(8:23), 2, 8), covariates=c('x1', 'x2', 'x3'))
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
	  stop('The dimensions of values and matrix must match.')
	}else if(!ncol(values)== nrow(values)*(length(covariates)+1)){
	  stop(paste0("The matrix values should have ",nrow(values)*(length(covariates)+1), " columns."))
	}
	
	sv <- extractValues(values, params)
	pn <- extractParams(params)
	
	if(length(sv) > nrow(values)*(length(covariates)+1)){
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
##' @param jacobian a list of formulas specifying the jacobian matrices of the drift/
##' state-transition, which can be missing
prep.formulaDynamics <- function(formula, startval, isContinuousTime=FALSE, jacobian){
  if(is.null(names(startval))){
    stop('startval must be a named vector')
  }
  x <- list(formula=formula, startval=startval, paramnames=names(startval), isContinuousTime=isContinuousTime)
  if (missing(jacobian)){
    autojcb=try(lapply(formula,autojacob,length(formula[[1]])))
    if (class(autojcb) == "try-error") {
      stop("Automatic differentiantion is not available. Please provide the jacobian functions.")
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
##' @param params.dyn the parameters matrix for the linear dynamics
##' @param values.dyn the values matrix for the linear dynamics
##' @param params.exo the parameters matrix for the effect of the covariates on the dynamic outcome (see details)
##' @param values.exo the values matrix for the effect of the covariates on the dynamic outcome (see details)
##' @param params.int the parameters vector for the intercept in the dynamic model (see details)
##' @param values.int the values vector for the intercept in the dynamic model (see details)
##' @param covariates the names or the index numbers of the covariates used for the dynamics
##' @param isContinuousTime logical. When TRUE, use a continuous time model.  When FALSE use a discrete time model. Either 'discrete' or 'continuous'.  Partial matching is used so 'c' or 'd' is sufficient. Capitalization is ignored.
##' 
##' @details
##' The dynamic outcome is the latent variable vector at the next time point in the discrete time case,
##' and the derivative of the latent variable vector at the current time point in the continuous time case.
prep.matrixDynamics <- function(params.dyn, values.dyn, params.exo, values.exo, params.int, values.int, 
                                covariates, isContinuousTime){
	# Handle numerous cases of missing or non-list arguments
	# General idea
	# If they give us a non-list argument, make it a one-element list
	# If they don't give us params, assume it's all fixed params
	# If they don't give us values, assume they don't want that part of the model
	if(!is.list(values.dyn)){
		values.dyn <- list(values.dyn)
	}
	if(missing(params.dyn)){
		params.dyn <- rep(list(matrix(0, nrow(values.dyn), ncol(values.dyn))), length(values.dyn))
	}
	if(!is.list(params.dyn)){
		params.dyn <- list(params.dyn)
	}
	if(missing(values.exo)){
		values.exo <- list()
		params.exo <- list()
	}
	if(missing(params.exo)){
		params.exo <- rep(list(matrix(0, nrow(values.exo[[1]]), ncol(values.exo[[1]]))), length(values.exo))
	}
	if(missing(values.int)){
		values.int <- list()
		params.int <- list()
	}
	if(missing(params.int)){
		params.int <- rep(list(matrix(0, nrow(values.int[[1]]), ncol(values.int[[1]]))), length(values.int))
	}
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


##' The translation function for initial conditions
##' Output a C function to set up initial conditions (i.e., intial state vector, initial error covariance matrix, and initial probabilities of being in each regime) and the starting values of the related parameters.
##' 
##' @param values.inistate a vector of the starting or fixed values of the initial state vector
##' @param params.inistate a vector of the parameter indices of the initial state vector. If an element is 0, the corresponding element is fixed at the value specified in the values vector; Otherwise, the corresponding element is to be estimated with the starting value specified in the values vector.
##' @param values.inicov a positive definite matrix of the starting or fixed values of the initial error covariance matrix. To ensure the matrix is positive definite in estimation, we apply LDL transformation to the matrix. Values are hence automatically adjusted for this purpose. If theorectically an element is of value 0, please adjust it to some small number (e.g., 0.000001).
##' @param params.inicov a matrix of the parameter indices of the initial error covariance matrix. If an element is 0, the corresponding element is fixed at the value specified in the values matrix; Otherwise, the corresponding element is to be estimated with the starting value specified in the values matrix.
##' @param values.regimep a vector of the starting or fixed values of the initial probabilities of being in each regime. By default, the initial probability of being in the first regime is fixed at 1.
##' @param params.regimep a vector of the parameter indices of the initial probabilities of being in each regime. If an element is 0, the corresponding element is fixed at the value specified in the values vector; Otherwise, the corresponding element is to be estimated with the starting value specified in the values vector.
prep.initial <- function(values.inistate, params.inistate, values.inicov, params.inicov, values.regimep=1, params.regimep=0){
	values.inistate <- preProcessValues(values.inistate)
	params.inistate <- preProcessParams(params.inistate)
	
	values.inicov <- preProcessValues(values.inicov)
	params.inicov <- preProcessParams(params.inicov)
	values.inicov.inv.ldl <- replaceDiagZero(values.inicov)
	values.inicov.inv.ldl <- reverseldl(values.inicov.inv.ldl)
	
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


matrix2formula <- function(x){
	if(!is.matrix(x)){
		stop("Dude! You have to give me a matrix. If you do, I'll give you a formula. Seriously.")
	}
	if(is.null(rownames(x))){
		rownames(x) <- paste0('y', 1:nrow(x))
	}
	if(is.null(colnames(x))){
		colnames(x) <- paste0('x', 1:ncol(x))
	}
	outcomes <- rownames(x)
	preds <- character(nrow(x))
	for(i in 1:nrow(x)){
		preds[i] <- paste(outcomes[i], paste(x[i,], colnames(x), sep='*', collapse=' + '), sep=' ~ ')
	}
	if(is.numeric(x)){
		preds <- gsub('1*', '', preds, fixed=TRUE)
		preds <- gsub(paste("0\\*(", paste(colnames(x), collapse = "|"), ") \\+ ", sep=""), '', preds)
		preds <- gsub(paste(" \\+ 0\\*(", paste(colnames(x), collapse = "|"), ")", sep=""), '', preds)
	} #else {
		#preds <- gsub('1*', '', preds, fixed=TRUE)
		#preds <- gsub(paste("0\\*(", paste(colnames(x), collapse = "|"), ") \\+ ", sep=""), '', preds)
		#preds <- gsub(paste(" \\+ 0\\*(", paste(colnames(x), collapse = "|"), ")", sep=""), '', preds)
	#}
	form <- lapply(preds, formula, env=.GlobalEnv)
	return(form)
}

## Examples
#meas <- prep.loadings(
#  map=list(
#    eta1=paste0('y', 1:3),
#    eta2=paste0('y', 4:6)),
#    params = c(1:4))
#matrix2formula(meas$values.load[[1]])
#matrix2formula(meas$params.load[[1]])
#lapply(matrix2formula(meas$values.load))

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


#------------------------------------------------------------------------------
prep.dP_dt <- "/**\n * The dP/dt function: depend on function_dF_dx, needs to be compiled on the user end\n * but user does not need to modify it or care about it.\n */\nvoid mathfunction_mat_to_vec(const gsl_matrix *mat, gsl_vector *vec){\n\tsize_t i,j;\n\tsize_t nx=mat->size1;\n\t/*convert matrix to vector*/\n\tfor(i=0; i<nx; i++){\n\t\tgsl_vector_set(vec,i,gsl_matrix_get(mat,i,i));\n\t\tfor (j=i+1;j<nx;j++){\n\t\t\tgsl_vector_set(vec,i+j+nx-1,gsl_matrix_get(mat,i,j));\n\t\t\t/*printf(\"%lu\",i+j+nx-1);}*/\n\t\t}\n\t}\n}\nvoid mathfunction_vec_to_mat(const gsl_vector *vec, gsl_matrix *mat){\n\tsize_t i,j;\n\tsize_t nx=mat->size1;\n\t/*convert vector to matrix*/\n\tfor(i=0; i<nx; i++){\n\t\tgsl_matrix_set(mat,i,i,gsl_vector_get(vec,i));\n\t\tfor (j=i+1;j<nx;j++){\n\t\t\tgsl_matrix_set(mat,i,j,gsl_vector_get(vec,i+j+nx-1));\n\t\t\tgsl_matrix_set(mat,j,i,gsl_vector_get(vec,i+j+nx-1));\n\t\t}\n\t}\n}\nvoid function_dP_dt(double t, size_t regime, const gsl_vector *p, double *param, size_t n_param, const gsl_vector *co_variate, gsl_vector *F_dP_dt){\n\t\n\tsize_t nx;\n\tnx = (size_t) sqrt(2*(double) p->size + .25) - .5;\n\tgsl_matrix *P_mat=gsl_matrix_calloc(nx,nx);\n\tmathfunction_vec_to_mat(p,P_mat);\n\tgsl_matrix *F_dx_dt_dx=gsl_matrix_calloc(nx,nx);\n\tfunction_dF_dx(t, regime, param, co_variate, F_dx_dt_dx);\n\tgsl_matrix *dFP=gsl_matrix_calloc(nx,nx);\n\tgsl_matrix *dP_dt=gsl_matrix_calloc(nx,nx);\n\tgsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, F_dx_dt_dx, P_mat, 0.0, dFP);\n\tgsl_matrix_transpose_memcpy(dP_dt, dFP);\n\tgsl_matrix_add(dP_dt, dFP);\n\tsize_t n_Q_vec=((1+nx)*nx)/2;\n\tgsl_vector *Q_vec=gsl_vector_calloc(n_Q_vec);\n\tsize_t i;\n\tfor(i=1;i<=n_Q_vec;i++){\n\t\t\tgsl_vector_set(Q_vec,n_Q_vec-i,param[n_param-i]);\n\t}\n\tgsl_matrix *Q_mat=gsl_matrix_calloc(nx,nx);\n\tmathfunction_vec_to_mat(Q_vec,Q_mat);\n\tgsl_matrix_add(dP_dt, Q_mat);\n\tmathfunction_mat_to_vec(dP_dt, F_dP_dt);\n\tgsl_matrix_free(P_mat);\n\tgsl_matrix_free(F_dx_dt_dx);\n\tgsl_matrix_free(dFP);\n\tgsl_matrix_free(dP_dt);\n\tgsl_vector_free(Q_vec);\n\tgsl_matrix_free(Q_mat);\n}\n"

#------------------------------------------------------------------------------
# Utility functions for writing GSL code

setGslMatrixElements <- function(values, params, name){
	ret <- ""
	numRow <- nrow(values)
	numCol <- ncol(values)
	for(j in 1:numCol){
		for(i in 1:numRow){
			if(params[i, j] > 0){
				ret <- paste(ret,
					'\tgsl_matrix_set(', name, ', ', i-1, ', ', j-1,
					', param[', params[i, j] - 1, ']);\n', sep='')
			} else if(values[i, j] != 0){
				ret <- paste(ret,
					'\tgsl_matrix_set(', name, ', ', i-1, ', ', j-1,
					', ', values[i, j], ');\n', sep='')
			}
		}
	}
	return(ret)
}

setGslVectorElements <- function(values, params, name){
  ret <- ""
  numLength <- length(values)
  for(i in 1:numLength){
      if(params[i] > 0){
        ret <- paste(ret,
                     '\tgsl_vector_set(', name, ', ', i-1,
                     ', param[', params[i] - 1, ']);\n', sep='')
      } else if(values[i] != 0){
        ret <- paste(ret,
                     '\tgsl_vector_set(', name, ', ', i-1,
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

