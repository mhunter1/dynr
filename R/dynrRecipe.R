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
           paramnum = "numeric"
         )
)

##DONE##
setClass(Class = "dynrMeasurement",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnum = "numeric",
           values = "matrix",
           params = "matrix"),
         contains = "dynrRecipe"
)

##DONE##
setClass(Class = "dynrDynamics",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnum = "numeric"
           ),
         contains = "dynrRecipe"
)

setClass(Class = "dynrDynamicsFormula",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnum = "numeric",
           formula = "formula",
           jacobian = "formula",
           isContinuousTime = "logical"
           ),
         contains = "dynrDynamics"
)

setClass(Class = "dynrDynamicsMatrix",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnum = "numeric",
           values.dyn = "matrix",
           params.dyn = "matrix",
           values.exo = "matrix",
           params.exo = "matrix",
           isContinuousTime = "logical"
           ),
         contains = "dynrDynamics"
)

##DONE##
setClass(Class = "dynrRegimes",
         representation = representation(
           c.string =  "character",
           startval = "numeric",
           paramnum = "numeric",
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
           paramnum = "numeric",
           values.inistate = "matrix",
           params.inistate = "matrix",
           values.inicov = "matrix",
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
           paramnum = "numeric",
           values.latent = "matrix",
           params.latent = "matrix",
           values.observed = "matrix",
           params.observed = "matrix"),
         contains = "dynrRecipe"
)

setClass(Class = "dynrTrans",
         representation = representation(
           c.string =  "character",
           startval = "numeric",#not sure if needed in dynrTrans
           paramnum = "numeric",#not sure if needed in dynrTrans
           tfun="function",
           inv.tfun="function",
           formula.trans="list"
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

setMethod("initialize", "dynrTrans",
          function(.Object, x){
            for(i in names(x)){
              slot(.Object, name=i, check = TRUE) <- x[[i]]
            }
            return(.Object)
          }
)

#------------------------------------------------------------------------------
# printex method definitions

setGeneric("printex", function(object, observed, latent, covariates, show=TRUE) { 
	return(standardGeneric("printex")) 
})



setMethod("printex", "dynrMeasurement",
	function(object, observed, latent, covariates, show=TRUE){
		lC <- .xtableMatrix(object$values, show)
		return(invisible(list(measurement=lC)))
	}
)


# not sure what to do here yet
setMethod("printex", "dynrDynamicsFormula",
	function(object, observed, latent, covariates, show=TRUE){
#		lx0 <- .xtableMatrix(object$values.inistate)
#		lP0 <- .xtableMatrix(object$values.inicov)
#		lr0 <- .xtableMatrix(object$values.regimep)
#		return(list(initial.state=lxo, initial.covariance=lP0, initial.probability=lr0))
		message('Sorry, mate! This part is still under development.')
	}
)

setMethod("printex", "dynrDynamicsMatrix",
	function(object, observed, latent, covariates, show=TRUE){
		lA <- .xtableMatrix(object$values.dyn, show)
		#TODO add something better for covariates
		lB <- ifelse(nrow(object$values.exo) != 0, .xtableMatrix(object$values.exo, show), "")
		return(invisible(list(dyn=lA, exo=lB)))
	}
)


setMethod("printex", "dynrRegimes",
	function(object, observed, latent, covariates, show=TRUE){
		lG <- ifelse(nrow(object$values) != 0, .xtableMatrix(object$values, show), "")
		return(invisible(list(regimes=lG)))
	}
)


setMethod("printex", "dynrInitial",
	function(object, observed, latent, covariates, show=TRUE){
		lx0 <- .xtableMatrix(object$values.inistate, show)
		lP0 <- .xtableMatrix(object$values.inicov, show)
		lr0 <- .xtableMatrix(object$values.regimep, show)
		return(invisible(list(initial.state=lx0, initial.covariance=lP0, initial.probability=lr0)))
	}
)


setMethod("printex", "dynrNoise",
	function(object, observed, latent, covariates, show=TRUE){
		lQ <- .xtableMatrix(object$values.latent, show)
		lR <- .xtableMatrix(object$values.observed, show)
		return(invisible(list(dynamic.noise=lQ, measurement.noise=lR)))
	}
)


.xtableMatrix <- function(m, show){
	x <- xtable::xtable(m, align=rep("", ncol(m)+1))
	out <- print(x, floating=FALSE, tabular.environment="bmatrix", 
		hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE,
		print.results=show)
	return(out)
}

#printex(dynrModel)
#printex(dynrRecipe)


#------------------------------------------------------------------------------
# writeCcode method definitions

setGeneric("writeCcode", function(object, show=TRUE) { 
	return(standardGeneric("writeCcode")) 
})



setMethod("writeCcode", "dynrMeasurement",
	function(object){
		values <- object$values
		params <- object$params
		ret <- "void function_measurement(size_t t, size_t regime, double *param, const gsl_vector *eta, const gsl_vector *co_variate, gsl_matrix *Ht, gsl_vector *y){\n\n"
		ret <- paste(ret, setGslMatrixElements(values, params, "Ht"), sep="\n")
		ret <- paste(ret, "\n\tgsl_blas_dgemv(CblasNoTrans, 1.0, Ht, eta, 0.0, y);\n")
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
	  
	  if (missing(jacob)){
	    autojcb=try(lapply(formula,autojacob,n[1]))
	    if (class(autojcb) == "try-error") {
	      stop("Automatic differentiantion is not available. Please provide the jacobian functions.")
	    }else{
	      row=lapply(autojcb,"[[","row")
	      col=lapply(autojcb,"[[","col")
	      rhsj=lapply(autojcb,"[[","rhsj")
	      jacob=lapply(autojcb,"[[","jacob")
	    }
	  }else{
	    fmlj=lapply(jacob,processFormula)
	    row=lapply(fmlj,function(x){lapply(x,"[[",1)})
	    col=lapply(fmlj,function(x){lapply(x,"[[",2)})
	    rhsj=lapply(fmlj,function(x){lapply(x,"[[",3)})
	  }
	  
	  
	  #TODO in the continuous case x is stacked at the end of param in function_dF_dx.
	  #TODO in the continuous case allow users to use d()
	  #TODO add covariate
	  
	  if (isContinuosTime){
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
	          rhs[[1]][[i]]=gsub(lhs[[1]][[j]],paste0("gsl_vector_get(xstart,",j-1,")"),rhs[[1]][[i]])
	        }
	        ret=paste(ret,paste0("\tgsl_vector_set(x_tend,",i-1,",",rhs[[1]][[i]],");"),sep="\n\t")    
	      }
	    }
	    
	    ret=paste0(ret,"\n\t}")
	    
	    #function_dF_dx
	    ret=paste0(ret,"\n\nvoid function_dF_dx(double t, size_t regime, double *param, const gsl_vector *co_variate, gsl_matrix *F_dx_dt_dx){")
	    if (nregime>1){
	      ret=paste(ret,"switch (regime) {",sep="\n\t")
	      for (r in 1:nregime){
	        ret=paste(ret,paste0("case ",r-1,":"),sep="\n\t")
	        for (i in 1:length(jacob[[r]])){
	          for (j in 1:length(lhs[[r]])){
	            rhsj[[r]][[i]]=gsub(lhs[[r]][[j]],paste0("gsl_vector_get(xstart,",j-1,")"),rhsj[[r]][[i]])
	          }
	          
	          ret=paste(ret,paste0("\tgsl_matrix_set(F_dx_dt_dx,",which(lhs[[r]]==row[[r]][[i]])-1,",",which(lhs[[r]]==col[[r]][[i]])-1,",",rhsj[[r]][[i]],");"),sep="\n\t")    
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

setMethod("writeCcode", "dynrDynamicsMatrix",
	function(object){
		isContinuousTime <- object$isContinuousTime
		params.dyn <- object$params.dyn
		values.dyn <- object$values.dyn
		params.exo <- object$params.exo
		values.exo <- object$values.exo
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
		
		# Create dynamics (state-transition or drift matrix) with covariate effects
		ret <- paste(dynHead,
			createGslMatrix(nrow(params.dyn), ncol(params.dyn), "Amatrix"),
			setGslMatrixElements(values=values.dyn, params=params.dyn, name="Amatrix"),
			blasMV(FALSE, "1.0", "Amatrix", inName, "0.0", outName),
			destroyGslMatrix("Amatrix"),
			sep="\n")
		
		if(nrow(params.exo) != 0){
			ret <- paste(ret,
				createGslMatrix(nrow(params.exo), ncol(params.exo), "Bmatrix"),
				setGslMatrixElements(values=values.exo, params=params.exo, name="Bmatrix"),
				blasMV(FALSE, "1.0", "Bmatrix", "co_variate", "1.0", outName),
				destroyGslMatrix("Bmatrix"),
				sep="\n")
		}
		
		ret <- paste(ret, "}\n\n", sep="\n")
		
		
		# Create jacobian function
		ret <- paste(ret,
			jacHead,
			setGslMatrixElements(values=values.dyn, params=params.dyn, name=jacName),
			"}\n\n",
			sep="\n")
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
		values.inicov <- object$values.inicov
		params.inicov <- object$params.inicov
		values.regimep <- object$values.regimep
		params.regimep <- object$params.regimep
		ret <- "void function_initial_condition(double *param, gsl_vector **co_variate, gsl_vector *pr_0, gsl_vector **eta_0, gsl_matrix **error_cov_0){\n"
		ret <- paste(ret, setGslVectorElements(values.regimep,params.regimep, "pr_0"), sep="\n")
		ret <- paste0(ret,"\tsize_t num_regime=pr_0->size;\n\tsize_t dim_latent_var=error_cov_0[0]->size1;\n\tsize_t num_sbj=(eta_0[0]->size)/(dim_latent_var);\n\tsize_t i,j;\n\tfor(j=0;j<num_regime;j++){\n\t\tfor(i=0;i<num_sbj;i++){\n")
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
		values.inicov <- replaceDiagZero(values.inicov)
		values.inicov <- reverseldl(values.inicov)
		ret <- paste(ret, setGslMatrixElements(values.inicov,params.inicov, "(error_cov_0)[j]"), sep="\t\t}\n")
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
		values.latent <- replaceDiagZero(object@values.latent)
		values.observed <- replaceDiagZero(object@values.observed)
		values.latent <- reverseldl(values.latent)
		values.observed <- reverseldl(values.observed)
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
            
            n=length(object@formula.trans)
            if (n>0){
              formula.trans=object$formula.trans
              fml=processFormula(formula.trans)
              lhs=lapply(fml,function(x){x[1]})
              rhs=lapply(fml,function(x){x[2]})
              
              for (i in 1:n){
                ret=paste(ret,paste0(lhs[i],"=",rhs[i]),sep="\n\t") 
              }
            }
            
            ret=paste0(ret,"\n\t}\n")              
            object@c.string <- ret
            return(object)
          }
)

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
	} else {
		numRow <- nrow(x)
		numCol <- ncol(x)
	}
	x <- tolower(c(x))
	sel <- pmatch(x, "freed", duplicates.ok=TRUE)
	x[sel %in% 1] <- 0
	x <- matrix(as.numeric(x), numRow, numCol)
	return(x)
}

preProcessParams <- function(x){
	x[is.na(x)] <- 'fixed'
	if(is.null(dim(x))){
		numRow <- length(x)
		numCol <- 1
	} else {
		numRow <- nrow(x)
		numCol <- ncol(x)
	}
	x <- tolower(c(x))
	sel <- pmatch(x, "fixed", duplicates.ok=TRUE)
	x[sel %in% 1] <- 0
	x <- matrix(as.numeric(x), numRow, numCol)
	return(x)
}

extractWhichParams <- function(p){
	p!=0 & !duplicated(p, MARGIN=0)
}

extractParams <- function(p){
	p[extractWhichParams(p)]
}

extractValues <- function(v, p){
	#if(is.list(v)){
	#	ret <- c()
	#	for(i in 1:length(v)){
	#		ret <- c(ret, extractValues(v[[i]], p[[i]]))
	#	}
	#	return(ret)
	#}
  
	v[extractWhichParams(p)]
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
##'
prep.loadings <- function(map, params, idvar){
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
	x <- prep.measurement(values=valuesMat, params=paramsMat)
	return(x)
}

# Examples
# Single factor model with one latent variable
#prep.loadings( list(eta1=paste0('y', 1:4)), 4:6)

# Two factor model with simple structure
#prep.loadings( list(eta1=paste0('y', 1:4), eta2=paste0('y', 5:7)), c(4:6, 1:2))

# Two factor model with repeated use of a free parameter
#prep.loadings( list(eta1=paste0('y', 1:4), eta2=paste0('y', 5:8)), c(4:6, 1:2, 4))

# Two factor model with a cross loading
#prep.loadings( list(eta1=paste0('y', 1:4), eta2=c('y5', 'y2', 'y6')), c(4:6, 1:2))


#--------------------------------------
# matrix input version

# values, and params are all MxN matrices
# a zero param is taken to be fixed.

prep.measurement <- function(values, params){
	values <- preProcessValues(values)
	params <- preProcessParams(params)
	return(new("dynrMeasurement", list(startval=extractValues(values, params), paramnum=extractParams(params), values=values, params=params)))
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
	sv <- c(extractValues(values.latent, params.latent), extractValues(values.observed, params.observed))
	pn <- c(extractParams(params.latent), extractParams(params.observed))
	sv <- extractValues(sv, pn)
	pn <- extractParams(pn)
	x <- list(startval=sv, paramnum=pn, values.latent=values.latent, values.observed=values.observed, params.latent=params.latent, params.observed=params.observed)
	return(new("dynrNoise", x))
}

# Examples
#prep.noise(values.latent=diag(c('Free', 1)), params.latent=diag(c('fixed', 3)),
#	values.observed=diag(1.5,1), params.observed=diag(4, 1))


reverseldl<-function(values){
  if (dim(values)[1]==1){
    return(log(values))
  }else{
    mat<-KFAS::ldl(values)
    diag(mat)<-log(diag(mat))
    return(mat)
  }
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

##' Recipe function for creating Regime Switching
##' 
##' @param values matrix giving the values. Should have Number of Regimes rows and Number of Regimes time Number of Covariates columns
##' @param param matrix of the same size as values giving the free parameters
##' 
##' Note that the ROW sums for the transition probability matrix must be one.
prep.regimes <- function(values, params, covariates){
	if(!missing(values)){
		values <- preProcessValues(values)
	}
	if(!missing(params)){
		params <- preProcessParams(params)
	}
	#TODO check matrix dimensions
	#TODO check that some form of identification is made
	if(missing(values)){
		values <- matrix(0, 0, 0)
		params <- matrix(0, 0, 0)
	}
	sv <- extractValues(values, params)
	pn <- extractParams(params)
	x <- list(startval=sv, paramnum=pn, values=values, params=params)
	return(new("dynrRegimes", x))
}

# Examples
# Regime-switching with no covariates (self-transition ID)
#b <- prep.regimes(values=matrix(0, 3, 3), params=matrix(c(0, 1, 2, 3, 0, 4, 5, 6, 0), 3, 3))
#
# Regime switching with no covariates (second regime ID)
#b <- prep.regimes(values=matrix(0, 3, 3), params=matrix(c(1, 2, 3, 0, 0, 0, 4, 5, 6), 3, 3))
#
# 2 regimes with three covariates
#b <- prep.regimes(values=matrix(c(0), 2, 8), params=matrix(c(8:23), 2, 8), covariates=c('x1', 'x2', 'x3'))

# B <- matrix(c(8:(8+24-1)), nr, (nc+1)*nr, byrow=TRUE)
#matrix(t(B), nrow=nr*nr, ncol=nc+1, byrow=TRUE)

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

##' The translation function for the dynamic functions 
##' 
##' @param formula a list of formulas specifying the drift or state-transition equations for the latent variables in continuous or discrete time, respectively
##' @param jacobian a list of formulas specifying the jacobian matrices of the drift/state-transition
##' @param isContinuousTime If True, the left hand side of the formulas represent the first-order derivatives of the specified variables; if False, the left hand side of the formulas represent the current state of the specified variable while the same variable on the righ hand side is its previous state.  
##' @param ... 
prep.nonlindynamics <- function(formula, isContinuousTime=FALSE, jacobian){
  x <- list(formula=formula, jacobian=jacobian, isContinuosTime=isContinuousTime)
  return(new("dynrDynamicsFormula", x))
}

##' Recipe function for creating Linear Dynamcis
##' 
##' @param params.dyn the parameters matrix for the linear dynamics
##' @param values.dyn the values matrix for the linear dynamics
##' @param params.exo the parameters matrix for the effect of the covariates on the dynamic outcome (see details)
##' @param values.exo the values matrix for the effect of the covariates on the dynamic outcome (see details)
##' @param covariates the names or the index numbers of the covariates used for the dynamics
##' @param isContinuousTime logical. When TRUE, use a continuous time model.  When FALSE use a discrete time model. Either 'discrete' or 'continuous'.  Partial matching is used so 'c' or 'd' is sufficient. Capitalization is ignored.
##' 
##' @details
##' The dynamic outcome is the latent variable vector at the next time point in the discrete time case,
##' and the derivative of the latent variable vector at the current time point in the continuous time case.
prep.linearDynamics <- function(params.dyn, values.dyn, params.exo, values.exo, covariates, isContinuousTime){
	#time <- checkAndProcessTimeArgument(time)
	values.dyn <- preProcessValues(values.dyn)
	params.dyn <- preProcessParams(params.dyn)
	if(!missing(values.exo)){
		values.exo <- preProcessValues(values.exo)
	}
	if(!missing(params.exo)){
		params.exo <- preProcessParams(params.exo)
	}
	
	if(missing(params.exo)){
		params.exo <- matrix(0, 0, 0)
		values.exo <- matrix(0, 0, 0)
	}
	sv <- c(extractValues(values.dyn, params.dyn), extractValues(values.exo, params.exo))
	pn <- c(extractParams(params.dyn), extractParams(params.exo))
	sv <- extractValues(sv, pn)
	pn <- extractParams(pn)

	x <- list(startval=sv, paramnum=pn, params.dyn=params.dyn, values.dyn=values.dyn, params.exo=params.exo, values.exo=values.exo, isContinuousTime=isContinuousTime)
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
##' @param values.regimep a vector of the starting or fixed values of the initial probalities of being in each regime. By default, the initial probability of being in the first regime is fixed at 1.
##' @param params.regimep a vector of the parameter indices of the initial probalities of being in each regime. If an element is 0, the corresponding element is fixed at the value specified in the values vector; Otherwise, the corresponding element is to be estimated with the starting value specified in the values vector.
prep.initial <- function(values.inistate, params.inistate, values.inicov, params.inicov, values.regimep=1, params.regimep=0){
	values.inistate <- preProcessValues(values.inistate)
	params.inistate <- preProcessParams(params.inistate)
	values.inicov <- preProcessValues(values.inicov)
	params.inicov <- preProcessParams(params.inicov)
	values.regimep <- preProcessValues(values.regimep)
	params.regimep <- preProcessParams(params.regimep)
	sv <- c(extractValues(values.inistate, params.inistate), extractValues(values.inicov, params.inicov), extractValues(values.regimep, params.regimep))
	pn <- c(extractParams(params.inistate), extractParams(params.inicov), extractParams(params.regimep))
	sv <- extractValues(sv, pn)
	pn <- extractParams(pn)
	x <- list(startval=sv, paramnum=pn, values.inistate=values.inistate, params.inistate=params.inistate, values.inicov=values.inicov, params.inicov=params.inicov, values.regimep=values.regimep, params.regimep=params.regimep)
	return(new("dynrInitial", x))
}

##' The translation function for transformation functions of free parameters that are not in covariance structures
##' Output a C function to set up transformation functions
##' 
##' @param formula.trans a list of formulae that transform free parameters in the model that are not in the covariance structures
##' @param formula.inv a list of formulae that inverse the transformation on the free parameters. If formula.inv is missing, point-wise inverse functions that are based on the uniroot function will be used to calculate starting values.
prep.tfun<-function(formula.trans,formula.inv){
  #input: formula.trans=list(a~exp(a),b~b^2)
  #input: formula.inv=list(a~log(a),b~sqrt(b))
  if (missing(formula.trans)){
    #TODO modify this part if needed
    x<-list()
  }else{
    if (missing(formula.inv)){
      #TODO when formula.inv is not specified, use automatic inverse functions
      inverse = function (f, lower = -100, upper = 100) {function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)$root} 
      eval(parse(text="inv.tf<-function(namedvec){return(c(a=inverse(exp)(namedvec[\"a\"]),b=inverse(function(x)x^2,0.0001,100)(namedvec[\"b\"])))}")) 
    }else{
      fml.str=formula2string(formula.inv)
      lhs=fml.str$lhs
      rhs=fml.str$rhs
      sub=paste0("namedvec[\"",lhs,"\"]")
      f.string<-"inv.tf<-function(namedvec){return(c("
      for (i in 1:length(lhs)){
        #TODO modify the sub pattern, and make sure "a" in "abs" will not be substituted
        if (i!=1){f.string<-paste0(f.string,",")}
        f.string<-paste0(f.string,lhs[i],"=",gsub(lhs[i],sub[i],rhs[i]))
      }
      f.string<-paste0(f.string,"))}")
      eval(parse(text=f.string))
    }
    fml.str=formula2string(formula.trans)
    lhs=fml.str$lhs
    rhs=fml.str$rhs
    sub=paste0("namedvec[\"",lhs,"\"]")
    f.string<-"tf<-function(namedvec){return(c("
    for (i in 1:length(lhs)){
      #TODO modify the sub pattern, and make sure "a" in "abs" will not be substituted
      if (i!=1){f.string<-paste0(f.string,",")}
      f.string<-paste0(f.string,lhs[i],"=",gsub(lhs[i],sub[i],rhs[i]))
    }
    f.string<-paste0(f.string,"))}")
    eval(parse(text=f.string))
    x <- list(tfun=tf,inv.tfun=inv.tf,formula.trans=formula.trans)
  }
  
  return(new("dynrTrans", x))
}
formula2string<-function(formula.list){
  tuple=lapply(formula.list,as.list)
  lhs=sapply(tuple,function(x){deparse(x[[2]])})
  rhs=sapply(tuple,function(x){deparse(x[[3]])})
  return(list(lhs=lhs,rhs=rhs))
}

#------------------------------------------------------------------------------
prep.dP_dt <- "/**\n * The dP/dt function: depend on function_dF_dx, needs to be compiled on the user end\n * but user does not need to modify it or care about it.\n */\nvoid mathfunction_mat_to_vec(const gsl_matrix *mat, gsl_vector *vec){\n\tsize_t i,j;\n\tsize_t nx=mat->size1;\n\t/*convert matrix to vector*/\n\tfor(i=0; i<nx; i++){\n\t\tgsl_vector_set(vec,i,gsl_matrix_get(mat,i,i));\n\t\tfor (j=i+1;j<nx;j++){\n\t\t\tgsl_vector_set(vec,i+j+nx-1,gsl_matrix_get(mat,i,j));\n\t\t\t/*printf(\"%lu\",i+j+nx-1);}*/\n\t\t}\n\t}\n}\nvoid mathfunction_vec_to_mat(const gsl_vector *vec, gsl_matrix *mat){\n\tsize_t i,j;\n\tsize_t nx=mat->size1;\n\t/*convert vector to matrix*/\n\tfor(i=0; i<nx; i++){\n\t\tgsl_matrix_set(mat,i,i,gsl_vector_get(vec,i));\n\t\tfor (j=i+1;j<nx;j++){\n\t\t\tgsl_matrix_set(mat,i,j,gsl_vector_get(vec,i+j+nx-1));\n\t\t\tgsl_matrix_set(mat,j,i,gsl_vector_get(vec,i+j+nx-1));\n\t\t}\n\t}\n}\nvoid function_dP_dt(double t, size_t regime, const gsl_vector *p, double *param, size_t n_param, const gsl_vector *co_variate, gsl_vector *F_dP_dt){\n\t\n\tsize_t nx;\n\tnx = (size_t) floor(sqrt(2*(double) p->size));\n\tgsl_matrix *P_mat=gsl_matrix_calloc(nx,nx);\n\tmathfunction_vec_to_mat(p,P_mat);\n\tgsl_matrix *F_dx_dt_dx=gsl_matrix_calloc(nx,nx);\n\tfunction_dF_dx(t, regime, param, co_variate, F_dx_dt_dx);\n\tgsl_matrix *dFP=gsl_matrix_calloc(nx,nx);\n\tgsl_matrix *dP_dt=gsl_matrix_calloc(nx,nx);\n\tgsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, F_dx_dt_dx, P_mat, 0.0, dFP);\n\tgsl_matrix_transpose_memcpy(dP_dt, dFP);\n\tgsl_matrix_add(dP_dt, dFP);\n\tsize_t n_Q_vec=(1+nx)*nx/2;\n\tgsl_vector *Q_vec=gsl_vector_calloc(n_Q_vec);\n\tsize_t i;\n\tfor(i=1;i<=n_Q_vec;i++){\n\t\t\tgsl_vector_set(Q_vec,n_Q_vec-i,param[n_param-i]);\n\t}\n\tgsl_matrix *Q_mat=gsl_matrix_calloc(nx,nx);\n\tmathfunction_vec_to_mat(Q_vec,Q_mat);\n\tgsl_matrix_add(dP_dt, Q_mat);\n\tmathfunction_mat_to_vec(dP_dt, F_dP_dt);\n\tgsl_matrix_free(P_mat);\n\tgsl_matrix_free(F_dx_dt_dx);\n\tgsl_matrix_free(dFP);\n\tgsl_matrix_free(dP_dt);\n\tgsl_vector_free(Q_vec);\n\tgsl_matrix_free(Q_mat);\n}\n"

#TODO change 
#	nx = (size_t) floor(sqrt(2*(double) p->size));
# to
#	nx = (size_t) sqrt(2*double p->size + .25) - .5
#
# 	size_t n_Q_vec=((1+nx)*nx)/2;



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

