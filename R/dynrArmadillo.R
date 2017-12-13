#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2017-10-30 08:50:19
# Filename: dynrArmadillo.R
# Purpose: Replicate much of the dynrRecipe methods for writing C code but for
#  armadillo code instead.
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Define Generic

setGeneric("writeArmadilloCode", 
           function(object, covariates, show=TRUE) { 
             return(standardGeneric("writeArmadilloCode")) 
           })

#------------------------------------------------------------------------------
# Define method for dynrMeasurement class

setMethod("writeArmadilloCode", "dynrMeasurement",
	function(object, covariates){
		ret <- ""
		object@c.string <- ret
		return(object)
	}
)

#------------------------------------------------------------------------------
# Define method for dynrMatrixDynamics class

setMethod("writeArmadilloCode", "dynrMatrixDynamics",
	function(object, covariates){
		ret <- ""
		object@c.string <- ret
		return(object)
	}
)

#------------------------------------------------------------------------------
# Define method for dynrFormulaDynamics class

setMethod("writeArmadilloCode", "dynrFormulaDynamics",
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
                  ret="void  function_dx_dt(double t, size_t regime, const vec *x, double *param, size_t n_param, const vec *co_variate, vec *F_dx_dt){"
                  
                  if (nregime>1){
                      ret=paste(ret,"switch (regime) {",sep="\n\t")
                      for (r in 1:nregime){
                          ret=paste(ret,paste0("\tcase ",r-1,":"),sep="\n\t")
                          for (i in 1:n[r]){
                              for (j in 1:length(lhs[[r]])){
                                  rhs[[r]][[i]]=gsub(paste0("\\<",lhs[[r]][[j]],"\\>"),paste0("(*x)(",j-1,")"),rhs[[r]][[i]])
                              }
                              #ret=paste(ret,paste0("\tgsl_vector_set(F_dx_dt,",i-1,",",rhs[[r]][[i]],");"),sep="\n\t")   
                              ret=paste(ret,paste0("\t(*F_dx_dt)(",i-1,") =",rhs[[r]][[i]],";"),sep="\n\t")  
                          }
                          ret=paste(ret,paste0("break;\n"),sep="\n\t")
                          
                      }
                      ret=paste(ret,paste0("\t}"),sep="\n\t")
                      
                  }else{
                      for (i in 1:n){
                          for (j in 1:length(lhs[[1]])){
                              rhs[[1]][[i]]=gsub(paste0("\\<",lhs[[1]][[j]],"\\>"),paste0("(*x)(",j-1,")"),rhs[[1]][[i]])
                          }
                          #ret=paste(ret,paste0("\tgsl_vector_set(F_dx_dt,",i-1,",",rhs[[1]][[i]],");"),sep="\n\t")  
                          ret=paste(ret,paste0("\tF_dx_dt(",i-1,")=",rhs[[1]][[i]],";"),sep="\n\t")            
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

#------------------------------------------------------------------------------
# Define method for dynrNoise class

# TODO This is probably the easiest armadillo code to write
#  It will be the most similar to the corresponding writeCcode method.

setMethod("writeArmadilloCode", "dynrNoise",
	function(object, covariates){
		ret <- ""
		object@c.string <- ret
		return(object)
	}
)

#------------------------------------------------------------------------------
# Define method for dynrInitial class

setMethod("writeArmadilloCode", "dynrInitial",
	function(object, covariates){
		ret <- ""
		object@c.string <- ret
		return(object)
	}
)

#------------------------------------------------------------------------------
# Define method for dynrRegime class

setMethod("writeArmadilloCode", "dynrRegimes",
	function(object, covariates){
		ret <- ""
		object@c.string <- ret
		return(object)
	}
)

#------------------------------------------------------------------------------
# End
