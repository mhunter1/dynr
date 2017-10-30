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
		ret <- ""
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
