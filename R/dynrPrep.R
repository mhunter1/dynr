# 
# dynr model CLASS
# 
setClass(Class =  "dynrModel",
         representation = representation(
           dynamics =  "dynrDynamics",
           measurement = "dynrMeasurement",
           noise = "dynrNoise",
           initial = "dynrInitial",
           regimes= "dynrRegimes",
           num_regime="integer",
           dim_latent_var="integer",
           infile="character",
           outfile="character",
           isContinuousTime="logical",
           verbose="logical",
           compileLib="logical",
           xstart="vector",
           ub="vector",
           lb="vector",
           options="list"
         ),
         prototype(
           num_regime=as.integer(1),
           isContinuousTime=TRUE,
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

setMethod("printex", "dynrModel",
	function(object, observed, latent, covariates, show=TRUE){
		meas <- printex(object$measurement, show=FALSE)
		dyn <- printex(object$dynamics, show=FALSE)
		reg <- printex(object$regimes, show=FALSE)
		noise <- printex(object$noise, show=FALSE)
		init <- printex(object$initial, show=FALSE)
		message(' :(  Dagnabbit. This part is not quite working yet.')
		measTex <- paste("The measurement model is given by\n\\begin{equation}\n",
			.xtableMatrix(matrix(observed, nrow=length(observed), ncol=1), show=FALSE),
			" = ",
			meas$measurement,
			.xtableMatrix(matrix(latent, nrow=length(latent), ncol=1), show=FALSE),
			" + \\vec{r}\n",
			"\\end{equation}\nwith\n",
			"\\begin{equation}\n\\text{Cov}(\\vec{r}) = ",
			noise$measurement.noise,
			"\\end{equation}\n", sep="")
		dynTex <- paste("The dynamic model is given by\n\\begin{equation}\n",
			ifelse(object$isContinuousTime, "\\frac{d}{dt} ", ""),
			.xtableMatrix(matrix(latent, nrow=length(latent), ncol=1), show=FALSE),
			ifelse(object$isContinuousTime, "", "_t"),
			" = ",
			dyn$dyn,
			.xtableMatrix(matrix(latent, nrow=length(latent), ncol=1), show=FALSE),
			ifelse(object$isContinuousTime, "", "_t"),
			" + ",
			"\\vec{q}\n",
			"\\end{equation}\nwith\n",
			"\\begin{equation}\n\\text{Cov}(\\vec{q}) = ",
			noise$dynamic.noise,
			"\\end{equation}\n", sep="")
		return(paste(measTex, dynTex, sep="\n"))
		#
		# make equations
		# y = C x + r with
		# Cov(r) = measurement.noise
		# Make a matrix of the names of the observed variables for y
		# Make a matreis of the names of the latent variables for x
		# C is the meas$measurement factor loadings
		#
		# x = dynamics(x) + q with
		# Cov(q) = dynamic.noise
	}
)


# modeling is what happens to recipes.
# The alpha version of this file just takes a bunch of recipes and puts
#  them together into a C file of the user's naming.
#

.logisticCFunction <- "/**\n * This function takes a double and gives back a double\n * It computes the logistic function (i.e. the inverse of the logit link function)\n * @param x, the double value e.g. a normally distributed number\n * @return logistic(x), the double value e.g. a number between 0 and 1\n */\ndouble mathfunction_logistic(const double x){\n\tdouble value = 1.0/(1.0 + exp(-x));\n\treturn value;\n}\n"

.softmaxCFunction <- "/**\n * This function takes a gsl_vector and modifies its second argument (another gsl_vector)\n * It computes the softmax function (e.g. for multinomial logistic regression)\n * @param x, vector of double values e.g. a vector of normally distributed numbers\n * @param result, softmax(x), e.g. a vector of numbers between 0 and 1 that sum to 1\n */\nvoid mathfunction_softmax(const gsl_vector *x, gsl_vector *result){\n\t/* Elementwise exponentiation */\n\tsize_t index=0;\n\tfor(index=0; index < x->size; index++){\n\t\tgsl_vector_set(result, index, exp(gsl_vector_get(x, index)));\n\t}\n\t\n\t/* Sum for the scaling coeficient */\n\tdouble scale = 0.0;\n\tfor(index=0; index < x->size; index++){\n\t\tscale += gsl_vector_get(result, index);\n\t}\n\t\n\t/* Multiply all elements of result by 1/scale */\n\tgsl_blas_dscal(1/scale, result);\n}\n"

.cfunctions <- paste(.logisticCFunction, .softmaxCFunction, sep="\n")

dynr.model <- function(dynamics, measurement, noise, initial, ..., infile=tempfile(),outfile="./demo/cooked"){
  #initiate a dynrModel object
  obj.dynrModel=new("dynrModel",list(infile=infile, outfile=outfile, dynamics=dynamics, measurement=measurement, noise=noise, initial=initial, ...))
  obj.dynrModel@dim_latent_var=dim(obj.dynrModel@noise@values.latent)[1]
  inputs <- list(dynamics=dynamics, measurement=measurement, noise=noise, initial=initial,...)
  obj.dynrModel@xstart<-unlist(sapply(inputs, slot, name='startval'))
  obj.dynrModel@ub<-rep(9999,length(obj.dynrModel@xstart))
  obj.dynrModel@lb<-rep(9999,length(obj.dynrModel@xstart))
  #write out the C script
  cparts <- sapply(inputs, slot, name='c.string')
  includes <- "#include <math.h>\n#include <gsl/gsl_matrix.h>\n#include <gsl/gsl_blas.h>\n"
  body <- paste(cparts, collapse="\n\n")
  if( length(grep("void function_regime_switch", body)) == 0 ){ # if regime-switching function isn't provided, fill in 1 regime model
    body <- paste(body, writeCcode(prep.regimes())$c.string, sep="\n\n")
  }
  glom <- paste(includes, body, prep.dP_dt, .cfunctions, sep="\n\n")
  cat(glom, file=obj.dynrModel@infile)
  
  return(obj.dynrModel)
  #modify the object slot, including starting values, etc.
}


