# Preping is what happens to recipes.
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


