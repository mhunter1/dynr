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
           transform="dynrTrans",
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

#TODO: Print RS model and initial condition
setMethod("printex", "dynrModel",
function(object, show=TRUE){
            model2<-PopBackModel(object, object$param.names)
            inlist <- list(model2$dynamics, model2$measurement, model2$noise, model2$initial, model2$regimes)
            outlist <- lapply(inlist, printex,show)
            

            #Gather dynamic model
            dynequ="\\begin{eqnarray}\n"
            if (class(model2$dynamics) == 'dynrDynamicsFormula'){
              exp1 <- .formulatoTex(outlist[[1]],model2)
              exp1 <- lapply(exp1,function(x){gsub('$','',x,fixed=TRUE)})
              for (j in 1:length(exp1)){
                a = NULL
                space=ifelse(j<length(exp1),"\\\\","")
                if (length(exp1)>1){a = paste0("\\text{Regime ",j,":}\\\\\n")}
                dynequ <- c(dynequ,paste0(a,exp1[j],space,"\n"))
              }#loop through regimes
            }else{

              for (j in 1:length(outlist[[1]]$dyn_tran)){
                a <- NULL
                space=ifelse(j<length(outlist[[1]]$dyn_tran),"\\\\","")
                if (length(outlist[[1]]$dyn_tran) > 1) {
                  a = paste0("\\text{Regime ",j,":}\\\\\n")}
                
                if ((model2$dynamics)$isContinuousTime){
                  b=paste0(.xtableMatrix(matrix(paste0("d",(model2$measurement)$state.names,"(t)"),ncol=1),F))
                  b1=paste0(.xtableMatrix(matrix(paste0((model2$measurement)$state.names,"(t)"),ncol=1),F))
                  b2="dt"
                }else{
                  exo=NULL
                  dint=NULL
                  if (length((model2$dynamics)$covariates) > 0){
                    exo=paste0("+",outlist[[1]]$dyn_exo[[j]],outlist[[1]]$dyn_exo.names)  
                  }
                  if (length((model2$dynamics)$values.int) > 0){
                    dint=paste0(outlist[[1]]$dyn_int[[j]],"+")  
                  }
                  b=.xtableMatrix(matrix(paste0((model2$measurement)$state.names,"(t+1)"),ncol=1),F)
                  b1=.xtableMatrix(matrix(paste0((model2$measurement)$state.names,"(t)"),ncol=1),F)
                  b2=NULL
                }
                dynequ<-c(dynequ,paste0(a,b,"=",dint,
                                          outlist[[1]]$dyn_tran[[j]],b1,b2,
                                          exo,
                                          "+ w\\\\",
                                          "w\\sim N(",
                                          .xtableMatrix(matrix(rep(0,length((model2$measurement)$state.names)),ncol=1),F),
                                          ",",outlist[[3]]$dynamic.noise[[1]],")",space,"\n")) #TODO: REPLACE 1 W j AFTER DYNRNOISE BECOMES LIST
              }#Loops through regimes
            }#If continous-time loop
            dynequ=c(dynequ,"\\end{eqnarray}")
            
            #Gather measurement model
            measequ="\\begin{eqnarray}\n"
            for (j in 1:length((model2$measurement)$values.load)){
              space=ifelse(j<length((model2$measurement)$values.load),"\\\\","")
              a <- NULL
              if (length((model2$measurement)$values.load) > 1) {
                a = paste0("\\text{Regime ",j,":}\\\\\n")}
              exom=NULL
              mint=NULL
              if (length((model2$measurement)$exo.names) > 0){
                exom=paste0("+",outlist[[2]]$meas_exo[[j]],outlist[[2]]$meas_exo.names)  
              }
              if (length((model2$measurement)$values.int) > 0){
                mint=paste0(outlist[[2]]$meas_int[[j]])  
              }
              measequ<-c(measequ,paste0(a,
                                        (model2$measurement)$obs.names,"=",mint,
                                        "+",outlist[[2]]$meas_loadings[[j]],
                                        b1,exom,
                                        "+ e\\\\",
                                        "e\\sim N(",
                                        .xtableMatrix(matrix(rep(0,length((model2$measurement)$obs.names)),ncol=1),F),
                                        ",",outlist[[3]]$measurement.noise[[1]],")",space,"\n")) #TODO: REPLACE 1 W j AFTER DYNRNOISE BECOMES LIST
            }
            measequ=c(measequ,"\\end{eqnarray}")
            
            cat("\nHere is the measurement model:\n")
            cat(measequ) 
            cat("\n\nHere is the dynamic model:\n")
            cat(dynequ)
            
          }
)


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
##' @param ... additional arguments specifying other dynrRecipe objects. Argument regimes is for 
##' a dynrRegimes object prepared with \code{\link{prep.regimes}} and argument transform is for 
##' a dynrTrans object prepared with \code{\link{prep.tfun}}.
##' @param infile a character string of the name of the input C script of model functions to be compiled 
##' for parameter estimation. The default is a temporary file name.
##' @param outfile a character string of the name of the output C script of model functions to be compiled 
##' for parameter estimation.
##' 
##' @examples
##' model <- dynr.model(dynamics=dynm, measurement=meas,noise=mdcov, initial=initial, 
##' regimes=regime, transform=trans,outfile="Recipe")
dynr.model <- function(dynamics, measurement, noise, initial, ..., infile=tempfile(), outfile){
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
    inputs$transform<-createRfun(prep.tfun(), param.data, 
                                 params.observed=inputs$noise$params.observed, params.latent=inputs$noise$params.latent, params.inicov=inputs$initial$params.inicov,
                                 values.observed=inputs$noise$values.observed.inv.ldl, values.latent=inputs$noise$values.latent.inv.ldl, values.inicov=inputs$initial$values.inicov.inv.ldl)
  }
  # paramName2Number on each recipe (this changes are the params* matrices to contain parameter numbers instead of names
  inputs <- sapply(inputs, paramName2Number, names=param.data$param.name)
  
  # writeCcode on each recipe
  inputs <- sapply(inputs, writeCcode) 
  all.values <- unlist(sapply(inputs, slot, name='startval'))
  unique.values <- extractValues(all.values, all.params)
  
  if(length(inputs$transform$formula.inv)>0){
    unique.values<-inputs$transform$inv.tfun(unique.values)
  }
  param.data$param.value=unique.values
  
  #initiate a dynrModel object
  obj.dynrModel <- new("dynrModel", c(list(infile=infile, outfile=outfile, param.names=as.character(param.data$param.name)), inputs))
  obj.dynrModel@dim_latent_var <- dim(inputs$noise$params.latent)[1]
  
  obj.dynrModel@xstart <- param.data$param.value
  obj.dynrModel@ub<-rep(9999,length(obj.dynrModel@xstart))
  obj.dynrModel@lb<-rep(9999,length(obj.dynrModel@xstart))
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
    glom <- paste(glom, prep.dP_dt, sep="\n\n")
  }
  cat(glom, file=obj.dynrModel@outfile)
  
  return(obj.dynrModel)
  #modify the object slot, including starting values, etc.
}


