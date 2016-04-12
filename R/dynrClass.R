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

