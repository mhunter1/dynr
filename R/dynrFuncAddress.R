#--------------------------------------------------
# .C2funcaddresses
# Purpose: 
# takes in a C file or C scripts 
# returns a list of addresses of the compiled model functions and maybe R functions for debug purposes
#------------------------------------------------
# Changed DLL name and directory to be user-specified and permanent
.C2funcaddress<-function(verbose,isContinuousTime, infile, outfile,compileLib){
	
	#-------Set some variables: This function may later be extended----------
	language <- "C"
	#-------Get the full name of the library----------
	if ( .Platform$OS.type == "windows" ) outfile <- gsub("\\\\", "/", outfile)
	libLFile  <- paste(outfile, .Platform$dynlib.ext, sep="")
	if (compileLib|(!file.exists(libLFile))){#when the compileLib flag is TRUE or when the libLFile does not exist
		#-------Check the input arguments----------------------------
		code <- readLines(infile)
		# ---- Write and compile the code ----
		
		#filename<- basename(tempfile())
		CompileCode(code, language, verbose, libLFile)
		#---- SET A FINALIZER TO PERFORM CLEANUP: register an R function to be called upon garbage collection of object or at the end of an R session---  
		#cleanup <- function(env) {
		#    if ( filename %in% names(getLoadedDLLs()) ) dyn.unload(libLFile)
		#    unlink(libLFile)
		#  }
		#reg.finalizer(environment(), cleanup, onexit=TRUE)
	}
	
	#-----dynamically load the library-------
	DLL <- dyn.load( libLFile )
	if (isContinuousTime){
		res <- list(f_measure=getNativeSymbolInfo("function_measurement", DLL)$address,
			f_dx_dt=getNativeSymbolInfo("function_dx_dt", DLL)$address,
			f_dF_dx=getNativeSymbolInfo("function_dF_dx", DLL)$address,
			f_dP_dt=getNativeSymbolInfo("function_dP_dt", DLL)$address,
			f_initial_condition=getNativeSymbolInfo("function_initial_condition", DLL)$address,
			f_regime_switch=getNativeSymbolInfo("function_regime_switch", DLL)$address,
			f_noise_cov=getNativeSymbolInfo("function_noise_cov", DLL)$address,
			f_transform=getNativeSymbolInfo("function_transform", DLL)$address)
	}else{
		res <- list(f_measure=getNativeSymbolInfo("function_measurement", DLL)$address,
			f_dynamic=getNativeSymbolInfo("function_dynam", DLL)$address,
			f_jacob_dynamic=getNativeSymbolInfo("function_jacob_dynam", DLL)$address,
			f_initial_condition=getNativeSymbolInfo("function_initial_condition", DLL)$address,
			f_regime_switch=getNativeSymbolInfo("function_regime_switch", DLL)$address,
			f_noise_cov=getNativeSymbolInfo("function_noise_cov", DLL)$address,
			f_transform=getNativeSymbolInfo("function_transform", DLL)$address)
	}
	return(list(address=res, libname=libLFile))
}

#--------------------------------------------------
# CompileCode: A function adapted from the compileCode function in the inline pacakge
# Purpose: compiles a C file to create a shared library
#------------------------------------------------

CompileCode <- function(code, language, verbose, libLFile) {
	wd <- getwd()
	on.exit(setwd(wd))
	## Prepare temp file names
	if ( .Platform$OS.type == "windows" ) {
		## windows files
		#outfile <- gsub("\\\\", "/", outfile)
		
		## windows gsl flags
		LIB_GSL <- Sys.getenv("LIB_GSL")
		LIB_GSL <- gsub("\\\\", "/", LIB_GSL) # replace "\" with "/"
		LIB_GSL <- gsub("\"", "", LIB_GSL) # remove "
		gsl_cflags <- sprintf( "-I\"%s/include\"", LIB_GSL)
		gsl_libs   <- sprintf( "-L\"%s/lib/%s\" -lgsl -lgslcblas", LIB_GSL, .Platform$r_arch)
	}else {
		## UNIX-alike build
		
		## Unix gsl flags
		gsl_cflags <- system( "gsl-config --cflags" , intern = TRUE )
		gsl_libs   <- system( "gsl-config --libs"   , intern = TRUE )
	}
	if (verbose) cat("Setting PKG_CPPFLAGS to", gsl_cflags, "\n")
	Sys.setenv(PKG_CPPFLAGS=gsl_cflags)
	if (verbose) cat("Setting PKG_LIBS to", gsl_libs, "\n")
	Sys.setenv(PKG_LIBS=gsl_libs)
	
	#libCFile  <- paste(outfile, ".EXT", sep="")
	libCFile <-sub(.Platform$dynlib.ext,".EXT",libLFile)
	extension <- switch(language, "C++"=".cpp", C=".c", Fortran=".f", F95=".f95",
						ObjectiveC=".m", "ObjectiveC++"=".mm")
	libCFile <- sub(".EXT$", extension, libCFile)
	#libLFile  <- paste(outfile, .Platform$dynlib.ext, sep="")
	
	## Write the code to the temp file for compilation
	write(code, libCFile)
	
	## Compile the code using the running version of R if several available
	#if ( file.exists(libLFile) ) {file.remove( libLFile )}
	setwd(dirname(libCFile))
	errfile <- paste( basename(libCFile), ".err.txt", sep = "" )
	cmd <- paste(R.home(component="bin"), "/R CMD SHLIB ", basename(libCFile), " ", gsl_libs, " ", " 2> ", errfile, sep="")
	if (verbose) cat("Compilation argument:\n", cmd, "\n")
	compiled <- system2(paste0(R.home(component="bin"), "/R"), args=c("CMD", "SHLIB", basename(libCFile)), stderr=errfile, stdout=verbose)
	errmsg <- readLines(errfile)
	unlink(errfile)
	if(length(errmsg) > 0){cat("May I present to you your error messages?\n")}
	writeLines(errmsg)
	setwd(wd)
	
	#### Error Messages
	if ( !file.exists(libLFile) | length(errmsg > 0) ) {
		cat("\nERROR(s) during compilation: source code errors or compiler configuration errors!\n")
		cat("\nProgram source:\n")
		codeWithLineNums <- paste(sprintf(fmt="%3d: ", 1:length(code)), code, sep="")
		writeLines(codeWithLineNums)
		stop("Compilation ERROR, function(s)/method(s) not created!")
	}
  #return( libLFile )
}

#------------------------------------------------------------------------------
# Check configuration

##' Check that dynr in configured properly
##' 
##' @param verbose logical.  Whether to print messages during/after checks
##' 
##' @details
##' The 'dynr' package requires additional set-up and configuration beyond
##' just installing the package.  In particular, it requires compiling C code
##' along with GSL to run (cook) models.  This function runs some basic checks
##' of the configuration.  We check that (1) R is on the PATH variable, (2)
##' Rtools exists and is on the PATH variable for Windows, (3) a C compiler
##' is available, and (4) GSL is available and on the PATH.
##' 
##' In general, see the 'Installation for Users' vignette for set-up and
##' configuration instructions.
##' 
##' @return No return value.
##' 
##' @examples
##' dynr.config()
dynr.config <- function(verbose=FALSE){
	genmsg <- paste0(
		"\nPlease read the 'Installation for Users' vignette at\n",
		"https://cran.r-project.org/web/packages/dynr/",
		"vignettes/InstallationForUsers.pdf",
		"\nor\n",
		"vignette(package='dynr', 'InstallationForUsers')\n")
	# Check that R is on the path
	noRmsg <- "R did not appear to be on the 'PATH' environment variable."
	# Check that Rtools exists and is on the path
	noRtoolsmsg <- "No Rtools found in 'PATH' environment variable."
	# Check for a C compiler
	noCmsg <- "No C compiler found."
	# Check for GSL
	noGSLmsg <- "LIB_GSL variable not found."
	if ( .Platform$OS.type == "windows" ) {
		path <- Sys.getenv("PATH")
		path <- gsub("\\\\", "/", path)
		if(grep(R.home(component="bin"), path) != 1){
			stop(paste0(noRmsg, genmsg))
		}
		if(grep('rtools', path, ignore.case=TRUE) != 1){
			stop(paste0(noRtoolsmsg, genmsg))
		}
		if(!checkForCompiler() ){
			stop(paste0(noCmsg, genmsg))
		}
		
		## windows gsl flags
		LIB_GSL <- Sys.getenv("LIB_GSL")
		LIB_GSL <- gsub("\\\\", "/", LIB_GSL) # replace "\" with "/"
		LIB_GSL <- gsub("\"", "", LIB_GSL) # remove "
		if(nchar(LIB_GSL) < 1){
			stop(paste0(noGSLmsg, genmsg))
		}
	}else {
		## UNIX-alike build
		
		if(!checkForCompiler() ){
			stop(paste0(noCmsg, genmsg))
		}
		
		## Unix gsl flags
		gsl_cflags <- system( "gsl-config --cflags" , intern = TRUE )
		gsl_libs   <- system( "gsl-config --libs"   , intern = TRUE )
		if(0){ # TODO fill with something!
			stop(paste0(noGSLmsg, genmsg))
		}
	}
	if(verbose){
		message("Configuration check complete.  Ready to rock and roll.")
	}
}


# Taken from Rcpp
checkForCompiler <- function(minVersion = package_version("4.6.0")){
	binaries <- c("g++", Sys.getenv("CXX", unset = ""),
		Sys.getenv("CXX1X", unset = ""))
	binpaths <- lapply(binaries, function(b) {
		if (b == "")
			NULL
		else Sys.which(b)
	})
	allgood <- FALSE
	rl <- lapply(binpaths, function(b) {
		if (is.null(b)) 
			return(NULL)
		con <- pipe(paste(b, "-v 2>&1"), "r")
		lines <- readLines(con)
		close(con)
		lines <- lines[grepl("^g.. version", lines)]
		if (length(lines) == 0)
			return(NULL)
		ver <- strsplit(lines, " ")[[1]][3]
		package_version(ver) >= minVersion
	})
	all(do.call(c, rl))
}
