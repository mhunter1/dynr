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
			f_noise_cov=getNativeSymbolInfo("function_noise_cov", DLL)$address,f_transform=getNativeSymbolInfo("function_transform", DLL)$address,
			f_transform=getNativeSymbolInfo("function_transform", DLL)$address#,
			#f_test=getNativeSymbolInfo("function_arma_hello_world", DLL)$addresses
			)
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

#--------------------------------------------------
# .C2funcaddressSAEM 
# Purpose: 
# takes in a C file or C scripts 
# returns a list of addresses of the compiled model functions and maybe R functions for debug purposes
# simialr with .C2funcaddress but compiles functions for SAEM (not used now)
#------------------------------------------------
# Changed DLL name and directory to be user-specified and permanent
.C2funcaddressSAEM<-function(verbose,isContinuousTime, infile, outfile,compileLib){
  
  #-------Set some variables: This function may later be extended----------
  language="C++"
  #-------Get the full name of the library----------
  if ( .Platform$OS.type == "windows" ) outfile <- gsub("\\\\", "/", outfile)
  libLFile  <- paste(outfile, .Platform$dynlib.ext, sep="")
  if (compileLib|(!file.exists(libLFile))){#when the compileLib flag is TRUE or when the libLFile does not exist
    #-------Check the input arguments----------------------------
    code<-readLines(infile)
    # ---- Write and compile the code ----
    
    #filename<- basename(tempfile())
    CompileCodeSAEM(code, language, verbose, libLFile)
    #---- SET A FINALIZER TO PERFORM CLEANUP: register an R function to be called upon garbage collection of object or at the end of an R session---  
    #cleanup <- function(env) {
    #    if ( filename %in% names(getLoadedDLLs()) ) dyn.unload(libLFile)
    #    unlink(libLFile)
    #  }
    #reg.finalizer(environment(), cleanup, onexit=TRUE)
    print('CompileCodeSAEM done')
    print(names(getLoadedDLLs()))
    print('--------------------')
  }
  
  #-----dynamically load the library-------
  DLL <- dyn.load( libLFile )
  if (isContinuousTime){
    res=list(
      f_test  = getNativeSymbolInfo("function_hello3", DLL)$address,
      f_dyn   = getNativeSymbolInfo("dynfunICM", DLL)$address,
      f_dfdx  = getNativeSymbolInfo("dfdxFreeICM", DLL)$address,
      f_dfdp  = getNativeSymbolInfo("dfdparFreeIC", DLL)$address,
      f_dfdx2 = getNativeSymbolInfo("dfdx2FreeIC", DLL)$address,
      f_dfdxdp= getNativeSymbolInfo("dfdxdpFreeIC", DLL)$address,
      f_dfdpdx= getNativeSymbolInfo("dfdpdxFreeIC", DLL)$address,
      f_dfdp2 = getNativeSymbolInfo("dfdpar2FreeIC", DLL)$address,
      f_setpars = getNativeSymbolInfo("setParsFreeICwb", DLL)$address)
  }else{
    # todo: Implement the time point regularization and corresponding model modification to support isContinuousTime = FALSE
    warning('SAEM process only supports isContinuousTime = TRUE now.')
    # res=list(
    # f_test  = getNativeSymbolInfo("function_hello3", DLL)$address,
    # f_dyn   = getNativeSymbolInfo("dynfunICM", DLL)$address,
    # f_dfdx  = getNativeSymbolInfo("dfdxFreeICM", DLL)$address,
    # f_dfdp  = getNativeSymbolInfo("dfdparFreeIC", DLL)$address,
    # f_dfdx2 = getNativeSymbolInfo("dfdx2FreeIC", DLL)$address,
    # f_dfdxdp= getNativeSymbolInfo("dfdxdpFreeIC", DLL)$address,
    # f_dfdpdx= getNativeSymbolInfo("dfdpdxFreeIC", DLL)$address,
    # f_dfdp2 = getNativeSymbolInfo("dfdpar2FreeIC", DLL)$address)   
  }
  print(res)
  .Call(res$f_test, as.double(seq(.1, .5, by = 0.1)))
  return(list(address=res, libname=libLFile))
}

#--------------------------------------------------
# CompileCodeSAEM: A function adapted from the compileCode function in the inline pacakge
# Purpose: compiles a C file to create a shared library (armadillo version, not used now)
#------------------------------------------------
CompileCodeSAEM <- function(code, language, verbose, libLFile) {
  wd = getwd()
  on.exit(setwd(wd))
  ## Prepare temp file names
  if ( .Platform$OS.type == "windows" ) {
    ## windows files
    #outfile <- gsub("\\\\", "/", outfile)

    ## windows gsl flags
	LIB_ARMADILLO <- Sys.getenv("LIB_ARMADILLO") 
	#print(LIB_ARMADILLO)
    gsl_cflags <- sprintf( "-I %s/include", LIB_ARMADILLO)
    gsl_libs   <- sprintf( "-L%s/lib_win64 -lblas -llapack ", LIB_ARMADILLO)
  }else {
    ## UNIX-alike build

    ## Unix gsl flags
    gsl_cflags <- sprintf("-larmadillo")
    #gsl_cflags <- system( "-larmadillo" , intern = TRUE )
    #gsl_libs   <- system( ""   , intern = TRUE )
    
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
  cmd <- paste(R.home(component="bin"), "/R CMD SHLIB ", basename(libCFile), " 2> ", errfile, sep="")
  if (verbose) cat("Compilation argument:\n", cmd, "\n")
  compiled <- system(cmd, intern=!verbose)
  errmsg <- readLines(errfile)
  unlink(errfile)
  writeLines(errmsg)
  setwd(wd)
  
  #### Error Messages
  if ( !file.exists(libLFile) ) {
    cat("\nERROR(s) during compilation: source code errors or compiler configuration errors!\n")
    cat("\nProgram source:\n")
    code <- strsplit(code, "\n")
    for (i in 1:length(code[[1]])) cat(format(i,width=3), ": ", code[[1]][i], "\n", sep="")
    stop( paste( "Compilation ERROR, function(s)/method(s) not created!", paste( errmsg , collapse = "\n" ) ) )
    }
  #return( libLFile )
}

#--------------------------------------------------
# CompileCodeSAEMRCpp: A function adapted from the compileCode function in the inline pacakge
# Purpose: compiles a C file to create a shared library (RcppArmadiilo version)
#------------------------------------------------
CompileCodeSAEMRCpp <- function(code, language, verbose, libLFile) {
  #browser()
  wd = getwd()
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
		
		
		arma_libs <- paste(Sys.getenv("SHLIB_OPENMP_CXXFLAGS"), Sys.getenv("LAPACK_LIBS"), Sys.getenv("BLAS_LIBS"), Sys.getenv("FLIBS"))
	}else {
		## UNIX-alike build
		
		## Unix gsl flags
		gsl_cflags <- system( "gsl-config --cflags" , intern = TRUE )
		gsl_libs   <- system( "gsl-config --libs"   , intern = TRUE )
		arma_libs <- ""
	}
	if (verbose) cat('Setting PKG_CPPFLAGS to', gsl_cflags, '-I"$(RCPP_LIBS)/Rcpp/include" -I"$(RCPP_LIBS)/RcppArmadillo/include" -std=c++11')
	Sys.setenv(PKG_CPPFLAGS=paste(gsl_cflags, '-I"$(RCPP_LIBS)/Rcpp/include" -I"$(RCPP_LIBS)/RcppArmadillo/include" -std=c++11'))
	if (verbose) cat("Setting PKG_LIBS to", gsl_libs, arma_libs, "-fopenmp", paste0("-L",Sys.getenv("LIB_ARMADILLO"), "/lib_win64 -lblas -llapack"), "\n")
	#Sys.setenv(PKG_LIBS=paste(gsl_libs, arma_libs, "-fopenmp"))
	Sys.setenv(PKG_LIBS = paste(gsl_libs, arma_libs, "-fopenmp", paste0("-L",Sys.getenv("LIB_ARMADILLO"), "/lib_win64 -lblas -llapack")))

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
  cmd <- paste(R.home(component="bin"), "/R CMD SHLIB ", basename(libCFile), " 2> ", errfile, sep="")
  if (verbose) cat("Compilation argument:\n", cmd, "\n")
  compiled <- system(cmd, intern=!verbose)
  errmsg <- readLines(errfile)
  unlink(errfile)
  writeLines(errmsg)
  setwd(wd)
  
  #browser()
  #### Error Messages
  if ( !file.exists(libLFile) ) {
    cat("\nERROR(s) during compilation: source code errors or compiler configuration errors!\n")
    cat("\nProgram source:\n")
    code <- strsplit(code, "\n")
    for (i in 1:length(code[[1]])) cat(format(i,width=3), ": ", code[[1]][i], "\n", sep="")
    stop( paste( "Compilation ERROR, function(s)/method(s) not created!", paste( errmsg , collapse = "\n" ) ) )
    }
  return( libLFile )
}
#--------------------------------------------------
# .C2funcaddressSAEMRcpp 
# Purpose: 
# takes in a C file or C scripts 
# returns a list of addresses of the compiled model functions and maybe R functions for debug purposes
# simialr with .C2funcaddress but compiles functions for SAEM
#------------------------------------------------
.C2funcaddressSAEMRcpp<-function(verbose,isContinuousTime, infile, outfile,compileLib){
  
  #-------Set some variables: This function may later be extended----------
  language="C++"
  #-------Get the full name of the library----------
  if ( .Platform$OS.type == "windows" ) outfile <- gsub("\\\\", "/", outfile)
  libLFile  <- paste(outfile, .Platform$dynlib.ext, sep="")
  if (compileLib|(!file.exists(libLFile))){#when the compileLib flag is TRUE or when the libLFile does not exist
    #-------Check the input arguments----------------------------
    code<-readLines(infile)
    # ---- Write and compile the code ----
    
    #filename<- basename(tempfile())
    CompileCodeSAEMRCpp(code, language, verbose, libLFile)
    #---- SET A FINALIZER TO PERFORM CLEANUP: register an R function to be called upon garbage collection of object or at the end of an R session---  
    #cleanup <- function(env) {
    #    if ( filename %in% names(getLoadedDLLs()) ) dyn.unload(libLFile)
    #    unlink(libLFile)
    #  }
    #reg.finalizer(environment(), cleanup, onexit=TRUE)
  }
  
  #-----dynamically load the library-------
  #browser()
  DLL <- dyn.load( libLFile )
  if (isContinuousTime){
    res=list(
      f_test  = getNativeSymbolInfo("hello_world", DLL)$address,
      f_dyn   = getNativeSymbolInfo("dynfunICM", DLL)$address,
      f_dfdx  = getNativeSymbolInfo("dfdxFreeICM", DLL)$address,
      f_dfdp  = getNativeSymbolInfo("dfdparFreeIC", DLL)$address,
      f_dfdx2 = getNativeSymbolInfo("dfdx2FreeIC", DLL)$address,
      f_dfdxdp= getNativeSymbolInfo("dfdxdpFreeIC", DLL)$address,
      f_dfdpdx= getNativeSymbolInfo("dfdpdxFreeIC", DLL)$address,
      f_dfdp2 = getNativeSymbolInfo("dfdpar2FreeIC", DLL)$address,
      f_setpars = getNativeSymbolInfo("setParsFreeICwb", DLL)$address)
  }else{
    # todo: Implement the time point regularization and corresponding model modification to support isContinuousTime = FALSE
    warning('SAEM process only supports isContinuousTime = TRUE now.')
    # res=list(
    # f_test  = getNativeSymbolInfo("function_hello3", DLL)$address,
    # f_dyn   = getNativeSymbolInfo("dynfunICM", DLL)$address,
    # f_dfdx  = getNativeSymbolInfo("dfdxFreeICM", DLL)$address,
    # f_dfdp  = getNativeSymbolInfo("dfdparFreeIC", DLL)$address,
    # f_dfdx2 = getNativeSymbolInfo("dfdx2FreeIC", DLL)$address,
    # f_dfdxdp= getNativeSymbolInfo("dfdxdpFreeIC", DLL)$address,
    # f_dfdpdx= getNativeSymbolInfo("dfdpdxFreeIC", DLL)$address,
    # f_dfdp2 = getNativeSymbolInfo("dfdpar2FreeIC", DLL)$address)   
  }
  
  #browser()
  return(list(address=res, libname=libLFile))
}
