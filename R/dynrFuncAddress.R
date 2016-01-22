#--------------------------------------------------
# dynr.funcaddresses
# Purpose: 
# takes in a C file or C scripts 
# returns a list of addresses of the compiled model functions and maybe R functions for debug purposes
#------------------------------------------------
# Changed DLL name and directory to be user-specified and permanent
dynr.funcaddress<-function(includes=character(), func_noise_cov=character(), verbose=TRUE,isContinuousTime, infile, outfile=tempfile()){

  #-------Set some variables: This function may later be extended----------
  language="C"
  #-------Check the input arguments----------------------------
  if(missing(infile)){
    #if(missing(func_noise_cov)){
    #  stop("The function of the noise covariance matrix is missing")
    #}
    #-------Generate the code-----------  
    code<-""#paste("#include <R.h>\n#include <Rdefines.h>\n","#include <R_ext/Error.h>\n", sep="")
    code<-paste(c(code,includes, ""), collapse="\n")
    code<-paste(c(code,func_noise_cov, ""), collapse="\n")
  }else{
    code<-readLines(infile)
  }
  # ---- Write and compile the code ----

  #filename<- basename(tempfile())
  libLFile <- CompileCode(code, language, verbose, outfile)
  #---- SET A FINALIZER TO PERFORM CLEANUP: register an R function to be called upon garbage collection of object or at the end of an R session---  
  #cleanup <- function(env) {
  #    if ( filename %in% names(getLoadedDLLs()) ) dyn.unload(libLFile)
  #    unlink(libLFile)
  #  }
  #reg.finalizer(environment(), cleanup, onexit=TRUE)
  #-----dynamically load the library-------
  DLL <- dyn.load( libLFile )  
  if (isContinuousTime==TRUE){
  res=list(f_measure=getNativeSymbolInfo("function_measurement", DLL)$address,
           f_dx_dt=getNativeSymbolInfo("function_dx_dt", DLL)$address,
           f_dF_dx=getNativeSymbolInfo("function_dF_dx", DLL)$address,
           f_dP_dt=getNativeSymbolInfo("function_dP_dt", DLL)$address,
           f_initial_condition=getNativeSymbolInfo("function_initial_condition", DLL)$address,
           f_regime_switch=getNativeSymbolInfo("function_regime_switch", DLL)$address,
           f_noise_cov=getNativeSymbolInfo("function_noise_cov", DLL)$address,
           f_transform=getNativeSymbolInfo("function_transform", DLL)$address)
  }else{
    res=list(f_measure=getNativeSymbolInfo("function_measurement", DLL)$address,
             f_dynamic=getNativeSymbolInfo("function_dynamic", DLL)$address,
             f_dF_dx=getNativeSymbolInfo("function_dF_dx", DLL)$address,
             f_dP_dt=null,
             f_initial_condition=getNativeSymbolInfo("function_initial_condition", DLL)$address,
             f_regime_switch=getNativeSymbolInfo("function_regime_switch", DLL)$address,
             f_noise_cov=getNativeSymbolInfo("function_noise_cov", DLL)$address,
             f_transform=getNativeSymbolInfo("function_transform", DLL)$address)    
  }
  return(res)
}    
#--------------------------------------------------
# CompileCode: A function adapted from the compileCode function in the inline pacakge
# Purpose: compiles a C file to create a shared library and returns its name.
#------------------------------------------------

CompileCode <- function(code, language, verbose, outfile) {
  wd = getwd()
  on.exit(setwd(wd))
  ## Prepare temp file names
  if ( .Platform$OS.type == "windows" ) {
    ## windows files
    outfile <- gsub("\\\\", "/", outfile)

    ## windows gsl flags
    LIB_GSL <- Sys.getenv("LIB_GSL")
    gsl_cflags <- sprintf( "-I%s/include", LIB_GSL )
    gsl_libs   <- sprintf( "-L%s/lib/%s -lgsl -lgslcblas", LIB_GSL, .Platform$r_arch)
  }
  else {
    ## UNIX-alike build

    ## Unix gsl flags
    gsl_cflags <- system( "gsl-config --cflags" , intern = TRUE )
    gsl_libs   <- system( "gsl-config --libs"   , intern = TRUE )
    
  }
  if (verbose) cat("Setting PKG_CPPFLAGS to", gsl_cflags, "\n")
  Sys.setenv(PKG_CPPFLAGS=gsl_cflags)
  if (verbose) cat("Setting PKG_LIBS to", gsl_libs, "\n")
  Sys.setenv(PKG_LIBS=gsl_libs)

  libCFile  <- paste(outfile, ".EXT", sep="")
  extension <- switch(language, "C++"=".cpp", C=".c", Fortran=".f", F95=".f95",
                      ObjectiveC=".m", "ObjectiveC++"=".mm")
  libCFile <- sub(".EXT$", extension, libCFile)
  libLFile  <- paste(outfile, .Platform$dynlib.ext, sep="")
  ## Write the code to the temp file for compilation
  write(code, libCFile)
  
  ## Compile the code using the running version of R if several available
  if ( file.exists(libLFile) ) file.remove( libLFile )
  
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
  return( libLFile )
}
