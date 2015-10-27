#--------------------------------------------------
# compileCode: A function adapted from the compileCode function in the inline pacakge
# Purpose: compiles a C file to create a shared library and returns its name.
#------------------------------------------------
compileCode <- function(f, code, language, verbose) {
  wd = getwd()
  on.exit(setwd(wd))
  ## Prepare temp file names
  if ( .Platform$OS.type == "windows" ) {
    ## windows files
    dir <- gsub("\\\\", "/", tempdir())
    libCFile  <- paste(dir, "/", f, ".EXT", sep="")
    libLFile  <- paste(dir, "/", f, ".dll", sep="")
    libLFile2 <- paste(dir, "/", f, ".dll", sep="")
    ## windows gsl flags
    LIB_GSL <- Sys.getenv("LIB_GSL")
    gsl_cflags <- sprintf( "-I%s/include", LIB_GSL )
    gsl_libs   <- sprintf( "-L%s/lib -lgsl -lgslcblas", LIB_GSL )
  }
  else {
    ## UNIX-alike build
    libCFile  <- paste(tempdir(), "/", f, ".EXT",               sep="")
    libLFile  <- paste(tempdir(), "/", f, .Platform$dynlib.ext, sep="")
    libLFile2 <- paste(tempdir(), "/", f, ".sl",                sep="")
    
    ## Unix gsl flags
    gsl_cflags <- system( "gsl-config --cflags" , intern = TRUE )
    gsl_libs   <- system( "gsl-config --libs"   , intern = TRUE )
    
  }
  extension <- switch(language, "C++"=".cpp", C=".c", Fortran=".f", F95=".f95",
                      ObjectiveC=".m", "ObjectiveC++"=".mm")
  libCFile <- sub(".EXT$", extension, libCFile)
  
  ## Write the code to the temp file for compilation
  write(code, libCFile)
  
  ## Compile the code using the running version of R if several available
  if ( file.exists(libLFile) ) file.remove( libLFile )
  if ( file.exists(libLFile2) ) file.remove( libLFile2 )
  
  setwd(dirname(libCFile))
  errfile <- paste( basename(libCFile), ".err.txt", sep = "" )
  cmd <- paste(R.home(component="bin"), "/R CMD SHLIB ", basename(libCFile), " 2> ", errfile, sep="")
  if (verbose) cat("Compilation argument:\n", cmd, "\n")
  compiled <- system(cmd, intern=!verbose)
  errmsg <- readLines( errfile )
  unlink( errfile )
  writeLines( errmsg )
  setwd(wd)
  
  if ( !file.exists(libLFile) && file.exists(libLFile2) ) libLFile <- libLFile2
  if ( !file.exists(libLFile) ) {
    cat("\nERROR(s) during compilation: source code errors or compiler configuration errors!\n")
    cat("\nProgram source:\n")
    code <- strsplit(code, "\n")
    for (i in 1:length(code[[1]])) cat(format(i,width=3), ": ", code[[1]][i], "\n", sep="")
    stop( paste( "Compilation ERROR, function(s)/method(s) not created!", paste( errmsg , collapse = "\n" ) ) )
  }
  return( libLFile )
}

#------Compile the Model script in C-------
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## CFunc is an S4 class derived from 'function'. This inheritance allows objects
## to behave exactly as functions do, but it provides a slot @code that keeps the
## source C or Fortran code used to create the inline call
setClass("CFunc",
         representation(
           code="character"
         ),
         contains="function"
)

setClass( "CFuncList", contains = "list" )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cfunction <- function(RArgs.List=character(), body=character(), includes=character(),
                      language="C",verbose=FALSE, convention=".C") {
  
  f <- basename(tempfile())
  
  if ( !is.list(Arg.List) ) {
    Arg.List <- list(Arg.List)
    CfileNameWoExt <- f
    names(body) <- f
  }
  if( length(Arg.List) != length(body) )
    stop("mismatch between the number of functions declared in 'Arg.List' and the number of function bodies provided in 'body'")
  
  types <- vector(mode="list", length=length(Arg.List))
  ## ---- GENERATE THE CODE ----
  
  ## ---- WRITE AND COMPILE THE CODE ----
  libLFile <- compileCode(f, code, language, verbose)
  
  ## SET A FINALIZER TO PERFORM CLEANUP
  #register an R function to be called at the end of an R session
  cleanup <- function(env) {
    if ( f %in% names(getLoadedDLLs()) ) dyn.unload(libLFile)
    unlink(libLFile)
  }
  reg.finalizer(environment(), cleanup, onexit=TRUE)
  
  res <- vector("list", length(Arg.List))
  names(res) <- CfileNameWoExt
  
  ## GENERATE R FUNCTIONS
  for ( i in seq_along(Arg.List) ) {
    ## Create new objects of class CFunc, each containing the code of ALL inline
    ## functions. This will be used to recompile the whole shared lib when needed
    res[[i]] <- new("CFunc", code = code)
    
    ## this is the skeleton of the function, the external call is added below using 'body'
    ## important here: all variables are kept in the local environment
    fn <- function(arg) {
      NULL
    }
    
    DLL <- dyn.load( libLFile )
    
    ## Modify the function formals to give the right argument list
    args <- formals(fn)[ rep(1, length(Arg.List[[i]])) ]
    names(args) <- names(Arg.List[[i]])
    formals(fn) <- args
    
    ## create .C/.Call function call that will be added to 'fn'
    if (convention == ".Call") {
      body <- quote( CONVENTION("EXTERNALNAME", ARG) )[ c(1:2, rep(3, length(Arg.List[[i]]))) ]
      for ( j in seq(along = Arg.List[[i]]) ) body[[j+2]] <- as.name(names(Arg.List[[i]])[j])
    }
    else {
      body <- quote( CONVENTION("EXTERNALNAME", as.logical(ARG), as.integer(ARG),
                                as.double(ARG), as.complex(ARG), as.character(ARG),
                                as.raw(ARG), as.double(ARG)) )[ c(1:2,types[[i]]+2) ]
      names(body) <- c( NA, "", names(Arg.List[[i]]) )
      for ( j in seq(along = Arg.List[[i]]) ) body[[j+2]][[2]] <- as.name(names(Arg.List[[i]])[j])
    }
    body[[1]] <- get(convention)
    body[[2]] <- getNativeSymbolInfo("function_noise_cov", DLL)$address
    ## update the body of 'fn'
    body(fn) <- body
    ## set fn as THE function in CFunc of res[[i]]
    res[[i]]@.Data <- fn
  }
  
  ## OUTPUT PROGRAM CODE IF DESIRED
  if ( verbose ) {
    cat("Program source:\n")
    lines <- strsplit(code, "\n")
    for ( i in 1:length(lines[[1]]) )
      cat(format(i,width=3), ": ", lines[[1]][i], "\n", sep="")
  }
  
  ## Remove unnecessary objects from the local environment
  remove(list = c("args", "body", "fn", "funCArg.List", "i", "includes", "j"))
  
  ## RETURN THE FUNCTION
  if (length(res) == 1 && names(res) == f) return( res[[1]] )
  else return( new( "CFuncList", res ) )
}

