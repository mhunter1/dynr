/*==============================================================================
# NLOPT CMake configuration file
#
# NLopt is a free/open-source library for nonlinear optimization, providing
# a common interface for a number of different free optimization routines
# available online as well as original implementations of various other
# algorithms
# WEBSITE: http://ab-initio.mit.edu/wiki/index.php/NLopt
# AUTHOR: Steven G. Johnson
#
# This config.cmake.h.in file was created to compile NLOPT with the CMAKE utility.
# Benoit Scherrer, 2010 CRL, Harvard Medical School
# Copyright (c) 2008-2009 Children's Hospital Boston
#
# Minor changes to the source was applied to make possible the compilation with
# Cmake under Linux/Win32
#============================================================================*/

/* Bugfix version number. */
#define BUGFIX_VERSION 2

/* Define to enable extra debugging code. */
#undef DEBUG

/* Define to 1 if you have the `BSDgettimeofday' function. */
#undef HAVE_BSDGETTIMEOFDAY


/* Define to 1 if you have the `m' library (-lm). */
#undef HAVE_LIBM

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#undef LT_OBJDIR

/* Major version number. */
#define MAJOR_VERSION 2

/* Minor version number. */
#define MINOR_VERSION 6

/* Name of package */
#undef PACKAGE

/* Define to the address where bug reports for this package should be sent. */
#undef PACKAGE_BUGREPORT

/* Define to the full name of this package. */
#undef PACKAGE_NAME

/* Define to the full name and version of this package. */
#undef PACKAGE_STRING

/* Define to the one symbol short name of this package. */
#undef PACKAGE_TARNAME

/* Define to the home page for this package. */
#undef PACKAGE_URL

/* Define to the version of this package. */
#undef PACKAGE_VERSION

/* replacement for broken HUGE_VAL macro, if needed */
#undef REPLACEMENT_HUGE_VAL

/* The size of `unsigned int', as computed by sizeof. */
//#define SIZEOF_UNSIGNED_INT @SIZEOF_UNSIGNED_INT@

/* The size of `unsigned long', as computed by sizeof. */
//#define SIZEOF_UNSIGNED_LONG @SIZEOF_UNSIGNED_LONG@

/* Define to 1 if you have the ANSI C header files. */
#undef STDC_HEADERS

/* Define to C thread-local keyword, or to nothing if this is not supported in
   your compiler. */
#define THREADLOCAL 


/* Version number of package */
#undef VERSION



/* Define to empty if `const' does not conform to ANSI C. */
#undef const

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
#undef inline
#endif
