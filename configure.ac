AC_INIT(dynr, 0.1, mhunter.ou@gmail.com)

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

## Use gsl-config to find arguments for compiler and linker flags
##
## Check for non-standard programs: gsl-config(1)
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
## If gsl-config was found, let’s use it
if test "${GSL_CONFIG}" != ""; then
# Use gsl-config for header and linker arguments (without BLAS which we get from R)
GSL_CFLAGS='${GSL_CONFIG} --cflags'
GSL_LIBS='${GSL_CONFIG} --libs-without-cblas'
else
AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi
# Now substitute these variables in src/Makevars.in to create src/Makevars
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)
AC_CONFIG_FILES([
  src/Makevars
])
AC_OUTPUT
