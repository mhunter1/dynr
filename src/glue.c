

#include <stdio.h>
#include <sys/types.h>
#include <errno.h>

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "estimation_nloptR.h"


static R_NativePrimitiveArgType main_R_t[] = {
    VECSXP, VECSXP
};

static R_CMethodDef cMethods[] = {
   {"dynrBackend", (DL_FUNC) &main_R, 2, main_R_t},
   {NULL, NULL, 0}
};


static R_CallMethodDef callMethods[] = {
	{"dynrBackend", (DL_FUNC) main_R, 2},
	{NULL, NULL, 0}
};


void R_init_dynr(DllInfo *info) {
	R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
}






