

/*#include <stdio.h>*/
#include <sys/types.h>
#include <errno.h>

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/*#include <armadillo>*/

#include "mainR.h"
#include "mainSAEM.h"

static R_NativePrimitiveArgType main_R_t[] = {
    VECSXP, VECSXP, LGLSXP, LGLSXP, LGLSXP, LGLSXP, LGLSXP
};


static R_CMethodDef cMethods[] = {
   {".BackendC", (DL_FUNC) &main_R, 7, main_R_t},
   {".BackendCS", (DL_FUNC) &main_SAEM, 7, main_R_t},
   {NULL, NULL, 0}
};


static R_CallMethodDef callMethods[] = {
	{".Backend", (DL_FUNC) main_R, 7},
	{".BackendS", (DL_FUNC) main_SAEM, 7},
	{NULL, NULL, 0}
};



void R_init_dynr(DllInfo *info) {
	R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
}






