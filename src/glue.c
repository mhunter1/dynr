
#include <stdio.h>
#include <sys/types.h>
#include <errno.h>

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static R_CallMethodDef callMethods[] = {
	{"main_R", (DL_FUNC) main_R, 2}
};

#ifdef  __cplusplus
extern "C" {
#endif

void R_init_dynr(DllInfo *info) {
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

void R_unload_dynr(DllInfo *) {
	/* keep this stub in case we need it */
}

#ifdef  __cplusplus
}
#endif

