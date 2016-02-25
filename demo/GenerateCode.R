#-------Generate the code-----------
#if(missing(func_noise_cov))
#stop("The function of the noise covariance matrix is missing")
#code<-readLines("demo/LinearSDE.c")
code<-"#include <math.h>\n#include <gsl/gsl_matrix.h>\n#include <gsl/gsl_blas.h>\n#include <stdio.h>"                                                                                                                                   
#   [6] ""
#   [7] "/**"
#   [8] " * The measurement function"
#   [9] " */"
#  [10] "void function_measurement(size_t t, size_t regime, double *param, const gsl_vector *eta, const gsl_vector *co_variate, gsl_matrix *Ht,gsl_vector *y){"
#  [11] ""
#  [12] "    gsl_matrix_set(Ht,0,0,1.0);"
#  [13] "    gsl_matrix_set(Ht,0,1,0.0);"
#  [14] "\t"
#  [15] "\tgsl_blas_dgemv(CblasNoTrans, 1.0, Ht, eta, 0.0, y);"
#  [16] ""
#  [17] "}"
if (isContinuousTime){
	#  [19] "/**\n * The dx/dt function\n */"
	#  [23] "void function_dx_dt(double t, size_t regime, const gsl_vector *x, double *param, size_t n_param, const gsl_vector *co_variate, gsl_vector *F_dx_dt){"
	#  [24] ""
	if (isLinear){
		#  [79] "\t/*size_t nrow=F_dx_dt_dx->size1, ncol=F_dx_dt_dx->size2,index_col,index_row;"
		
		#  [26] "\tsize_t nrow=F_dx_dt->size, ncol=x->size;"
		#  [27] "\t/*size_t index_col,index_row;*/"
		#  [28] "\t"
		#  [29] "\tgsl_matrix *A=gsl_matrix_calloc(nrow,ncol);"
		#  [30] "\tgsl_matrix_set(A,0,0,0.0);"
		#  [31] "\tgsl_matrix_set(A,0,1,1.0);"
		#  [32] "\tgsl_matrix_set(A,1,0,param[0]);"
		#  [33] "\tgsl_matrix_set(A,1,1,param[1]);"
		function_gsl_matrix_set(Pattern, StartVal, Fit, MatrixName="A")		
		#  [40] "\tgsl_blas_dgemv(CblasNoTrans, 1.0, A, x, 0.0, F_dx_dt);"
		#  [41] "\tgsl_matrix_free(A);"
	}else{#NonLinear
		#XXXXXXXXX
	}
	#  [42] "\t"
	if (CoVariate){
		if (isLinear){
			#  [44] "\t\tncol=co_variate->size;"
			#  [45] "\t\tgsl_matrix *B=gsl_matrix_calloc(nrow,ncol);"
			code<-paste(c(code,function_gsl_matrix_set(Pattern, StartVal, Fit, MatrixName="B"),collapse="\n")	
			#  [52] "\t\tgsl_blas_dgemv(CblasNoTrans, 1.0, B, co_variate, 1.0, F_dx_dt);"
			#  [53] "\t\tgsl_matrix_free(B);"
		}else{#NonLinear
			#XXXXXXXXX
		}
	}
		
	#  [57] "}"
	
	"/**\n * The dF/dx function\n * The partial derivative of the jacobian of the DE function with respect to the variable x\n * @param param includes at the end the current state estimates in the same order as the states following the model parameters\n */"
	#  [66] "void function_dF_dx(double t, size_t regime, double *param, const gsl_vector *co_variate, gsl_matrix *F_dx_dt_dx){"
	#  [67] ""
	#  [68] ""
	#  [69] "    /*Supply the Jacobian matrix for the ODEs"
	#  [70] "      ODE functions go down the rows; latent states go across columns*/"
	#  [71] "\t"
	#  [72] "\tgsl_matrix_set(F_dx_dt_dx,0,0,0.0);"
	#  [73] "\tgsl_matrix_set(F_dx_dt_dx,0,1,1.0);"
	#  [74] "\tgsl_matrix_set(F_dx_dt_dx,1,0,param[0]);"
	#  [75] "\tgsl_matrix_set(F_dx_dt_dx,1,1,param[1]);"
	#  [76] "\t"
	if (isLinear){
		#  [79] "\t/*size_t nrow=F_dx_dt_dx->size1, ncol=F_dx_dt_dx->size2,index_col,index_row;"
		function_gsl_matrix_set(Pattern, StartVal, Fit, MatrixName="F_dx_dt_dx")
		
	}else{#NonLinear
		
	}
	#  [88] "}"
	#  [89] ""
	dP_dt="/**\n * The dP/dt function: depend on function_dF_dx, needs to be compiled on the user end\n * but user does not need to modify it or care about it.\n */\nvoid mathfunction_mat_to_vec(const gsl_matrix *mat, gsl_vector *vec){\n\tsize_t i,j;\n\tsize_t nx=mat->size1;\n\t/*convert matrix to vector*/\n\tfor(i=0; i<nx; i++){\n\t\tgsl_vector_set(vec,i,gsl_matrix_get(mat,i,i));\n\t\tfor (j=i+1;j<nx;j++){\n\t\t\tgsl_vector_set(vec,i+j+nx-1,gsl_matrix_get(mat,i,j));\n\t\t\t/*printf(\"%lu\",i+j+nx-1);}*/\n\t\t}\n\t}\n}\nvoid mathfunction_vec_to_mat(const gsl_vector *vec, gsl_matrix *mat){\n\tsize_t i,j;\n\tsize_t nx=mat->size1;\n\t/*convert vector to matrix*/\n\tfor(i=0; i<nx; i++){\n\t\tgsl_matrix_set(mat,i,i,gsl_vector_get(vec,i));\n\t\tfor (j=i+1;j<nx;j++){\n\t\t\tgsl_matrix_set(mat,i,j,gsl_vector_get(vec,i+j+nx-1));\n\t\t\tgsl_matrix_set(mat,j,i,gsl_vector_get(vec,i+j+nx-1));\n\t\t}\n\t}\n}\nvoid function_dP_dt(double t, size_t regime, const gsl_vector *p, double *param, size_t n_param, const gsl_vector *co_variate, gsl_vector *F_dP_dt){\n\t\n\tsize_t nx;\n\tnx = (size_t) floor(sqrt(2*(double) p->size));\n\tgsl_matrix *P_mat=gsl_matrix_calloc(nx,nx);\n\tmathfunction_vec_to_mat(p,P_mat);\n\tgsl_matrix *F_dx_dt_dx=gsl_matrix_calloc(nx,nx);\n\tfunction_dF_dx(t, regime, param, co_variate, F_dx_dt_dx);\n\tgsl_matrix *dFP=gsl_matrix_calloc(nx,nx);\n\tgsl_matrix *dP_dt=gsl_matrix_calloc(nx,nx);\n\tgsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, F_dx_dt_dx, P_mat, 0.0, dFP);\n\tgsl_matrix_transpose_memcpy(dP_dt, dFP);\n\tgsl_matrix_add(dP_dt, dFP);\n\tsize_t n_Q_vec=(1+nx)*nx/2;\n\tgsl_vector *Q_vec=gsl_vector_calloc(n_Q_vec);\n\tsize_t i;\n\tfor(i=1;i<=n_Q_vec;i++){\n\t\t\tgsl_vector_set(Q_vec,n_Q_vec-i,param[n_param-i]);\n\t}\n\tgsl_matrix *Q_mat=gsl_matrix_calloc(nx,nx);\n\tmathfunction_vec_to_mat(Q_vec,Q_mat);\n\tgsl_matrix_add(dP_dt, Q_mat);\n\tmathfunction_mat_to_vec(dP_dt, F_dP_dt);\n\tgsl_matrix_free(P_mat);\n\tgsl_matrix_free(F_dx_dt_dx);\n\tgsl_matrix_free(dFP);\n\tgsl_matrix_free(dP_dt);\n\tgsl_vector_free(Q_vec);\n\tgsl_matrix_free(Q_mat);\n}\n"
	code<-paste(c(code,dP_dt), collapse="\n")
}else{#DiscreteTime
	#XXXXXXXXX
}

#  [90] "/**"
#  [91] " * Set the initial condition"
#  [92] " */"
#  [93] ""
#  [94] "void function_initial_condition(double *param, gsl_vector **co_variate, gsl_vector *pr_0, gsl_vector **eta_0, gsl_matrix **error_cov_0){"
#  [95] ""
#  [96] "    gsl_vector_set(pr_0,0,1);"
#  [97] ""
#  [98] "    size_t num_regime=pr_0->size;"
#  [99] "    size_t dim_latent_var=error_cov_0[0]->size1;"
# [100] "    size_t num_sbj=(eta_0[0]->size)/(dim_latent_var);"
# [101] "    /*printf(\"%lu %lu %lu\\n\",num_regime,dim_latent_var,num_sbj);*/"
# [102] "    "
# [103] "\tsize_t i,j;"
# [104] ""
# [105] "    for(j=0;j<num_regime;j++){"
# [106] "        for(i=0;i<num_sbj;i++){"
# [107] "            /*printf(\"%lu %lu\\n\",i,j);*/"
# [108] "\t\t\tgsl_vector_set((eta_0)[j],i*dim_latent_var,param[4]);"
# [109] "\t\t\tgsl_vector_set((eta_0)[j],i*dim_latent_var+1,1.0);"
# [110] "\t\t\t"
# [111] "\t\t\t/*in R, use code like below to produce the lines */"
# [112] "\t\t\t/*for(index_col=0; index_col<dim_latent_var; index_col++){"
# [113] "\t\t\t\t\tgsl_vector_set((eta_0)[j],i*dim_latent_var+index_col,param[index_col]);"
# [114] "\t\t\t}*/"
# [115] "\t\t\t"
# [116] "        }/*statevar_1_p1 statevar_2_p1 statevar_1_p2 statevar_2_p2 ..., eta_0[] with a length of num_sbj*dim_latent_var*/"
# [117] ""
# [118] "\t\t"
# [119] "\t\tgsl_matrix_set((error_cov_0)[j],0, 0, log(1.0));"
# [120] "\t\tgsl_matrix_set((error_cov_0)[j],1, 1, log(1.0));"
# [121] "\t\t"
function_gsl_matrix_set(Pattern, StartVal, Fit, MatrixName="error_cov_0)[j]")
# [139] ""
# [140] "    }"
# [141] ""
# [142] "}"
# [143] ""
# [144] "/**"
# [145] " * Set the regime-switch transition probability matrix"
# [146] " */"
# [147] ""
# [148] "void function_regime_switch(size_t t, size_t type, double *param, const gsl_vector *co_variate, gsl_matrix *regime_switch_mat){"
# [149] "     /*if there is only one regime*/"
# [150] "\tgsl_matrix_set_identity(regime_switch_mat);"
# [151] "\t"
# [152] "}"
# [153] ""
# [154] "/**"
# [155] " * Set the noise covariance matrix"
# [156] " *They are to be LDL' transformed."
# [157] " *e.g., [a b"
# [158] "         b c]"
# [159] " *-->LDL', L=[1 0;b 1], D=diag(a,c)"
# [160] " */"
# [161] ""
# [162] "void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov){"
# [163] ""
# [164] "    gsl_matrix_set(eta_noise_cov,0,0,-10.0);"
# [165] "\tgsl_matrix_set(eta_noise_cov,1,1,param[2]);"
# [166] "\tgsl_matrix_set(y_noise_cov,0,0, param[3]);"
# [167] "\t/*in R use for loops to generate the lines above for symmetric and diagonal matrices*/"
# [168] ""
# [169] "}"
# [170] ""
# [171] "/**"
# [172] " * This function modifies function parameters, excluding parameters in noise_cov matrices, so that it satisfies the model constraint."
# [174] " */"
# [175] "void function_transform(double *param){"
# [176] "}"
# [177] ""
# [178] ""

function_gsl_matrix_set(Pattern,StartVal,Fit,MatrixName){
	code=""
	if (Pattern=="Sym"){
		for (index_col in 1:ncol){
			if (fit){
				code<-paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ",index_row, index_col,param[",index_col*(index_col+1)+index_col,"]);"),""), collapse="\n")
			}else{
				code<-paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ",index_row, index_col,",StartVal[index_col*(index_col+1)+index_col],");"),""), collapse="\n")
			}
			for (index_row in (index_col+1):nrow){
				if (fit){
					code<-paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ",index_row, index_col,param[",index_row*(index_col+1)+index_col,"]);"),""), collapse="\n")
				}else{
					code<-paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ",index_row, index_col,",StartVal[index_row*(index_col+1)+index_col],");"),""), collapse="\n")
					code<-paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ",index_col, index_row,",StartVal[index_row*(index_col+1)+index_col],");"),""), collapse="\n")
				}
			}
		}
	}else if (Pattern=="Diag"){
		for (index_col in 1:ncol){
			if (fit){
				code<-paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ",index_row, index_col,param[",index_col,"]);"),""), collapse="\n")
			}else{
				code<-paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ",index_row, index_col,",StartVal[index_col],");"),""), collapse="\n")
			}
		}
	}else{#Full
		for (index_col in 1:ncol){
			for (index_row in 1:nrow){
				if (fit){
					code<-paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ",index_row, index_col,param[",index_row*ncol+index_col,"]);"),""), collapse="\n")
				}else{
					code<-paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ",index_row, index_col,",StartVal[index_row*ncol+index_col],");"),""), collapse="\n")
				}
	
			}
		}
	}
	return(code)
}


