/********************************************************
 * Estimation Routines*
 * @Author: Lu Ou.                       *
 * @created Summer 2014
 * Purpose:    *
 * Usage:                                        *
 * References: A                                   *
 * File formats:                                 *
 * Restrictions:*
 * Revision history:*
 * Error handling:*
 * Notes:*
 ********************************************************/
/**
 * This file implements the brekfis: b* regime-switch extended kim filter.
 */

#include "brekfis.h"
#include "ekf.h"
#include "data_structure.h"
#include "math_function.h"
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <time.h>
#include "print_function.h"

/**
 * This method implements a one-step brekfis
 * @param y the observations
 * @param total_time the number of total time points
 * @param config the configuration of the model
 * @param init the initial values for some parameters
 * @param param the model and user-defined function parameters
 * @return log-likelihood
 */
double brekfis(gsl_vector ** y, gsl_vector **co_variate, size_t total_time, double *y_time, const ParamConfig *config, ParamInit *init, Param *param){
	int DEBUG_BREKFIS = 0; /*0=false/no; 1=true/yes*/
    size_t t, regime_j, regime_k, sbj;
    double neg_log_p, p, log_like=0, innov_determinant;

    size_t col_index;
    double sum_overj;
    size_t type;

    /**************initialization*****************************************************************/
     /** Input for EKF **/
    gsl_vector **eta_j_t=(gsl_vector **)malloc(config->num_regime*sizeof(gsl_vector *));
    for(regime_j=0; regime_j<config->num_regime; regime_j++){
        eta_j_t[regime_j]=gsl_vector_alloc(config->dim_latent_var);
    }

    gsl_matrix **error_cov_j_t=(gsl_matrix **)malloc(config->num_regime*sizeof(gsl_matrix *));
    for(regime_j=0; regime_j<config->num_regime; regime_j++){
        error_cov_j_t[regime_j]=gsl_matrix_alloc(config->dim_latent_var, config->dim_latent_var);

    }

    /** output for EKF **/
    gsl_vector ***eta_jk_t_plus_1=(gsl_vector ***)malloc(config->num_regime*sizeof(gsl_vector **));
    for(regime_j=0; regime_j<config->num_regime; regime_j++){
        eta_jk_t_plus_1[regime_j]=(gsl_vector **)malloc(config->num_regime*sizeof(gsl_vector *));
    }
    for(regime_j=0; regime_j<config->num_regime; regime_j++){
        for(regime_k=0; regime_k<config->num_regime; regime_k++)
            eta_jk_t_plus_1[regime_j][regime_k]=gsl_vector_calloc(config->dim_latent_var);
    }
    gsl_matrix ***error_cov_jk_t_plus_1=(gsl_matrix ***)malloc(config->num_regime*sizeof(gsl_matrix **));
    for(regime_j=0; regime_j<config->num_regime; regime_j++){
        error_cov_jk_t_plus_1[regime_j]=(gsl_matrix **)malloc(config->num_regime*sizeof(gsl_matrix *));
    }
    for(regime_j=0; regime_j<config->num_regime; regime_j++){
        for(regime_k=0; regime_k<config->num_regime; regime_k++)
            error_cov_jk_t_plus_1[regime_j][regime_k]=gsl_matrix_calloc(config->dim_latent_var, config->dim_latent_var);
    }
    gsl_vector ***innov_v=(gsl_vector ***)malloc(config->num_regime*sizeof(gsl_vector **));
    for(regime_j=0; regime_j<config->num_regime; regime_j++){
        innov_v[regime_j]=(gsl_vector **)malloc(config->num_regime*sizeof(gsl_vector *));
    }
    for(regime_j=0; regime_j<config->num_regime; regime_j++){
        for(regime_k=0; regime_k<config->num_regime; regime_k++)
            innov_v[regime_j][regime_k]=gsl_vector_calloc(config->dim_obs_var);
    }
    gsl_matrix ***residual_cov=(gsl_matrix ***)malloc(config->num_regime*sizeof(gsl_matrix **));
    for(regime_j=0; regime_j<config->num_regime; regime_j++){
        residual_cov[regime_j]=(gsl_matrix **)malloc(config->num_regime*sizeof(gsl_matrix *));
    }
    for(regime_j=0; regime_j<config->num_regime; regime_j++){
        for(regime_k=0; regime_k<config->num_regime; regime_k++)
            residual_cov[regime_j][regime_k]=gsl_matrix_calloc(config->dim_obs_var, config->dim_obs_var);/*may be the inverse of the residual/innov cov*/
    }

    /** input for hamilton filter **/
    /* handling missing data */
    gsl_vector *cp_y_t=gsl_vector_alloc(y[0]->size);
    gsl_vector *y_non_miss=gsl_vector_alloc(y[0]->size);
    size_t miss_case;
	
    gsl_vector *pr_t=gsl_vector_alloc(config->num_regime);

    double tran_prob_jk;

    /** output for hamilton filter **/
    gsl_matrix *like_jk=gsl_matrix_alloc(config->num_regime, config->num_regime);
    /*gsl_vector *pr_t_plus_1=gsl_vector_alloc(config->num_regime);*/

    /** input for collapse_process**/
    gsl_vector *diff_eta_vec=gsl_vector_alloc(config->dim_latent_var);
    gsl_matrix *diff_eta=gsl_matrix_alloc(config->dim_latent_var, 1);
    gsl_matrix *modif_p=gsl_matrix_alloc(config->dim_latent_var, config->dim_latent_var);


    /********************************************************************************/


        /*for(index=0; index<config->dim_latent_var; index++){
            fprintf(h_file, " %lf", gsl_vector_get(eta_j_t[0], index));
        }
        fprintf(h_file,"\n");*/
        /*FILE *eta_file=fopen("results/eta_t_ada.txt","w");*/
        /*FILE *pr_file=fopen("regimeprob.txt","w");*/
	
	
	/********************************************************************************/
	for(sbj=0; sbj < config->num_sbj; sbj++){
		for(t=(config->index_sbj)[sbj]; t < (config->index_sbj)[sbj+1]; t++){
			
			gsl_vector_memcpy(cp_y_t, y[t]);
			miss_case = find_miss_data(cp_y_t, y_non_miss); /* 0 - no miss, 1 - part miss, 2 - all miss*/
			
			/** step 1: call cda ekalman filter for each possible regime switch **/
			for(regime_j=0; regime_j < config->num_regime; regime_j++){/*from regime j*/
				
				/**set the regime switch matrix**/
				if (t==(config->index_sbj)[sbj]){
					gsl_matrix_set_identity(param->regime_switch_mat);
					gsl_vector_memcpy(pr_t, init->pr_0[sbj]);
					if(DEBUG_BREKFIS){
						MYPRINT("initial regime probabilities:\n");
						print_vector(pr_t);
						MYPRINT("\n");
					}
				}else{
					type=1;
					config->func_regime_switch(t, type, param->func_param, co_variate[t], param->regime_switch_mat);
				}
				
				config->func_noise_cov(t, regime_j, param->func_param, param->y_noise_cov, param->eta_noise_cov);
				model_constraint_par(config, param);
				
				if(DEBUG_BREKFIS){
					MYPRINT("sbj %lu at time %lu in regime %lu:\n",sbj,t,regime_j);
					MYPRINT("\n");
					MYPRINT("regime_switch_matrix:\n");
					print_matrix(param->regime_switch_mat);
					MYPRINT("\n");
					MYPRINT("parameters:\n");
					print_array(param->func_param,config->num_func_param);
					MYPRINT("\n");
					MYPRINT("measurement error:\n");
					print_matrix(param->y_noise_cov);
					MYPRINT("\n");
					MYPRINT("process noise: \n");
					print_matrix(param->eta_noise_cov);
					MYPRINT("\n");
				}
				
				
				for(regime_k=0; regime_k<config->num_regime; regime_k++){/*to regime k*/
					
					if (t==(config->index_sbj)[sbj]){
						for(col_index=0; col_index<config->dim_latent_var; col_index++){
							gsl_vector_set(eta_j_t[regime_j], col_index, gsl_vector_get((init->eta_0)[regime_j], config->dim_latent_var*sbj+col_index));
						}
						gsl_matrix_memcpy(error_cov_j_t[regime_j], (init->error_cov_0)[regime_j]);
						
						/*MYPRINT("eta_S_at_a_previous_time_point:\n");
						print_vector(eta_j_t[regime_j]);
						MYPRINT("\n");
						MYPRINT("error_cov_at_a_previous_time_point:\n");
						print_matrix(error_cov_j_t[regime_j]);
						MYPRINT("\n");*/
						
						innov_determinant = ext_kalmanfilter_updateonly(t, regime_k,
							eta_j_t[regime_j], error_cov_j_t[regime_j],
							y[t],co_variate[t],y_time,
							param->eta_noise_cov, param->y_noise_cov,
							param->func_param,
							config->func_measure,
							eta_jk_t_plus_1[regime_j][regime_k],
							error_cov_jk_t_plus_1[regime_j][regime_k],
							innov_v[regime_j][regime_k],
							residual_cov[regime_j][regime_k]);/*inverse*/
						
						/*MYPRINT("From regime %lu to regime %lu:\n",regime_j,regime_k);
						MYPRINT("\n");
						MYPRINT("eta_jk:\n");
						print_vector(eta_jk_t_plus_1[regime_j][regime_k]);
						MYPRINT("\n");
						MYPRINT("error_cov_jk:\n");
						print_matrix(error_cov_jk_t_plus_1[regime_j][regime_k]);
						MYPRINT("\n");
						MYPRINT("innov_vector:\n");
						print_vector(innov_v[regime_j][regime_k]);
						MYPRINT("\n");
						MYPRINT("inverse of the residual covariance:\n");
						print_matrix(residual_cov[regime_j][regime_k]);
						MYPRINT("\n");*/
						
					} else {
						/*MYPRINT("eta_S_at_a_previous_time_point:\n");
						print_vector(eta_j_t[regime_j]);
						MYPRINT("\n");
						MYPRINT("error_cov_at_a_previous_time_point:\n");
						print_matrix(error_cov_j_t[regime_j]);
						MYPRINT("\n");*/
						
						innov_determinant=ext_kalmanfilter(t, regime_k,
							eta_j_t[regime_j], error_cov_j_t[regime_j],
							y[t],co_variate[t],y_time,
							param->eta_noise_cov, param->y_noise_cov,
							param->func_param,config->num_func_param,
							config->isContinuousTime,
							config->func_measure,
							config->func_dx_dt,
							config->func_dP_dt,
							config->func_dF_dx,
							config->func_dynam,
							config->func_jacob_dynam,
							eta_jk_t_plus_1[regime_j][regime_k],
							error_cov_jk_t_plus_1[regime_j][regime_k],
							innov_v[regime_j][regime_k],
							residual_cov[regime_j][regime_k]);/*inverse*/
						
						/*MYPRINT("From regime %lu to regime %lu:\n",regime_j,regime_k);
						MYPRINT("\n");
						MYPRINT("eta_jk:\n");
						print_vector(eta_jk_t_plus_1[regime_j][regime_k]);
						MYPRINT("\n");
						MYPRINT("error_cov_jk:\n");
						print_matrix(error_cov_jk_t_plus_1[regime_j][regime_k]);
						MYPRINT("\n");
						MYPRINT("innov_vector:\n");
						print_vector(innov_v[regime_j][regime_k]);
						MYPRINT("\n");
						MYPRINT("inverse of the residual covariance:\n");
						print_matrix(residual_cov[regime_j][regime_k]);
						MYPRINT("\n");*/
					}
					
					/*for(col_index=0; col_index<config->dim_latent_var; col_index++){
						fprintf(eta_file, " %lf", gsl_vector_get(eta_jk_t_plus_1[0][0], col_index));
					}
					fprintf(eta_file,"\n");*/
					
					/*change, for random effect estimation*/
					/*if (t==(config->index_sbj)[sbj+1]-1){
						gsl_vector_set((init->eta_0)[regime_j],config->dim_latent_var*sbj+2,gsl_vector_get(eta_jk_t_plus_1[regime_j][regime_k], 2));
					}*/
					
					/** step 2: call hamilton filter to compute the probability of moving one step ahead **/
					
					/** Step 2.1: compute transition probability matrix, Pr(S_{t-1} = j,S_{t} = k|Y_{t-1}) given the pr_t_1 **/
					tran_prob_jk=gsl_vector_get(pr_t, regime_j)*gsl_matrix_get(param->regime_switch_mat, regime_j, regime_k);
					
					/*MYPRINT("prob_regime:\n");
					print_vector(pr_t);
					MYPRINT("\n");*/
					
					/** Step 2.2: compute log value of function f(.), i.e., prediction error decomposition function **/
					neg_log_p=mathfunction_negloglike_multivariate_normal_invcov(innov_v[regime_j][regime_k], residual_cov[regime_j][regime_k], y_non_miss, innov_determinant);
					
					/** compare the p with the (0.0001) and get the bigger one. We do not like probability that is too small. :)**/
					double numNotMissingVars = mathfunction_sum_vector(y_non_miss);
					double tooSmallNumber = numNotMissingVars < 30 ? pow(1e-10, numNotMissingVars):1e-300;
					double tryP = exp(-neg_log_p);
					p = ( isfinite(tryP) && (tryP > tooSmallNumber) ) ? tryP:tooSmallNumber;
					
					if(DEBUG_BREKFIS){
						MYPRINT("likelihood f(y_it|S_it=k,S_i,t-1=j,Y_i,t-1):\n");
						MYPRINT("before exponential %lf\n",neg_log_p);
						MYPRINT("original %lf\n",exp(-neg_log_p));
						MYPRINT("adjusted %lf\n",p);
					}
					
					/*p=exp(-neg_log_p)*tran_prob_jk;*/
					gsl_matrix_set(like_jk, regime_j, regime_k, p*tran_prob_jk);
					
				}/*end of from regime j*/
			}/*end of to regime k*/
			/*Still inside the subject and time loops*/
			
			/** Step 2.3: update transit probability Pr(S_{t-1} = j,S_{t} = k|Y_t) given Pr(S_{t-1} = j,S_{t} = k|Y_{t-1})**/
			if (config->isnegloglikeweightedbyT){
                log_like+=log(mathfunction_matrix_normalize(like_jk))/((config->index_sbj)[sbj+1]-(config->index_sbj)[sbj]);
            }else{
                log_like+=log(mathfunction_matrix_normalize(like_jk));
            }
            /*like_jk scaled, Pr(S_{t-1} = j,S_{t} = k|Y_{t}, like=sum*/
            /*individual ÐLL is now divided by each individualÕs number of occasions. This helps speed convergence and ensure that each individualÕs data get weighted equally regardless of T_i.*/

			if(DEBUG_BREKFIS){
	            MYPRINT("Pr(S_{t-1} = j,S_{t} = k|Y_{t}):\n");
	            print_matrix(like_jk);
	            MYPRINT("\n");
	
	            MYPRINT("negative log likelihood at time %lu:\n",t);
	            MYPRINT("%lf",-log_like);
	            MYPRINT("\n");
	
	            MYPRINT("Pr(S_it=k|Y_it) original and adjusted:\n");
			}

            for(regime_k=0; regime_k<config->num_regime; regime_k++){
                    sum_overj=0;

            	    for(regime_j=0; regime_j<config->num_regime; regime_j++){
            	    	    sum_overj+=gsl_matrix_get(like_jk, regime_j, regime_k);
            	    }

            	    /** step 2.4: sum transition probability to obtain pr_t**/
            	    gsl_vector_set(pr_t, regime_k, sum_overj);/*write pr_t_plus_1 into pr_t*/

            	    /*MYPRINT("%lf ",sum_overj);*/

            	}/*end of k*/

            /** step 2.4.1: check whether there is zero probability. If so, a small amount of value is added. Again we do not like too small and zero probability **/
			double tooSmallRegimeNumber = config->num_regime < 30 ? pow(1e-10, config->num_regime):1e-300;
			if(gsl_vector_min(pr_t) < tooSmallRegimeNumber){
				gsl_vector_add_constant(pr_t, tooSmallRegimeNumber);
				mathfunction_vector_normalize(pr_t);
			}
			

	    /*MYPRINT("\n");
	    print_vector(pr_t);
	    MYPRINT("\n");*/
			
            for(regime_k=0; regime_k<config->num_regime; regime_k++){
                    gsl_vector_set_zero(eta_j_t[regime_k]);/*here, corresponds to eta_k_t in the paper*/
                    gsl_matrix_set_zero(error_cov_j_t[regime_k]);/*here, corresponds to error_cov_k_t in the paper*/

            	    for(regime_j=0; regime_j<config->num_regime; regime_j++){
            	    	    /** step 3: call collapse process **/
            	    	    /** step 3.1: collapse the latent variable to get eta_k **/


            	    	    gsl_blas_daxpy(gsl_matrix_get(like_jk, regime_j, regime_k), eta_jk_t_plus_1[regime_j][regime_k], eta_j_t[regime_k]); /*sum over j through loop*/

            	    	    /*if(regime_k==1){
            	    	    MYPRINT("Here!");
            	    	    print_vector(eta_j_t[regime_k]);}*/
            	    }
            	}/*end of k*/
				
            for(regime_k=0; regime_k<config->num_regime; regime_k++){
            	    gsl_vector_scale(eta_j_t[regime_k], 1.0/gsl_vector_get(pr_t,regime_k));

            	    /*MYPRINT("eta collapsed estimate in regime %lu:\n",regime_k);
            	    print_vector(eta_j_t[regime_k]);
            	    MYPRINT("\n");*/

            	    /** step 3.2: collapse the covariance matrix to get error_cov_k **/
            	    for(regime_j=0; regime_j<config->num_regime; regime_j++){
						
						mathfunction_collapse(eta_j_t[regime_k], eta_jk_t_plus_1[regime_j][regime_k], 
						error_cov_jk_t_plus_1[regime_j][regime_k], gsl_matrix_get(like_jk,regime_j, regime_k), error_cov_j_t[regime_k],
						diff_eta_vec, diff_eta, modif_p);

                    }/*end of j*/
                    gsl_matrix_scale(error_cov_j_t[regime_k], 1.0/gsl_vector_get(pr_t,regime_k));

                    /*MYPRINT("error_cov collapsed estimate in regime %lu:\n",regime_k);
            	    print_matrix(error_cov_j_t[regime_k]);
            	    MYPRINT("\n");*/
            }/*end of k*/



	   /*fprintf(pr_file,"%lu %lu %lf %lf\n",sbj,t,gsl_vector_get(pr_t,0),gsl_vector_get(pr_t,1));*/
	   /*MYPRINT("Pr_St|t:\n");
           print_vector(pr_t);
           MYPRINT("\n");*/

        }/*end of t*/

       /*if (sbj==2){exit(0);}*/
         /*fprintf(h_file, "%d", t+1);*/
    }/*end of sbj*/
	
	/*fclose(h_file);*/
	/*fclose(eta_file);*/
	/*fclose(pr_file);*/
	/****************************** free allocated space ***************************/
	for(regime_j=0; regime_j<config->num_regime; regime_j++)
		gsl_vector_free(eta_j_t[regime_j]);
	free(eta_j_t);
	for(regime_j=0; regime_j<config->num_regime; regime_j++)
		gsl_matrix_free(error_cov_j_t[regime_j]);
	free(error_cov_j_t);
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
		for(regime_k=0; regime_k<config->num_regime; regime_k++)
		gsl_vector_free(eta_jk_t_plus_1[regime_j][regime_k]);
	}
	for(regime_j=0; regime_j<config->num_regime; regime_j++)
		free(eta_jk_t_plus_1[regime_j]);
	free(eta_jk_t_plus_1);
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
		for(regime_k=0; regime_k<config->num_regime; regime_k++)
			gsl_matrix_free(error_cov_jk_t_plus_1[regime_j][regime_k]);
	}
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
		free(error_cov_jk_t_plus_1[regime_j]);
	}
	free(error_cov_jk_t_plus_1);
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
		for(regime_k=0; regime_k<config->num_regime; regime_k++)
			gsl_vector_free(innov_v[regime_j][regime_k]);
	}
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
		free(innov_v[regime_j]);
	}
	free(innov_v);
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
		for(regime_k=0; regime_k<config->num_regime; regime_k++)
			gsl_matrix_free(residual_cov[regime_j][regime_k]);
	}
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
		free(residual_cov[regime_j]);
	}
	free(residual_cov);
	
	gsl_vector_free(pr_t);
	gsl_matrix_free(like_jk);
	gsl_vector_free(cp_y_t);
	gsl_vector_free(y_non_miss);
	gsl_vector_free(diff_eta_vec);
	gsl_matrix_free(diff_eta);
	gsl_matrix_free(modif_p);
	
	return(-log_like);
}



/**
 * This function modifies some of the parameters so that it satisfies the model constraint.
 */
void model_constraint_par(const ParamConfig *pc, Param *par){
    /*double min_p=1e-4;*/
    double v=0;
    size_t ri, ci;

    /*gsl_vector *temp=gsl_vector_alloc(pc->num_regime);*/
    /** each element of regime switch matrix must be non-negative and the row sum is equal to 1 **/
    /*for(ri=0; ri<pc->num_regime; ri++){
        gsl_matrix_get_row(temp, par->regime_switch_mat, ri);
        gsl_vector_set(temp, pc->num_regime-1, 0);
        for(ci=0; ci<pc->num_regime; ci++){
            v=gsl_vector_get(temp, ci);
            if(v<min_p){
                gsl_vector_get(temp, ci);
            }
        }
        mathfunction_normalize_log_vector(temp);
        for(ci=0; ci<pc->num_regime; ci++){
            v=gsl_vector_get(temp, ci);
            if(v<min_p){
                gsl_vector_set(temp, ci, min_p);
            }
        }
        mathfunction_vector_normalize(temp);
        gsl_matrix_set_row(par->regime_switch_mat, ri, temp);
    }
    gsl_vector_free(temp);*/

     /** covariance of process noise must be positive definite: LDL' decomposition applied**/
    gsl_matrix *D;
    gsl_matrix *L;

    L=gsl_matrix_calloc(pc->dim_latent_var, pc->dim_latent_var);
    D=gsl_matrix_calloc(pc->dim_latent_var, pc->dim_latent_var);
    gsl_matrix *temp_eta_noise=gsl_matrix_calloc(pc->dim_latent_var, pc->dim_latent_var);
    gsl_matrix_set_zero(D);
    gsl_matrix_memcpy(L, par->eta_noise_cov);/*set the lower diagnoal*/
    for(ri=0; ri<pc->dim_latent_var; ri++){
        gsl_matrix_set(L,ri,ri,1);/*set diagonals to 1*/
        v=gsl_matrix_get(par->eta_noise_cov, ri, ri);
        v=exp(v);
        gsl_matrix_set(D, ri, ri, v);/*D = diag(exp(V11)~exp(Vnn))*/
        for(ci=ri+1; ci<pc->dim_latent_var; ci++){
            gsl_matrix_set(L, ri, ci, 0);/*set upper diagonal to zero*/
        }
    }
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, L, D, 0.0, temp_eta_noise);/*temp_eta_noise=LD*/
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, temp_eta_noise, L, 0.0, par->eta_noise_cov);/*eta_noise_cov=LDL'*/

    gsl_matrix_free(temp_eta_noise);
    gsl_matrix_free(L);
    gsl_matrix_free(D);

     /** covariance of measurement noise must be positive definite: LDL' decomposition applied **/
    L=gsl_matrix_calloc(pc->dim_obs_var, pc->dim_obs_var);
    D=gsl_matrix_calloc(pc->dim_obs_var, pc->dim_obs_var);
    gsl_matrix *temp_y_noise=gsl_matrix_calloc(pc->dim_obs_var, pc->dim_obs_var);
    gsl_matrix_set_zero(D);
    gsl_matrix_memcpy(L, par->y_noise_cov);/*set the lower diagnoal*/
    for(ri=0; ri<pc->dim_obs_var; ri++){
        gsl_matrix_set(L,ri,ri,1);/*set diagonals to 1*/
        v=gsl_matrix_get(par->y_noise_cov, ri, ri);
        v=exp(v);
        gsl_matrix_set(D, ri, ri, v);/*D = diag(exp(V11)~exp(Vnn))*/
        for(ci=ri+1; ci<pc->dim_obs_var; ci++){
            gsl_matrix_set(L, ri, ci, 0);/*set upper diagonal to zero*/
        }
    }
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, L, D, 0.0, temp_y_noise);/*temp_eta_noise=LD*/
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, temp_y_noise, L, 0.0, par->y_noise_cov);/*eta_noise_cov=LDL'*/

    gsl_matrix_free(temp_y_noise);
    gsl_matrix_free(L);
    gsl_matrix_free(D);
    


}
/**
 * This function modifies some of the parameters in the initial condition so that it satisfies the model constraint.
 */
void model_constraint_init(const ParamConfig *pc, ParamInit *pi){
    double v=0;
    size_t ri,ci;

    /** each initial distribution is non-negative and the sum is equal to 1 **/
    /*double min_p=1e-4;
    gsl_vector_set(pi->pr_0, pc->num_regime-1, 0);

    mathfunction_normalize_log_vector(pi->pr_0);

    for(ci=0; ci<pc->num_regime; ci++){
        v=gsl_vector_get(pi->pr_0, ci);
        v=v<min_p?min_p:v;
        gsl_vector_set(pi->pr_0, ci, v);
    }
    mathfunction_vector_normalize(pi->pr_0*/

    /** covariance of the initial condition must be positive definite: LDL' decomposition applied **/
    size_t regime_j;
    gsl_matrix *L=gsl_matrix_calloc(pc->dim_latent_var, pc->dim_latent_var);
    gsl_matrix *D=gsl_matrix_calloc(pc->dim_latent_var, pc->dim_latent_var);
    gsl_matrix *temp_error_cov_0=gsl_matrix_calloc(pc->dim_latent_var, pc->dim_latent_var);
    
    for (regime_j=0;regime_j<pc->num_regime;regime_j++){
        gsl_matrix_set_zero(D);
        /*set the lower diagnoal*/
        gsl_matrix_memcpy(L, pi->error_cov_0[regime_j]);
        for(ri=0; ri<pc->dim_latent_var; ri++){
            /*set diagonals to 1*/
            gsl_matrix_set(L,ri,ri,1);
            /*set upper diagonal to zero*/
            for(ci=ri+1; ci<pc->dim_latent_var; ci++){
            gsl_matrix_set(L, ri, ci, 0);
            }
            /*D = diag(exp(V11)~exp(Vnn))*/
            v=gsl_matrix_get(pi->error_cov_0[regime_j], ri, ri);
            v=exp(v);            
            gsl_matrix_set(D, ri, ri, v);
           
        }
        /*temp_error_cov_0=LD*/
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, L, D, 0.0, temp_error_cov_0);
        /*eta_noise_cov=LDL'*/
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, temp_error_cov_0, L, 0.0, pi->error_cov_0[regime_j]);
    }
    
    gsl_matrix_free(temp_error_cov_0);
    gsl_matrix_free(L);
    gsl_matrix_free(D);




}
/****************************Extended Kim Filter************************/
/**
* This function implements the extended Kim Filter
* Parameters/Input *
* *
* **>>Parameters/Input Pointers<<**
* y -- the data
* y_time -- continous time points in the data
* co_variate -- covariates
* config -- model configuration,
* param -- parameters
*
* eta_regime_j_t -- eta^k_it|t*
* error_cov_regime_j_t -- error_cov^k_it|t*
* *
* Output*
* *
* **>>Output via using pointers<<**
* eta_regime_j_t -- eta^k_it|t*
* error_cov_regime_j_t -- error_cov^k_it|t*
* pr_t_given_t_minus_1 -- Pr(S_it=k|Y_i,t-1) *
* pr_t -- Pr(S_it=k|Y_it) *
* eta_regime_jk_pred -- eta^regime_jk_it|t-1 *
* error_cov_regime_jk_pred -- error_cov^regime_jk_it|t-1 *
* eta_t -- filtered state estimate
* error_cov_t -- filtered error covariance estimate
* eta_pred_t -- predicted state estimate
* error_cov_pred_t -- predicted error covariance estimate
* innov_v_t -- innovation vector
* residual_cov_t -- inverse of the residual covariance
**/

double EKimFilter(gsl_vector ** y, gsl_vector **co_variate, double *y_time, const ParamConfig *config, ParamInit *init, Param *param,
    gsl_vector ***eta_regime_j_t, gsl_matrix ***error_cov_regime_j_t,
	gsl_vector ****eta_regime_jk_pred, gsl_matrix ****error_cov_regime_jk_pred,
    gsl_vector **pr_t, gsl_vector **pr_t_given_t_minus_1,
	gsl_vector **eta_t, gsl_matrix **error_cov_t, 
	gsl_vector **eta_pred_t, gsl_matrix **error_cov_pred_t,
	gsl_vector **innov_v_t, gsl_matrix **residual_cov_t){


    /************** initialization *****************************************************************/
    size_t t, index_sbj_t, regime_j, regime_k, sbj;
    double neg_log_p,p, log_like=0, innov_determinant;

    size_t col_index;
    double sum_overj;
    size_t type;

    /* handling missing data */
    gsl_vector *cp_y_t=gsl_vector_alloc(y[0]->size);
    gsl_vector *y_non_miss=gsl_vector_alloc(y[0]->size);
    size_t miss_case;
	
	/** output of extended Kalman filter **/
    /*eta^regime_jk_it|t -- eta_regime_jk_t_plus_1 -- filtered regime specific state estimate */
    gsl_vector ****eta_regime_jk_t_plus_1=(gsl_vector ****)malloc(config->total_obs*sizeof(gsl_vector ***));
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	eta_regime_jk_t_plus_1[index_sbj_t]=(gsl_vector ***)malloc(config->num_regime*sizeof(gsl_vector *));
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    eta_regime_jk_t_plus_1[index_sbj_t][regime_j]=(gsl_vector **)malloc(config->num_regime*sizeof(gsl_vector *));
	}
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    for(regime_k=0; regime_k<config->num_regime; regime_k++)
	    eta_regime_jk_t_plus_1[index_sbj_t][regime_j][regime_k]=gsl_vector_calloc(config->dim_latent_var);
	}
    }
    
	/*error_cov^regime_jk_it|t -- error_cov_regime_jk_t_plus_1 -- filtered regime specific error covariance estimate */
    gsl_matrix ****error_cov_regime_jk_t_plus_1=(gsl_matrix ****)malloc(config->total_obs*sizeof(gsl_matrix ***));
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	error_cov_regime_jk_t_plus_1[index_sbj_t]=(gsl_matrix ***)malloc(config->num_regime*sizeof(gsl_matrix **));
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    error_cov_regime_jk_t_plus_1[index_sbj_t][regime_j]=(gsl_matrix **)malloc(config->num_regime*sizeof(gsl_matrix *));
	}
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    for(regime_k=0; regime_k<config->num_regime; regime_k++)
		error_cov_regime_jk_t_plus_1[index_sbj_t][regime_j][regime_k]=gsl_matrix_calloc(config->dim_latent_var, config->dim_latent_var);
	}
    }
	
    /*output of filter: innovation vector*/
    gsl_vector ****innov_v=(gsl_vector ****)malloc(config->total_obs*sizeof(gsl_vector ***));
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	innov_v[index_sbj_t]=(gsl_vector ***)malloc(config->num_regime*sizeof(gsl_vector *));
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    innov_v[index_sbj_t][regime_j]=(gsl_vector **)malloc(config->num_regime*sizeof(gsl_vector *));
	}
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    for(regime_k=0; regime_k<config->num_regime; regime_k++)
	    innov_v[index_sbj_t][regime_j][regime_k]=gsl_vector_calloc(config->dim_obs_var);
	}
    }
	
    /*output of filter: inverse_residual_cov*/
    gsl_matrix ****inv_residual_cov=(gsl_matrix ****)malloc(config->total_obs*sizeof(gsl_matrix ***));
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	inv_residual_cov[index_sbj_t]=(gsl_matrix ***)malloc(config->num_regime*sizeof(gsl_matrix **));
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    inv_residual_cov[index_sbj_t][regime_j]=(gsl_matrix **)malloc(config->num_regime*sizeof(gsl_matrix *));
	}
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    for(regime_k=0; regime_k<config->num_regime; regime_k++)
		inv_residual_cov[index_sbj_t][regime_j][regime_k]=gsl_matrix_calloc(config->dim_obs_var, config->dim_obs_var);
	}
    }
	
    /*output of filter: residual_cov*/
    gsl_matrix ****residual_cov=(gsl_matrix ****)malloc(config->total_obs*sizeof(gsl_matrix ***));
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
		residual_cov[index_sbj_t]=(gsl_matrix ***)malloc(config->num_regime*sizeof(gsl_matrix **));
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    residual_cov[index_sbj_t][regime_j]=(gsl_matrix **)malloc(config->num_regime*sizeof(gsl_matrix *));
	}
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    for(regime_k=0; regime_k<config->num_regime; regime_k++)
		residual_cov[index_sbj_t][regime_j][regime_k]=gsl_matrix_calloc(config->dim_obs_var, config->dim_obs_var);
	}
    }
	

    gsl_vector ***eta_pred_regime_t=(gsl_vector ***)malloc(config->total_obs*sizeof(gsl_vector **));
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	eta_pred_regime_t[index_sbj_t]=(gsl_vector **)malloc(config->num_regime*sizeof(gsl_vector *));
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    eta_pred_regime_t[index_sbj_t][regime_j]=gsl_vector_calloc(config->dim_latent_var);
	}
    }

    gsl_matrix ***error_cov_pred_regime_t=(gsl_matrix ***)malloc(config->total_obs*sizeof(gsl_matrix **));
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	error_cov_pred_regime_t[index_sbj_t]=(gsl_matrix **)malloc(config->num_regime*sizeof(gsl_matrix *));
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    error_cov_pred_regime_t[index_sbj_t][regime_j]=gsl_matrix_calloc(config->dim_latent_var, config->dim_latent_var);
	}
    }

    gsl_vector ***innov_v_regime_t=(gsl_vector ***)malloc(config->total_obs*sizeof(gsl_vector **));
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	innov_v_regime_t[index_sbj_t]=(gsl_vector **)malloc(config->num_regime*sizeof(gsl_vector *));
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    innov_v_regime_t[index_sbj_t][regime_j]=gsl_vector_calloc(config->dim_latent_var);
	}
    }

    gsl_matrix ***residual_cov_regime_t=(gsl_matrix ***)malloc(config->total_obs*sizeof(gsl_matrix **));
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	residual_cov_regime_t[index_sbj_t]=(gsl_matrix **)malloc(config->num_regime*sizeof(gsl_matrix *));
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    residual_cov_regime_t[index_sbj_t][regime_j]=gsl_matrix_calloc(config->dim_latent_var, config->dim_latent_var);
	}
    }
	
	
    /** output for hamilton filter **/
	gsl_matrix *tran_prob_jk = gsl_matrix_alloc(config->num_regime, config->num_regime);/*given t-1*/
    gsl_matrix *like_jk = gsl_matrix_alloc(config->num_regime, config->num_regime);/*given t*/


    /** input for collapse_process**/
    gsl_vector *diff_eta_vec=gsl_vector_alloc(config->dim_latent_var);
    gsl_matrix *diff_eta=gsl_matrix_alloc(config->dim_latent_var, 1);
    gsl_matrix *modif_p=gsl_matrix_alloc(config->dim_latent_var, config->dim_latent_var);


    /********************************************************************************/


        /*for(index=0; index<config->dim_latent_var; index++){
            fprintf(h_file, " %lf", gsl_vector_get(eta_regime_j_t[t][0], index));
        }
        fprintf(h_file,"\n");*/
        /*FILE *eta_file=fopen("eta_t.txt","w");*/
        /*FILE *pr_file=fopen("regimeprob.txt","w");*/


    for(sbj=0; sbj<config->num_sbj; sbj++){

   /********************************************************************************/

        for(t=(config->index_sbj)[sbj]; t<(config->index_sbj)[sbj+1]; t++){
			gsl_vector_memcpy(cp_y_t, y[t]);
			miss_case=find_miss_data(cp_y_t, y_non_miss); /* 0 - no miss, 1 - part miss, 2 - all miss*/
			
        /** step 1: call cda ekalman filter for each possible regime switch **/

            for(regime_j=0; regime_j<config->num_regime; regime_j++){/*from regime j*/

            	/**set the regime switch matrix**/
            	if (t==(config->index_sbj)[sbj]){
            	    gsl_matrix_set_identity(param->regime_switch_mat);
                    gsl_vector_memcpy(pr_t[t], init->pr_0[sbj]);
                }else{
                    type=1;
                    config->func_regime_switch(t, type, param->func_param, co_variate[t], param->regime_switch_mat);
                }

                config->func_noise_cov(t, regime_j, param->func_param, param->y_noise_cov, param->eta_noise_cov);
                model_constraint_par(config, param);

                /*MYPRINT("sbj %lu at time %lu in regime %lu:\n",sbj,t,regime_j);
                MYPRINT("\n");
                MYPRINT("regime_switch_matrix:\n");
                print_matrix(param->regime_switch_mat);
                MYPRINT("\n");
                MYPRINT("parameters:\n");
                print_array(param->func_param,config->num_func_param);
                MYPRINT("\n");
                MYPRINT("measurement error:\n");
                print_matrix(param->y_noise_cov);
                MYPRINT("\n");
                MYPRINT("process noise: \n");
                print_matrix(param->eta_noise_cov);
                MYPRINT("\n");*/


                for(regime_k=0; regime_k<config->num_regime; regime_k++){/*to regime k*/

                    if (t==(config->index_sbj)[sbj]){
                        for(col_index=0; col_index<config->dim_latent_var; col_index++){
            	        gsl_vector_set(eta_regime_j_t[t][regime_j], col_index, gsl_vector_get((init->eta_0)[regime_j], config->dim_latent_var*sbj+col_index));
            	        }
            	        gsl_matrix_memcpy(error_cov_regime_j_t[t][regime_j], (init->error_cov_0)[regime_j]);
            	        
						/*MYPRINT("eta_S_at_a_previous_time_point:\n");
            	        print_vector(eta_regime_j_t[t][regime_j]);
            	        MYPRINT("\n");
            	        MYPRINT("error_cov_at_a_previous_time_point:\n");
            	        print_matrix(error_cov_regime_j_t[t][regime_j]);
            	        MYPRINT("\n");*/


                    innov_determinant=ext_kalmanfilter_updateonly_smoother(t, regime_k,
                        eta_regime_j_t[t][regime_j], error_cov_regime_j_t[t][regime_j],
                        y[t],co_variate[t],y_time,
                        param->eta_noise_cov, param->y_noise_cov,
                        param->func_param,
                        config->func_measure,
                        eta_regime_jk_pred[t][regime_j][regime_k], error_cov_regime_jk_pred[t][regime_j][regime_k],
                        eta_regime_jk_t_plus_1[t][regime_j][regime_k], error_cov_regime_jk_t_plus_1[t][regime_j][regime_k],
                        innov_v[t][regime_j][regime_k], inv_residual_cov[t][regime_j][regime_k], residual_cov[t][regime_j][regime_k]);/*inverse*/

                        /*MYPRINT("From regime %lu to regime %lu:\n",regime_j,regime_k);
                        MYPRINT("\n");
                        MYPRINT("eta_jk_pred:\n");
                        print_vector(eta_regime_jk_pred[t][regime_j][regime_k]);
                        MYPRINT("\n");
                        MYPRINT("error_cov_jk_pred:\n");
                        print_matrix(error_cov_regime_jk_pred[t][regime_j][regime_k]);
                        MYPRINT("\n");
                        MYPRINT("eta_jk:\n");
                        print_vector(eta_regime_jk_t_plus_1[t][regime_j][regime_k]);
                        MYPRINT("\n");
                        MYPRINT("error_cov_jk:\n");
                        print_matrix(error_cov_regime_jk_t_plus_1[t][regime_j][regime_k]);
                        MYPRINT("\n");
                        MYPRINT("innov_v[t]ector:\n");
                        print_vector(innov_v[t][regime_j][regime_k]);
                        MYPRINT("\n");
                        MYPRINT("inverse of the residual covariance:\n");
                        print_matrix(inv_residual_cov[t][regime_j][regime_k]);
                        MYPRINT("\n");*/


                   }else{
					   
            	        /*MYPRINT("eta_S_at_a_previous_time_point:\n");
            	        print_vector(eta_regime_j_t[t-1][regime_j]);
            	        MYPRINT("\n");
            	        MYPRINT("error_cov_at_a_previous_time_point:\n");
            	        print_matrix(error_cov_regime_j_t[t-1][regime_j]);
            	        MYPRINT("\n");*/

                    innov_determinant=ext_kalmanfilter_smoother(t, regime_k,
                        eta_regime_j_t[t-1][regime_j], error_cov_regime_j_t[t-1][regime_j],
                        y[t],co_variate[t],y_time,
                        param->eta_noise_cov, param->y_noise_cov,
                        param->func_param,config->num_func_param,
						config->isContinuousTime,
                        config->func_measure,
                        config->func_dx_dt,
                        config->func_dP_dt,
						config->func_dF_dx,
                        config->func_dynam,
						config->func_jacob_dynam,
                        eta_regime_jk_pred[t][regime_j][regime_k], error_cov_regime_jk_pred[t][regime_j][regime_k],
                        eta_regime_jk_t_plus_1[t][regime_j][regime_k], error_cov_regime_jk_t_plus_1[t][regime_j][regime_k], 
						innov_v[t][regime_j][regime_k], inv_residual_cov[t][regime_j][regime_k], residual_cov[t][regime_j][regime_k]);/*inverse*/

                        /*MYPRINT("From regime %lu to regime %lu:\n",regime_j,regime_k);
                        MYPRINT("\n");
                        MYPRINT("eta_jk_pred:\n");
                        print_vector(eta_regime_jk_pred[t][regime_j][regime_k]);
                        MYPRINT("\n");
                        MYPRINT("error_cov_jk_pred:\n");
                        print_matrix(error_cov_regime_jk_pred[t][regime_j][regime_k]);
                        MYPRINT("\n");
                        MYPRINT("eta_jk:\n");
                        print_vector(eta_regime_jk_t_plus_1[t][regime_j][regime_k]);
                        MYPRINT("\n");
                        MYPRINT("error_cov_jk:\n");
                        print_matrix(error_cov_regime_jk_t_plus_1[t][regime_j][regime_k]);
                        MYPRINT("\n");
                        MYPRINT("innov_vector:\n");
                        print_vector(innov_v[t][regime_j][regime_k]);
                        MYPRINT("\n");
                        MYPRINT("inverse of the residual covariance:\n");
                        print_matrix(inv_residual_cov[t][regime_j][regime_k]);
                        MYPRINT("\n");
                        MYPRINT("the residual covariance:\n");
                        print_matrix(residual_cov[t][regime_j][regime_k]);
                        MYPRINT("\n");*/


                   }

                   /*for(col_index=0; col_index<config->dim_latent_var; col_index++){
                        fprintf(eta_file, " %lf", gsl_vector_get(eta_regime_jk_t_plus_1[t][0][0], col_index));
                    }
                   fprintf(eta_file,"\n");*/

                   /*change, for random effect estimation*/
                   /*if (t==(config->index_sbj)[sbj+1]-1){
                        gsl_vector_set((init->eta_0)[regime_j],config->dim_latent_var*sbj+2,gsl_vector_get(eta_regime_jk_t_plus_1[t][regime_j][regime_k], 2));
                   }*/

				   
                   /** step 2: call hamilton filter to compute the probability of moving one step ahead **/

                   /** Step 2.1: compute transition probability matrix, Pr(S_{t-1} = j,S_{t} = k|Y_{t-1}) given the pr_t_1 **/
                   if (t==(config->index_sbj)[sbj]){
                   	gsl_matrix_set(tran_prob_jk, regime_j, regime_k, gsl_vector_get(pr_t[t], regime_j)*gsl_matrix_get(param->regime_switch_mat, regime_j, regime_k));
                   }else{
                   	gsl_matrix_set(tran_prob_jk, regime_j, regime_k, gsl_vector_get(pr_t[t-1], regime_j)*gsl_matrix_get(param->regime_switch_mat, regime_j, regime_k));
                   }
				   
                   gsl_vector_set(pr_t_given_t_minus_1[t], regime_k, 
				   		gsl_vector_get(pr_t_given_t_minus_1[t],regime_k) + gsl_matrix_get(tran_prob_jk, regime_j, regime_k));
                   
				   /*MYPRINT("prob_regime:\n");
                   print_vector(pr_t[t-1]);
                   MYPRINT("\n");*/

                   /** Step 2.2: compute log value of function f(.), i.e., prediction error decomposition function **/
                   neg_log_p=mathfunction_negloglike_multivariate_normal_invcov(innov_v[t][regime_j][regime_k], inv_residual_cov[t][regime_j][regime_k], y_non_miss, innov_determinant);

					double numNotMissingVars = mathfunction_sum_vector(y_non_miss);
					double tooSmallNumber = numNotMissingVars < 30 ? pow(1e-10, numNotMissingVars):1e-300;
					double tryP = exp(-neg_log_p);
					p = ( isfinite(tryP) && (tryP > tooSmallNumber) ) ? tryP:tooSmallNumber;

                   /*MYPRINT("likelihood f(y_it|S_it=k,S_i,t-1=j,Y_i,t-1):\n");
                   MYPRINT("oringinal %lf\n",exp(-neg_log_p));
                   MYPRINT("adjusted %lf\n",p);*/

                   /** compare the p with the (0.0001) and get the bigger one. We do not like probability that is too small. :)**/
                   /*p=exp(-neg_log_p)*tran_prob_jk;*/
                   gsl_matrix_set(like_jk, regime_j, regime_k, p*gsl_matrix_get(tran_prob_jk, regime_j, regime_k));
			   	   
				}/*end of to regime k*/
            }/*end of from regime j*/	 
		 
			
            /** Step 2.3: update transit probability Pr(S_{t-1} = j,S_{t} = k|Y_t) given Pr(S_{t-1} = j,S_{t} = k|Y_{t-1})**/
            log_like+=log(mathfunction_matrix_normalize(like_jk));/*like_jk scaled, Pr(S_{t-1} = j,S_{t} = k|Y_{t}, like=sum*/

            /*MYPRINT("Pr(S_{t-1} = j,S_{t} = k|Y_{t}):\n");
            print_matrix(like_jk);
            MYPRINT("\n");

            MYPRINT("negative log likelihood at time %lu:\n",t);
            MYPRINT("%lf",-log_like);
            MYPRINT("\n");

            MYPRINT("Pr(S_it=k|Y_it) oringinal and adjusted:\n");*/

            for(regime_k=0; regime_k<config->num_regime; regime_k++){
                    sum_overj=0;

            	    for(regime_j=0; regime_j<config->num_regime; regime_j++){
            	    	    sum_overj+=gsl_matrix_get(like_jk, regime_j, regime_k);
            	    }

            	    /** step 2.4: sum transition probability to obtain pr_t[t]**/
            	    gsl_vector_set(pr_t[t], regime_k, sum_overj);/*pr_t_plus_1*/

            	    /*MYPRINT("%lf ",sum_overj);*/
					
            }/*end of k*/
			
            /** step 2.4.1: check whether there is zero probability. If so, a small amount of value is added. Again we do not like too small and zero probability **/
			double tooSmallRegimeNumber = config->num_regime < 30 ? pow(1e-10, config->num_regime):1e-300;
	    	
			if(gsl_vector_min(pr_t[t]) < tooSmallRegimeNumber){
	        	gsl_vector_add_constant(pr_t[t], tooSmallRegimeNumber);
	        	mathfunction_vector_normalize(pr_t[t]);
	    		}	
				
				/*TODO same for pr_t_given_t_minus_1*/
			
				/* miss_case!=0; When there is missingness*/
				/*if (t!=(config->index_sbj)[sbj]){
				gsl_vector_memcpy(pr_t[t],pr_t[t-1]);
				}
				gsl_vector_memcpy(pr_t_given_t_minus_1[t],pr_t[t]);*/
			

			/*MYPRINT("\n");
			MYPRINT("miss_case: %lu\n",miss_case);
			MYPRINT("log_like: %f\n",log_like);
			print_vector(pr_t[t]);
			print_matrix(like_jk);
			MYPRINT("\n");*/
			
	        /** step 3: call collapse process **/
	        /** step 3.1: collapse the latent variable to get eta_k **/
				
			/* obtain eta_t and error_cov_t*/
			/** Other optional outputs of the Kalman Filter **/
			
	        for(regime_k=0; regime_k<config->num_regime; regime_k++){
				
                if (t==(config->index_sbj)[sbj]){
        	        
					gsl_vector_set_zero(eta_regime_j_t[t][regime_k]);
					gsl_matrix_set_zero(error_cov_regime_j_t[t][regime_k]);
                
				}

	            	    for(regime_j=0; regime_j<config->num_regime; regime_j++){

							gsl_blas_daxpy(gsl_matrix_get(like_jk, regime_j, regime_k), 
								eta_regime_jk_t_plus_1[t][regime_j][regime_k], 
								eta_regime_j_t[t][regime_k]); /*sum over j through loop*/
					   		   
							gsl_blas_daxpy(gsl_matrix_get(tran_prob_jk, regime_j, regime_k), 
								eta_regime_jk_pred[t][regime_j][regime_k], 
								eta_pred_regime_t[t][regime_k]); 
					   		
							gsl_blas_daxpy(gsl_matrix_get(tran_prob_jk, regime_j, regime_k), 
								innov_v[t][regime_j][regime_k], 
								innov_v_regime_t[t][regime_k]); 

	            	    	    /*if(regime_k==1){
	            	    	    MYPRINT("Here!");
	            	    	    print_vector(eta_regime_j_t[t][regime_k]);}*/
	            	    }
	        
			}/*end of k*/
			
            for(regime_k=0; regime_k<config->num_regime; regime_k++){
				
            	    gsl_vector_scale(eta_regime_j_t[t][regime_k], 1.0/gsl_vector_get(pr_t[t],regime_k));
					
					gsl_vector_scale(eta_pred_regime_t[t][regime_k], 1.0/gsl_vector_get(pr_t_given_t_minus_1[t],regime_k));
					
					gsl_vector_scale(innov_v_regime_t[t][regime_k], 1.0/gsl_vector_get(pr_t_given_t_minus_1[t],regime_k));

            	    /*MYPRINT("eta collapsed estimate in regime %lu:\n",regime_k);
            	    print_vector(eta_regime_j_t[t][regime_k]);
            	    MYPRINT("\n");*/

            	    /** step 3.2: collapse the covariance matrix to get error_cov_k **/
            	    for(regime_j=0; regime_j<config->num_regime; regime_j++){
						
						mathfunction_collapse(eta_regime_j_t[t][regime_k], eta_regime_jk_t_plus_1[t][regime_j][regime_k], 
						error_cov_regime_jk_t_plus_1[t][regime_j][regime_k], gsl_matrix_get(like_jk,regime_j, regime_k), 
						error_cov_regime_j_t[t][regime_k],
						diff_eta_vec, diff_eta, modif_p);
						
						mathfunction_collapse(eta_pred_regime_t[t][regime_k], eta_regime_jk_pred[t][regime_j][regime_k], 
						error_cov_regime_jk_pred[t][regime_j][regime_k], gsl_matrix_get(like_jk,regime_j, regime_k), 
						error_cov_pred_regime_t[t][regime_k],
						diff_eta_vec, diff_eta, modif_p);
						
						mathfunction_collapse(innov_v_regime_t[t][regime_k], innov_v[t][regime_j][regime_k], 
						residual_cov[t][regime_j][regime_k], gsl_matrix_get(like_jk,regime_j, regime_k), 
						residual_cov_regime_t[t][regime_k],
						diff_eta_vec, diff_eta, modif_p);

                    }/*end of j*/
                    
					gsl_matrix_scale(error_cov_regime_j_t[t][regime_k], 1.0/gsl_vector_get(pr_t[t],regime_k));
					
					gsl_matrix_scale(error_cov_pred_regime_t[t][regime_k], 1.0/gsl_vector_get(pr_t_given_t_minus_1[t],regime_k));
					
					gsl_matrix_scale(residual_cov_regime_t[t][regime_k], 1.0/gsl_vector_get(pr_t_given_t_minus_1[t],regime_k));

                    /*MYPRINT("error_cov collapsed estimate in regime %lu:\n",regime_k);
            	    print_matrix(error_cov_regime_j_t[t][regime_k]);
            	    MYPRINT("\n");*/
            }/*end of k*/
			
			/*gsl_vector_set_zero(eta_t[t]);*/
			/*gsl_matrix_set_zero(error_cov_t[t]);*/
    	    
			for(regime_k=0; regime_k<config->num_regime; regime_k++){
 
    	    	    gsl_blas_daxpy(gsl_vector_get(pr_t[t], regime_k), eta_regime_j_t[t][regime_k], eta_t[t]); /*sum over k through loop*/
					
    	    	    gsl_blas_daxpy(gsl_vector_get(pr_t_given_t_minus_1[t], regime_k), eta_pred_regime_t[t][regime_k], eta_pred_t[t]); /*sum over k through loop*/
					
    	    	    gsl_blas_daxpy(gsl_vector_get(pr_t_given_t_minus_1[t], regime_k), innov_v_regime_t[t][regime_k], innov_v_t[t]); /*sum over k through loop*/
					
					/*if(regime_k==1){
    	    	    MYPRINT("Here!");
    	    	    print_vector(eta_t[t]);}*/
    	    }
    	    
			for(regime_k=0; regime_k<config->num_regime; regime_k++){
					
					mathfunction_collapse(eta_t[t], eta_regime_j_t[t][regime_k], 
					error_cov_regime_j_t[t][regime_k], gsl_vector_get(pr_t[t], regime_k), error_cov_t[t],
					diff_eta_vec, diff_eta, modif_p);
					
					mathfunction_collapse(eta_pred_t[t], eta_pred_regime_t[t][regime_k], 
					error_cov_pred_regime_t[t][regime_k], gsl_vector_get(pr_t_given_t_minus_1[t], regime_k), error_cov_pred_t[t],
					diff_eta_vec, diff_eta, modif_p);
					
					mathfunction_collapse(innov_v_t[t], innov_v_regime_t[t][regime_k], 
					residual_cov_regime_t[t][regime_k], gsl_vector_get(pr_t_given_t_minus_1[t], regime_k), residual_cov_t[t],
					diff_eta_vec, diff_eta, modif_p);

			}
	

	   /*fprintf(pr_file,"%lu %lu %lf %lf\n",sbj,t,gsl_vector_get(pr_t[t],0),gsl_vector_get(pr_t[t],1));*/
	   	   /*MYPRINT("Pr_St|t:\n");
           print_vector(pr_t[t]);
           MYPRINT("\n");*/

        }/*end of t*/

       /*if (sbj==2){exit(0);}*/
         /*fprintf(h_file, "%d", t+1);*/
		
    }/*end of sbj*/

    /*fclose(h_file);*/
    /*fclose(eta_file);*/
    /*fclose(pr_file);*/
    /****************************** free allocated space ***************************/

    gsl_matrix_free(tran_prob_jk);
    gsl_matrix_free(like_jk);
	gsl_vector_free(cp_y_t);
	gsl_vector_free(y_non_miss);	


    gsl_vector_free(diff_eta_vec);
    gsl_matrix_free(diff_eta);
    gsl_matrix_free(modif_p);

	/** output of extended Kalman filter **/
    /*eta^regime_jk_it|t*/
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            for(regime_k=0; regime_k<config->num_regime; regime_k++)
            gsl_vector_free(eta_regime_jk_t_plus_1[index_sbj_t][regime_j][regime_k]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            free(eta_regime_jk_t_plus_1[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        free(eta_regime_jk_t_plus_1[index_sbj_t]);
    }
    free(eta_regime_jk_t_plus_1);

    /*error_cov^regime_jk_it|t*/
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            for(regime_k=0; regime_k<config->num_regime; regime_k++)
                gsl_matrix_free(error_cov_regime_jk_t_plus_1[index_sbj_t][regime_j][regime_k]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            free(error_cov_regime_jk_t_plus_1[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        free(error_cov_regime_jk_t_plus_1[index_sbj_t]);
    }
    free(error_cov_regime_jk_t_plus_1);
	
    /*output of filter: innovation vector*/
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            for(regime_k=0; regime_k<config->num_regime; regime_k++)
            gsl_vector_free(innov_v[index_sbj_t][regime_j][regime_k]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            free(innov_v[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        free(innov_v[index_sbj_t]);
    }
    free(innov_v);

    /*output of filter: inverse_residual_cov*/
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            for(regime_k=0; regime_k<config->num_regime; regime_k++)
                gsl_matrix_free(inv_residual_cov[index_sbj_t][regime_j][regime_k]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            free(inv_residual_cov[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        free(inv_residual_cov[index_sbj_t]);
    }
    free(inv_residual_cov);
	
    /*output of filter: residual_cov*/
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            for(regime_k=0; regime_k<config->num_regime; regime_k++)
                gsl_matrix_free(residual_cov[index_sbj_t][regime_j][regime_k]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            free(residual_cov[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        free(residual_cov[index_sbj_t]);
    }
    free(residual_cov);
	
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            gsl_vector_free(eta_pred_regime_t[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        free(eta_pred_regime_t[index_sbj_t]);
    }
    free(eta_pred_regime_t);

    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            gsl_matrix_free(error_cov_pred_regime_t[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        free(error_cov_pred_regime_t[index_sbj_t]);
    }
    free(error_cov_pred_regime_t);
	
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            gsl_vector_free(innov_v_regime_t[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        free(innov_v_regime_t[index_sbj_t]);
    }
    free(innov_v_regime_t);

    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            gsl_matrix_free(residual_cov_regime_t[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        free(residual_cov_regime_t[index_sbj_t]);
    }
    free(residual_cov_regime_t);

    return(-log_like);
}



/****************************Extended Kim Smoother************************/
/* *
* This function implements the extended Kim Smoother
* Parameters/Input *
* *
* **>>Parameters/Input Pointers<<**
* y_time -- continous time points in the data
* co_variate -- covariates
* config -- model configuration,
* param -- parameters
* pr_t_given_t_minus_1 -- Pr(S_it=k|Y_i,t-1) *
* pr_t -- Pr(S_it=k|Y_it) *
* eta_regime_jk_pred -- eta^regime_jk_it|t-1 *
* error_cov_regime_jk_pred -- error_cov^regime_jk_it|t-1 *
* eta_regime_j_t -- eta^k_it|t*
* error_cov_regime_j_t -- error_cov^k_it|t*
* *
* Output*
* *
* **>>Output via using pointers<<**
* transprob_T -- Pr(S_i,t+1=h, S_it=k|Y_iT)*
* pr_T -- Pr(S_it=k|Y_iT)*
* eta_smooth -- eta_it|T *
* error_cov_smooth -- error_cov_it|T *
* *
* */

void EKimSmoother(double *y_time, gsl_vector **co_variate, const ParamConfig *config, const Param *param,
    gsl_vector **pr_t_given_t_minus_1, gsl_vector **pr_t, 
	gsl_vector ****eta_regime_jk_pred, gsl_matrix ****error_cov_regime_jk_pred,
	gsl_vector ***eta_regime_j_t, gsl_matrix ***error_cov_regime_j_t,
	gsl_vector **eta_smooth, gsl_matrix **error_cov_smooth,
	gsl_vector **pr_T, gsl_vector ***transprob_T){

    /**initialization**/
    gsl_vector *temp_diff_eta_vec=gsl_vector_alloc(config->dim_latent_var);
    gsl_matrix *temp_diff_eta=gsl_matrix_alloc(config->dim_latent_var, 1);
    gsl_matrix *temp_modif_p=gsl_matrix_alloc(config->dim_latent_var, config->dim_latent_var);
    /*Pr[S_i,t+1=regime_k|Y_iT]*/
    gsl_vector *p_next_regime_T=gsl_vector_alloc(config->num_regime);
    gsl_matrix *Jacob_dyn_x=gsl_matrix_calloc(config->dim_latent_var, config->dim_latent_var);
    size_t sbj, t, index_sbj_t, regime_j,regime_k;
    /*size_t i;
      double params_aug[config->num_func_param+config->dim_latent_var];
        for (i=0;i<config->num_func_param;i++)
            params_aug[i]=param->func_param[i];*/
    gsl_matrix *P_tilde_regime_jk=gsl_matrix_alloc(config->dim_latent_var,config->dim_latent_var);
    gsl_matrix *pb=gsl_matrix_alloc(config->dim_latent_var,config->dim_latent_var);
    gsl_matrix *inv_P_jk_pred=gsl_matrix_alloc(config->dim_latent_var,config->dim_latent_var);
    double sum_overk;

    gsl_matrix *eta_regime_jk_T=gsl_matrix_alloc(config->dim_latent_var,1);
    gsl_vector *eta_regime_jk_T_vec=gsl_vector_alloc(config->dim_latent_var);
    gsl_matrix *error_cov_regime_jk_T=gsl_matrix_alloc(config->dim_latent_var,config->dim_latent_var);
    gsl_matrix *temp_diff_P=gsl_matrix_alloc(config->dim_latent_var,config->dim_latent_var);

    /*eta^k_it|T*/
    gsl_vector ***eta_regime_j_smooth=(gsl_vector ***)malloc(config->total_obs*sizeof(gsl_vector **));
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	eta_regime_j_smooth[index_sbj_t]=(gsl_vector **)malloc(config->num_regime*sizeof(gsl_vector *));
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    eta_regime_j_smooth[index_sbj_t][regime_j]=gsl_vector_calloc(config->dim_latent_var);
	}
    }
    /*error_cov^k_it|T*/
    gsl_matrix ***error_cov_regime_j_smooth=(gsl_matrix ***)malloc(config->total_obs*sizeof(gsl_matrix **));
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	error_cov_regime_j_smooth[index_sbj_t]=(gsl_matrix **)malloc(config->num_regime*sizeof(gsl_matrix *));
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
	for(regime_j=0; regime_j<config->num_regime; regime_j++){
	    error_cov_regime_j_smooth[index_sbj_t][regime_j]=gsl_matrix_calloc(config->dim_latent_var,config->dim_latent_var);
	}
    }

    for(sbj=0; sbj<config->num_sbj; sbj++){/*start of the sbj loop*/

        /**set eta_regime_j_smooth and error_cov_regime_j_smooth at time T to filtered estimates at time T**/

        t=(config->index_sbj)[sbj+1]-1;
        gsl_vector_memcpy(pr_T[t],pr_t[t]);

        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            gsl_vector_memcpy(eta_regime_j_smooth[t][regime_j],eta_regime_j_t[t][regime_j]);
            gsl_blas_daxpy(gsl_vector_get(pr_T[t],regime_j), eta_regime_j_smooth[t][regime_j], eta_smooth[t]);/*eta_smooth[t]=eta_smooth[t]+gsl_vector_get(pr_T[t],regime_j)*eta_regime_j_smooth[t][regime_j]*/


            gsl_matrix_memcpy(error_cov_regime_j_smooth[t][regime_j], error_cov_regime_j_t[t][regime_j]);
			
			mathfunction_collapse(eta_smooth[t], eta_regime_j_smooth[t][regime_j], 
			error_cov_regime_j_smooth[t][regime_j], gsl_vector_get(pr_T[t],regime_j), error_cov_smooth[t],
			temp_diff_eta_vec, temp_diff_eta, temp_modif_p);

        }

        /**start iteration**/

        for(t=(config->index_sbj)[sbj+1]-1;t-->(config->index_sbj)[sbj];){/*Notice the decreasing for loop for size_t t*/
            gsl_vector_memcpy(p_next_regime_T,pr_T[t+1]);

            /**set the regime switch matrix**/
            config->func_regime_switch(t+1, 1, param->func_param, co_variate[t+1], param->regime_switch_mat);/*type=1*/

            for(regime_j=0; regime_j<config->num_regime; regime_j++){/*from regime regime_j*/
                sum_overk=0;
		 	   	
				/*MYPRINT("pr_t_given_t_minus_1[t+1]:\n");
		        print_vector(pr_t_given_t_minus_1[t+1]);
		        MYPRINT("\n");*/
				
                for(regime_k=0; regime_k<config->num_regime; regime_k++){/*to regime k*/

                    /*MYPRINT("sbj %lu t %lu from regime %lu to regime %lu\n",sbj,t,regime_j,regime_k);
                    MYPRINT("\n");*/


					/*Cf. Chow & Zhang Equation A.9*/
					/*Check for division zero (near) zero*/
                    /*Pr[S_i,t+1=regime_k, S_it=regime_j|Y_iT]*/
                    gsl_vector_set(transprob_T[t][regime_j],regime_k, gsl_vector_get(p_next_regime_T,regime_k)*gsl_vector_get(pr_t[t],regime_j)*gsl_matrix_get(param->regime_switch_mat,regime_j,regime_k)/gsl_vector_get(pr_t_given_t_minus_1[t+1],regime_k));
                    

					/*Pr[S_it=j|Y_iT] sum over k*/
					/*Cf. Denominator of A.9*/
                    sum_overk+=gsl_vector_get(transprob_T[t][regime_j],regime_k);
					/*print this sum_overk , i.e. check denominator not to near zero*/


                    /*Jacobian matrix of the dynamic function*/
                    /*Notice that the parameters input into function_dF_dx and function_dP_dt*/
                    /*for (i=0;i<config->dim_latent_var;i++)
                        params_aug[config->num_func_param+i]=gsl_vector_get(eta_regime_j_t[t][regime_k],i);*/
                    config->func_jacob_dynam(y_time[t],y_time[t+1],regime_k,eta_regime_j_t[t][regime_k],param->func_param,config->num_func_param, co_variate[t],config->func_dF_dx, Jacob_dyn_x);

                    /*P_tilde_regime_jk=error_cov_regime_j_t[t][regime_j] %*% Jacob_dyn_x %*% inv(error_cov_regime_jk_pred[t+1][regime_j][regime_k])*/

                    gsl_matrix_set_zero(P_tilde_regime_jk);
                    gsl_matrix_set_zero(pb);
                    gsl_matrix_set_zero(inv_P_jk_pred);
                    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, error_cov_regime_j_t[t][regime_j], Jacob_dyn_x, 0.0, pb); /* compute P*B'*/

                    /*MYPRINT("sbj %lu t-1 %lu \n",sbj,t);
                    print_matrix(error_cov_regime_jk_pred[t+1][regime_j][regime_k]);
                    MYPRINT("\n");*/
					
                    mathfunction_inv_matrix(error_cov_regime_jk_pred[t+1][regime_j][regime_k], inv_P_jk_pred);/*obtain the inverse matrix*/
					
					/*print_matrix(inv_P_jk_pred);
                    MYPRINT("\n");*/
					
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, pb, inv_P_jk_pred, 0.0, P_tilde_regime_jk); /* compute P*B*Pjk^{-1}*/
                    
					/*print_matrix(P_tilde_regime_jk);
                    MYPRINT("\n");*/

                    /*eta_regime_jk_T*/
                    /* compute eta_regime_j_t[t+1][regime_k] - eta_regime_jk_pred[t+1][regime_j][regime_k]*/
                    gsl_vector_memcpy(temp_diff_eta_vec, eta_regime_j_t[t+1][regime_k]);
                    gsl_vector_sub(temp_diff_eta_vec, eta_regime_jk_pred[t+1][regime_j][regime_k]);
                    gsl_matrix_set_col(temp_diff_eta, 0, temp_diff_eta_vec);
                    /*eta_regime_j_t[t][regime_j]+P_tilde_regime_jk%*%temp_diff_eta*/
                    gsl_matrix_set_col(eta_regime_jk_T,0, eta_regime_j_t[t][regime_j]);
                    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, P_tilde_regime_jk, temp_diff_eta, 1, eta_regime_jk_T);
                    gsl_matrix_get_col(eta_regime_jk_T_vec, eta_regime_jk_T, 0);


                    /*error_cov_regime_jk_T*/
                    /* compute error_cov_regime_j_t[t+1][regime_k] - error_cov_regime_jk_pred[t+1][regime_j][regime_k]*/
                    gsl_matrix_memcpy(temp_diff_P, error_cov_regime_j_t[t+1][regime_k]);
                    gsl_matrix_sub(temp_diff_P, error_cov_regime_jk_pred[t+1][regime_j][regime_k]);
                    /*error_cov_regime_j_t[t][regime_j]+P_tilde_regime_jk%*%temp_modif_p%*%P_tilde_regime_jk*/
                    gsl_matrix_memcpy(error_cov_regime_jk_T,error_cov_regime_j_t[t][regime_j]);
                    gsl_matrix_set_zero(temp_modif_p);
                    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, P_tilde_regime_jk, temp_diff_P, 0.0, temp_modif_p);
                    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, temp_modif_p,P_tilde_regime_jk, 1, error_cov_regime_jk_T);



                    /*eta_regime_j_smooth[t][regime_j] = eta_regime_j_smooth[t][regime_j] + gsl_vector_get(transprob_T[t][regime_j],regime_k)*eta_regime_jk_T*/
                    gsl_blas_daxpy(gsl_vector_get(transprob_T[t][regime_j],regime_k), eta_regime_jk_T_vec, eta_regime_j_smooth[t][regime_j]);;

                    /*error_cov_regime_j_smooth[t][regime_j]*/
					mathfunction_collapse(eta_regime_j_smooth[t][regime_j], eta_regime_jk_T_vec, 
					error_cov_regime_jk_T, gsl_vector_get(transprob_T[t][regime_j],regime_k), error_cov_regime_j_smooth[t][regime_j],
					temp_diff_eta_vec, temp_diff_eta, temp_modif_p);

                    /*if(t==998){
                    print_matrix(Jacob_dyn_x);
                    MYPRINT("\n");
                    print_matrix(P_tilde_regime_jk);
                    MYPRINT("\n");
                    print_matrix(eta_regime_jk_T);
                    MYPRINT("\n");
                    print_matrix(error_cov_regime_jk_T);
                    MYPRINT("\n");
                    }*/

               }/*end of to regime k*/


               if (sum_overk>1){
                   gsl_vector_scale(transprob_T[t][regime_j],1/sum_overk);
                   gsl_vector_set(pr_T[t], regime_j, 1.0);
               }else{
                   gsl_vector_set(pr_T[t], regime_j, sum_overk);
               }

               gsl_vector_scale(eta_regime_j_smooth[t][regime_j],1/sum_overk);
               gsl_matrix_scale(error_cov_regime_j_smooth[t][regime_j],1/sum_overk);
   			
			/*MYPRINT("\n");
   			MYPRINT("sum_overk: %f\n",sum_overk);*/

            }/*end of from regime regime_j*/
   			
			/*MYPRINT("pr_T[t]:\n");
   			print_vector(pr_T[t]);
			MYPRINT("\n");*/
			/*MYPRINT("eta_regime_j_smooth[t]:\n");
   			print_vector(eta_regime_j_smooth[t]);
			MYPRINT("\n");
			MYPRINT("error_cov_regime_j_smooth[t]:\n");
   			print_matrix(error_cov_regime_j_smooth[t]);
   			MYPRINT("\n");*/

            /*Final collapse over j*/
            for(regime_j=0; regime_j<config->num_regime; regime_j++){

                gsl_blas_daxpy(gsl_vector_get(pr_T[t],regime_j), eta_regime_j_smooth[t][regime_j], eta_smooth[t]);/*eta_smooth[t]=eta_smooth[t]+gsl_vector_get(pr_T[t],regime_j)*eta_regime_j_smooth[t][regime_j]*/

				mathfunction_collapse(eta_smooth[t], eta_regime_j_smooth[t][regime_j], 
				error_cov_regime_j_smooth[t][regime_j], gsl_vector_get(pr_T[t],regime_j), error_cov_smooth[t],
				temp_diff_eta_vec, temp_diff_eta, temp_modif_p);

            }
			
			/*MYPRINT("eta_smooth[t]:\n");
   			print_vector(eta_smooth[t]);
			MYPRINT("\n");
			MYPRINT("error_cov_smooth[t]:\n");
   			print_matrix(error_cov_smooth[t]);
			MYPRINT("\n");*/
			
        }/*end of the t loop*/
    }/*end of the sbj loop*/


    /**free allocated space**/
    gsl_vector_free(temp_diff_eta_vec);
    gsl_matrix_free(temp_diff_eta);
    gsl_matrix_free(temp_modif_p);
    /*Pr[S_i,t+1=regime_k|Y_iT]*/
    gsl_vector_free(p_next_regime_T);
    gsl_matrix_free(Jacob_dyn_x);

    gsl_matrix_free(P_tilde_regime_jk);
    gsl_matrix_free(pb);
    gsl_matrix_free(inv_P_jk_pred);

    gsl_matrix_free(eta_regime_jk_T);
    gsl_vector_free(eta_regime_jk_T_vec);
    gsl_matrix_free(error_cov_regime_jk_T);
    gsl_matrix_free(temp_diff_P);
	
    /*output of smooth: eta^k_it|T*/
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            gsl_vector_free(eta_regime_j_smooth[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        free(eta_regime_j_smooth[index_sbj_t]);
    }
    free(eta_regime_j_smooth);


    /*output of smooth: error_cov^k_it|T*/
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<config->num_regime; regime_j++){
            gsl_matrix_free(error_cov_regime_j_smooth[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<config->total_obs;index_sbj_t++){
        free(error_cov_regime_j_smooth[index_sbj_t]);
    }
    free(error_cov_regime_j_smooth);
	

}/*end of function EKimSmoother*/



