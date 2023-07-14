#include <time.h>
//#include <stdlib.h>
//#include<math.h>
  
//#include <armadillo>
  
//using namespace arma;

//#include "structure_prototype.h"
//#include "supplementary_function.h"
//#include "converted_function.h"
  
  
  
// Step 3 in the MainUseThins.m
void saem_estimation(C_INFDS &InfDS, C_INFDS0 &InfDS0, arma::mat upperb, arma::mat lowerb, arma::mat x1, char *filenamePar, char *filenameSE, char *filenameconv, char *filenamebhat, char *filenamebhat2, int kk, int trueInit, int batch, int seed, int freeIC, struct C_OUTPUT &output, int observedFlag){
  //Rprintf("in MainUseThis\n");
  
  arma::mat sgnTH, L, QQ, D, mscore2, OMEGAb, infoMat, minfoMat, meanb,tpOld, score, Covscore;
  arma::vec mscore;
  int k, stage, gmm, MAXGIB, setScaleb, noIncrease, yesMean, switchFlag, useMultN, GIB, stop, isBlock1Only, redFlag, convFlag, k2;
  double bAccept, ss, ttt, ssmin;
  int prev_stage;
  int is_meanb = 0; //Not use. To get rid of later. Old isPar.
  time_t timer;
  //int i, j;
  
  //C_OUTPUT output;
  //--for writing output files--
    //int i, j, fitInit;
  //FILE *p_filenamePar, *p_filenameSE, *p_filenameconv, *p_filenamebhat, *p_filenamebhat2;
  //----
    
	//InfDS.y0.print("y0 at the start of saem_estimation");
	
    //freeIC = 1;
    timer = time(NULL);
    
    //don't use this to avoid warning
	//arma_rng::set_seed(seed);
	

	//InfDS.Xtild_12.reset();
	output.avebAccept = 0;
	isBlock1Only = 0;	
	switchFlag = 0;
	//upperb = "5,5,5"; 
	//lowerb = "-5,-5,-5";
	
	// ------- 
	yesMean = 0;
	
	//InfDS.par = join_vert(x1, InfDS.par);
	//InfDS.par.print("par");

	InfDS.Nx = InfDS.NxState;
	InfDS.G = eye(InfDS.Nx, InfDS.Nx);
	//InfDS.Nbeta = 5;
	//InfDS.Ntheta = 3;
	InfDS.alp = 1;
	MAXGIB = InfDS.MAXGIB;
	
	//Rprintf("checkpoint M32\n");	
	//to be checked
	/*
	InfDS.odefn    = @dynfunIC;      % Function for ODE solver 
	InfDS.dfdx = @dfdxFreeIC;
	InfDS.dfdxdp = @dfdxdpFreeIC;
	InfDS.dfdpdx = @dfdpdxFreeIC;
	InfDS.dfdp = @dfdparFreeIC;
	InfDS.dfdp2 = @dfdpar2FreeIC;
	InfDS.dfdx2 = @dfdx2FreeIC;
	*/
	
	
	//New set parameter: 
	//dXtildAll contains dXtild for all time points, which is a field with totalT points. Each cell is a cube (dimension the same as dXtild).
	//d2XtildAll contains d2Xtild for all time points, which is a field with totalT points. Each cell is a cube (dimension the same as d2Xtild).
	//Rprintf("InfDS.tspan.n_cols = %d\n", InfDS.tspan.n_cols);
	InfDS.dXtildAll.set_size(InfDS.tspan.n_cols);
	InfDS.d2XtildAll.set_size(InfDS.tspan.n_cols);
	
	InfDS.sy = arma::zeros<arma::mat>(InfDS.par.n_elem, 1);
	InfDS.ES = arma::zeros<arma::mat>(InfDS.par.n_elem, InfDS.par.n_elem);
	InfDS.EI = arma::zeros<arma::mat>(InfDS.par.n_elem,InfDS.par.n_elem);
	sgnTH = arma::zeros<arma::mat>(InfDS.par.n_elem,InfDS.maxIterStage1);
	InfDS.Iy = arma::zeros<arma::mat>(InfDS.par.n_elem,InfDS.par.n_elem);
	InfDS.thetatild = arma::zeros<arma::mat>(InfDS.par.n_elem,1);
	InfDS.sytild = arma::zeros<arma::mat>(InfDS.par.n_elem,1); 
	InfDS.EStild = arma::zeros<arma::mat>(InfDS.par.n_elem,InfDS.par.n_elem);
	InfDS.EItild = arma::zeros<arma::mat>(InfDS.par.n_elem,InfDS.par.n_elem);
	InfDS.Iytild = arma::zeros<arma::mat>(InfDS.par.n_elem,InfDS.par.n_elem);
	

	//Rprintf("InfDS.par.n_elem = %d \n", InfDS.par.n_elem);	

	/*
	InfDS.lowBound = arma::ones<arma::mat>(InfDS.par.n_elem,1);
	InfDS.lowBound.fill(10e-8);
	InfDS.upBound = arma::ones<arma::mat>(InfDS.par.n_elem,1);
	InfDS.upBound.fill(10);
	*/
	/*
	InfDS.lowBound(span(8,16), span::all).fill(-10);
	InfDS.upBound(span(8,16), span::all).fill(10);
	InfDS.lowBound.print("InfDS.lowBound");
	InfDS.upBound.print("InfDS.upBound");
	*/
	
	//Rprintf("checkpoint M63\n");	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%self change, with caution and ask
	//InfDS.lowBound(span(8,9), span::all).fill(10e-8);
	//InfDS.lowBound(span(10, 15), span::all).fill(log(10e-8));
	//InfDS.lowBound(16) = -2; //After transformation, covariance between -5 and 5

	//InfDS.upBound(span(8,9), span::all).fill(10);
	//InfDS.upBound(span(10, 15), span::all).fill(2);
	
	//Rprintf("checkpoint M73\n");	

	k2 = 1;
	k = 1; 
	stage = 1; 
	gmm = 1;
	bAccept = 0; // to avoid uninitialized warning
	
	setScaleb = InfDS.setScaleb;
	//InfDS.bAdaptParams = ".5, 2.5, .5";
	ssmin = 100; 
	noIncrease = 0;
	//freeIC = 1;
	//InfDS.scaleb = 10; //Used in drawbGenera6_opt3.m to determine whether to apply scaling constant on drawb.
	//InfDS.KKO = 20; //Used in SAEM. Only starts to evaluate whether to transition to stage 2 after KKO iterations.
	stop = 0;
    //Rprintf("check point 2 MAXITER %d freeIC %d\n", InfDS.MAXITER, freeIC);
	
	//InfDS.par.print("InfDS.par");
	//Rprintf("checkpoint M92 entering the k loop\n");	
	int getdxFlag = 1;
	if (InfDS.Nbeta > 0) getdxFlag=1;
	InfDS.meanb = InfDS0.trueb; //TEMP For debugging - REMOVE LATER
	InfDS.b = InfDS0.trueb;  //TEMP For debugging
	InfDS.useb = InfDS.meanb;
	InfDS = getXtildIC3(0, getdxFlag ,freeIC, InfDS); //%Get Xtilde and dXtilddthetaf before saem
	C_INFDS InfDS_meanb = InfDS;
	InfDS.Xtild_meanb = InfDS.Xtild;
	InfDS.dXtildthetafAll_meanb = InfDS.dXtildthetafAll;
	InfDS.dXtildthetafAll2_meanb = InfDS.dXtildthetafAll2;
	int Npar = InfDS.Nmu + InfDS.Ny + InfDS.NLambda + InfDS.Nbeta + InfDS.Nbpar;
	
	while (k <= InfDS.MAXITER && stop == 0){

		//disp 'iteration';
		Rprintf("k = %d\n",k);
		
//if(k >= 1){
			// in the first iteration we adopt the parameters from dynr interface
			//if(k > 1)
			InfDS.fp.setParsFreeICwb(InfDS); //qqqq	
			Rprintf("The matrices at the beginning of iteration %d: \n",k);
			InfDS.par.print("InfDS.par");
			InfDS.Sigmab.print("Sigmab");
			InfDS.Sigmae.print("InfDS.Sigmae");
					
		if (stage==2 && switchFlag==0){
			switchFlag = 1; 
		}

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
		Rprintf("Observed flag = %d\n", observedFlag);
		unsigned int k_stop;
		if (observedFlag) k_stop = kk;
		else k_stop = gmm;
		yesMean = 1;
		//if (stage==1 && bAccept < .1){ //Not looking too good. Pump up no. of chains
		if (((observedFlag == 0 || k <= 5) || bAccept < .001 || bAccept > .99)){ //TEMP; CHANGE BACK TO FOLLOWING AFTER DEBUGGING
//		if ((observedFlag == 0 && k <= 5) || bAccept < .001 || bAccept > .99){ 
            useMultN = 1; // Lu modified, 04-12-13,5;
		        MAXGIB=10;
		}
		//else if(stage==1 && k == 3){
		//SMC commented out 7/5/22
		//else if(observedFlag==0 && k > 5 && k <= InfDS.KKO){
    //        useMultN = 1; MAXGIB=20;
		//
		else {
			useMultN = 0; 
			MAXGIB = InfDS.MAXGIB;
		}
		
		if (InfDS.Nbpar ==0 && InfDS.Nbeta ==0){
		  useMultN = 0;
		  MAXGIB = 1;
		}

		mscore = arma::zeros<arma::mat>(InfDS.par.n_elem, 1);
		mscore2 = arma::zeros<arma::mat>(InfDS.par.n_elem, InfDS.par.n_elem);
		minfoMat = arma::zeros<arma::mat>(InfDS.par.n_elem, InfDS.par.n_elem);  
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//%% Gibbs sampler
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
		//Since beta changes, get Xtild and InfDS.InfDS.dxstardthetafAll and InfDS.dxstardthetafAll2
	
		//InfDS.y0.print("y0 before setting useb"); //already odd
		if (InfDS.Nbpar > 0 || InfDS.Nbeta > 0){
		Rprintf("Set useb to meanb in saem_estimation.h\n");
		//Rprintf("End of getXtildIC3 with meanb\n");
		//if(k == 1){
		//	InfDS.Xtild_p1 = InfDS.Xtild;
		//	InfDS.dXtild_p1 = InfDS.dXtild;
		//	InfDS.d2Xtild_p1 = InfDS.d2Xtild;
    //
		//}
		//else if (stage == 2 && InfDS.Xtild_12.is_empty()){			
		//	//copy Xtild, dXtild, and d2Xtild before transitting to stage 2
		//	InfDS.Xtild_12 = InfDS.Xtild;
		//	InfDS.dXtild_12 = InfDS.dXtild;
		//	InfDS.d2Xtild_12 = InfDS.d2Xtild;
		//	//InfDS.d2Xtild.print("InfDS.d2Xtild");
		//	
		//}
    Rprintf("\nGet covariance of proposal distribution of b");
		PropSigb(InfDS);  
    Rprintf("\nEvaluate likelihood function under old b");
		tpOld = ekfContinuous10(InfDS.Nsubj, InfDS.N, InfDS.Ny, InfDS.Nx, InfDS.Nb, InfDS.NxState, InfDS.Lambda, InfDS.totalT, InfDS.Sigmae, InfDS.Sigmab, InfDS.mu, InfDS.b, InfDS.allT, InfDS.Xtild, InfDS.Y); //%get density of full conditional distribution of b 
		
		InfDS.bacc = arma::zeros<arma::mat>(InfDS.Nsubj,1);	
		meanb = arma::zeros<arma::mat>(InfDS.Nsubj, InfDS.Nb); //Holds running MCMC means for b
		}
		
		Rprintf("[DEBUG] MAXGIB = %d \n", MAXGIB);
 
		for(GIB = 1; GIB <= MAXGIB; GIB++){
			  //Rprintf("GIB = %d\n", GIB);
				//Rprintf("checkpoint enter drowbGeneral6_opt3\n");	
				if (InfDS.Nbpar > 0 || InfDS.Nbeta > 0){
				//temporarility commnet out drawbGeneral6_opt3 function
				//drawbGeneral6_opt3(is_meanb, InfDS, yesMean, meanb, upperb, lowerb, useMultN, tpOld, freeIC, isBlock1Only, setScaleb, bAccept, MAXGIB);
				//Rprintf("checkpoint leave drowbGeneral6_opt3\n");	
				arma::mat cor_b = cor(InfDS.b,InfDS0.trueb);
				cor_b.print("Correlation between b and trueb");
				Rprintf("\nRange of true b estimates: [%lf, %lf]\n", InfDS0.trueb.min(), InfDS0.trueb.max());
				Rprintf("\nRange of b estimates: [%lf, %lf]\n", InfDS.b.min(), InfDS.b.max());
				}else{
				  InfDS.useb = InfDS.b;
				  InfDS = getXtildIC3(0, getdxFlag ,freeIC, InfDS); //%Get Xtilde and dXtilddthetaf before saem
				}
				
				
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			//Rprintf("checkpoint enter getScoreInfoY_tobs_opt\n");	
			getScoreInfoY_tobs_opt(InfDS, stage, k, freeIC, score, infoMat);
			//Rprintf("checkpoint leave getScoreInfoY_tobs_opt\n");		
			
			//score.print("Score");
			//infoMat.print("InfoMat");
			
			mscore = mscore + (1.0/MAXGIB)*score;
			mscore2 = mscore2 + (1.0/MAXGIB)*(score*score.t());
			minfoMat = minfoMat + (1.0/MAXGIB)*infoMat;
		} //end of Gibbs sampler loop
		
		if (InfDS.Nbpar > 0 || InfDS.Nbeta > 0){
		  //TEMP DEBUG - UNCOMMENT THIS WHOLE BLOCK LATER
		  //----------------------------------------------
		  //InfDS.meanb = meanb; 
		  //InfDS.useb = InfDS.meanb;
		  ////Rprintf("Start of getXtildIC3 with meanb\n");
		  //InfDS_meanb = getXtildIC3(0, 1 ,freeIC, InfDS); //%Get updated Xtilde and dXtilddthetaf
		  //InfDS.Xtild_meanb = InfDS_meanb.Xtild;
		  //InfDS.dXtildthetafAll_meanb = InfDS_meanb.dXtildthetafAll;
		  //InfDS.dXtildthetafAll2_meanb = InfDS_meanb.dXtildthetafAll2;
		  //TEMP DEBUG - END OF UNCOMMENT THIS WHOLE BLOCK LATER
		  //----------------------------------------------
		arma::mat cor_b = cor(InfDS.meanb,InfDS0.trueb);
		cor_b.print("Correlation between meanb and trueb");
		Rprintf("\nRange of true b estimates: [%lf, %lf]\n", InfDS0.trueb.min(), InfDS0.trueb.max());
		Rprintf("\nRange of meanb estimates: [%lf, %lf]\n", InfDS.meanb.min(), InfDS.meanb.max());
		Rprintf("\nCurrent b acceptance rate = %6f\n",bAccept);
		output.avebAccept += bAccept;
		}
		
		Covscore = mscore2 - mscore*mscore.t();
		if (Npar > 0){
		saem(InfDS, gmm, k, stage, redFlag, convFlag, noIncrease, stop, ssmin, ss, sgnTH, mscore, mscore2, minfoMat, Covscore,k_stop);
		
		Rprintf("\nStage = %5d, iteration = %5d\n",stage,k);
		//temporarily printing out messages
		if(1 || prev_stage != stage){
			//Rprintf("length of InfDS.par: %d\n", InfDS.par.n_elem);
			InfDS.par.print("InfDS.par (free)");
			InfDS.thetatild.print("InfDS.thetatild (free)");
			//InfDS.par(span(0,7), span::all).print("InfDS.par(1:8)");
			//exp(InfDS.par(span(8,11), span::all)).t().print("InfDS.par(9:12)");
			//Rprintf("Averaging:\n");
			//InfDS.thetatild(span(0,7), span::all).print("InfDS.thetatild(1:8)");
			//exp(InfDS.thetatild(span(8,11), span::all)).t().print("InfDS.par(9:12)");
			//D = diagmat(exp(InfDS.thetatild(span(12,14), span::all)));
			//L = "1 0 0;	0 1 0;0 0 1";
			//L(2,1) = InfDS.thetatild(16);
			//QQ = L*D*L.t();
			
			Rprintf("ss = %lf, InfDS.errtrol = %lf, InfDS.errtrol1 = %lf\n", ss, InfDS.errtrol, InfDS.errtrol1);
		}
		}

		//%%%%%%%%%%
		k = k+1;
		prev_stage = stage;
		if (stage == 2){ 
			k2 = k2 +1;
		}
    
	} //end of scoring iteration loop

	ttt = difftime(time(NULL), timer);
	if( convFlag == 1)
		Rprintf("\n\nThe estimation converged. There are totally %5d iterations in SAEM. Total running time is %5f seconds\n", k - 1, ttt);
	else
		Rprintf("\n\nThe estimation did not converge. There are totally %5d iterations in SAEM. Total running time is %5f seconds\n", k - 1, ttt);
	
	if (Npar > 0) {InfDS.par = InfDS.thetatild;	}
	

	Rprintf("(4) Wrap up estimation and write out results\n");

  //InfDS.Iytild.print("InfDS.Iytild");
	//InfDS.thetatild.print("InfDS.thetatild");

	/* setting outputs*/
	output.convFlag = convFlag;
	output.nIterStage1 = (k - 1) - (k2-1);
	output.nIterStage2 = k2 - 1;
	output.ss = ss;
	output.avebAccept = (output.avebAccept)/(k-1);
	
	output.Iytild = InfDS.Iytild;
	output.thetatild = InfDS.thetatild;
	
	/*
	output.Iytild = (double *)malloc((InfDS.par.n_elem * InfDS.par.n_elem + 1)* sizeof(double));
	output.thetatild = (double *)malloc((InfDS.par.n_elem + 1)* sizeof(double));
	for(j = 0; j < InfDS.par.n_elem; j++){
		output.thetatild[j] = InfDS.thetatild(j);
		for(i = 0; i < InfDS.par.n_elem;i++){
			output.Iytild[j*InfDS.par.n_elem + i] = InfDS.Iytild(i, j);
		}
	}
	*/


	return;
}
