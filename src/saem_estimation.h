#include <time.h>
//#include <stdlib.h>
//#include<math.h>

//#include <armadillo>

//using namespace arma;

//#include "structure_prototype.h"
//#include "supplementary_function.h"
//#include "converted_function.h"

// Step 3 in the MainUseThins.m
void saem_estimation(C_INFDS &InfDS, C_INFDS0 &InfDS0, arma::mat upperb, arma::mat lowerb, arma::mat x1, char *filenamePar, char *filenameSE, char *filenameconv, char *filenamebhat, char *filenamebhat2, int kk, int trueInit, int batch, int seed, int freeIC){
	//printf("in MainUseThis\n");
	arma::mat sgnTH, meanb, L, QQ, D, mscore2, OMEGAb, infoMat, minfoMat, tpOld, score, Covscore;
	arma::vec mscore;
	int k, stage, gmm, MAXGIB, setScaleb, noIncrease, isPar, yesMean, switchFlag, useMultN, GIB, STARTGIB, stop, isBlock1Only, redFlag, convFlag;
	double bAccept, ss, ttt, ssmin;
	int prev_stage;
	time_t timer;
	//--for writing output files--
	//int i, j, fitInit;
	//FILE *p_filenamePar, *p_filenameSE, *p_filenameconv, *p_filenamebhat, *p_filenamebhat2;
	//----
	
	//freeIC = 1;
	timer = time(NULL);
	
	//rma_rng::set_seed_random(); 
	arma_rng::set_seed(seed);

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
	InfDS.alp=1;
	MAXGIB = InfDS.MAXGIB;
	
	//printf("checkpoint M32\n");	
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
	
	//InfDS.y0 = arma::zeros<arma::mat> (InfDS.Nsubj, 2);
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

	printf("InfDS.par.n_elem = %d \n", InfDS.par.n_elem);	

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
	
	//printf("checkpoint M63\n");	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%self change, with caution and ask
	//InfDS.lowBound(span(8,9), span::all).fill(10e-8);
	//InfDS.lowBound(span(10, 15), span::all).fill(log(10e-8));
	//InfDS.lowBound(16) = -2; //After transformation, covariance between -5 and 5

	//InfDS.upBound(span(8,9), span::all).fill(10);
	//InfDS.upBound(span(10, 15), span::all).fill(2);
	
	//printf("checkpoint M73\n");	

	k = 1; 
	stage = 1; 
	gmm = 1;
	//MAXGIB=1;
	InfDS.errtrol1 = 1.5; //Stage 1 error tolerance
	InfDS.errtrol = 1.5; //Stage 2 error tolerance	


	
	setScaleb = 0;
	//InfDS.bAdaptParams = ".5, 2.5, .5";
	ssmin = 100; 
	noIncrease = 0;
	//freeIC = 1;
	InfDS.scaleb = 1; //Used in drawbGenera6_opt3.m to determine whether to apply scaling constant on drawb.
	//InfDS.KKO = 20; //Used in SAEM. Only starts to evaluate whether to transition to stage 2 after KKO iterations.
	stop = 0;
    //printf("check point 2 MAXITER %d freeIC %d\n", InfDS.MAXITER, freeIC);
	
	//InfDS.par.print("InfDS.par");
	//printf("checkpoint M92 entering the k loop\n");	
	while (k <= InfDS.MAXITER && stop == 0){

		//disp 'iteration';
		//k
		printf("k = %d\n",k);

		isPar = 0;
		
		//setPartsFreeICwb shoudl NOT be called
		//setParsFreeICwb(InfDS); //qqqq	
		//printf("checkpoint M101 setParsFreeICwb\n");	
		
		if (stage==2 && switchFlag==0){
			yesMean= 1;
			switchFlag = 1; 
			meanb = arma::zeros<arma::mat>(InfDS.Nsubj,InfDS.Nb); //Keeps b estimates averaged across Gibbs runs
		}

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

    
		if (stage==1 && bAccept < .1){ //Not looking too good. Pump up no. of chains
			//%InfDS.bAdaptParams(1) = .5;% Lu modified, 04-12-13 .1;
            //%InfDS.bAdaptParams(2) = 2.5;% Lu modified, 04-12-13 10.5;
            //%InfDS.bAdaptParams(3) = .5;
            useMultN = 1; // Lu modified, 04-12-13,5;
		}
		else if(stage==1 && k == 3){
            useMultN = 1; 
			//MAXGIB = 30;
		}
		else if (stage==1 && k <= 4 && k != 3){
            useMultN = 0; 
			//MAXGIB = 1;
		}
		else {
			useMultN = 0; 
			//MAXGIB = 10;
		}

		mscore = arma::zeros<arma::mat>(InfDS.par.n_elem, 1);
		mscore2 = arma::zeros<arma::mat>(InfDS.par.n_elem, InfDS.par.n_elem);
		minfoMat = arma::zeros<arma::mat>(InfDS.par.n_elem, InfDS.par.n_elem);  
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//%% Gibbs sampler
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
		//Since beta changes, get Xtild and InfDS.InfDS.dxstardthetafAll and InfDS.dxstardthetafAll2
	
		
		//printf("checkpoint M142\n");	
		//InfDS.par.print("InfDS.par");
		isPar = (InfDS.Nx == InfDS.NxState) ? 0 : 1;
		InfDS = getXtildIC3(isPar, 1 ,freeIC, InfDS); //%Get updated Xtilde
		//printf("checkpoint M145 getXtildIC3\n");	
		InfDS.Xtild(span(0,1), span(0,1), span(0,1)).print("InfDS.Xtild");

		PropSigb(InfDS);  //covariance of proposal distribution of b
		//printf("checkpoint M147PropSigb\n");	
		InfDS.par.print("InfDS.par");

		tpOld = ekfContinuous10(InfDS.Nsubj, InfDS.N, InfDS.Ny, InfDS.Nx, InfDS.Nb, InfDS.NxState, InfDS.Lambda, InfDS.totalT, InfDS.Sigmae, InfDS.Sigmab, InfDS.mu, InfDS.b, InfDS.allT, InfDS.Xtild, InfDS.Y); //%get density of full conditional distribution of b 
		
		InfDS.bacc = arma::zeros<arma::mat>(InfDS.Nsubj,1);	
	
		//printf("checkpoint M153 enter GIB loop\n");	
		//MAXGIB = 5;
		printf("[DEBUG] MAXGIB = %d \n", MAXGIB);

		for(GIB = 1; GIB <= MAXGIB; GIB++){
			if (stage == 2){
				yesMean = 1;
				STARTGIB = STARTGIB+1;
			}

			if (k >= 4){        
				//printf("checkpoint enter drowbGeneral6_opt3\n");	
				drawbGeneral6_opt3(isPar, InfDS, meanb, yesMean, upperb, lowerb, useMultN, tpOld, freeIC, isBlock1Only, setScaleb, bAccept);
				//printf("checkpoint leave drowbGeneral6_opt3\n");	
	                        //InfDS.par.print("InfDS.par");
			}
        
			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			//printf("checkpoint enter getScoreInfoY_tobs_opt\n");
			getScoreInfoY_tobs_opt(InfDS, stage, k, freeIC, score, infoMat);			
			//printf("checkpoint leave getScoreInfoY_tobs_opt\n");
			printf("GIB=%d\n", GIB);
			score.print("score");
        
			//InfDS.par.print("InfDS.par");
			
			mscore = mscore + (1.0/MAXGIB)*score;
			mscore2 = mscore2 + (1.0/MAXGIB)*(score*score.t());
			minfoMat = minfoMat + (1.0/MAXGIB)*infoMat;
        
		} //end of Gibbs sampler loop
  
			
		//printf("checkpoint M181 leave GIB loop\n");
		//InfDS.par.print("InfDS.par");
	
		Covscore = mscore2 - mscore*mscore.t();
	
	
		//printf("checkpoint enter saem\n");
		saem(InfDS, gmm, stage, redFlag, convFlag, noIncrease, stop, ssmin, ss, sgnTH, mscore, mscore2, minfoMat, Covscore);
		//printf("checkpoint leave saem\n");
    
		//printf("checkpoint M145 saem\n");	
	
    
		printf("\nStage = %5d, iteration = %5d\n",stage,k);
		printf("\nCurrent b acceptance rate = %6f\n",bAccept);
		//printf("\nRange of InfDS0.trueb = %6f, %6f",min(InfDS0.trueb),max(InfDS0.trueb));
		//printf("\nRange of bhat = %6f, %6f\n",(double)min(InfDS.b),(double)max(InfDS.b));
		//corr(InfDS.b(:,1:size(InfDS0.trueb,2)),InfDS0.trueb)

		//temporarily printing out messages
		if(1 || prev_stage != stage){
			InfDS.par(span(0,9), span::all).print("InfDS.par(1:10)");
			//exp(InfDS.par(span(10,12), span::all)).t().print("InfDS.par(11:13)");
			///InfDS.Sigmab.print("Sigmab");
			printf("Averaging:\n");
			InfDS.thetatild(span(0,9),span::all).print("InfDS.thetatild(1:10)");
			//exp(InfDS.thetatild(span(10,12), span::all)).t().print("exp(InfDS.thetatild)");
			//D = diagmat(exp(InfDS.thetatild(span(12,14), span::all)));
			L = "1 0 0;	0 1 0;0 0 1";
			//L(2,1) = InfDS.thetatild(16);
			//QQ = L*D*L.t();
		}

		//%%%%%%%%%%
		k = k+1;
		prev_stage = stage;

		/*remove later*/
		//break;
    
	} //end of scoring iteration loop

	ttt = difftime(time(NULL), timer);
	if( convFlag == 1)
		printf("\nThe estimation converged. There are totally %5d iterations. Total running time is %5f seconds\n", k, ttt);
	else
		printf("\nThe estimation did not converge. There are totally %5d iterations. Total running time is %5f seconds\n", k, ttt);
	

	meanb = meanb/STARTGIB;
	InfDS.par = InfDS.thetatild;	
	

	printf("(4) Wrap up estimation and write out results\n");


	//arma::mat dgdpar;

	//If using transformation functions
	//dgdpar = eye(InfDS.par.n_rows, InfDS.par.n_rows);
	//dgdpar(span(10,12), span(10,12)) = diagmat(exp(InfDS.par(span(10,12),0)));

/*
	printf("(41) Wrap up estimation and write out results\n");

	//Columns -- par, rows --transformation function
	dgdpar(13, 13) = exp(InfDS.par(13));
	dgdpar(14, 14) = exp(InfDS.par(14));

	dgdpar(15,14) = InfDS.par(16)*exp(InfDS.par(14));//diff(f1,par2)
	dgdpar(15,16) = exp(InfDS.par(14));//%diff(f1,par4)
	dgdpar(16,14) = InfDS.par(16)*InfDS.par(16)*exp(InfDS.par(14));//diff(f2,par2)
	dgdpar(16,15) = exp(InfDS.par(15));//%diff(f2,par3)
	dgdpar(16,16) = 2*InfDS.par(16)*exp(InfDS.par(14)); //diff(f2,par4)
	
	printf("(42) Wrap up estimation and write out results\n");
	
	arma::mat SE;

	//SE = sqrt(diagvec(dgdpar/InfDS.Iytild*dgdpar.t()));
        SE = sqrt(diagvec(dgdpar/InfDS.Iytild*dgdpar.t()));
	InfDS.par(10) = exp(InfDS.par(10));
	InfDS.par(11) = exp(InfDS.par(11));
	InfDS.par(12) = exp(InfDS.par(12));
	//to be converted
	//InfDS.par(span(13,16)) = [QQ(0,0), QQ(1,1), QQ(2,2), QQ(1,2)];
	InfDS.par(13) = QQ(0,0);
	InfDS.par(14) = QQ(1,1);
	InfDS.par(15) = QQ(2,2);
	InfDS.par(16) = QQ(1,2);
	//InfDS.lastb = InfDS.b; 
	//InfDS.meanb = meanb;
	//InfDS.convFlag = convFlag;
	//InfDS.ss = ss;

	printf("(43) Wrap up estimation and write out results\n");
*/
/*
	fitInit = 1; //Fit models with freely estimated IC.
	//dlmwrite(filenamePar,[trueInit fitInit batch kk reshape(InfDS.par,1,length(InfDS.par))],'-append');
	p_filenamePar = fopen(filenamePar, "a");
	fprintf(p_filenamePar, "%d,%d,%d,%d,%d", seed, trueInit, fitInit, batch, kk);
	for (i = 0; i < InfDS.par.n_elem; i++)
		fprintf(p_filenamePar, ",%lf", InfDS.par(i));
	fprintf(p_filenamePar, "\n");
	
	//dlmwrite(filenameSE,[trueInit fitInit batch kk reshape(InfDS.SE,1,length(InfDS.par))],'-append');
	p_filenameSE = fopen(filenameSE, "a");
	fprintf(p_filenameSE, "%d,%d,%d,%d,%d", seed, trueInit, fitInit, batch, kk);
	for (i = 0; i < InfDS.par.n_elem; i++)
		fprintf(p_filenameSE, ",%lf", SE(i));
	fprintf(p_filenameSE, "\n");
	
	//dlmwrite(filenameconv,[trueInit fitInit batch kk InfDS.convFlag InfDS.ss],'-append');
	p_filenameconv = fopen(filenameconv, "a");
	fprintf(p_filenameconv, "%d,%d,%d,%d,%d,%d,%lf\n", seed, trueInit, fitInit, batch, kk, convFlag, ss);
	
	
	//dlmwrite(filenamebhat,[trueInit fitInit batch kk reshape(InfDS.lastb,1,InfDS.Nsubj*InfDS.Nb)],'-append'); %Reshape by column
	p_filenamebhat = fopen(filenamebhat, "a");
	fprintf(p_filenamebhat, "%d,%d,%d,%d,%d", seed, trueInit, fitInit, batch, kk);
	for (j = 0; j < InfDS.b.n_cols; j++){
		for (i = 0; i < InfDS.b.n_rows; i++)
			fprintf(p_filenamebhat, ",%lf", InfDS.b(i,j));
	}
	fprintf(p_filenamebhat, "\n");
	
	//dlmwrite(filenamebhat2,[trueInit fitInit batch kk reshape(InfDS.meanb,1,InfDS.Nsubj*InfDS.Nb)],'-append');%Reshape by column
	p_filenamebhat2 = fopen(filenamebhat2, "a");
	fprintf(p_filenamebhat2, "%d,%d,%d,%d,%d",seed, trueInit, fitInit, batch, kk);
	for (j = 0; j < meanb.n_cols; j++){
		for (i = 0; i < meanb.n_rows; i++)
			fprintf(p_filenamebhat2, ",%lf", meanb(i,j));
	}
	fprintf(p_filenamebhat2, "\n");
*/



	return;
}
