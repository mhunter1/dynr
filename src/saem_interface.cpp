#include <stdlib.h>
#include<math.h>

#include <armadillo>
using namespace arma;
#include "structure_prototype.h"
#include "supplementary_function.h"
#include "converted_function.h"

#include "saem_estimation.h"
#define ERROR 0.0000001



extern "C" void saem_interface(int, int, int, int, int, int, int, int, int, int, int, int, double, double **, double **, double **, double **, double, double *, double **, double, double, int, int, int, double, double, double, double, double *, int, double *, double *, double  *,double *, double **, double **, double **, int, double **, int *, double *, double **, double **, double **, double **, double **, double **, double *, int, double **, double **);

void saem_interface(int seed, int freeIC, int Nsubj, int NxState, int Ny, int Nu, int Ntheta, int Nbeta, int totalT, int NLambda, int Nmu, int Nb, double delt, double **U1, double **b, double **H, double **Z, double maxT, double *allT, double **y0 , double lb, double ub, int MAXGIB, int MAXITER, int maxIterStage1, double gainpara, double gainparb, double gainpara1, double gainparb1, double *bAdaptParams, int Nbpar, double *mu, double *tspan, double *lower_bound, double *upper_bound, double **Lambda, double **dmudparMu, double **dmudparMu2, int num_time, double **Y, int *tobs, double *timeDiscrete, double **Sigmab, double **Sigmae, double **dSigmabdb, double **dSigmabdb2, double **dLambdparLamb, double **dLambdparLamb2, double *par_value, int KKO, double **dSigmaede, double **dSigmaede2){


	//printf("check point 0\n");	
	int i, j, Npar, tobs_pointer;
	C_INFDS InfDS;
	C_INFDS0 InfDS0;
	arma::mat upperb, lowerb, x1;

	InfDS.NxState = NxState;
	InfDS.Ny = Ny;
	InfDS.Ntheta = Ntheta;
	InfDS.Nb = Nb;
	InfDS.Nmu = Nmu;
	InfDS.NLambda = NLambda;
	InfDS.Nbpar = Nbpar;
	InfDS.Nu = Nu;
	InfDS.Nbeta = Nbeta;
	InfDS.Nbetax = Nbeta;
	InfDS.Nx = NxState; // in the saem process, Nx remains as NxState; it only becaomes NxState + Nbeta during initial value estimation
	InfDS.Nsubj = Nsubj;
	InfDS.totalT = totalT;
	InfDS.delt = delt;
	InfDS.maxT = maxT;
	InfDS.KKO = KKO;
	
	/*constant*/
	InfDS.N = 1;
	InfDS.omega = 61.685028;
	
	
	printf("Seed = %d\n", seed);
	arma_rng::set_seed(seed);
	//printf("check point 1 Nb %d Nsubj %d KKO %d\n", Nb, Nsubj, KKO);	

	Npar = Ntheta + NxState + Nmu + NLambda + Ny + Nbpar;
	//Npar = NxState + Nmu + NLambda + Ny + Nbpar;
	//Npar = Ntheta + NLambda + Ny + Nbpar;
	
	InfDS.U1.set_size(Nsubj, Nu);
	for(i = 0; i < Nsubj; i++){
		for(j = 0; j < Nu; j++){
			InfDS.U1(i, j) = U1[i][j];
		}
	}
	//InfDS.U1.print("InfDS.U1");

	InfDS.b.set_size(Nsubj, Nb);
	for(i = 0; i < Nsubj; i++){
		for(j = 0; j < Nb; j++){
			InfDS.b(i, j) = b[i][j];
			//InfDS.b(i, j) = .1;
			//temporarily set it as 0.1
		}
	}
	//InfDS.b.print("InfDS.b");


	InfDS.H.set_size(Nsubj* Ntheta, InfDS.Nbetax);
	for(i = 0; i < Nsubj* Ntheta; i++){
		for(j = 0; j < InfDS.Nbetax; j++){
			InfDS.H(i, j) = H[i][j];
		}
	}
	//InfDS.H.print("InfDS.H");

	InfDS.Z.set_size(Ntheta, Nb);
	for(i = 0; i < Ntheta; i++){
		for(j = 0; j < Nb; j++){
			InfDS.Z(i, j) = Z[i][j];
		}
	}
	//InfDS.Z.print("InfDS.Z");


	
	InfDS.allT.set_size(1, Nsubj);
	for(i = 0; i < Nsubj; i++){
		InfDS.allT(0,i) = allT[i];
	}
	//InfDS.allT.print("InfDS.allT");

	InfDS.Tfilter.set_size(Nsubj, totalT);
	for(i = 0; i < Nsubj; i++){
		for(j = 0; j < totalT; j++){
			InfDS.Tfilter(i, j) = 1;
		}
	}
	//InfDS.Tfilter.print("InfDS.Tfilter");

	InfDS.lens.set_size(1, Nsubj);
	for(i = 0; i < Nsubj; i++){
		InfDS.lens(0,i) = num_time;
	}
	//InfDS.lens.print("InfDS.lens");

	InfDS.tspan.set_size(1, num_time);
	for(i = 0; i < num_time; i++){
		InfDS.tspan(0,i) = tspan[i];
	}
	//InfDS.tspan.print("tspan");


	// start from here to be modeified
	InfDS.y0.set_size(Nsubj, NxState);
	for(i = 0; i < Nsubj; i++){
		for(j = 0; j < NxState; j++){
			InfDS.y0(i, j) = y0[i][j];
		}
	} 
	//InfDS.y0.zeros();
	//InfDS.y0.print("InfDS.y0");


	InfDS.timeDiscrete.set_size(Nsubj,1);
	tobs_pointer = 0;
	for(i = 0; i < Nsubj; i++){
		InfDS.timeDiscrete(i).set_size(totalT, 1);
		InfDS.timeDiscrete(i).zeros();
		for(j = 0; j < InfDS.allT(i); j++){
			InfDS.timeDiscrete(i)(j,0) = timeDiscrete[tobs_pointer];
			tobs_pointer++;
		}
	}
	//InfDS.timeDiscrete(0).t().print("timeDiscrete(0)");
	/*
	for(i = 0; i < 1; j++){
		printf("InfDS.timeDiscrete{%d} = [", i);
		for(j = 0; j < InfDS.allT(i); j++){
			if(j > 0)
				printf(";");
			printf("%lf",InfDS.timeDiscrete(i)(j,0));
		}
		printf("]");
	
	}
	*/



	InfDS.tobs.set_size(Nsubj,1);
	tobs_pointer = 0;	
	for(i = 0; i < Nsubj; i++){
		InfDS.tobs(i).set_size(InfDS.allT(i), 1);
		//InfDS.tobs(i).col(0) = span_vec(1, InfDS.allT(i), 1);
		for(j = 0; j < InfDS.allT(i); j++){
			InfDS.tobs(i)(j,0) = tobs[tobs_pointer];
			tobs_pointer++;
		}	
	}
	//InfDS.tobs(49).print("InfDS.tobs(199)");
	
	/*
	for(i = 0; i < 1; j++){
		printf("InfDS.tobs{%d} = [", i);
		for(j = 0; j < InfDS.allT(i); j++){
			if(j > 0)
				printf(";");
			printf("%d",InfDS.tobs(i)(j,0));
		}
		printf("]");
	
	}
	*/
	
	InfDS.Y.set_size(Nsubj,1);
	for(i = 0; i < Nsubj; i++){
	    InfDS.Y(i).set_size(Ny, InfDS.allT(i));
	}
	for(int y = 0; y < Ny; y++){
		tobs_pointer = 0;
		for(i = 0; i < Nsubj; i++){
			for(j = 0; j < InfDS.allT(i); j++){
				InfDS.Y(i)(y,j) = Y[y][tobs_pointer];
				tobs_pointer++;
			}
		}
	}
	
	/*
	for(i = 0; i < 1; j++){
		printf("InfDS.Y{%d} = [", i);
		for(int y = 0; y < Ny; y++){
			for(j = 0; j < InfDS.allT(i); j++){
				if(j > 0)
					printf(",");
				printf("%d",InfDS.Y(i)(y,j));
			}
			if(y != Ny-1)
				printf(";\n");
		}
		printf("];\n");
	
	}
	*/
	
	

	InfDS.mu.set_size(Nmu,1);
	for(i = 0; i < Nmu; i++){
		InfDS.mu(i) = mu[i];
	}
	//InfDS.mu.print("InfDS.mu");
	
	
	InfDS.Lambda.set_size(Ny,NxState);
	for(i = 0; i < Ny; i++){
		for(j = 0; j < NxState; j++){
			InfDS.Lambda(i,j) = Lambda[i][j];
		}
	}
	//InfDS.Lambda.print("InfDS.Lambda");
	
	
	

	InfDS.dXtildthetafAll.set_size(Nsubj,1);
	for(i = 0; i < Nsubj; i++){
		InfDS.dXtildthetafAll(i).set_size(Ntheta, NxState*InfDS.allT(0,i));
		InfDS.dXtildthetafAll(i).zeros();
	}
	//InfDS.dXtildthetafAll(0).print("InfDS.dXtildthetafAll(0)");;

	InfDS.dXtildthetafAll2.set_size(Nsubj,1);
	for(i = 0; i < Nsubj; i++){
		InfDS.dXtildthetafAll2(i).set_size(Ntheta*NxState*InfDS.allT(0,i), Ntheta);
		InfDS.dXtildthetafAll2(i).zeros();
	}
	//InfDS.dXtildthetafAll2(0).print("InfDS.dXtildthetafAll2(0)");
	
	
	//the following part needs to be generalized 201-211
	InfDS.Sigmab = zeros(Nb, Nb);
	for(i = 0; i < Nb; i++){
		for(j = 0; j < Nb; j++){
			InfDS.Sigmab(i,j) = Sigmab[i][j];
		}
	}
	//InfDS.Sigmab.print("InfDS.Sigmab");
	
	InfDS.dSigmabdb = zeros(InfDS.Nbpar, Nb*Nb);
	for(i = 0; i < InfDS.Nbpar; i++){
		for(j = 0; j < Nb*Nb; j++){
			InfDS.dSigmabdb(i,j) = dSigmabdb[i][j];
		}
	}
	//InfDS.dSigmabdb.print("InfDS.dSigmabdb");
	
	InfDS.dSigmabdb2 = zeros(InfDS.Nbpar*Nb*Nb, InfDS.Nbpar);
	for(i = 0; i < InfDS.Nbpar*Nb*Nb; i++){
		for(j = 0; j < InfDS.Nbpar; j++){
			InfDS.dSigmabdb2(i,j) = dSigmabdb2[i][j];
		}
	}
	//InfDS.dSigmabdb2.print("InfDS.dSigmabdb2");
	
	InfDS.dLambdparLamb = zeros(NLambda, Ny * InfDS.Nx);
	for(i = 0; i < NLambda; i++){
		for(j = 0; j < Ny * InfDS.Nx; j++){
			InfDS.dLambdparLamb(i,j) = dLambdparLamb[i][j];
		}
	}
	//InfDS.dLambdparLamb.print("InfDS.dLambdparLamb");
	
	
	InfDS.dLambdparLamb2 = zeros(NLambda* Ny * InfDS.Nx, NLambda);
	for(i = 0; i < NLambda* Ny * InfDS.Nx; i++){
		for(j = 0; j < NLambda; j++){
			InfDS.dLambdparLamb2(i,j) = dLambdparLamb2[i][j];
		}
	}
	//InfDS.dLambdparLamb2.print("InfDS.dLambdparLamb2");
	
	InfDS.Sigmae = zeros(Ny, Ny);
	for(i = 0; i < Ny; i++){
		for(j = 0; j < Ny; j++){
			InfDS.Sigmae(i, j) = Sigmae[i][j];
		}
	}
	//InfDS.Sigmae.print("InfDS.Sigmae");
	
	
	InfDS.dSigmaede = zeros(Ny, Ny*Ny);
	for(i = 0; i < Ny; i++){
		for(j = 0; j < Ny*Ny; j++){
			InfDS.dSigmaede(i, j) = dSigmaede[i][j];
		}
	}
	//InfDS.dSigmaede.print("InfDS.dSigmaede");
	
	InfDS.dSigmaede2 = zeros(Ny*Ny*Ny, Ny);	
	for(i = 0; i < Ny*Ny*Ny; i++){
		for(j = 0; j < Ny; j++){
			InfDS.dSigmaede2(i, j) = dSigmaede2[i][j];
		}
	}
	//InfDS.dSigmaede2.print("InfDS.dSigmaede2");
	
	
	

	//to be obtained from prep.measurements
	InfDS.dmudparMu.set_size(Ny, Ny);
	for(i = 0; i < Ny; i++){
		for(j = 0; j < Ny; j++){
			InfDS.dmudparMu(i,j) = dmudparMu[i][j];
		}
	}
	//InfDS.dmudparMu.print("dmudparMu");
	
	InfDS.dmudparMu2.set_size(Ny * Ny, Ny);
	for(i = 0; i < Ny*Ny; i++){
		for(j = 0; j < Ny; j++){
			InfDS.dmudparMu2(i,j) = dmudparMu2[i][j];
		}
	}
	//InfDS.dmudparMu2.print("dmudparMu2");
	
	

	InfDS.sy = zeros(Npar, 1);
	InfDS.EI = zeros(Npar, 1);
	InfDS.ES = zeros(Npar, 1);
	InfDS.Iy = zeros(Npar, 1);

	InfDS.sytild= zeros(Npar, 1);
	InfDS.EItild = zeros(Npar, Npar);
	InfDS.EStild = zeros(Npar, Npar);
	InfDS.Iytild = zeros(Npar, Npar);
	

	InfDS.lowBound.set_size(Npar, 1);
	InfDS.upBound.set_size(Npar, 1);
	InfDS.par.set_size(Npar, 1);
	for(i = 0; i < Npar; i++){
		InfDS.lowBound(i) = lower_bound[i];
		InfDS.upBound(i) =upper_bound[i];
		InfDS.par(i) = par_value[i];
	}
	InfDS.par.print("InfDS.par");
	InfDS.lowBound.print("InfDS.lowBound");
	InfDS.upBound.print("InfDS.upBound");

	
	InfDS.thetatild = zeros(Npar, 1);
	InfDS.Xtild.set_size(InfDS.Nx, InfDS.Nsubj, InfDS.totalT);
	//InfDS.Xtild.zeros();

	//SAEM Control Parameters
	InfDS.MAXGIB = MAXGIB;
	InfDS.MAXITER = MAXITER;
	InfDS.maxIterStage1 = maxIterStage1;
	InfDS.gainpara = gainpara;
	InfDS.gainparb = gainparb;
	InfDS.gainpara1 = gainpara1;
	InfDS.gainparb1 = gainparb1;	
	InfDS.bAdaptParams.set_size(1,3);
	for(i = 0; i < 3; i++){
		InfDS.bAdaptParams(i) = bAdaptParams[i];
	}


	
	lowerb.set_size(1, Nb);
	lowerb.fill(lb);
	upperb.set_size(1, Nb);
	upperb.fill(ub);

	InfDS.Sigmab.set_size(Nb, Nb);
	for(i = 0; i < Nb; i++){
		for(j = 0; j < Nb; j++){
			InfDS.Sigmab(i,j) = Sigmab[i][j];
		}

	}
	//InfDS.Sigmab.print("sigmab");
	
	
	

	char filenamePar[64] = "./Results/TrueInitparG1.txt";
	char filenameSE[64] = "./Results/TrueInitSEG1.txt";
	char filenameconv[64] = "./Results/TrueInitconvG1.txt";
	char filenamebhat[64] = "./Results/TrueInitbhatG1.txt";
	char filenamebhat2[64] = "./Results/TrueInitbhat2G1.txt";
	int kk = 1;
	int trueInit = 1;
	int batch = 1;
	//int seed = 1;
	printf("calling...\n");
	saem_estimation(InfDS, InfDS0, upperb, lowerb, x1, filenamePar, filenameSE, filenameconv, filenamebhat, filenamebhat2, kk, trueInit, batch, seed, freeIC);
	return;
}
