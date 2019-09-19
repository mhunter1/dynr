#include <stdlib.h>
#include<math.h>

#include <armadillo>
using namespace arma;
#include "structure_prototype.h"
#include "supplementary_function.h"
#include "converted_function.h"

#include "MainUseThis.h"
#define ERROR 0.0000001

extern "C" void interface(int, int, int, int, int, int, int, int, int, int, int, double, double **, double **, double **, double **, double, double *, double **, double, double, int, int, int, double, double, double, double, double *, int, double *, double *, double  *,double *, double *, double **, double **);

void interface(int seed, int Nsubj, int NxState, int Ny, int Nu, int Ntheta, int Nbeta, int totalT, int NLambda, int Nmu, int Nb, double delt, double **U1, double **b, double **H, double **Z, double maxT, double *allT, double **y0 , double lb, double ub, int MAXGIB, int MAXITER, int maxIterStage1, double gainpara, double gainparb, double gainpara1, double gainparb1, double *bAdaptParams, int Nbpar, double *mu, double *tspan, double *lower_bound, double *upper_bound, double *Lambda, double **dmudparMu, double **dmudparMu2){
	int i, j, Npar;
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
	InfDS.Nx = NxState + Nbeta;
	InfDS.Nsubj = Nsubj;
	InfDS.totalT = totalT;
	InfDS.delt = delt;
	
	/*constant*/
	InfDS.N = 1;
	InfDS.omega = 61.685028;
	
	
	
	Npar = Ntheta + NxState + Nmu + NLambda + Ny + Nbpar;

	InfDS.U1.set_size(Nsubj, Nu);
	for(i = 0; i < Nsubj; i++){
		for(j = 0; j < Nu; j++){
			InfDS.U1(i, j) = U1[i][j];
		}
	}
	//InfDS.U1.print("InfDS.U1");

	InfDS.b.set_size(Nsubj, Nb);
	for(i = 0; i < Nsubj; i++){
		for(j = 0; j < Nu; j++){
			InfDS.b(i, j) = b[i][j];
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


	
	/*Need to be generalized*/
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

	/*Need to be generalized*/
	InfDS.lens.set_size(1, Nsubj);
	for(i = 0; i < Nsubj; i++){
		InfDS.lens(0,i) = 431;
	}

	InfDS.tspan.set_size(1, totalT);
	for(i = 0; i < Nsubj; i++){
		InfDS.tspan(0,i) = tspan[i];
	}


	// start from here
	InfDS.y0.set_size(Nsubj, InfDS.Nx);
	InfDS.y0.zeros();
	for(i = 0; i < Nsubj; i++){
		for(j = 0; j < NxState; j++){
			InfDS.y0(i, j) = y0[i][j];
		}
	}
	InfDS.y0.print("InfDS.y0");


	InfDS.timeDiscrete.set_size(Nsubj,1);
	for(i = 0; i < Nsubj; i++){
		InfDS.timeDiscrete(i).set_size(totalT, 1);
		InfDS.timeDiscrete(i).zeros();
	}

	InfDS.tobs.set_size(Nsubj,1);
	for(i = 0; i < Nsubj; i++){
		InfDS.tobs(i).set_size(totalT, 1);
		InfDS.tobs(i).col(0) = span_vec(1, totalT, 1);
	}
	InfDS.tobs(0).print("InfDS.tobs(0)");

	InfDS.mu.set_size(Nmu,1);
	for(i = 0; i < Nmu; i++){
		InfDS.mu(i) = mu[i];
	}
	InfDS.mu.print("InfDS.tobs(0)");
	
	
	InfDS.Lambda.set_size(Ny,NxState);
	for(i = 0; i < Ny; i++){
		for(j = 0; j < NxState; j++){
			InfDS.Lambda(i,j) = Lambda[i][j];
		}
	}
	InfDS.Lambda.print("InfDS.Lambda");
	
	//todo obeservations
	InfDS.Y.set_size(200,1);
	InfDS.Y(0).set_size(3, 300);

	InfDS.dXtildthetafAll.set_size(Nsubj,1);
	for(i = 0; i < Nsubj; i++){
		InfDS.dXtildthetafAll(i).set_size(Ntheta, NxState*InfDS.allT(0,i));
		InfDS.dXtildthetafAll(i).zeros();
	}
	InfDS.dXtildthetafAll(0).print("InfDS.dXtildthetafAll(0)");;

	InfDS.dXtildthetafAll2.set_size(Nsubj,1);
	for(i = 0; i < Nsubj; i++){
		InfDS.dXtildthetafAll2(i).set_size(Ntheta, NxState*InfDS.allT(0,i));
		InfDS.dXtildthetafAll2(i).zeros();
	}
	InfDS.dXtildthetafAll2(0).print("InfDS.dXtildthetafAll2(0)");
	
	
	InfDS.Sigmab = zeros(Nb, Nb);
	InfDS.dSigmabdb = zeros(InfDS.Nbpar, Nb*Nb);
	InfDS.dSigmabdb2 = zeros(InfDS.Nbpar*Nb*Nb, InfDS.Nbpar);
	
	InfDS.dLambdparLamb = zeros(NLambda, Ny * InfDS.Nx);
	
	InfDS.dSigmaede2 = zeros(Ny*Ny*Ny, Ny);	
	InfDS.dSigmaede = zeros(Ny, Ny*Ny);
	

	//to be obtained from prep.measurements
	InfDS.dmudparMu.set_size(Ny, Ny);
	for(i = 0; i < Ny; i++){
		for(j = 0; j < Ny; j++){
			InfDS.dmudparMu(i,j) = dmudparMu[i][j];
		}
	}
	InfDS.dmudparMu.print("dmudparMu");
	
	InfDS.dmudparMu2.set_size(Ny * Ny, Ny);
	for(i = 0; i < Ny*Ny; i++){
		for(j = 0; j < Ny; j++){
			InfDS.dmudparMu2(i,j) = dmudparMu2[i][j];
		}
	}
	InfDS.dmudparMu2.print("dmudparMu2");
	
	
	InfDS.par = zeros(Npar, 1);
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
	for(i = 0; i < Npar; i++){
		InfDS.lowBound(i) = lower_bound[i];
		InfDS.upBound(i) =upper_bound[i]
	}
	

	
	InfDS.thetatild = zeros(Npar, 1);
	InfDS.Xtild.set_size(InfDS.Nx, InfDS.Nsubj, InfDS.totalT);
	InfDS.Xtild.zeros();

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
	

	char filenamePar[64] = "./Results/TrueInitparG1.txt";
	char filenameSE[64] = "./Results/TrueInitSEG1.txt";
	char filenameconv[64] = "./Results/TrueInitconvG1.txt";
	char filenamebhat[64] = "./Results/TrueInitbhatG1.txt";
	char filenamebhat2[64] = "./Results/TrueInitbhat2G1.txt";
	int kk = 1;
	int trueInit = 1;
	int batch = 1;

	MainUseThis(InfDS, InfDS0, upperb, lowerb, x1, filenamePar, filenameSE, filenameconv, filenamebhat, filenamebhat2, kk, trueInit, batch, seed);
	return;
}