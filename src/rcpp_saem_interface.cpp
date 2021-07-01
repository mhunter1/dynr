// [[Rcpp::depends(RcppArmadillo)]]
// The above lined is required for RcppArmadillo

// required headers
#include <RcppArmadillo.h>
using namespace std;
using namespace arma;
using namespace Rcpp;



#include "structure_prototype.h"
#include "supplementary_function.h"
#include "converted_function.h"

#include "saem_estimation.h"






// [[Rcpp::export]]
Rcpp::List rcpp_saem_interface(Rcpp::List model_sexp, Rcpp::List data_sexp, bool weight_flag, bool debug_flag, bool optimization_flag, bool hessian_flag, bool verbose){
	struct C_OUTPUT output;
	C_INFDS InfDS;
	C_INFDS0 InfDS0;

	
		
	/*constant*/
	InfDS.N = 1;
	InfDS.omega = 61.685028;
	
	/*scalar inputs*/
	//int
	int seed = as<int>(model_sexp["seed"]);
	int num_time = as<int>(model_sexp["num_time"]);
	int isfreeIC = as<int>(model_sexp["freeIC"]);
	InfDS.NxState = as<int>(model_sexp["dim_latent_var"]);
	InfDS.Nx = InfDS.NxState; // in the saem process, Nx remains as NxState; it only becaomes NxState + Nbeta during initial value estimation
	InfDS.Ny = as<int>(model_sexp["dim_obs_var"]);
	InfDS.Ntheta = as<int>(model_sexp["num_theta"]);
	InfDS.Nb = as<int>(model_sexp["num_random"]);
	InfDS.Nmu = as<int>(model_sexp["num_mu"]);
	InfDS.NLambda = as<int>(model_sexp["num_lambda"]);
	InfDS.Nbpar = as<int>(model_sexp["num_bpar"]);
	InfDS.Nu = as<int>(model_sexp["dim_co_variate"]);
	InfDS.Nbeta = as<int>(model_sexp["num_beta"]);
	InfDS.Nbetax = InfDS.Nbeta;
	InfDS.Nsubj = as<int>(model_sexp["num_sbj"]);
	InfDS.totalT = as<int>(model_sexp["total_t"]);
	InfDS.KKO = as<int>(model_sexp["KKO"]);
	InfDS.MAXGIB = as<int>(model_sexp["MAXGIB"]);
	InfDS.MAXITER = as<int>(model_sexp["MAXITER"]);
	InfDS.maxIterStage1 = as<int>(model_sexp["maxIterStage1"]);
	InfDS.KKO = as<int>(model_sexp["KKO"]);
	InfDS.scaleb = as<double>(model_sexp["scaleb"]);
	//double
	InfDS.delt = as<double>(model_sexp["delt"]);
	InfDS.maxT = as<double>(model_sexp["max_t"]);
	InfDS.gainpara = as<double>(model_sexp["gainpara"]);
	InfDS.gainparb = as<double>(model_sexp["gainparb"]);
	InfDS.gainpara1 = as<double>(model_sexp["gainpara1"]);
	InfDS.gainparb1 = as<double>(model_sexp["gainparb1"]);
	
	double lb = as<double>(model_sexp["random.lb"]);
	mat lowerb;
	lowerb.set_size(1, InfDS.Nb);
	lowerb.fill(lb);
	
	
	double ub = as<double>(model_sexp["random.ub"]);
	mat upperb;
	upperb.set_size(1, InfDS.Nb);
	upperb.fill(ub);
	
	
	
	
	//print out the values
	Rprintf("flags %d %d %d %d %d\n", int(weight_flag), int(debug_flag), int(optimization_flag), int(hessian_flag), int(verbose));
	Rprintf("scaleb = %lf\n", InfDS.scaleb);
	Rprintf("seed = %d\n", seed);
	Rprintf("num_time = %d\n", num_time);
	Rprintf("isfreeIC = %d\n", isfreeIC);
	Rprintf("InfDS.NxState = %d\n", InfDS.NxState);
	Rprintf("InfDS.Ny = %d\n", InfDS.Ny);
	Rprintf("InfDS.Nb = %d\n", InfDS.Nb);
	Rprintf("InfDS.Nmu = %d\n", InfDS.Nmu);
	Rprintf("InfDS.NLambda = %d\n", InfDS.NLambda);
	Rprintf("InfDS.Nbpar = %d\n", InfDS.Nbpar);
	Rprintf("InfDS.Nu = %d\n", InfDS.Nu);
	Rprintf("InfDS.Nbeta = %d\n", InfDS.Nbeta);
	Rprintf("InfDS.Nsubj = %d\n", InfDS.Nsubj);
	Rprintf("InfDS.NxState = %d\n", InfDS.NxState);
	Rprintf("InfDS.totalT = %d\n", InfDS.totalT);
	Rprintf("InfDS.KKO = %d\n", InfDS.KKO);
	Rprintf("InfDS.delt = %lf\n", InfDS.delt);
	Rprintf("InfDS.maxT = %lf\n", InfDS.maxT);
	upperb.print("upperb");
	lowerb.print("lowerb");
	Rprintf("InfDS.MAXGIB = %d\n", InfDS.MAXGIB);
	Rprintf("InfDS.MAXITER = %d\n", InfDS.MAXITER);
	Rprintf("InfDS.maxIterStage1 = %d\n", InfDS.maxIterStage1);
	
	Rprintf("[SAEM Parameters] %lf %lf %lf %lf\n", InfDS.gainpara, InfDS.gainparb, InfDS.gainpara1, InfDS.gainparb1);
	
	

	if(InfDS.Nb > 0){

		InfDS.Sigmab = as<mat>(model_sexp["sigmab"]);
		InfDS.Sigmab.print("InfDS.Sigmab");
	
	
		InfDS.dSigmabdb = as<mat>(model_sexp["dSigmabdb"]);
		InfDS.dSigmabdb.print("InfDS.dSigmabdb");
		
		InfDS.dSigmabdb2 = as<mat>(model_sexp["dSigmabdb2"]);
		InfDS.dSigmabdb2.print("InfDS.dSigmabdb2");
	}
	
	if(InfDS.Ny > 0){
		InfDS.Sigmae = as<mat>(model_sexp["sigmae"]);
		InfDS.Sigmae.print("InfDS.Sigmae");
		
		InfDS.dSigmaede = as<mat>(model_sexp["dSigmaede"]);
		InfDS.dSigmaede.print("InfDS.dSigmaede");
		
		InfDS.dSigmaede2 = as<mat>(model_sexp["dSigmaede2"]);
		InfDS.dSigmaede2.print("InfDS.dSigmaede2");
	}
	

	InfDS.bAdaptParams = as<mat>(model_sexp["bAdaptParams"]);
	InfDS.bAdaptParams.print("InfDS.bAdaptParams");
	

	if (InfDS.NxState > 0 && InfDS.Ny > 0){
		InfDS.Lambda = as<mat>(model_sexp["lambda"]);
		InfDS.Lambda.print("InfDS.Lambda");
	}
	
	if (InfDS.Nb > 0 && InfDS.Nsubj > 0){
		InfDS.b = as<mat>(model_sexp["b"]);
		InfDS.b.print("InfDS.b");
	}
	
	if (InfDS.Ntheta > 0 && InfDS.Nsubj > 0 && InfDS.Nbetax > 0){
		InfDS.H = as<mat>(model_sexp["H"]);
		InfDS.H.print("InfDS.H");
	}
	
	/*-----*/
	if (InfDS.Nb > 0 && InfDS.Ntheta > 0){
		InfDS.Z = as<mat>(model_sexp["Z"]);
		InfDS.Z.print("InfDS.Z");
	}
	
	
	if (InfDS.Nsubj > 0){
		//InfDS.allT.set_size(1, InfDS.Nsubj);
		InfDS.allT = as<rowvec>(model_sexp["allT"]);
		InfDS.allT.print("InfDS.allT");
	}
	
	
	if (InfDS.totalT > 0 ){
		InfDS.tspan = as<rowvec>(model_sexp["tspan"]);
		InfDS.tspan.print("InfDS.tspan");
	}
	
	if (InfDS.Nmu){
		InfDS.mu = as<mat>(model_sexp["mu"]);
		InfDS.mu.print("InfDS.mu");
	}
	
	/*-----*/
	
	if(InfDS.Nmu > 0 || InfDS.Ny > 0 || InfDS.NLambda > 0 || InfDS.Nbeta > 0 || InfDS.Nbpar > 0){
		int Npar = InfDS.Nmu + InfDS.Ny + InfDS.NLambda + InfDS.Nbeta + InfDS.Nbpar;
		printf("Npar = %d\n", Npar);
		
		InfDS.lowBound.set_size(1, Npar);
		InfDS.upBound.set_size(1, Npar);
		InfDS.par.set_size(1, Npar);
		InfDS.lowBound = as<arma::mat>(model_sexp["lower_bound"]);
		InfDS.upBound = as<arma::mat>(model_sexp["upper_bound"]);
		InfDS.par = as<arma::mat>(model_sexp["par_value"]);

		InfDS.lowBound.print("InfDS.lower_bound");
		InfDS.upBound.print("InfDS.upper_bound");
		InfDS.par.print("InfDS.par_value");
		
	}
	
	if (InfDS.Ny > 0){
		InfDS.dmudparMu = as<mat>(model_sexp["dmudparMu"]);
		InfDS.dmudparMu.print("InfDS.dmudparMu");
		
		InfDS.dmudparMu2 = as<mat>(model_sexp["dmudparMu2"]);
		InfDS.dmudparMu2.print("InfDS.dmudparMu2");
	}
	
	if (InfDS.Ny > 0 && InfDS.NxState > 0 && InfDS.NLambda > 0){
		InfDS.dLambdparLamb = as<mat>(model_sexp["dLambdparLamb"]);
		InfDS.dLambdparLamb.print("InfDS.dLambdparLamb");
		
		InfDS.dLambdparLamb2 = as<mat>(model_sexp["dLambdparLamb2"]);
		InfDS.dLambdparLamb2.print("InfDS.dLambdparLamb2");
	}
	
	if (InfDS.Nsubj > 0 && InfDS.NxState > 0){
		InfDS.y0 = as<mat>(model_sexp["y0"]);
		//InfDS.y0.print("InfDS.y0");
	}

	
	if (InfDS.Nu > 0){
		Rcpp::List U1_sexp = data_sexp["covariates"];
		char str_name[64];
		
		for(int u = 0;u < InfDS.Nu; u++){
			sprintf(str_name, "covar%lu", (long unsigned int) u+1);
			//Rcpp::NumericVector x = U1_sexp[str_name];
			//InfDS.U1.col(u) = as<Col>U1_sexp[str_name];
			NumericVector temp = U1_sexp[str_name];
			if(u == 0)
				InfDS.U1.set_size(temp.length(), InfDS.Nu);
			
			InfDS.U1.col(u)= as<arma::vec>(temp);
		}
		
		//InfDS.U1.print("InfDS.U1");
    }else{
        InfDS.U1.set_size(0, 0);
    }
	

	//Y
	//todo:: remove Y from model list
	if (InfDS.Ny > 0){
		
		Rcpp::List Y_sexp = data_sexp["observed"];
		char str_name[64];
		
		InfDS.Y.set_size(InfDS.Nsubj);
		for(int i = 0; i < InfDS.Nsubj; i++){
			InfDS.Y(i).set_size(InfDS.Ny, InfDS.allT(i));
			//Rprintf("%d %d\n", i, InfDS.allT(i));
		}
		for(int u = 0;u < InfDS.Ny; u++){
			sprintf(str_name, "obs%lu", (long unsigned int) u+1);
			NumericVector temp = Y_sexp[str_name];
			
			int tobs_pointer = 0;
			for(int i = 0; i < InfDS.Nsubj; i++){
				for(int j = 0; j < InfDS.allT(i); j++){
					InfDS.Y(i)(u,j) = temp[tobs_pointer];
					tobs_pointer++;
				}
			}
		}
    }else{
        InfDS.Y.set_size(0);
    }
	InfDS.Y[0](span(0,2), span(0,9)).print("InfDS.Y[0]");
	InfDS.Y[1](span(0,2), span(0,9)).print("InfDS.Y[1]");
	InfDS.Y[199](span(0,2), span(0,9)).print("InfDS.Y[199]");

	//timeDiscrete & tobs
	//NumericVector timeDiscrete_vec = model_sexp["time_"];
	NumericVector timeDiscrete_vec = data_sexp["time"];
	NumericVector tobs_vec = model_sexp["tobs"];
	//Rcout << "The value of v : " << timeDiscrete_vec << "\n";
	InfDS.timeDiscrete.set_size(InfDS.Nsubj);
	InfDS.tobs.set_size(InfDS.Nsubj);
	int tobs_pointer = 0;
	for(int i = 0; i < InfDS.Nsubj; i++){
		InfDS.timeDiscrete(i).set_size(InfDS.allT(i), 1);
		InfDS.tobs(i).set_size(InfDS.allT(i), 1);
		for(int j = 0; j < InfDS.allT(i); j++){
			InfDS.timeDiscrete(i)(j,0) = timeDiscrete_vec[tobs_pointer];
			InfDS.tobs(i)(j,0) = tobs_vec[tobs_pointer];
			tobs_pointer++;
		}
	}
	InfDS.timeDiscrete[0](span(0,9), span::all).print("InfDS.timeDiscrete[0]");
	InfDS.timeDiscrete[199](span(0,9), span::all).print("InfDS.timeDiscrete[199]");
	InfDS.tobs[0](span(0,9), span::all).print("InfDS.tobs[0]");
	InfDS.tobs[199](span(0,9) , span::all).print("InfDS.tobs[199]");
	
	
	Rcpp::List func_addr_sexp = model_sexp["func_address"];
	
	//C_FUNC_POINTER func_addr;
	//func_addr.p_dynfunICM = func_addr_sexp["f_dyn"];
	//func_addr.p_test = func_addr_sexp["f_test"];
	
	*(void **) (&InfDS.fp.dynfunICM) = R_ExternalPtrAddr(func_addr_sexp["f_dyn"]);
	*(void **) (&InfDS.fp.dfdxFreeICM) = R_ExternalPtrAddr(func_addr_sexp["f_dfdx"]);
	*(void **) (&InfDS.fp.dfdparFreeIC) = R_ExternalPtrAddr(func_addr_sexp["f_dfdp"]);
	*(void **) (&InfDS.fp.dfdx2FreeIC) = R_ExternalPtrAddr(func_addr_sexp["f_dfdx2"]);
	*(void **) (&InfDS.fp.dfdxdpFreeIC) = R_ExternalPtrAddr(func_addr_sexp["f_dfdxdp"]);
	*(void **) (&InfDS.fp.dfdpdxFreeIC) = R_ExternalPtrAddr(func_addr_sexp["f_dfdpdx"]);
	*(void **) (&InfDS.fp.dfdpar2FreeIC) = R_ExternalPtrAddr(func_addr_sexp["f_dfdp2"]);
	*(void **) (&InfDS.fp.setParsFreeICwb) = R_ExternalPtrAddr(func_addr_sexp["f_setpars"]);
	*(void **) (&InfDS.fp.test) = R_ExternalPtrAddr(func_addr_sexp["f_test"]);


	//arma::mat x = InfDS.fp.test();
	//x.print("output of hello world");
	
	// assign output arbitrary values and return it back for testing
	output.convFlag =  output.nIterStage1 = output.nIterStage2 = 824;
	output.ss = 0.9;
	output.avebAccept = 0.78;
	output.Iytild.set_size(5,5);
	output.Iytild.fill(0.09);
	output.thetatild.set_size(5);
	output.thetatild.fill(0.7);
	

	char filenamePar[64] = "./Results/TrueInitparG1.txt";
	char filenameSE[64] = "./Results/TrueInitSEG1.txt";
	char filenameconv[64] = "./Results/TrueInitconvG1.txt";
	char filenamebhat[64] = "./Results/TrueInitbhatG1.txt";
	char filenamebhat2[64] = "./Results/TrueInitbhat2G1.txt";
	int kk = 1;
	int trueInit = 1;
	int batch = 1;

	arma::mat x1;
	
	//InfDS.trueb = as<rowvec>(model_sexp["trueb"]);
	//InfDS.trueb.print("InfDS.trueb");
	
	
	saem_estimation(InfDS, InfDS0, upperb, lowerb, x1, filenamePar, filenameSE, filenameconv, filenamebhat, filenamebhat2, kk, trueInit, batch, seed, isfreeIC, output);

	return Rcpp::List::create(Rcpp::Named("convFlag") = output.convFlag, 
	                          Rcpp::Named("nIterStage1") = output.nIterStage1, 
	                          Rcpp::Named("nIterStage2") = output.nIterStage2,  
							  Rcpp::Named("ss") = output.ss,
							  Rcpp::Named("avebAccept") = output.avebAccept,
							  Rcpp::Named("Iytild") = output.Iytild,
							  Rcpp::Named("Iytild") = output.thetatild);
}