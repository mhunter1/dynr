
arma::mat dynfunICM(const int isPar, const arma::mat &xin, arma::vec &i, const int t, const int isStart, struct C_INFDS &InfDS){

	//local parameters
	arma::mat y, r, thetaf;
	
	// if i is empty, traverse all vectors
	if(i.is_empty()){
		i = span_vec(1, InfDS.Nsubj, 1);
	}	

	//input parameters
	y = xin;
	r.set_size(InfDS.Nx, int(i.n_elem));
	r.clear();
	r = y;

	thetaf=calculateTheta(isPar, y, i,InfDS);

	if (isStart==1){
		r.zeros();
		int row, s;
		for (s = 0; s < int(i.n_elem); s++){
			for (row = 0; row < InfDS.NxState; row++)
				r(row, s) = thetaf(row +1, s);

			if (isPar == 1){
				for (row = InfDS.NxState; row < InfDS.NxState + InfDS.Nbeta; row++){
 					r(row, s) = y(row, s);
 				}
 			}
		}

	}
	else{
		r.zeros();
		int row, s;
		for (s = 0; s < int(i.n_elem); s++){
			r(0, s)= exp(thetaf(0,s)) * (InfDS.par(6,s) - y(0, s)) + InfDS.par(2,s) * (InfDS.par(7,s) - y(1, s));
			r(1, s)= exp(thetaf(1,s)) * (InfDS.par(7,s) - y(1, s)) + InfDS.par(5,s) * (InfDS.par(6,s) - y(0, s));
			if (isPar == 1){
				for (row = InfDS.NxState; row < InfDS.NxState + InfDS.Nbeta; row++){
					r(row, s) = 0;
				}
			}
		}
	}
	return r;
}

//----------------
arma::cube dfdxFreeICM(const int isPar, arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){
	arma::mat Hi, b, bb, thetaf, rowNum, y;
	arma::cube r;

	// if i is empty, traverse all vectors
	if(i.is_empty()){
		i = span_vec(1, InfDS.Nsubj, 1);
	}

	y = xin;
	r = arma::zeros<arma::cube>(InfDS.Nx, InfDS.Nx, y.n_cols) ;
	thetaf = arma::zeros<arma::mat>(InfDS.Ntheta, y.n_cols) ;

	thetaf=calculateTheta(isPar, y, i,InfDS);

 	if (isStart==1){
		if (isPar == 0)
			; //Undefined case in dfdxParIC
		else{
			int s;
			for (s = 0; s < int(y.n_cols); s++){
				r.slice(s)(2, 2)=1 ;
				r.slice(s)(3, 3)=1 ;
				r.slice(s)(4, 4)=1 ;
				r.slice(s)(5, 5)=1 ;
				r.slice(s)(6, 6)=1 ;
			} 
		}
	}
 	else{
		int s;

		for (s = 0; s < int(y.n_cols); s++){
			r.slice(s)(0,0) = -exp(thetaf(0,s));
			r.slice(s)(1,0) = -InfDS.par(2,s);
			r.slice(s)(0,1) = -InfDS.par(5,s);
			r.slice(s)(1,1) = -exp(thetaf(1,s)); 
			if(isPar == 1){ 
				r.slice(s) = r.slice(s).t();
			}
		}
	}
	return r;
}

//----------------
arma::cube dfdparFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){
	arma::mat y, thetaf;
	arma::cube r;

	// if i is empty, traverse all vectors
	if(i.is_empty()){
		i = span_vec(1, InfDS.Nsubj, 1);
	}
	y = xin;
	r = arma::zeros<arma::cube>(InfDS.Ntheta, InfDS.NxState, y.n_cols);

	thetaf=calculateTheta(0, y, i,InfDS);

 	int s;
	for (s = 0; s < int(y.n_cols); s++){		
		r.slice(s)(0,0) = exp(thetaf(0,s)) * (InfDS.par(6,s) - y(0,s));
		r.slice(s)(1,1) = exp(thetaf(1,s)) * (InfDS.par(7,s) - y(1,s)); 
	}
	return r;
}

//----------------
arma::cube dfdx2FreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){
	// i and t are dummy variables
	arma::mat thetaf, y;
	arma::cube r;

	// if i is empty, traverse all vectors
	if(i.is_empty()){
		i = span_vec(1, InfDS.Nsubj, 1);
	}

	y = xin ;  
	r = arma::zeros<arma::cube>(InfDS.Nx * InfDS.Nx, InfDS.Nx, y.n_cols) ;
	thetaf = arma::zeros<arma::mat>(InfDS.Ntheta, y.n_cols) ;
 	int s;
	for (s = 0; s < int(y.n_cols); s++){		 
	}
	return r;
}

//----------------
arma::cube dfdxdpFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){
	// i and t are dummy variables
	arma::mat y, thetaf;
	arma::cube r;

	// if i is empty, traverse all vectors
	if(i.is_empty()){
		i = span_vec(1, InfDS.Nsubj, 1);
	}

	y = xin ;  
	r = arma::zeros<arma::cube>(InfDS.Nx * InfDS.Nx, InfDS.Ntheta, y.n_cols);

	thetaf=calculateTheta(0, y, i,InfDS);

 	int s;
	for (s = 0; s < int(y.n_cols); s++){		
		r.slice(s)(0,0) = -exp(thetaf(0,s));
		r.slice(s)(3,1) = -exp(thetaf(1,s)); 
	}
	return r;
}

//----------------
arma::cube dfdpdxFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){
	// i and t are dummy variables
	arma::mat  y, thetaf;
	arma::cube r;

	// if i is empty, traverse all vectors
	if(i.is_empty()){
		i = span_vec(1, InfDS.Nsubj, 1);
	}

	y = xin ;  
	r = arma::zeros<arma::cube>(InfDS.Nx*InfDS.Ntheta, InfDS.Nx, y.n_cols) ;

	thetaf=calculateTheta(0, y, i,InfDS);

 	int s;
	for (s = 0; s < int(y.n_cols); s++){		
		r.slice(s)(0,0) = -exp(thetaf(0,s));
		r.slice(s)(3,1) = -exp(thetaf(1,s)); 
	}
	return r;
}

//----------------
arma::cube dfdpar2FreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){
	// i and t are dummy variables
	arma::mat y, thetaf;
	arma::cube r;

	// if i is empty, traverse all vectors
	if(i.is_empty()){
		i = span_vec(1, InfDS.Nsubj, 1);
	}

	y = xin ;  
	r = arma::zeros<arma::cube>(InfDS.NxState*InfDS.Ntheta, InfDS.Ntheta, y.n_cols) ;


	thetaf=calculateTheta(0, y, i,InfDS);

 	int s;
	for (s = 0; s < int(y.n_cols); s++){
		;
		r.slice(s)(0,0) = exp(thetaf(0,s)) * (InfDS.par(6,s) - y(0,s));
		r.slice(s)(3,1) = exp(thetaf(1,s)) * (InfDS.par(7,s) - y(1,s)); 
	}
	return r;
}

void setParsFreeICwb(C_INFDS &InfDS){
	int startM, startL, starte, startb;
	arma::mat D, L;
	
	InfDS.par = real(InfDS.par);
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//% Model specification: modify based on fitted model
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
	InfDS.mu = arma::zeros<arma::mat> (InfDS.Ny, 1);
	if (InfDS.Nmu > 0){
		startM = InfDS.Nbeta+1;
		InfDS.mu = InfDS.par(span(startM - 1, startM + 1), span::all);
		InfDS.dmudparMu = eye(InfDS.Ny, InfDS.Ny); //dmu/dtheta_mu
		InfDS.dmudparMu2 = arma::zeros<arma::mat>(InfDS.Ny*InfDS.Ny, InfDS.Ny);
	}
	
	if (InfDS.NLambda > 0){
		startL = InfDS.Nbeta+ InfDS.Nmu+1;
		InfDS.Lambda= "1, 0; 0, 0; 0,0";
		InfDS.Lambda(1,0) = InfDS.par(startL-1);
		InfDS.Lambda(2,0) = InfDS.par(startL);
		
		InfDS.dLambdparLamb = arma::zeros<arma::mat>(InfDS.NLambda, InfDS.Ny*InfDS.Nx);
		InfDS.dLambdparLamb(0,1) = 1;
		InfDS.dLambdparLamb(1,2) = 1;
		InfDS.dLambdparLamb2 = arma::zeros<arma::mat>(InfDS.Nx*InfDS.Ny*InfDS.NLambda, InfDS.NLambda);
	}

	
	starte = InfDS.Nbeta + InfDS.Nmu + InfDS.NLambda + 1;
	InfDS.Sigmae = diagmat(exp(InfDS.par(span(starte-1, starte+1), span::all)));
	InfDS.dSigmaede = arma::zeros<arma::mat>(3, 9);
	InfDS.dSigmaede(0, 0) = exp(InfDS.par(starte-1));	
	InfDS.dSigmaede(1, 4) = exp(InfDS.par(starte));
	InfDS.dSigmaede(2, 8) = exp(InfDS.par(starte+1));

                                 
	InfDS.dSigmaede2 = arma::zeros<arma::mat>(InfDS.Ny*InfDS.Ny*InfDS.Ny, InfDS.Ny);
	InfDS.dSigmaede2(0, 0) = exp(InfDS.par(starte-1));
	InfDS.dSigmaede2(13, 1) = exp(InfDS.par(starte));
	InfDS.dSigmaede2(26, 2) = exp(InfDS.par(starte+1));
	

	startb = InfDS.Nbeta + InfDS.Nmu + InfDS.NLambda + InfDS.Ny + 1; 
	InfDS.dSigmabdb(0, 0) = exp(InfDS.par(startb-1));
	InfDS.dSigmabdb2(0, 0) = exp(InfDS.par(startb-1));
	
	/*
	if (InfDS.Nb > 0){        
		startb = InfDS.Nbeta + InfDS.Nmu + InfDS.NLambda + InfDS.Ny + 1;    


		if (exp(InfDS.par(startb-1))==0){ 
			InfDS.par(startb-1) = log(10e-7);
		}
		if (exp(InfDS.par(startb))==0){
			InfDS.par(startb) = log(10e-7);
		}
		
		D = diagmat(exp(InfDS.par(span(startb-1,startb+InfDS.Nbpar-3), span::all)));
		L = diagmat(arma::ones<arma::vec>(3)); 
		L(2, 1) = InfDS.par(startb+2);
		
		InfDS.Sigmab = L * D * L.t();
		InfDS.dSigmabdb = arma::zeros<arma::mat>(InfDS.Nbpar,InfDS.Nb * InfDS.Nb);
		InfDS.dSigmabdb(0,0) = exp(InfDS.par(startb - 1));
		InfDS.dSigmabdb(1,4) =  exp(InfDS.par(startb));
		InfDS.dSigmabdb(1,5) = InfDS.par(startb+2) * exp(InfDS.par(startb));
		InfDS.dSigmabdb(1,7) = InfDS.par(startb+2) * exp(InfDS.par(startb));
		InfDS.dSigmabdb(1,8) = InfDS.par(startb+2) * InfDS.par(startb+2) * exp(InfDS.par(startb));
		
		InfDS.dSigmabdb(2,8) = exp(InfDS.par(startb+1));
		InfDS.dSigmabdb(3,5) = exp(InfDS.par(startb));
		InfDS.dSigmabdb(3,7) = exp(InfDS.par(startb));
		InfDS.dSigmabdb(3,8) = 2 * InfDS.par(startb+2) * exp(InfDS.par(startb));


		InfDS.dSigmabdb2 = arma::zeros<arma::mat>(InfDS.Nbpar * InfDS.Nb * InfDS.Nb, InfDS.Nbpar);
		InfDS.dSigmabdb2(0,0) = exp(InfDS.par(startb - 1));
		InfDS.dSigmabdb2(17,1) = exp(InfDS.par(startb));
		InfDS.dSigmabdb2(21,1) = InfDS.par(startb + 2) * exp(InfDS.par(startb));
		InfDS.dSigmabdb2(21,3) = exp(InfDS.par(startb));
		InfDS.dSigmabdb2(23,1) = exp(InfDS.par(startb));
		InfDS.dSigmabdb2(29,1) = InfDS.par(startb + 2) * exp(InfDS.par(startb));
		InfDS.dSigmabdb2(29,3) = exp(InfDS.par(startb));
		InfDS.dSigmabdb2(31,1) = exp(InfDS.par(startb));
		InfDS.dSigmabdb2(33,1) = InfDS.par(startb + 2) * InfDS.par(startb+2) * exp(InfDS.par(startb));
		InfDS.dSigmabdb2(33,3) = 2 * (InfDS.par(startb + 2)) * exp(InfDS.par(startb));
		InfDS.dSigmabdb2(34,2) = exp(InfDS.par(startb + 1));
		InfDS.dSigmabdb2(35,1) = 2 * (InfDS.par(startb + 2))*exp(InfDS.par(startb));
		InfDS.dSigmabdb2(35,3) = 2 * exp(InfDS.par(startb));
		
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	} //End of Nb > 0
	*/
	
	return;
}
   