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
			r(0, s)= y(1, s);
			r(1, s)= -61.68503 * y(0, s) + thetaf(0,s) * (1 - pow(y(0, s), 2)) * y(1, s);
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
			r.slice(s)(1,0) = 1;
			r.slice(s)(0,1) = -((thetaf(0,s)) * (2 * y(0,s)) * y(1,s) + 61.68503);
			r.slice(s)(1,1) = (thetaf(0,s)) * (1 - pow(y(0,s), 2)); 
			if(isPar == 1){
				r.slice(s)(2,1) = (1 - pow(y(0,s), 2)) * y(1,s);
				r.slice(s)(3,1) = InfDS.U1(s,0) * (1 - pow(y(0,s), 2)) * y(1,s);
				r.slice(s)(4,1) = InfDS.U1(s,1) * (1 - pow(y(0,s), 2)) * y(1,s); 
				r.slice(s) = r.slice(s).t();
			}
		}
	}
	return r;
}

//----------------
arma::cube dfdparFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){
	arma::mat y ;
	arma::cube r;

	// if i is empty, traverse all vectors
	if(i.is_empty()){
		i = span_vec(1, InfDS.Nsubj, 1);
	}
	y = xin;
	r = arma::zeros<arma::cube>(InfDS.Ntheta, InfDS.NxState, y.n_cols);
 	int s;
	for (s = 0; s < int(y.n_cols); s++){		
		r.slice(s)(0,1) = (1 - pow(y(0,s), 2)) * y(1,s); 
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

	thetaf=calculateTheta(0, y, i,InfDS);

 	int s;
	for (s = 0; s < int(y.n_cols); s++){		
		r.slice(s)(2,0) = -(thetaf(0,s) * 2 * y(1,s));
		r.slice(s)(2,1) = -(thetaf(0,s) * (2 * y(0,s)));
		r.slice(s)(3,0) = -(thetaf(0,s) * (2 * y(0,s))); 
	}
	return r;
}

//----------------
arma::cube dfdxdpFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){
	// i and t are dummy variables
	arma::mat y;
	arma::cube r;

	// if i is empty, traverse all vectors
	if(i.is_empty()){
		i = span_vec(1, InfDS.Nsubj, 1);
	}

	y = xin ;  
	r = arma::zeros<arma::cube>(InfDS.Nx * InfDS.Nx, InfDS.Ntheta, y.n_cols);
 	int s;
	for (s = 0; s < int(y.n_cols); s++){		
		r.slice(s)(2,0) = -((2 * y(0,s)) * y(1,s));
		r.slice(s)(3,0) = (1 - pow(y(0,s), 2)); 
	}
	return r;
}

//----------------
arma::cube dfdpdxFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){
	// i and t are dummy variables
	arma::mat  y;
	arma::cube r;

	// if i is empty, traverse all vectors
	if(i.is_empty()){
		i = span_vec(1, InfDS.Nsubj, 1);
	}

	y = xin ;  
	r = arma::zeros<arma::cube>(InfDS.Nx*InfDS.Ntheta, InfDS.Nx, y.n_cols) ;
 	int s;
	for (s = 0; s < int(y.n_cols); s++){		
		r.slice(s)(3,0) = -(2 * y(0,s) * y(1,s));
		r.slice(s)(3,1) = (1 - pow(y(0,s), 2)); 
	}
	return r;
}

//----------------
arma::cube dfdpar2FreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){
	// i and t are dummy variables
	arma::mat y ;
	arma::cube r;

	// if i is empty, traverse all vectors
	if(i.is_empty()){
		i = span_vec(1, InfDS.Nsubj, 1);
	}

	y = xin ;  
	r = arma::zeros<arma::cube>(InfDS.NxState*InfDS.Ntheta, InfDS.Ntheta, y.n_cols) ;

 	int s;
	for (s = 0; s < int(y.n_cols); s++){
		; 
	}
	return r;
}
     