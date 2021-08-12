#include <time.h>
#include <math.h>
#include <Rcpp.h>

// follow the function of matlab
double quantile(arma::vec X, double p){
	double threshold;
	double devision;
	int id;
	if (X.is_sorted("ascend") == false)
		X = sort(X, "ascend");
	
	threshold = (0.5)/X.n_elem;
	if(p <= threshold)
		return X(0);	//first number: 0.5/n_elem-qunatile
	else if(p >= 1 - threshold)
		return X(X.n_elem - 1);	// last number: (n-0.5)/n_elem
	
	// else case
	devision = (p - threshold) * X.n_elem;
	id = int(devision + 0.0000001);
	return X(id) + (X(id+1) - X(id)) * (devision - id);
}

//Reminder: judge whether X and Y are numbers
double round2(double x, double y){
	return round(x/y)*y;
}

//Reminder: the space needs to be created by caller
arma::vec randsample(arma::vec &source, int length, int size){
	int i, rand_id;
	arma::vec result;
	//srand (time(NULL));
	result.set_size(size);
	for (i = 0; i < size; i++){
		rand_id = rand() % length;
		result(i) = source(rand_id);
	}
	
	return result;
}

//This function return a column vector [start, start + step, start+ step * 2, ..., end]^T
//use span(...).t() to get the rowvector

arma::vec span_vec(int start, int end, int step){
	arma::vec x;
	int i;

	x.set_size((end - start)/step + 1);
	x[0] = start;
	for(i = 1; i < (end - start)/step + 1; i++)
		x[i]= x[i-1] + step;
	return x;
}


arma::mat rowProjection(arma::mat data, arma::mat index){
	//index is a row vector
	int i;
	arma::mat output;
	output.reset();
	//data.print("data =");
	//index.print("index = ");
	for (i = 0; i < (int)index.n_cols; i++){
		output.insert_rows(output.n_rows, data.row(index(i)-1));
	}
	//output.print("output");
	return output;
}

arma::mat colProjection(arma::mat data, arma::mat index){
	//index is a row vector
	int i;
	arma::mat output;
	output.reset();
	//data.print("data =");
	//index.print("index = ");
	for (i = 0; i < (int)index.n_cols; i++){
		output.insert_cols(output.n_cols, data.col(index(i)-1));
	}
	//output.print("output");
	return output;
}

arma::mat modProjection(arma::mat data, int divisor, int remainder){
	// equal to data(mod(vec, divisor)==1)
	// vec = [1 2 3 4 5]; vec(mod(vec, 3)==2) = [2 5]
	int i;

	data = data.row(0);
	arma::mat output;
	output.reset();
	for (i = 0; i < (int)data.n_cols; i++){
		if (int(data(i)) % divisor == remainder){
		  
			output.insert_cols(i/divisor, data.col(i));
		}
	}
	//output.print("data(output)");
	return output;
  
}

int maxElement(arma::imat integer_matrix){
	int i, j;
	int max = -2144783647;
	for (i = 0; i < (int)integer_matrix.n_rows; i++){
		for(j = 0; j < (int)integer_matrix.n_cols;j++){
			if (integer_matrix(i,j) > max)
				max = integer_matrix(i,j);
		}
	}
	return max;
}


arma::mat calculateTheta(const int isPar, const arma::mat &y, arma::vec &i, struct C_INFDS &InfDS){
	
	if(isPar == 1)
		printf("[HJ] Should NOT happen! isPar = 1");
	
	arma::mat b, betai, temp, Hi, thetaf, par;
	
	thetaf = arma::zeros<arma::mat>(InfDS.Ntheta, y.n_cols);
	par = InfDS.par;

        //printf("thetaf %d %d\n", thetaf.n_rows, thetaf.n_cols);
    if (isPar==0){
        //printf("InfDS.Nbeta %d\n", InfDS.Nbeta);
        temp = span_vec(1, InfDS.Nbeta, 1).t();//Nbeta is the number of fixed effects parameters
        par = rowProjection(par, temp);
       // par.print("par");
    }
    
	b = arma::zeros<arma::mat>(thetaf.n_rows, thetaf.n_cols) ;	
	betai = arma::zeros<arma::mat>(InfDS.Ntheta, 1);	// beta
    //printf("i n elem %d\n", i.n_elem);
    for (int ii = 0; ii < int(i.n_elem); ii++){
        //printf("ii = %d\n", ii);
        int current_i = int(i(ii));   
        
		Hi = InfDS.H.rows(1+(current_i-1)*InfDS.Ntheta-1, current_i*InfDS.Ntheta-1);   
                //Hi.print("Hi"); 
        
		if (isPar == 1)            
            betai = y.col(ii).rows(InfDS.NxState,InfDS.NxState + Hi.n_cols-1);
        else
            betai = reshape(par,InfDS.Nbeta, 1);
  		
	        //betai.print("betai");
		if (InfDS.Nb > 0)
			b.col(ii) = reshape(InfDS.b.row(current_i-1),InfDS.Nb,1);
                //b.col(ii).print("b col i");
	
		thetaf.col(ii) = Hi * betai + InfDS.Z * b.col(ii);
    }
    //printf("loop done");
    //thetaf.print("thetaf");
	
	return thetaf;
}


arma::vec ekfContinuous10(int Nsubj, const int N, const int Ny, const int Nx, const int Nb, const int NxState,const arma::mat Lambda, int totalT, const arma::mat Sigmae, const arma::mat Sigmab, const arma::mat mu, const arma::mat b, const arma::mat allT, const arma::cube Xtild, arma::field<arma::mat> Y){

	int i, tt, j;
	arma::mat Jytmp, chiYX, iSigmae, indexY, tp, yPredtmp;



	tp.set_size(Nsubj*N, 1);
	tp.zeros();

	chiYX.set_size(Nsubj, totalT);
	chiYX.zeros();


	iSigmae = inv(Sigmae); 
	//iSigmae.print("iSigmae");

	for (tt=1; tt<=totalT; tt++){
		for (i=1; i<=Nsubj; i++){
			if ( tt<=(int)allT(i-1)){
				indexY.set_size(1, Ny);
				indexY.zeros();
				for (j = 1; j <= Ny; j++){
					indexY(0,j-1) = j;
					if(std::isnan(Y(i-1, 0)(j-1,tt-1))){
					//if(std::isnan(Y(i-1,j-1,tt-1))){	
						indexY(0,j-1) = 0;
						//printf("%d %d %d\n",i-1,j -1,tt-1);
					}
				}

				indexY = sort(indexY);
				for (j = Ny; j > 0; j--){
					if (indexY(0, j-1) == 0)
						indexY.shed_col(j-1);
				}
				//mexPrintf("enter vdpmeas %d %d %d\n", Xtild.n_rows, Xtild.n_cols, Xtild.n_slices);
				//Xtild.slice(tt-1).col(i-1).print("Xtild.slice(tt-1).col(i-1)");
				VDPMeas(Xtild.slice(tt-1).col(i-1), Ny, Nx, NxState, Lambda, mu, &yPredtmp, &Jytmp);
				//printf("exit vdpmeas\n");
				
				arma::mat inovtmp, result;
				//inovtmp = rowProjection(Y.slice(tt-1).row(i-1).t(),indexY)- rowProjection(yPredtmp,indexY);
				inovtmp = rowProjection(Y(i-1, 0).col(tt-1), indexY)- rowProjection(yPredtmp,indexY);
				inovtmp.insert_cols(1,arma::zeros(indexY.n_cols, indexY.n_cols-1));

				result = (inovtmp.t()*colProjection(rowProjection(iSigmae, indexY), indexY)*inovtmp);
				chiYX(i-1,tt-1) = result(0,0);
			}
		}
	}
	//printf("middle of ekf\n");
	
	if (Nb>0){
		int i;
		arma::mat iSigmab;
		//Sigmab.print("Sigmab");
		iSigmab = inv(Sigmab);


		tp = b.row(1-1) * iSigmab;
		for (i = 2; i <= (int)b.n_rows; i++){
			tp = join_cols(tp, (b.row(i-1) * iSigmab));
			
		}
		tp = -0.5 *(sum(tp % b, 1)+ sum(chiYX, 1));
	}
	return tp ;
}

C_INFDS getXtildIC3(const int isPar, const int getDxFlag, const int freeIC, struct C_INFDS &InfDS){

	arma::mat dXtildPrev0, d2XtildPrev0, tcount, delt, tspan, XtildPrev;
	arma::cube fullX, xk1, xk2, xk21, tpsan, dXtildPrev, d2XtildPrev, dXtild, dXstar_t, d2Xstar_t, d2Xtild, dfdxATk21,d2vecdk1,first;
	arma::mat k1, k21, k2, Xtild_t; 
	arma::cube dk1, dk2, dk21, dk1dtheta, dk2dtheta, dvecdfdxATk21dtheta, dvecdfdxATk21dk21, dfdxATXtildprev, dvecdfdXtilddtheta, dveck1dtheta, dveck1dthetadXtild, dveck2dtheta, dveck2dthetadxATk21, dvecdfdXtilddXtild, second, third, d2vecdk2;
	arma::vec indext, currentt;
	arma::vec tindex, dt;
	
	static arma::vec empty_vec = span_vec(1, InfDS.Nsubj,1);
	int T, i, Nsubj;
	
	//printf("execution 1\n");
	//return InfDS;
	
	if(isPar == 1)
		Rprintf("[HJ] Should NOT happen! isPar = 1");
	
	InfDS.Xtild = arma::zeros<arma::cube>(InfDS.Nx,InfDS.Nsubj,InfDS.totalT);
	if (getDxFlag ==1){
		InfDS.dXstarAll.set_size(InfDS.Nsubj);
		InfDS.dXstarAll2.set_size(InfDS.Nsubj);
		InfDS.dXtildthetafAll.set_size(InfDS.Nsubj);
		InfDS.dXtildthetafAll2.set_size(InfDS.Nsubj);
	}  
	
	tspan = InfDS.tspan;
	fullX = arma::zeros<arma::cube>(InfDS.Nx,InfDS.Nsubj,InfDS.tspan.n_cols);
	//Rprintf("execution 1.1\n");
	
	T = InfDS.tspan.n_cols;
	xk1 = arma::zeros<arma::cube>(InfDS.Nx,InfDS.Nsubj,T);
	xk2 = arma::zeros<arma::cube>(InfDS.Nx,InfDS.Nsubj,T);
	tindex = InfDS.tspan.t();
	//tindex.print("tindex");
	//Rprintf("execution 1.2\n");
	
	delt.set_size(1,1);
	delt(0,0) = InfDS.delt;
	dt = repmat(delt, T, 1);
	Nsubj = InfDS.Nsubj;
	//tcount = arma::ones<arma::mat>(Nsubj,1);
	//Rprintf("execution 1.3\n");

	for (i = 0; i < InfDS.Nsubj;i++){
		if (getDxFlag ==1){
			d2Xtild = arma::zeros<arma::cube>(InfDS.Ntheta*InfDS.Nx, InfDS.Ntheta,1); //Running sum of second derivatives
			InfDS.dXstarAll(i) = arma::zeros<arma::mat>(InfDS.Ntheta, InfDS.Nx*InfDS.allT(i));
			InfDS.dXstarAll2(i) = arma::zeros<arma::mat>(InfDS.Ntheta*InfDS.Nx*InfDS.allT(i), InfDS.Ntheta);
			InfDS.dXtildthetafAll(i) = arma::zeros<arma::mat>(InfDS.Ntheta, InfDS.Nx*InfDS.allT(i));
			InfDS.dXtildthetafAll2(i) = arma::zeros<arma::mat>(InfDS.Ntheta*InfDS.Nx*InfDS.allT(i),InfDS.Ntheta);
		}
    
	}
	//printf("execution 1.4~~~ %d\n", isPar);
	//return InfDS;
	/*
	if(InfDS.NxState == InfDS.Nx)
		isPar = 0;
	else
		isPar = 1;
	*/
	
	//Rprintf("execution 1.4~~~\n");

	
	//odd setting
	//XtildPrev = InfDS.fp.dynfunICM(isPar, trans(InfDS.y0), empty_vec, 0, 1, InfDS);
	//change according to symiin's testing in July 2020
	XtildPrev = trans(InfDS.y0);

	dXtildPrev0 = arma::zeros<arma::mat>(InfDS.Ntheta,InfDS.Nx); 

	d2XtildPrev0 = arma::zeros<arma::mat>(InfDS.Nx*InfDS.Ntheta, InfDS.Ntheta); 

	
	
	/*
	if (freeIC==1){
		//printf("execution 1.4.1\n");
		//trans(InfDS.y0).print("ddd");
		XtildPrev = InfDS.fp.dynfunICM(isPar, trans(InfDS.y0), empty_vec, 0, 1, InfDS);	
		//InfDS.par.print("InfDS.par getXtildIC3 dynfunICM");

		// ****** size of dXtildPrev0 should be determined dynamically
		//dXtildPrev0 =  "0 0; 1 0; 0 1"; 
		dXtildPrev0 = arma::zeros<arma::mat>(InfDS.Ntheta,InfDS.Nx);
		d2XtildPrev0 = arma::zeros<arma::mat>(InfDS.Nx*InfDS.Ntheta,InfDS.Ntheta);
		//printf("execution 1.4.1\n");
	}
	else{
		//printf("execution 1.4.2\n");
		XtildPrev = InfDS.fp.dynfunICM(isPar, trans(InfDS.y0), empty_vec, 0, 1, InfDS);
		dXtildPrev0 = arma::zeros<arma::mat>(InfDS.Ntheta,InfDS.Nx); 
		d2XtildPrev0 = arma::zeros<arma::mat>(InfDS.Nx*InfDS.Ntheta, InfDS.Ntheta); 
		//printf("execution 1.4.2\n");
	}*/

	
	//InfDS.par.print("InfDS.par getXtildIC3 293");
	if (getDxFlag ==1){
		// Note that: reshape CANNOT convert a matrix into a cube in armadillo
		dXtildPrev.set_size(InfDS.Ntheta, InfDS.Nx * Nsubj, 1);
		dXtildPrev.slice(0) = repmat(dXtildPrev0, 1, Nsubj);
		dXtildPrev = reshape(dXtildPrev, InfDS.Ntheta, InfDS.Nx, Nsubj);
		dXtild = dXtildPrev;
		
		d2XtildPrev.set_size(InfDS.Nx*InfDS.Ntheta, InfDS.Ntheta * Nsubj, 1);
		d2XtildPrev.slice(0) = repmat(d2XtildPrev0, 1, Nsubj);
		d2XtildPrev = reshape(d2XtildPrev, InfDS.Ntheta*InfDS.Nx, InfDS.Ntheta, InfDS.Nsubj);
		d2Xtild = d2XtildPrev;

		dXstar_t = dXtild;
		d2Xstar_t = d2Xtild;
	}

	//printf("execution 2\n");
	
	InfDS.Xtild.slice(0) = XtildPrev;
	fullX.slice(0) = XtildPrev;
	
	if (getDxFlag == 1){ 	
		for( i = 0;i < InfDS.Nsubj; i++){
			InfDS.dXtildthetafAll(i) = arma::zeros<arma::mat>(InfDS.Ntheta, InfDS.Nx*InfDS.allT(i));
			InfDS.dXtildthetafAll(i).cols(0, InfDS.Nx - 1) = dXtildPrev.slice(i);
			InfDS.dXtildthetafAll2(i) = arma::zeros<arma::mat>(InfDS.Ntheta*InfDS.Nx*InfDS.allT(i),InfDS.Ntheta);
			InfDS.dXtildthetafAll2(i).rows(0,InfDS.Ntheta*InfDS.Nx - 1) = d2XtildPrev.slice(i);
		}
	}
    
	//InfDS.par.print("InfDS.par getXtildIC3 325");
	//Rprintf("execution 3\n");
	//Do interpolation at small, equal intervals, as opposed to at observed intervals to avoid numerical problems
	for (int t = 1; t < T; t++){
		//printf("t = %d T= %d\n",t,T);

		//ODE solver	
		/*
		if (freeIC == 0)
			k1 = InfDS.fp.dynfunICM(1, XtildPrev, empty_vec, tindex(t), 0, InfDS);
		else
			k1 = InfDS.fp.dynfunICM(0, XtildPrev, empty_vec, tindex(t), 0, InfDS);
		*/
		
		k1 = InfDS.fp.dynfunICM(isPar, XtildPrev, empty_vec, tindex(t), 0, InfDS);
		k21= XtildPrev+dt(t)*k1;
		//[to be checked] case 1 goes for isPar = 1, case 2 goes for isPar = 0
		/*
		if (freeIC == 0)
			k2 = InfDS.fp.dynfunICM(1, k21, empty_vec, tindex(t)+dt(t), 0, InfDS);
		else
			k2 = InfDS.fp.dynfunICM(0, k21, empty_vec, tindex(t)+dt(t), 0, InfDS);
		*/
		k2 = InfDS.fp.dynfunICM(isPar, k21, empty_vec, tindex(t)+dt(t), 0, InfDS);
		Xtild_t = XtildPrev+dt(t)*(k1+k2)/2;
		//if(t == 1)
		//	Xtild_t.print("Xtild_t @ t = 1");
		xk1.slice(t) = k1;
		xk2.slice(t) = k2;
		
		//printf("execution 4\n");
		if (getDxFlag==1 && t > 1){
			
			dk1dtheta = InfDS.fp.dfdparFreeIC(XtildPrev, empty_vec, tindex(t), 0, InfDS);
			dk2dtheta = InfDS.fp.dfdparFreeIC(k21, empty_vec, tindex(t)+dt(t), 0, InfDS);
			//printf("execution 4.1\n");
			
			//Note that dfdx has x in rows, f in columns, unlike the Jacobian function
			//used in ekf, specified in InfDS.Jdyn
			dfdxATXtildprev = InfDS.fp.dfdxFreeICM(0, XtildPrev, empty_vec, tindex(t), 0, InfDS);
		    dk1 = arma::zeros<arma::cube>(InfDS.Ntheta, InfDS.Nx, InfDS.Nsubj);
			for (i = 0; i < Nsubj; i++){
				//Rprintf("i %d\n",i); 
				dk1.slice(i) = dXtildPrev.slice(i)*dfdxATXtildprev.slice(i) + dk1dtheta.slice(i);
			}
			dk21 = dXtildPrev + dt(t)*dk1; 
			
			
			
			//dfdxATk21 = feval(InfDS.dfdx,k21,InfDS,[],tindex(t)); //df/dx evaluated at k21
			// Does freeIC need to be the freeIC in getXtildIC3?
			dfdxATk21 = InfDS.fp.dfdxFreeICM(0, k21, empty_vec, tindex(t), 0, InfDS);
			dk2 = arma::zeros<arma::cube>(InfDS.Ntheta, InfDS.Nx, InfDS.Nsubj);
			for (i = 0; i < Nsubj; i++){
				dk2.slice(i) = dk21.slice(i)*dfdxATk21.slice(i) + dk2dtheta.slice(i);
			}
			
			dXstar_t = dt(t)/2*(dk1 + dk2);
			dXtild = dXtildPrev + dXstar_t;
			//printf("t = %d\n", t);
			//dXtild.slice(0).rows(0, InfDS.Ntheta-1).print("dXtild"); //correct here
			//printf("execution 4.2\n");
			
			dvecdfdXtilddtheta = InfDS.fp.dfdxdpFreeIC(XtildPrev, empty_vec, tindex(t), 0, InfDS);
			//printf("execution 4.2.1\n");
			dvecdfdXtilddXtild = InfDS.fp.dfdx2FreeIC(XtildPrev, empty_vec, tindex(t), 0, InfDS);
			//printf("execution 4.2.2\n");
			dveck1dthetadXtild = InfDS.fp.dfdpdxFreeIC(XtildPrev, empty_vec, tindex(t), 0, InfDS);
			//printf("execution 4.2.3\n");
			dveck1dtheta = InfDS.fp.dfdpar2FreeIC(XtildPrev, empty_vec, tindex(t), 0, InfDS);
			//printf("execution 4.2.4\n");
			
			d2vecdk1= arma::zeros<arma::cube>(InfDS.Nx*InfDS.Ntheta, InfDS.Ntheta, InfDS.Nsubj);
			first= arma::zeros<arma::cube>(InfDS.Nx*InfDS.Ntheta, InfDS.Ntheta, InfDS.Nsubj);
			
			//Rprintf("execution 5\n");
			for (i = 0; i < InfDS.Nsubj; i++){
				arma:: mat eye1, eye2;
				d2vecdk1.slice(i) = 
				(kron(eye1.eye(InfDS.Nx, InfDS.Nx),dXtildPrev.slice(i))*(dvecdfdXtilddtheta.slice(i)+dvecdfdXtilddXtild.slice(i)*trans(dXtildPrev.slice(i))))+ 
				trans(kron(dfdxATXtildprev.slice(i),eye2.eye(InfDS.Ntheta, InfDS.Ntheta)))* d2XtildPrev.slice(i)+ 
				dveck1dtheta.slice(i) + 
				dveck1dthetadXtild.slice(i)*trans(dXtildPrev.slice(i));
				first.slice(i) = (d2XtildPrev.slice(i) + dt(t)*d2vecdk1.slice(i));
			}
			//printf("execution 5.1\n");

			dvecdfdxATk21dtheta = InfDS.fp.dfdxdpFreeIC(k21, empty_vec, tindex(t)+dt(t), 0, InfDS);
			dvecdfdxATk21dk21 = InfDS.fp.dfdx2FreeIC(k21, empty_vec, tindex(t)+dt(t), 0, InfDS);
			second = arma::zeros<arma::cube>(pow(InfDS.Nx, 2), InfDS.Ntheta, InfDS.Nsubj);
			//printf("execution 5.2\n");
			
			for (i = 0; i < Nsubj; i++){
				second.slice(i) = dvecdfdxATk21dtheta.slice(i) + dvecdfdxATk21dk21.slice(i)*trans(dk21.slice(i)); //Second term, dvecfxdxatk21_2
			}
			//printf("execution 5.3\n");
			
			dveck2dtheta = InfDS.fp.dfdpar2FreeIC(k21, empty_vec, tindex(t)+dt(t), 0, InfDS);
			dveck2dthetadxATk21 = InfDS.fp.dfdpdxFreeIC(k21, empty_vec, tindex(t)+dt(t), 0, InfDS);
			//printf("execution 5.4\n");
			
			third = arma::zeros<arma::cube>(InfDS.Nx*InfDS.Ntheta, InfDS.Ntheta, InfDS.Nsubj);
			d2vecdk2 = arma::zeros<arma::cube>(InfDS.Nx*InfDS.Ntheta, InfDS.Ntheta, InfDS.Nsubj);
			for (i = 0; i < Nsubj; i++){
				third.slice(i)= dveck2dtheta.slice(i)  + dveck2dthetadxATk21.slice(i)*trans(dk21.slice(i)); //third term, dvecd2k2dtheta_2
				
				arma:: mat eye1, eye2;
				
				d2vecdk2.slice(i) = 
				(kron(dfdxATk21.slice(i),eye1.eye(InfDS.Ntheta, InfDS.Ntheta))).t()* first.slice(i) + 
				kron(eye2.eye(InfDS.Nx, InfDS.Nx), dk21.slice(i))*(second.slice(i)) +
				third.slice(i);
			}
			//printf("execution 5.5\n");
			
			d2Xstar_t = dt(t)/2*(d2vecdk1 + d2vecdk2);
			d2Xtild = d2XtildPrev + d2Xstar_t;
			dXtildPrev = dXtild;
			d2XtildPrev = d2Xtild;
			
			
			
			//printf("t = %d\n", t);
			//dXtildPrev.slice(0).cols(0, InfDS.Nx-1).print("dXtildPrev"); 
			//correct here
		}
		
		
		//Rprintf("execution 6\n");
	
		//printf("t = %d\n", t);
		
		XtildPrev = Xtild_t;
		//printf("t = %d\n", t);
		//XtildPrev.print("XtildPrev");
		
		
		for (i = 0; i < InfDS.Nsubj; i++){
			
			//This needs to be fixed. No need to sum again in getScoreInfo
			indext = span_vec(0,InfDS.allT(i)-1,1);
			currentt = indext(InfDS.timeDiscrete(i)==tspan(t));
			
			
			if(i== 0 && t == 0){
				InfDS.timeDiscrete(0).t().print("*timeDiscrete(0)");
				printf("t = %d tspan(t) = %lf\n", t, tspan(t));
			}
			
			fullX.slice(t).col(i) = Xtild_t.col(i);
			
			for(int j = 0;  j < int(currentt.n_elem);j++){
				if( int(currentt(j)) == 1){
					if (getDxFlag ==1){
						InfDS.dXtildthetafAll(i).cols((j)*InfDS.Nx, (j+1)*InfDS.Nx - 1) = dXtild.slice(i);
						InfDS.dXtildthetafAll2(i).rows((j)*InfDS.Ntheta*InfDS.Nx, (j+1)*InfDS.Ntheta*InfDS.Nx - 1) = d2Xtild.slice(i);
						InfDS.dXstarAll(i).cols((j)*InfDS.Nx,(j+1)*InfDS.Nx - 1) = dXstar_t.slice(i);
						InfDS.dXstarAll2(i).rows((j)*InfDS.Ntheta*InfDS.Nx, (j+1)*InfDS.Ntheta*InfDS.Nx - 1) = d2Xstar_t.slice(i);
						
						/*// correct here
						if (i == 49){
							//currentt.print("currentt");
							printf("i = %d, j = %d t = %d\n",i,j,t);
							dXtild.slice(i).print("dXtild.slice(i)");
						}*/
					}
				}
				//printf("end of loop\n");
			}
			//tcount(i) = tcount(i) + 1;
			//printf("end of loop2\n");
			
		}
		//printf("end of loop3\n");
	}
	//Rprintf("execution 6.2\n");
	//correct here
	//InfDS.dXtildthetafAll(0).cols((0)*InfDS.Nx, (0+1)*InfDS.Nx - 1).print("InfDS.dXtildthetafAll(0) time 0");
	//InfDS.dXtildthetafAll(0).cols((199)*InfDS.Nx, (199+1)*InfDS.Nx - 1).print("InfDS.dXtildthetafAll(0) time 199");

	
	for (i = 0; i < InfDS.Nsubj; i++){
		for(int j = 0; j < InfDS.allT(i); j++){
			InfDS.Xtild.slice(j).col(i) = fullX.slice(int(InfDS.tobs(i)(j) - 1)).col(i);
		}
	}
	
	//Rprintf("execution 7\n");
	return InfDS;
}

void VDPMeas(arma::vec x, int Ny, int Nx, int NxState, arma::mat InfDS_Lambda, arma::mat mu, arma::mat *yPred, arma::mat *Jy){
  
  // return arma::fmat& yPred, arma::fmat& Jy
  //arma::mat yPred, Jy;
  //mexPrintf("\n\ncheck point VD\n\n");
  /*  main function starts from here */ 
  // Line 1: Lambda = zeros(Ny, Nx);
  arma::mat Lambda;	
  Lambda.set_size(Ny, Nx) ;
  Lambda.zeros();
  
  // Line 2: Lambda(1:Ny,1:NxState) = InfDS_Lambda;
  Lambda.submat(0, 0, Ny-1, NxState-1)= InfDS_Lambda;
  //Lambda.print("Lambda=");
  
  
  
  //Line 3: yPred = Lambda*x + kron(ones(1,size(x,2)),mu);
  arma::mat temp(1, x.n_cols);
  temp.ones();
  *yPred = Lambda * x + kron(temp,mu);
  //(*yPred).print("y_Pred=");
  
  // Line 4: Jy = reshape(kron(ones(1,size(x,2)), Lambda),Ny, Nx, size(x,2));
  *Jy = reshape(kron(temp, Lambda), Ny, Nx, x.n_cols);
  //(*Jy).print("Jy = \n");
  
  /*  main function ends from here */ 
  return;
  //return List::create(Named("yPred") = yPred, Named("Jy") = Jy);
  
}

void getScoreInfoY_tobs_opt(struct C_INFDS &InfDS, int stage, int iter, int freeIC, arma::mat &score, arma::mat &infoMat){
	
	arma::mat SIndex, SIndex2, LIndex, dlikdBetaAll, dlikdBetaAll2, dlikdmuAll, dlikdmuAll2, dlikdLambdaAll, dlikdLambdaAll2, dlikdMudLambda, dlikdLambdMu, dlikdMudBeta, dlikdBetadMu, dlikdLambdadBeta, dlikdBetadLambda, dlikdEpsilonAll, dlikdEpsilonAll2, dlikdbAll, dlikdbAll2, ivSigmab, Lambda, Zt, dXtildthetaf, dXtildthetaf2, vvv, aa, b, ivSigmae2, dlik, dlik2, Z, SIndex3, a;
	arma::cube X(0,0,0);
	int mulp2, mulp, isPar, t, T, curpos2, curpos1, startBeta, endBeta, startMu, endMu, startLamb, endLamb;
	
	
	
	int	i, ii, jj;
	int count_not_finite, count_finite;
	//unsigned int ui, uj;
	
	

	arma::mat indexY, indexYt, dlik_, dlik2_, InfDS_dmudparMu_, ivSigmae2_, Lambda_, Zt_, InfDS_dSigmaede_, InfDS_dSigmaede2_, SIndex2_, indexY2, InfDS_dLambdparLamb_, InfDS_Sigmae_;
	
	//printf("**** %d\n",InfDS.Ny);
	//span_vec(1, InfDS.Ny*InfDS.Ny, 1).t().print("?");
	
	SIndex = reshape(span_vec(1, InfDS.Ny*InfDS.Ny, 1).t(), InfDS.Ny, InfDS.Ny);
	//SIndex.print("SIndex");
	SIndex2 = reshape(span_vec(1, InfDS.Ny*InfDS.Ny*InfDS.Ny, 1).t(), InfDS.Ny, InfDS.Ny*InfDS.Ny);
	//SIndex2.print("SIndex2");
	LIndex = reshape(span_vec(1, InfDS.Nx*InfDS.Ny, 1).t(),InfDS.Ny,InfDS.Nx);
	//LIndex.print("LIndex");

	if (InfDS.Nbeta > 0){
		////printf("*1\n");
		dlikdBetaAll = arma::zeros<arma::mat>(InfDS.Nbeta, 1);
		////printf("*2\n");
		dlikdBetaAll2 = arma::zeros<arma::mat>(InfDS.Nbeta, InfDS.Nbeta);
		////printf("*3\n");
	}

	if (InfDS.Nmu > 0){
		dlikdmuAll = arma::zeros<arma::mat>(InfDS.Nmu, 1);
		dlikdmuAll2 = arma::zeros<arma::mat>(InfDS.Nmu, InfDS.Nmu);
	}

	if (InfDS.NLambda > 0){
		dlikdLambdaAll = arma::zeros<arma::mat>(InfDS.NLambda, 1);
		dlikdLambdaAll2 = arma::zeros<arma::mat>(InfDS.NLambda, InfDS.NLambda);
	}
	

	if (InfDS.Nmu > 0 && InfDS.NLambda > 0){
	  dlikdMudLambda = arma::zeros<arma::mat>(InfDS.Nmu, InfDS.NLambda); 
	  dlikdLambdMu = arma::zeros<arma::mat>(InfDS.NLambda, InfDS.Nmu);
	}
	
	
	if (InfDS.Nbeta > 0 && InfDS.NLambda > 0){
	  dlikdMudBeta = arma::zeros<arma::mat>(InfDS.Nmu, InfDS.Nbeta); 
	  dlikdBetadMu = arma::zeros<arma::mat>(InfDS.Nbeta, InfDS.Nmu);
	}
	
	if (InfDS.Nbeta > 0 && InfDS.NLambda > 0){
	  dlikdLambdadBeta = arma::zeros<arma::mat>(InfDS.NLambda,InfDS.Nbeta); 
	  dlikdBetadLambda = arma::zeros<arma::mat>(InfDS.Nbeta, InfDS.NLambda);
	}


	dlikdEpsilonAll = arma::zeros<arma::mat>(InfDS.Ny,1);
	dlikdEpsilonAll2 = arma::zeros<arma::mat>(InfDS.Ny,InfDS.Ny);


	if (InfDS.Nb > 0){
		dlikdbAll = arma::zeros<arma::mat>(InfDS.Nbpar, 1);
		dlikdbAll2 = arma::zeros<arma::mat>(InfDS.Nbpar, InfDS.Nbpar);
		ivSigmab = inv(InfDS.Sigmab);
	}

	//InfDS.par.print("InfDS.par");
	infoMat = arma::zeros<arma::mat>(InfDS.par.n_rows, InfDS.par.n_rows);
	ivSigmae2 = inv(diagmat(InfDS.Sigmae));


    if (InfDS.Nbetax > 0){
        Lambda = arma::zeros<arma::mat>(InfDS.Ny,InfDS.Nx);
        Lambda(span::all, span(0, InfDS.NxState-1)) = InfDS.Lambda(span::all, span::all);
	}
    else {
        Lambda = InfDS.Lambda;
    }
	
	//printf("*\n");
	isPar = (InfDS.Nx == InfDS.NxState) ? 0 : 1;
	//Rprintf("isPar %d enter getXtildIC3\n", isPar);
	if (InfDS.Nbeta > 0){
		InfDS = getXtildIC3(isPar, 1 ,freeIC, InfDS);
	}
	else{
		InfDS = getXtildIC3(isPar, 0 ,freeIC, InfDS);
	}
	//Rprintf("isPar %d leave getXtildIC3\n", isPar);
	mulp2 = 1;

    
	//InfDS.dXtildthetafAll(0)(span::all,span(0, 9)).print("InfDS.dXtildthetafAll(0)");
	for (i = 0; i < InfDS.Nsubj; i++){
		T = InfDS.allT(i); //InfDS.allT(i)=length(InfDS.timeDiscrete{i}) only at observed
		X = InfDS.Xtild(span::all, span(i,i), span(0,T-1));
		X = reshape(X, InfDS.Nx, T, 1);
		mulp = 1; 

		Z = InfDS.Y(i)(span::all,span(0, T-1)) - (kron(arma::ones<arma::mat>(1,T),InfDS.mu) + Lambda* X.slice(0));
		
		
		for (t = 0; t < T; t++){

			Zt = Z(span::all,span(t, t));
			
			//[HJ note] different original indexY contains indices that satisfies the if; now indexY contains Ny zeros and ones indicating the corresponding index satisfies the index or not
			indexYt = ones(1, InfDS.Ny);	
			//indexY = span_vec(0, InfDS.Ny-1, 1);
			count_finite = count_not_finite = 0;
			for(ii = 0; ii < InfDS.Ny; ii++){
				if( is_finite(InfDS.Y(i)(ii,t)) && is_finite(Zt(ii))){
					count_finite += 1;
				}
				else{
					count_not_finite += 1;
					indexYt(i) = 0;
				}
			}
			indexY.set_size(1, count_finite);
			jj = 0;
			for(ii = 0; ii < InfDS.Ny; ii++){
				if( is_finite(InfDS.Y(i)(ii,t)) && is_finite(Zt(ii))){
					indexY(jj) = ii + 1;	// store the index in C start from 0)
					jj += 1;
				}
			}
			 
			
			
			if(count_not_finite < InfDS.Ny){ // Not all indicators are missing
			
				//repeated used results
				InfDS_dmudparMu_ = colProjection(rowProjection(InfDS.dmudparMu, indexY),indexY);
				ivSigmae2_ = colProjection(rowProjection(ivSigmae2, indexY),indexY);
				Lambda_ = rowProjection(Lambda, indexY);
				Zt_ = rowProjection(Zt, indexY);
				InfDS_Sigmae_ = colProjection(rowProjection(InfDS.Sigmae, indexY), indexY);
				

				if (InfDS.Nbeta > 0){
					dXtildthetaf = InfDS.dXtildthetafAll(i);
					dXtildthetaf2 = InfDS.dXtildthetafAll2(i);
					
					/*
					if(i == 0 && t == 0){
						InfDS.H(span((i)*InfDS.Ntheta, ((i+1)*InfDS.Ntheta)-1), span::all).t().print("H 706");
						dXtildthetaf(span::all, span((t)*InfDS.Nx, (t+1)*InfDS.Nx-1)).print("dXtildthetaf");
						Lambda_.print("Lambda_");
						ivSigmae2_.print("ivSigmae2_");
						Zt_.print("Zt_");
					}
					*/
					
				   
				    
					dlik = InfDS.H(span((i)*InfDS.Ntheta, ((i+1)*InfDS.Ntheta)-1), span::all).t() * dXtildthetaf(span::all, span((t)*InfDS.Nx, (t+1)*InfDS.Nx-1)) * Lambda_.t() * ivSigmae2_ * Zt_;
					
					
					
					/*
					if(i == 0){
						printf("t = %d\n",t);
						dXtildthetaf(span::all, span((t)*InfDS.Nx, (t+1)*InfDS.Nx-1)).print("dXtildthetaf");
						dlik.print("dlik");
					}
					*/
					
					
					vvv=kron(eye(InfDS.Nx, InfDS.Nx),InfDS.H(span((i)*InfDS.Ntheta,(i+1)*InfDS.Ntheta-1),span::all).t());
					
					dlik2 = kron(Zt_.t() * ivSigmae2_ * Lambda_, eye(InfDS.Nbeta, InfDS.Nbeta)) * vvv * dXtildthetaf2(span((t)*InfDS.Ntheta*InfDS.Nx,(t+1)*InfDS.Ntheta*InfDS.Nx-1),span::all) * InfDS.H(span((i)*InfDS.Ntheta,((i+1)*InfDS.Ntheta-1)),span::all)-InfDS.H(span((i)*InfDS.Ntheta,(i+1)*InfDS.Ntheta-1),span::all).t()*dXtildthetaf(span::all,span((t)*InfDS.Nx, (t+1)*InfDS.Nx - 1))*Lambda_.t()*ivSigmae2_*Lambda_*dXtildthetaf(span::all,span((t)*InfDS.Nx, (t+1)*InfDS.Nx-1)).t()*InfDS.H(span((i)*InfDS.Ntheta,(i+1)*InfDS.Ntheta-1),span::all);
					
					//dlik_: replace the not finite numbers to 0 in dlik;
					dlik_ = dlik;
					dlik_.elem( find_nonfinite(dlik_) ).zeros(); 
					
					//dlik2_: replace the not finite numbers to 0 in dlik2;
					dlik2_ = dlik2;
					dlik2_.elem( find_nonfinite(dlik2_) ).zeros(); 
					
					dlikdBetaAll = dlikdBetaAll+ dlik_/mulp;
					dlikdBetaAll2 = dlikdBetaAll2+ dlik2_/mulp;
					
					//Rprintf("checkpoint 776\n");
					
					if (dlik.has_nan() || dlik2.has_nan()){
						printf("i = %d t = %d is Nan\n", i, t);
						
					}
					if (dlik.has_inf() || dlik2.has_inf()){
						printf("i = %d t = %d is inf\n", i, t);
						
					}
					
				}//End of Nbeta > 0 loop            
				
				aa =  reshape(span_vec(1, InfDS.Ny*InfDS.Ny, 1).t(),InfDS.Ny,InfDS.Ny); 
				SIndex3 = reshape(colProjection(rowProjection(aa, indexY),indexY),indexY.n_cols,1);


				
				if (InfDS.Nmu > 0){
					dlik = InfDS_dmudparMu_ * ivSigmae2_ * Zt_;
					
					//kron(Zt_.t()*ivSigmae2_,eye(indexY.n_cols, indexY.n_cols)).print("eval shoud be ? *3");
					//colProjection(InfDS.dmudparMu2,indexY).print("eval 22"); 
					//(InfDS_dmudparMu_ * ivSigmae2_* InfDS_dmudparMu_.t()).print("eval 33");
					dlik2 = kron(Zt_.t()*ivSigmae2_,eye(indexY.n_cols, indexY.n_cols))* colProjection(InfDS.dmudparMu2,indexY) - InfDS_dmudparMu_ * ivSigmae2_* InfDS_dmudparMu_.t();
					
					for(ii = 0; ii < InfDS.Ny; ii++){
						if(indexYt(ii)){
							dlikdmuAll(ii) = dlikdmuAll(ii) + dlik(ii)/mulp;	// check dlik is  1-by-n or n-by-1
							dlikdmuAll2(ii,ii) = dlikdmuAll2(ii,ii)+  dlik2(ii, ii)/mulp;							
						}
					
					}	
				}
				//Rprintf("checkpoint 807\n");
				
				a = reshape(colProjection(rowProjection(SIndex, indexY),indexY), 1, InfDS.Ny*InfDS.Ny);
				
				InfDS_dSigmaede_ = colProjection(rowProjection(InfDS.dSigmaede,indexY),a);
				
				SIndex2_ = reshape(colProjection(rowProjection(SIndex2, indexY),a), 1, InfDS.Ny*InfDS.Ny*InfDS.Ny);
				
				InfDS_dSigmaede2_ = colProjection(rowProjection(InfDS.dSigmaede2, SIndex2_),indexY);
			
				
				dlik =.5*colProjection(rowProjection(InfDS.dSigmaede, indexY),a) * (kron(ivSigmae2_, ivSigmae2_) * reshape(Zt_*Zt_.t() - InfDS_Sigmae_, InfDS.Ny*InfDS.Ny, 1));
				
				//Rprintf("checkpoint 820\n");
				
				dlik2 = .5*kron(reshape(ivSigmae2_*(Zt_*Zt_.t())* ivSigmae2_, 1, InfDS.Ny*InfDS.Ny), eye(InfDS.Ny, InfDS.Ny))*InfDS_dSigmaede2_ - InfDS_dSigmaede_*kron(ivSigmae2_,ivSigmae2_*(Zt_*Zt_.t())*ivSigmae2_)*InfDS_dSigmaede_.t() - .5*kron(reshape(ivSigmae2_, 1, InfDS.Ny*InfDS.Ny),eye(InfDS.Ny, InfDS.Ny))*InfDS_dSigmaede2_+ .5*InfDS_dSigmaede_*kron(ivSigmae2_,ivSigmae2_)*InfDS_dSigmaede_.t();
				
				
				/*
				if(i == 0){
					printf("t = %d\n",t);
					dlik.print("dlik");
					dlik2.print("dlik2");
				}
				*/
				
				for (ii = 0; ii < (int)indexY.n_cols; ii++){
					dlikdEpsilonAll(indexY(ii)-1) = dlikdEpsilonAll(indexY(ii)-1)+ dlik(indexY(ii)-1)/mulp;
					dlikdEpsilonAll2(indexY(ii)-1, indexY(ii)-1) = dlikdEpsilonAll2(indexY(ii)-1,indexY(ii)-1) + dlik2(indexY(ii)-1, indexY(ii)-1)/mulp;
				}

				
				if(InfDS.NLambda > 0){
					jj = 0;
					for(ii = 0; ii < (int)indexY.n_cols; ii++){
						if(indexY(ii) >1){
							jj += 1;
						}
					}
					indexY2.set_size(jj);
					jj = 0;
					for(ii = 0; ii < (int)indexY.n_cols; ii++){
						if(indexY(ii) >1){
							indexY2(jj) = indexY(ii);		
							jj += 1; 
						}
					}
					
					arma::mat X_ = X(span::all, span(t, t), span::all);
					a = reshape(rowProjection(LIndex,indexY), 1, InfDS.Ny*InfDS.Nx);
					InfDS_dLambdparLamb_ = colProjection(rowProjection(InfDS.dLambdparLamb, (indexY2-1).t()), a);
					
					dlik = (InfDS_dLambdparLamb_*reshape(ivSigmae2_*Zt_*X_.t(),InfDS.Ny*InfDS.Nx,1));
					dlik2 = (-InfDS_dLambdparLamb_*(kron(X_*X_.t(), ivSigmae2_))*InfDS_dLambdparLamb_.t());
					
					//ask symiin
					//Don't need this part because InfDS.dLambdaparLambd2 = 0. If not, add it
					//back to dlik2
					//Lambda = length(indexY(indexY>1));
					//a2 = reshape(LIndex2(indexY2-1,:),1,Ny*InfDS.Nx*NLambda);			//+kron(reshape(ivSigmae2_*Zt_*Xh_ekf(:,i,t)',1,Ny*InfDS.Nx),eye(NLambda))*InfDS.dLambdaparLambd2(a2,indexY2-1)
					
					for(ii = 0; ii < (int)indexY2.n_rows; ii ++){
						jj = indexY2(ii);
						dlikdLambdaAll(jj - 2) =  dlikdLambdaAll(jj-2) + dlik(jj-2)/mulp;
						dlikdLambdaAll2(jj - 2, jj - 2) =dlikdLambdaAll2(jj - 2, jj - 2)+ dlik2(jj - 2, jj - 2)/mulp;
					}
				}
			}//end of if then check if all Nys are missing at time t.
		}//end of t loop
		
		if (InfDS.Nb > 0){
			b = reshape(InfDS.b.row(i),InfDS.Nb,1);
			dlik = .5 * InfDS.dSigmabdb * (kron(ivSigmab, ivSigmab)*reshape(b*b.t()-InfDS.Sigmab, InfDS.Nb*InfDS.Nb, 1));

			dlik2 = .5 * kron(reshape(ivSigmab*b*b.t()*ivSigmab,1,InfDS.Nb*InfDS.Nb), eye(InfDS.Nbpar, InfDS.Nbpar)) * InfDS.dSigmabdb2 - InfDS.dSigmabdb * (kron(ivSigmab, ivSigmab*b*b.t()*ivSigmab)) * InfDS.dSigmabdb.t() - .5 * kron(reshape(ivSigmab,1,InfDS.Nb*InfDS.Nb), eye(InfDS.Nbpar, InfDS.Nbpar)) * InfDS.dSigmabdb2 + .5 * InfDS.dSigmabdb * kron(ivSigmab,ivSigmab) * InfDS.dSigmabdb.t();
			
			dlikdbAll = dlikdbAll+ dlik/mulp2;
			dlikdbAll2 = dlikdbAll2+ dlik2/mulp2;
			
		}    
	} //end of subject loop
	//Rprintf("end of subject loop\n");

	//score.set_size(0,0);
	curpos2 = 0;
	curpos1 = 0;
	if (InfDS.Nbeta > 0){
		//dlikdBetaAll.print("dlikdBetaAll");
		score = dlikdBetaAll;
		curpos1 = 1; 
		curpos2 = curpos2 + InfDS.Nbeta;
		infoMat(span(curpos1 - 1,curpos2-1), span(curpos1 - 1,curpos2-1)) = dlikdBetaAll2;
		startBeta = curpos1; 
		endBeta = curpos2;
	}
	//Rprintf("beta\n");
	
	if (InfDS.Nmu > 0){
		//dlikdmuAll.print("dlikdmuAll");
		score = join_cols(score, dlikdmuAll);
		//score.print("score");
		curpos1 = curpos2 + 1; 
		curpos2 = curpos2 + InfDS.Nmu;
		//printf("%d %d\n", curpos1, curpos2);
		//infoMat.print("size of infomat");
		infoMat(span(curpos1 - 1,curpos2-1), span(curpos1 - 1,curpos2-1)) = dlikdmuAll2;
		//infoMat.print("infoMat");
		startMu = curpos1; 
		endMu = curpos2;   
	}
	//Rprintf("mu\n");

	if (InfDS.NLambda > 0){
		//dlikdLambdaAll.print("dlikdLambdaAll");
		score = join_cols(score, dlikdLambdaAll);
		curpos1 = curpos2 + 1; 
		curpos2 = curpos2 + InfDS.NLambda;
		infoMat(span(curpos1 - 1, curpos2-1), span(curpos1 - 1, curpos2-1)) = dlikdLambdaAll2;
		startLamb = curpos1; 
		endLamb = curpos2;  
	}
	//Rprintf("LAMBDA\n");
	
	//Measurement error variances
	//dlikdEpsilonAll.print("dlikdEpsilonAll");
	score = join_cols(score, dlikdEpsilonAll);
	curpos1 = curpos2+1; 
	curpos2 = curpos2+InfDS.Ny;
	infoMat(span(curpos1 - 1,curpos2-1), span(curpos1 - 1,curpos2-1)) = dlikdEpsilonAll2;

	if (InfDS.Nb > 0){
		//dlikdbAll.print("dlikdbAll");
		score = join_cols(score, dlikdbAll);
		curpos1 = curpos2 + 1; 
		curpos2 = curpos2 + InfDS.Nbpar;
		infoMat(span(curpos1 - 1, curpos2-1), span(curpos1 - 1,curpos2-1)) = dlikdbAll2;
	}

	
	//Rprintf("b\n");
	if (InfDS.Nmu > 0 && InfDS.NLambda > 0){
		 infoMat(span(startMu-1, endMu-1), span(startLamb-1, endLamb-1)) = dlikdMudLambda;
		 infoMat(span(startLamb-1, endLamb-1), span(startMu-1, endMu-1)) = dlikdLambdMu;
	}
	
	//Rprintf("mu lamdba\n");

	if (InfDS.Nmu > 0 && InfDS.Nbeta > 0){
		 infoMat(span(startMu-1, endMu-1), span(startBeta-1, endBeta-1)) = dlikdMudBeta;
		 infoMat(span(startBeta-1, endBeta-1), span(startMu-1, endMu-1 )) = dlikdBetadMu;
	}
	
	//Rprintf("mu beta\n");

	if (InfDS.NLambda > 0 && InfDS.Nbeta){
		 infoMat(span(startLamb-1, endLamb-1),span(startBeta-1, endBeta-1)) = dlikdLambdadBeta;
		 infoMat(span(startBeta-1, endBeta-1),span(startLamb-1, endLamb-1)) = dlikdBetadLambda;
	}
	
	//score.print("score in function");
	//infoMat.print("infoMat in fucntion");
	//printf("Lamdba Nbeta\n");
	infoMat = -1*infoMat;
	return;
}



void saem(struct C_INFDS &InfDS, int &gmm, int &stage, int &redFlag, int &convFlag, int &noIncrease, int &stop, double &ssmin, double &ss, arma::mat &sgnTH, arma::vec &mscore, arma::mat &mscore2, arma::mat &minfoMat, arma::mat &Covscore){


	// variable declaration
	arma::vec thetam, dc2;
	arma::mat EI, ES, Iy, R, sy;
	double gainpara, gainparb, gmm1, t, gain;
	int flag;
	
	ss = -1; //Here ss is a multiplicative constant (t in Zhu et al.) to make sure that the information matrix is positive definite. 
	stop = 0; 


	thetam = reshape(InfDS.par,InfDS.par.n_elem,1);
	convFlag = 0;

	/*
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Gain constant for stage 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	*/
    gainpara = InfDS.gainpara;
    gainparb = InfDS.gainparb;
 
    
	if (stage == 2){ 
		gainpara = InfDS.gainpara1; 
		gainparb = InfDS.gainparb1;
	}

	if (stage < 2 && gmm <= 5)
		gmm1 = 1;
	else
		gmm1 = gmm;
	
	gain = gainparb/(pow(gmm1, gainpara) + gainparb - 1);  
	//Rprintf("gain = %lf gainpara= %lf gainparab= %lf\n", gain, gainpara, gainparb);


	redFlag = 0; 
	sy = InfDS.sy; 
	ES = InfDS.ES; 
	EI = InfDS.EI; 
	Iy = InfDS.Iy;

	/*
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%"M"-step
	%Here ss is a multiplicative constant (t in Zhu et al.) to make sure that the
	%information matrix is positive definite.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	*/	
	//printf("Line 1015: InfDS.ES(12,12) = %lf\n",InfDS.ES(11,11));
 
	InfDS.sy = InfDS.sy + gain*(mscore - InfDS.sy);
	InfDS.ES = InfDS.ES + gain*(mscore2 - InfDS.ES);
	InfDS.EI = InfDS.EI + gain*(minfoMat - InfDS.EI);
	t = 1.0;
	InfDS.Iy = -t*InfDS.ES + (InfDS.sy*InfDS.sy.t()) + InfDS.EI;
	//Rprintf("checkpoint 947\n");
	//printf("Line 1022: InfDS.ES(12,12) = %lf\n",InfDS.ES(11,11));

	
	flag = chol(R, InfDS.Iy);
	if(!flag){
		//printf("\nIy is not positive definite. ss = %5f\n",ss);
		Rprintf("Serious error!\n");
	}
	
	while (!flag && t >= 0.4){ //decomposition fails and t >= 0.4
		t = 0.5 * t;
		InfDS.Iy = -t*InfDS.ES + (InfDS.sy*InfDS.sy.t()) + InfDS.EI;
		flag = chol(R, InfDS.Iy);
	}
	//Rprintf("checkpoint enter 955\n");
 
	
	flag = chol(R, InfDS.Iy);
	if (!flag){
		redFlag=1;
		printf("redFlag = %d @ Line 1112\n", redFlag);
	}
	

	
	
	unsigned int i, j;
	arma::sp_mat Iy_s(InfDS.Iy.n_rows, InfDS.Iy.n_cols);
	for(i = 0; i < InfDS.Iy.n_rows;i++){
		for(j = 0; j < InfDS.Iy.n_cols;j++)
			Iy_s(i, j) = InfDS.Iy(i, j);
	}
	//Rprintf("checkpoint enter 969\n");	
 
	//Iy_s.print("Iy_s"); //the same
	//mscore.print("mscore"); //the same
	// If spsolve can not find solution, try solve. If solve also doesn't work, set is as zero
	if(!spsolve(dc2, Iy_s, mscore, "lapack")){
		if(!solve(dc2, mat(Iy_s), mscore)){
			dc2= zeros(mscore.n_elem);
		}
	}
	//dc2.print("dc2"); //the same

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if (!redFlag){
		thetam = thetam + gain * dc2;
	}

	//Rprintf("checkpoint enter 1090, Stage %d\n", stage);
	
	if (stage > 0){
		/*
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%Stopping rule for stage 1
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		*/
		//Rprintf("checkpoint 986\n");
		for(i = 0; i < dc2.n_elem; i++){
			if(abs(dc2(i)) < .0001)
				dc2(i) = 0;
		}
		if (stage==1){
			//Rprintf("sgnTH %d %d dc2 %d\n",sgnTH.n_rows, sgnTH.n_cols, dc2.n_elem);
			sgnTH(span::all, gmm - 1) = gain * sign(dc2); 
			//Rprintf("sgnTH %d %d dc2 %d\n",sgnTH.n_rows, sgnTH.n_cols, dc2.n_elem);
			if (gmm >= InfDS.KKO){    
				ss = norm(mean(sgnTH(span::all,span(gmm - InfDS.KKO, gmm - 1)),1));
				Rprintf("\nStage 1 Error tolerance at convergence = %6f\n", ss); 
	   
				if (ss <= InfDS.errtrol1 || gmm == InfDS.maxIterStage1){
					stage = 2; 
					gmm = 1;
				}
			}
		}
		//Rprintf("checkpoint enter 1003\n");
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		unsigned int q;
		for (q = 0; q < InfDS.par.n_elem; q++){
			if (!std::isnan(InfDS.lowBound(q))){
				if (thetam(q) < InfDS.lowBound(q)){
					thetam(q) = InfDS.par(q);
					redFlag = 1;
					Rprintf("redFlag = %d @ Line 1064\n", redFlag);
					InfDS.sy(q) = sy(q);
					InfDS.ES(q, q) = ES(q, q);
					InfDS.EI(q, q) = EI(q, q);
					InfDS.Iy (q, q) = Iy(q, q); 
				}
			}
			
			if (!std::isnan(InfDS.upBound(q))){
				if (thetam(q) > InfDS.upBound(q)){
					thetam(q) = InfDS.par(q);
					redFlag = 1;
					Rprintf("redFlag = %d @ Line 1064\n", redFlag);
					InfDS.sy(q) = sy(q);
					InfDS.ES(q, q) = ES(q, q);
					InfDS.EI(q, q) = EI(q, q);
					InfDS.Iy (q, q) = Iy(q, q); 
				}
			}
			
			if (std::isnan(thetam(q))){
				thetam(q) = InfDS.par(q);
				redFlag=1;
				Rprintf("redFlag = %d @ Line 1074\n", redFlag);
				InfDS.sy(q) = sy(q);
				InfDS.ES(q, q) = ES(q, q);
				InfDS.EI(q, q) = EI(q, q);
				InfDS.Iy(q, q) = Iy(q, q);  
			}
		}
		//Rprintf("checkpoint enter 1027\n");
		//printf("Line 1022: InfDS.ES(12,12) = %lf\n",InfDS.ES(11,11));

		//if ((~any(any(any(isnan(InfDS.Xtild))))==0 || redFlag==1) && InfDS.Nbeta > 0)
		if ((InfDS.Xtild.has_nan() || redFlag==1) && InfDS.Nbeta > 0){
			thetam(span(0, InfDS.Nbeta-1), 0) = InfDS.par(span(0, InfDS.Nbeta - 1), 0);
		}

		/*
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%Off-line averaging procedure
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		*/
		if (stage==2){ 
			//careful
			if ( max(max(InfDS.thetatild)) == 0 && min(min(InfDS.thetatild)) == 0){
				//Rprintf("checkpoint enter 1042\n");
				InfDS.sytild = InfDS.sy; 
				InfDS.EStild = InfDS.ES; 
				InfDS.EItild = InfDS.EI;
				InfDS.Iytild = InfDS.Iy;
				InfDS.thetatild = thetam;
				//Rprintf("checkpoint leave 1042\n");
			}
			else{
				//Rprintf("checkpoint enter 1051\n");
				InfDS.sytild = InfDS.sytild + (InfDS.sy-InfDS.sytild)/gmm;
				InfDS.EStild = InfDS.EStild + (InfDS.ES-InfDS.EStild)/gmm;
				InfDS.EItild = InfDS.EItild + (InfDS.EI-InfDS.EItild)/gmm;
				InfDS.Iytild = InfDS.Iytild + (InfDS.Iy-InfDS.Iytild)/gmm;
				InfDS.thetatild = InfDS.thetatild + (thetam-InfDS.thetatild)/gmm;
				//Rprintf("checkpoint enter 1051\n");
				
				//ss = abs(InfDS.sytild.t()*inv(InfDS.Iytild)*InfDS.sytild + trace(solve(InfDS.Iytild,Covscore))/gmm);
				arma::mat temp = InfDS.sytild.t()*inv(InfDS.Iytild)*InfDS.sytild;
				ss = abs(temp(0,0) + trace(solve(InfDS.Iytild,Covscore))/gmm);
				
				/*
				if (ss > ssmin)
				   noIncrease = noIncrease + 1;
				else{
				   ssmin = ss;
				   noIncrease = 0;
				}
				*/
				if (ss <= ssmin){
					ssmin = ss;
				}
				
				//Rprintf("checkpoint enter 1069\n");
	   
				if (ss < InfDS.errtrol || gmm == InfDS.MAXITER){
					if (gmm == InfDS.MAXITER)
						convFlag=0;
					else 
						convFlag=1;
					
					stop = 1;
				}
			}
		}
	}

	//printf("checkpoint enter 1083\n");

	gmm = gmm+1;
	InfDS.par = thetam;
	return;
}

void PropSigb(struct C_INFDS &InfDS){
	arma::mat temp, dfdbt, OMEGAb, iSigmae, Sigmae;
	arma::field<arma::mat> dXtilddthetafAll;
	
	int i, T, t;
	
	
	Sigmae = InfDS.Sigmae;
	iSigmae= eye(size(Sigmae)) / Sigmae;
	iSigmae.replace(datum::nan, 0);

	//%dxstardthetafAll = InfDS.dXstarAll; 
	dXtilddthetafAll = InfDS.dXtildthetafAll;
	OMEGAb.set_size(InfDS.Nb * InfDS.Nsubj, InfDS.Nb * InfDS.Nsubj);
	OMEGAb.zeros();

	for(i = 0; i < InfDS.Nsubj; i++){
    
		T = InfDS.timeDiscrete(i).n_elem;//%(InfDS.lens(i));
		temp.zeros(size(InfDS.Sigmab));
		
		//dt = InfDS.Deltat{i};%InfDS.fulldt{i};
		for(t = 0; t < T-1; t++){
            dfdbt = InfDS.Z.t()* dXtilddthetafAll(i)(span::all, span(t*InfDS.NxState, (t+1)*InfDS.NxState-1))*InfDS.Lambda.t();
            temp = temp + dfdbt * iSigmae * dfdbt.t(); 
			/*
			if(i == 0){
				printf("t = %d\n",t);
				dfdbt.print("dfdbt");
			}*/
		}//end of t loop
        	//inv(inv(InfDS.Sigmab) + temp).print("omegab");
		//printf("i = %d\n",i);
		//temp.print("temp");
		OMEGAb(span(i*InfDS.Nb,(i+1)*InfDS.Nb-1), span(i*InfDS.Nb,(i+1)*InfDS.Nb-1)) = inv(inv(InfDS.Sigmab) + temp);    
	}

	InfDS.OMEGAb = OMEGAb;
	
	return;

}




void drawbGeneral6_opt3(const int isPar, struct C_INFDS &InfDS, arma::mat &meanb, int yesMean, arma::mat upperb, arma::mat lowerb, int useMultN, arma::mat &tpOld, int freeIC, int isBlock1Only, int setScaleb, double &bAccept){
	arma::mat iSigmae, oldb, s, propden_new, propden_new1, tpNew1, OMEGAb, OMEGAi, R, MUb, normtmp, avgScalingb, cOMEGAb, bTemp, bdtmp, cOMEGAb0, tempi;
	arma::vec indexKept, cindexKept, filteredAvgScalingb, idx, tp, tp1, propden_old, tpNew;
	arma::ucolvec tpidx;
	arma::cube newb1, xtild;;
	double low1, high1, by1, scaleb;
	int Nb, q, Nkept, Ntmp, T, t;
	bool r;
	struct C_INFDS InfDS1;
	int i, j, jj;
	unsigned int ui, uj;
	
	arma::field<arma::mat> dXtilddthetafAll, dXtilddthetafAll2;
	
	//printf("execution point 0\n");
	if(isPar == 1)
		printf("[HJ] Should NOT happen! isPar = 1");
	
	
	iSigmae = inv(InfDS.Sigmae);
	//iSigmae.print("iSigmae");

	if (useMultN==1)
		InfDS.N = 4;
	else
		InfDS.N=1;


	if (InfDS.bAdaptParams.n_elem == 0){
		low1 = 1;
		high1 = 2;
		by1 = .1;
	}
	else{
		low1 = InfDS.bAdaptParams(0);
		high1 = InfDS.bAdaptParams(1);
		by1 = InfDS.bAdaptParams(2);
	}
	


	if (isBlock1Only)
		Nb = 1;
	else
		Nb = InfDS.Nb;


	
	MUb = reshape(InfDS.b(span::all, span(0, Nb -1)).t(), Nb * InfDS.Nsubj, 1);
	oldb = InfDS.b; 
	bTemp = zeros<arma::mat>(InfDS.Nsubj, Nb);
	s = zeros<arma::mat>(InfDS.Nsubj, InfDS.N);


	propden_new = zeros<arma::mat>(InfDS.Nsubj, 1);
	propden_new1 = zeros<arma::mat>(InfDS.Nsubj, InfDS.N);
	propden_old = zeros<arma::vec>(InfDS.Nsubj);
	tpNew1 = zeros<arma::mat>(InfDS.Nsubj, InfDS.N);
	newb1 = zeros<arma::cube>(InfDS.Nsubj, Nb, InfDS.N);



	xtild = InfDS.Xtild;
	dXtilddthetafAll = InfDS.dXtildthetafAll;
	dXtilddthetafAll2 = InfDS.dXtildthetafAll2;
	//printf("execution point 1\n");


	
	arma::mat rand_result(InfDS.Nsubj* InfDS.N,1);
	//old setting (uniformly random from low1-high1 (continious)
	rand_result = low1 + randu(InfDS.Nsubj* InfDS.N,1) * (high1-low1);
	// derandomized
	//for(i = 0; i < rand_result.n_rows; i++)
	//	rand_result(i, 0) = low1;
	
		
	//rand_result.print("rand_result");
	//printf("var %lf\n", var(rand_result));
	s(span::all,span::all) = reshape(rand_result, InfDS.Nsubj, InfDS.N); 
	for(ui = 0; ui < s.n_rows; ui++)
		for(uj = 0; uj < s.n_cols; uj++)
			s(ui, uj) = sqrt(s(ui, uj));

/////
	OMEGAb = InfDS.OMEGAb;
	scaleb = InfDS.scaleb;
	//printf("execution point 2\n");

	for (i = 0; i < InfDS.Nsubj; i++){
		//mexPrintf("i = %d\n", i);
		//OMEGAi = scaleb*OMEGAb(1+(i-1)*Nb:i*Nb,1+(i-1)*Nb:i*Nb);
		OMEGAi = scaleb*OMEGAb(span(i*Nb, (i+1)*Nb - 1), span(i*Nb, (i+1)*Nb - 1));
	        //printf("OMEGAi %d %d\n", OMEGAi.n_rows, OMEGAi.n_cols);
		//OMEGAi.print("OMEGAi");
	
		//If the decomposition fails:
		//chol(R,X) resets R (to be size zero) and returns a bool set to false (exception is not thrown)
		//[cOMEGAb0,r] = chol(OMEGAi);
		//OMEGAi.ones();
		//OMEGAi.print("OMEGAi (1)");
		r = chol(cOMEGAb0,OMEGAi);
		if (!r){
			//OMEGAi.print("OMEGAi");
			cOMEGAb0 = diagmat(diagvec(sqrt(OMEGAi)));
			OMEGAi.print("OMEGAi (1)");
			//mexPrintf("first chol\n");
			//cOMEGAb0.print("cOMEGAb0");
		}
		//cOMEGAb0.print("cOMEGAb0");
		
		//printf("cOMEGAb0 %d %d\n", cOMEGAb0.n_rows, cOMEGAb0.n_cols);
	      
    
		for (q = 0; q < InfDS.N; q++){
			cOMEGAb = s(i,q)*cOMEGAb0;
			/*
			if(q < 3 && i  < 10){
				//cOMEGAb.print("cOMEGAb");
				MUb(span(i*Nb, (i+1)*Nb - 1), span::all).print("MUb(span(i*Nb, (i+1)*Nb - 1), span::all)");
			}
			*/
			normtmp.set_size(Nb,1);
			normtmp.randn();
			//normtmp.ones(); // de randomized
			tempi = (MUb(span(i*Nb, (i+1)*Nb - 1), span::all) + cOMEGAb*normtmp).t();
			//printf("q %d  tempi %d %d \n", q, tempi.n_rows, tempi.n_cols);

			//tempi = (MUb((1+(i-1)*Nb):(i*Nb)) + cOMEGAb*normtmp)';
			//need to write a function: round2 [done test]
			for (uj = 0; uj < tempi.n_elem; uj++)
				tempi(uj) = round2(tempi(uj), 10e-4);
			//tempi.print("tempi");
			
			//tempi.print("tempi");
			//if(i  == 10){
			//	tempi.print("tempi");
			//	mexPrintf("min %lf max %lf\n", max(tempi), min(tempi));
			//}
			//printf("* newb1 %d %d %d\n", newb1.n_rows, newb1.n_cols, newb1.n_slices);
			//printf("* MUb %d %d\n", MUb.n_rows, MUb.n_cols);
			//printf("* InfDS.N %d\n", InfDS.N);
			for (jj = 0; jj < Nb; jj++){
				//printf("i %d jj %d q %d\n", i, jj,q);
				if ( (max(tempi(jj)) < upperb(jj) && min(tempi(jj)) > lowerb(jj))|| std::isnan(upperb(jj))){	
				//if ( (max(tempi(jj)) < upperb(jj) + ERROR && min(tempi(jj)) > lowerb(jj) - ERROR )|| 
				//	  std::isnan(upperb(jj)) ){
					//printf("tempi");
					newb1(i,jj,q) = tempi(jj);
					//printf(" **\n");
				}
				else{
					//printf("MUb");
					//oldline
					//newb1(i,jj,q) = MUb((jj+(i-1)*Nb));
					newb1(i,jj,q) = MUb((jj+(i)*Nb));
					//printf(" **\n");
					//newb1(i,jj,q) = round2(MUb((jj+(i-1)*Nb)), 10e-4);
				}
				//tempi.print("tempi(jj)");
				
				
			}
			//mexPrintf("loop end\n");
			//propden_new1(i, q) = -0.5* sum(sum(normtmp % normtmp, 1), 0) - sum(sum( log(diagvec(cOMEGAb)), 1), 0);
			propden_new1(i, q) = -0.5* sum(sum(normtmp % normtmp, 1))- sum(sum( log(diagvec(cOMEGAb)), 1));
			//printf("propden %d %d %lf\n", i, q, propden_new1(i, q));
		}
	}
	//newb1(span(0,2),span(0,2), span(0,0)).print("newb1");
	//propden_new1(span(0, 9), span(0,0)).t().print("propden_new1");
	//printf("execution point 3\n");
	
	//InfDS.dXtildthetafAll(0)(span::all, span(0, 9)).print("zero");
	//InfDS.dXtildthetafAll(0)(span::all, span(3*InfDS.NxState, (3 + 1)*InfDS.NxState -1)).print("before");
	for (q = 0; q < InfDS.N; q++){   
		InfDS1 = InfDS;
		InfDS1.N = 1;
		InfDS1.b(span::all, span(0,Nb-1)) = reshape(newb1.slice(q)(span::all,span(0, Nb-1)),InfDS.Nsubj, Nb);
		//mexPrintf("execution point 3.1\n");
		// to be examined: isPar 0: no need
		//C_INFDS getXtildIC3(const int isPar, const int getDxFlag, const int freeIC, struct C_INFDS &InfDS)
		InfDS1 = getXtildIC3(0, 0 ,freeIC,InfDS1);
		//mexPrintf("execution point 3.2\n");
		//arma::vec ekfContinuous10(int Nsubj, const int N, const int Ny, const int Nx, const int Nb, const int NxState,const arma::mat Lambda, int totalT, const arma::mat Sigmae, const arma::mat Sigmab, const arma::mat mu, const arma::mat b, const arma::mat allT, const arma::cube Xtild, arma::cube Y)
		tpNew = ekfContinuous10(InfDS1.Nsubj, InfDS1.N, InfDS1.Ny, InfDS1.Nx, InfDS1.Nb, InfDS1.NxState, InfDS1.Lambda, InfDS1.totalT, InfDS1.Sigmae, InfDS1.Sigmab, InfDS1.mu, InfDS1.b, InfDS1.allT, InfDS1.Xtild, InfDS1.Y);		
		tpNew1(span::all, q) = tpNew;	
	}
	//mexPrintf("InfDS1.b = %lf\n",InfDS1.b );
	//InfDS1.b.submat(0,0,2,2).print("InfDS1.b");
	//InfDS.dXtildthetafAll(0)(span::all, span(3*InfDS.NxState, (3 + 1)*InfDS.NxState -1)).print("after");
	//printf("execution point 4 %d %d\n", tpNew1.n_rows, tpNew1.n_cols);

	
	tpNew = max(tpNew1,1);
	

	
	//arma::uvec max_index_temp = index_max(tpNew1,1);
	//tpidx = ones<arma::mat>(tpNew.n_elem, 1);
	tpidx = index_max(tpNew1, 1);
	//tpidx.print("tpidx");
	
	
	for (i = 0; i < InfDS.Nsubj; i++){
		//mexPrintf("i = %d\n",i);
		bTemp.row(i) = reshape(newb1.slice(tpidx(i))(span(i,i), span(0, Nb-1)), 1, Nb);
		propden_new(i) = propden_new1(i,tpidx(i));
		//bTemp.row(i).print("bTemp");
	}
	//bTemp.print("bTemp");
	//mean(bTemp).print("mean btemp");
	//var(bTemp).print("var btemp");
	

	InfDS.b(span::all,span(0,Nb-1)) = bTemp;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Ntmp = InfDS.N;
	InfDS.N = 1;
	//InfDS = getXtildIC3(InfDS, 1, freeIC);
	// to be examined: isPar
	//mexPrintf("enter getXtildIC3\n");
	InfDS = getXtildIC3(0, 1 ,freeIC, InfDS);
	//mexPrintf("leave getXtildIC3\n");
	InfDS.N = Ntmp;
	//InfDS.dXtildthetafAll(0)(span::all, span(3*InfDS.NxState, (3 + 1)*InfDS.NxState -1)).print("after2");

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//Proposal for oldb
	//printf("execution point 5\n");

	for (i = 0; i < InfDS.Nsubj; i++){
		T = InfDS.timeDiscrete(i).n_elem;
		arma::mat temp = zeros<arma::mat>(Nb,Nb);
		//iSigmae.print("iSigmae");
		for (t = 0; t < T-1; t++){
		                                        
			arma::mat dfdbt;
			//dfdbt = InfDS.Z(:,1:Nb)'*InfDS.dXtildthetafAll{i}(:,(1+(t-1)*InfDS.NxState):(t*InfDS.NxState))*InfDS.Lambda';
			dfdbt = InfDS.Z(span::all,span(0,Nb-1)).t()*InfDS.dXtildthetafAll(i)(span::all, span(t*InfDS.NxState, (t + 1)*InfDS.NxState -1))*InfDS.Lambda.t();
			//InfDS.dXtildthetafAll(i)(span::all, span(t*InfDS.NxState, (t + 1)*InfDS.NxState -1)).print("dXtildthetafall");
			//dfdbt.print("dfdbt");
			temp = temp + dfdbt * iSigmae * dfdbt.t();
					
		}
		//mexPrintf("Z %lf dX %lf Lambda %lf\n",InfDS.Z(1,1), InfDS.dXtildthetafAll(1)(1,1),InfDS.Lambda(1,1));
			
		//OMEGAb(1+(i-1)*Nb:i*Nb,1+(i-1)*Nb:i*Nb) = inv(inv(InfDS.Sigmab(1:Nb,1:Nb)) + temp);    
		OMEGAb(span(i*Nb,(i+1)*Nb-1),span(i*Nb,(i+1)*Nb-1)) = inv(inv(InfDS.Sigmab(span(0,Nb-1), span(0,Nb-1))) + temp);  
		//OMEGAb(span(i*Nb,(i+1)*Nb-1),span(i*Nb,(i+1)*Nb-1)).print("OMEGAb");
	}

	//printf("execution point 6\n");
	
	bdtmp = oldb(span::all,span(0,Nb-1))-InfDS.b(span::all,span(0,Nb-1));
	//bdtmp.print("bdtmp");
	avgScalingb = zeros<arma::mat>(InfDS.Nsubj,1);

	for (i = 0; i < InfDS.Nsubj; i++){    
		//mexPrintf("i = %d\n",i);
		OMEGAi = s(i,tpidx(i)) * InfDS.scaleb * (OMEGAb(span(i*Nb,(i+1)*Nb-1),span(i*Nb,(i+1)*Nb-1))); 
		//OMEGAi.print("OMEGAi");
		//printf("s(i,tpidx(i)) = %lf InfDS.scaleb = %lf\n", s(i,tpidx(i)), InfDS.scaleb);
		//scanf("%d");
		avgScalingb(i) =s(i,tpidx(i));
		
		//[cOMEGAb,r] = chol(OMEGAi);
		r = chol(cOMEGAb,OMEGAi);
		//OMEGAi.print("OMEGAi (2)");
		if (!r){
			//OMEGAi.print("OMEGAi");
			cOMEGAb = diagmat(diagvec(sqrt(OMEGAi)));
			printf("second chol\n");
			OMEGAi.print("OMEGAi (2)");
			//cOMEGAb.print("cOMEGAb");
		}
		//cOMEGAb.print("cOMEGAb");
		
		//printf("inv error happens here:");
		normtmp = reshape(bdtmp(span(i,i),span::all),1,Nb) * inv(cOMEGAb); 
		//printf("exit");
		//normtmp.print("normtmp");		
		// to be check
		propden_old(i) = -0.5* sum(sum(normtmp % normtmp, 1)) - sum(sum(log(diagvec(cOMEGAb)), 1));
		//if (i <= 5)
		//	mexPrintf("propden_old(%d) = %lf\n",i, propden_old(i));

	}
	
	//propden_old(span(0,4)).print("propden_old");
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//printf("execution point 7\n");

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
	//Evaluate MH ratio and which set to keep%
	//totalb = totalb+InfDS.Nsubj;
	arma::vec index = span_vec(1, InfDS.Nsubj, 1);
	tp = 0 + randu(InfDS.Nsubj, 1) * (1-0); //sample from unif(0,1) (should be recovered)
	//tp = 0 + ones(InfDS.Nsubj,1) *(1-0);	//de randomized *****
	tp1.set_size(tp.n_elem);
	for (i = 0; i < (int)tp1.n_elem; i++){
		tp1(i)= fmin(1.0, exp(tpNew(i) + propden_old(i) - tpOld(i) - propden_new(i)));
	}
	
	
	
	
	//to be converted: the following three lines
	//idx.zeros(size(tp));
	//idx = (tp(i)<=tp1(i));//&&~isnan(tpNew); //(& : and)
	//idx.print("idx");
	indexKept.set_size(size(tpNew));
	indexKept.zeros();
	cindexKept.set_size(size(tpNew));
	cindexKept.zeros();
	Nkept = 0;
	for (i = 0; i < (int)tp.n_elem; i++){
		//mexPrintf("i %d %lf %lf\n",i, tp(i), tp1(i));
		//if (tp(i) <= tp1(i)+ERROR && tpNew( span(i, i) ).is_finite()){
		if (tp(i) <= tp1(i) && tpNew( span(i, i) ).is_finite()){
			indexKept(i) = 1;
			Nkept++;
		}
		else
			cindexKept(i) = 1;
	}


	bAccept = double(Nkept)/InfDS.Nsubj;
	

	InfDS.bacc = InfDS.bacc + indexKept;
	//InfDS.bacc.t().print("bacc");
	//printf("execution point 8 Nkept = %d\n", Nkept);

	if (Nkept  > 0){    
		//tpOld(indexKept) = tpNew(indexKept);
		for (i = 0; i < (int)indexKept.n_elem; i++){
			if(indexKept(i)){
				//mexPrintf("indexKept %d\n", i + 1);
				tpOld(i) = tpNew(i);
			}
		}
		if (setScaleb==1){
			//avgScalingb(indexKept);
			filteredAvgScalingb.set_size(Nkept);
			i = 0; 
			j = 0;
			while(i < Nkept && j < (int)indexKept.n_elem){
				if(indexKept(i)){
					filteredAvgScalingb(i) = avgScalingb(j);
					i++;
				}
				j++;
			}
			
			//to be convert: quantile
			InfDS.scaleb = InfDS.scaleb*quantile(filteredAvgScalingb, 1 - InfDS.setAccept);
			Rprintf("quantile(filteredAvgScalingb, 1 - InfDS.setAccept) = %lf\nbAccept = %lf\n", quantile(filteredAvgScalingb, 1 - InfDS.setAccept), bAccept);
			//%mean(avgScalingb);
		}
		//for i = indexKept
		for (i = 0; i < (int)indexKept.n_elem; i++){
			if(indexKept(i)){
				//i = indexKept(i);
				InfDS.OMEGAb(span(i*Nb,(i+1)*Nb-1), span(i*Nb,(i+1)*Nb-1)) = OMEGAb(span(i*Nb, (i+1)*Nb-1),span(i*Nb, (i+1)*Nb-1));
			}
		}
		//end
	}
	if (Nkept < InfDS.Nsubj){
		//InfDS.b(cindexKept,span(0,Nb-1)) = oldb(cindexKept,span(0,Nb-1));    
		//InfDS.Xtild(span::all,cindexKept,span::all) = xtild(span::all,cindexKept,span:all);
		for (i = 0; i < (int)cindexKept.n_elem; i++){
			//i = cindexKept(i);
			if(cindexKept(i)){
				InfDS.b(i,span(0,Nb-1)) = oldb(i,span(0,Nb-1));
				InfDS.Xtild(span::all,span(i,i),span::all) = xtild(span::all,span(i,i),span::all);
			}
		}
	
		//for i = cindexKept
		for (i = 0; i < (int)cindexKept.n_elem; i++){
			//i = cindexKept(i);
			if(cindexKept(i)){
				InfDS.dXtildthetafAll(i) = dXtilddthetafAll(i);
				InfDS.dXtildthetafAll2(i) = dXtilddthetafAll2(i);  
			}
		}
			   
		//end 
	}
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//printf("execution point 9\n");
	
	//printf("meanb %d*%d, InfDS.b %d*%d\n", meanb.n_rows, meanb.n_cols, InfDS.b.n_rows, InfDS.b.n_cols);
	
	if (yesMean ==1){
		meanb = meanb + InfDS.b;
	}

	if (useMultN==1){
		InfDS.N = 1;
	}
	
	/*Debug Printing*/
	//InfDS.b.submat(0,0,2,2).print("InfDS.b");
	//mexPrintf("low1=%lf high1=%lf by1=%lf\nisBlock1Only=%d InfDS.scaleb=%lf\nbAccept=%lf\n", low1, high1, by1,isBlock1Only, InfDS.scaleb,bAccept);
	
	//printf("execution point 10 end\n");

	return;
}




