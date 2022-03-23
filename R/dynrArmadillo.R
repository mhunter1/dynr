#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2017-10-30 08:50:19
# Filename: dynrArmadillo.R
# Purpose: Replicate much of the dynrRecipe methods for writing C code but for
#  armadillo code instead.
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Define Generic

setGeneric("writeArmadilloCode", 
	function(object, covariates, param.names=list(), dmudparMu, dmudparMu2, dLambdparLamb, dLambdparLamb2, dSigmaede, dSigmaede2, dSigmabdb, dSigmabdb2, Sigmab, known.vars, show=TRUE) { 
		return(standardGeneric("writeArmadilloCode")) 
	})


#------------------------------------------------------------------------------
# Define method for dynrMeasurement class

setMethod("writeArmadilloCode", "dynrMeasurement",
	function(object, covariates, param.names, dmudparMu, dmudparMu2, dLambdparLamb, dLambdparLamb2){
		#browser()
		ret <- ""
		
		

		#mu
		if(length(object$params.int) > 0){
			ret = paste0(ret, "\n\tif (InfDS.Nmu > 0){\n\t\tstartM = InfDS.Nbeta+1;\n")
			index = nrow(object$params.int[[1]]) - 2
			if(index < 0)
				ret = paste0(ret, "\t\tInfDS.mu = InfDS.par(span(startM - 1, startM -", -index, "), span::all);\n")
			else
				ret = paste0(ret, "\t\tInfDS.mu = InfDS.par(span(startM - 1, startM +", index, "), span::all);\n")
			ret = paste0(ret, "\t\tInfDS.dmudparMu = arma::zeros<arma::mat>(InfDS.Ny, InfDS.Ny); \n\t\tInfDS.dmudparMu2 = arma::zeros<arma::mat>(InfDS.Ny*InfDS.Ny, InfDS.Ny);\n")
			
			if(nrow(dmudparMu) > 0 && ncol(dmudparMu) > 0){
				for(i in 1: nrow(dmudparMu)){
					for(j in 1: ncol(dmudparMu)){
						if(dmudparMu[i,j] != 0)
							ret = paste0(ret, "\t\tInfDS.dmudparMu(", i-1, ", ", j-1, ") =", dmudparMu[i, j], ';\n')
					}
				}
			}
			if(nrow(dmudparMu2) > 0 && ncol(dmudparMu2) > 0){
				for(i in 1: nrow(dmudparMu2)){
					for(j in 1: ncol(dmudparMu2)){
						if(dmudparMu2[i,j] != 0)
							ret = paste0(ret, "\t\tInfDS.dmudparMu2(", i-1, ", ", j-1, ") =", dmudparMu2[i, j], ';\n')
					}
				}
			}
			ret = paste(ret, "\n\t}\n")
		}
		
		
		#lambda
		if(length(object$values.load) > 0){
			ret = paste(ret, "\tif (InfDS.NLambda > 0){\n\t\tstartL = InfDS.Nbeta+ InfDS.Nmu+1;\n\t\tInfDS.Lambda=\"")
			Lambda = object$values.load[[1]]
			if(nrow(Lambda) > 0 && ncol(Lambda) > 0){
				for(i in 1: nrow(Lambda)){
					for(j in 1: ncol(Lambda)){
						if(j == 1)
							ret = paste0(ret, Lambda[i, j])
						else
							ret = paste0(ret, ',', Lambda[i, j])
					}
					if(i != nrow(Lambda))
						ret = paste0(ret, ';')
				}
				ret = paste0(ret, '\";\n')
				
				Lambda = object$params.load[[1]]
				startL = min(Lambda[Lambda>0]) + 1
			
				for(i in 1: nrow(Lambda)){
					for(j in 1: ncol(Lambda)){
						if(Lambda[i, j] > 0){
							#print(param.names[Lambda[i, j]])
							if(Lambda[i, j] - startL < 0)
								ret = paste0(ret, "\t\tInfDS.Lambda(", i-1, ", ", j-1, ") = InfDS.par(startL - ", -(Lambda[i, j] - startL), ');\n')
							else
								ret = paste0(ret, "\t\tInfDS.Lambda(", i-1, ", ", j-1, ") = InfDS.par(startL + ", Lambda[i, j] - startL, ');\n')
						}
					}
					
				}
			}
			
			ret = paste(ret, "\t\tInfDS.dLambdparLamb = arma::zeros<arma::mat>(InfDS.NLambda, InfDS.Ny*InfDS.Nx);\n")
			if(nrow(dLambdparLamb) > 0 && ncol(dLambdparLamb) > 0){
				for(i in 1: nrow(dLambdparLamb)){
					for(j in 1: ncol(dLambdparLamb)){
						if(dLambdparLamb[i,j] != 0)
							ret = paste0(ret, "\t\tInfDS.dLambdparLamb(", i-1, ", ", j-1, ") =", dLambdparLamb[i, j], ';\n')
					}
				}
			}
			
			ret = paste(ret, "\t\tInfDS.dLambdparLamb2 = arma::zeros<arma::mat>(InfDS.Nx*InfDS.Ny*InfDS.NLambda, InfDS.NLambda);\n")
			if(nrow(dLambdparLamb2) > 0 && ncol(dLambdparLamb2) > 0){
				for(i in 1: nrow(dLambdparLamb2)){
					for(j in 1: ncol(dLambdparLamb2)){
						if(dLambdparLamb2[i,j] != 0)
							ret = paste0(ret, "\t\tInfDS.dLambdparLamb2(", i-1, ", ", j-1, ") =", dmudparMu2[i, j], ';\n')
					}
				}
			}
			
			ret = paste(ret, "\t}")
		}
		
		object@c.string <- ret
		return(object)
	}
)


#------------------------------------------------------------------------------
# Define method for dynrTrans class

setMethod("writeArmadilloCode", "dynrTrans",
	function(object, covariates){
		ret <- ""
		object@c.string <- ret
		return(object)
	}
)


#------------------------------------------------------------------------------
# Define method for dynrMatrixDynamics class

setMethod("writeArmadilloCode", "dynrDynamicsMatrix",
	function(object, covariates){
		ret <- ""
		object@c.string <- ret
		return(object)
	}
)


#------------------------------------------------------------------------------
# Define method for dynrFormulaDynamics class

setMethod("writeArmadilloCode", "dynrDynamicsFormula",
	#function(object, covariates, param.names, dmudparMu, dmudparMu2, dLambdparLamb, dLambdparLamb2,dSigmaede, dSigmaede2, dSigmabdb, dSigmabdb2){
	function(object, covariates, param.names, dSigmabdb, dSigmabdb2, Sigmab, known.vars){

		formula <- object$formulaOriginal
		formula2 <- object$formula2
		#transfer them as list
		jacob <- list(vectorizeMatrix(object$jacobianOriginal[[1]], byrow= FALSE))
		dfdtheta<- list(vectorizeMatrix(object$dfdtheta[[1]],byrow=FALSE))
		dfdx2<- list(vectorizeMatrix(object$dfdx2[[1]], byrow=TRUE))
		dfdxdtheta<- list(vectorizeMatrix(object$dfdxdtheta[[1]], byrow=TRUE))
		dfdthetadx<- list(vectorizeMatrix(object$dfdthetadx[[1]], byrow=TRUE))
		dfdtheta2<- list(vectorizeMatrix(object$dfdtheta2[[1]], byrow=TRUE))
		state.names <- object$state.names
		theta.names <- object$theta.names
		#intercept.names <- object$intercept.names
		random.names <- object$random.names
		theta.formula <-object$theta.formula
		nregime=length(formula)
		n=sapply(formula,length) #NxState+Nbeta (Nx)
		n.theta = length(object$theta.names)
		n.state = length(state.names) #NxState
		covariate.names <- covariates
		beta.names<- object$beta.names
		
		
		#Parse dyfun
		fml=lapply(formula,processFormula)
		lhs=lapply(fml,function(x){lapply(x,"[[",1)})
		rhs=lapply(fml,function(x){lapply(x,"[[",2)})
		
		
		
		#Parse jacobian (dfdx)
		fmlj=lapply(jacob,processFormula)
		row=lapply(fmlj,function(x){lapply(x,"[[",1)})
		col=lapply(fmlj,function(x){lapply(x,"[[",2)})
		rhsj=lapply(fmlj,function(x){lapply(x,"[[",3)})
		
		#Parse dfdtheta
		fmlp=lapply(dfdtheta,processFormula)
		row=lapply(fmlp,function(x){lapply(x,"[[",1)})
		col=lapply(fmlp,function(x){lapply(x,"[[",2)})
		rhsp=lapply(fmlp,function(x){lapply(x,"[[",3)})
		
		#Parse dfdx2
		fmlx2=lapply(dfdx2,processFormula)
		row=lapply(fmlx2,function(x){lapply(x,"[[",1)})
		col=lapply(fmlx2,function(x){lapply(x,"[[",2)})
		rhsx2=lapply(fmlx2,function(x){lapply(x,"[[",4)})
		
		#Parse dfdxdtheta
		fmlxp=lapply(dfdxdtheta,processFormula)
		row=lapply(fmlxp,function(x){lapply(x,"[[",1)})
		col=lapply(fmlxp,function(x){lapply(x,"[[",2)})
		rhsxp=lapply(fmlxp,function(x){lapply(x,"[[",4)})
		
		#Parse dfdthetadx
		fmlpx=lapply(dfdthetadx,processFormula)
		row=lapply(fmlpx,function(x){lapply(x,"[[",1)})
		col=lapply(fmlpx,function(x){lapply(x,"[[",2)})
		rhspx=lapply(fmlpx,function(x){lapply(x,"[[",4)})
		
		#Parse dfdtheta2
		fmlp2=lapply(dfdtheta2,processFormula)
		row=lapply(fmlp2,function(x){lapply(x,"[[",1)})
		col=lapply(fmlp2,function(x){lapply(x,"[[",2)})
		rhsp2=lapply(fmlp2,function(x){lapply(x,"[[",4)})
		
		#browser()
		#Parse theta.formula
		fmlt = processFormula(theta.formula)
		for (i in 1:length(theta.formula)){
			lhst = fmlt[[i]][1]
			rhst = fmlt[[i]][2]
			
			
			# Replace theta_i in formula with thetaf
			# - thetaf is calculated by calculateTheta()
			# - the variable name of thetaf is from the LHS of theta.formula
			# - in jacobian (dfdx), LHS of theta.formula is already replaced by RHS of theta.formula in rhsj (to get correct differentiation), thus, we replace the RHS of theta formula by thetaf
			rhs <- list(lapply(rhs[[1]], function(x){gsub(paste0(lhst),paste0("thetaf(",i-1,",s)"),x, fixed = TRUE)}))
			rhsj <- list(lapply(rhsj[[1]], function(x){gsub(paste0(lhst),paste0("thetaf(",i-1,",s)"),x, fixed = TRUE)}))
			rhsp <- list(lapply(rhsp[[1]], function(x){gsub(paste0(lhst),paste0("thetaf(",i-1,",s)"),x, fixed = TRUE)}))
			rhsx2 <- list(lapply(rhsx2[[1]], function(x){gsub(paste0(lhst),paste0("thetaf(",i-1,",s)"),x, fixed = TRUE)}))
			rhsxp <- list(lapply(rhsxp[[1]], function(x){gsub(paste0(lhst),paste0("thetaf(",i-1,",s)"),x, fixed = TRUE)}))
			rhspx <- list(lapply(rhspx[[1]], function(x){gsub(paste0(lhst),paste0("thetaf(",i-1,",s)"),x, fixed = TRUE)}))
			rhsp2 <- list(lapply(rhsp2[[1]], function(x){gsub(paste0(lhst),paste0("thetaf(",i-1,",s)"),x, fixed = TRUE)}))
		}
		
		#browser()
		#Replace variables (to be estimated) in formula and differentiations with InfDS.par
		# - Hui-Ju: Here, I assume that variables in InfDS.par follows the order in model@param.names, 
		# - Need to examine whether it works well for OSC model
		# - param.names: model@param.names
		# - replace it with the descreasing order of variable name length
		#print(param.names)
		index <- 1:length(param.names)
		repalce_order <- param.names[order(nchar(param.names), param.names, decreasing= TRUE)]
		for (i in 1:length(param.names)){
			var.name <- repalce_order[i]
			pattern <- var.name
			ind <- index[param.names == var.name]
			rhs  <- list(lapply(rhs[[1]], function(x){gsub(pattern, paste0("InfDS.par(",ind-1,",s)"),x, fixed = TRUE)}))
			rhsj <- list(lapply(rhsj[[1]], function(x){gsub(pattern, paste0("InfDS.par(",ind-1,",s)"),x, fixed = TRUE)}))
			rhsp <- list(lapply(rhsp[[1]], function(x){gsub(pattern, paste0("InfDS.par(",ind-1,",s)"),x, fixed = TRUE)}))
			rhsx2 <- list(lapply(rhsx2[[1]], function(x){gsub(pattern, paste0("InfDS.par(",ind-1,",s)"),x, fixed = TRUE)}))
			rhsxp <- list(lapply(rhsxp[[1]], function(x){gsub(pattern, paste0("InfDS.par(",ind-1,",s)"),x, fixed = TRUE)}))
			rhspx <- list(lapply(rhspx[[1]], function(x){gsub(pattern, paste0("InfDS.par(",ind-1,",s)"),x, fixed = TRUE)}))
			rhsp2 <- list(lapply(rhsp2[[1]], function(x){gsub(pattern, paste0("InfDS.par(",ind-1,",s)"),x, fixed = TRUE)}))
		}
		
		# Replace the covariate to corresponding variables in SAEM (i.e., InfDS.U1)
        for (i in 1:length(covariate.names)){
            selected <- covariate.names[i]
            #get <- paste0("covariate(s,", which(covariate.names == selected)-1,")")
			get <- paste0("InfDS.U1(s,", which(covariate.names == selected)-1,")")
            rhs <- lapply(gsub(paste0("\\<",selected,"\\>"), get, rhs), function(x){eval(parse(text=x))})
            rhsj <- lapply(gsub(paste0("\\<",selected,"\\>"), get, rhsj), function(x){eval(parse(text=x))})
			rhsp <- lapply(gsub(paste0("\\<",selected,"\\>"), get, rhsp), function(x){eval(parse(text=x))})
			rhsx2 <- lapply(gsub(paste0("\\<",selected,"\\>"), get, rhsx2), function(x){eval(parse(text=x))})
			rhsxp <- lapply(gsub(paste0("\\<",selected,"\\>"), get, rhsxp), function(x){eval(parse(text=x))})
			rhspx <- lapply(gsub(paste0("\\<",selected,"\\>"), get, rhspx), function(x){eval(parse(text=x))})
			rhsp2 <- lapply(gsub(paste0("\\<",selected,"\\>"), get, rhsp2), function(x){eval(parse(text=x))})
		}
		
		
		
		# Obtain isStateVariables[i]
		# - isStateVariables[i] indicates whether the i-th equation is for state variable (value = 1) or beta (value = 0), variables of covariates are stored in beta.names
		isStateVariables <- 1:length(lhs[[1]])
		for (i in 1:length(lhs[[1]])){
		    isStateVariables[[i]] = TRUE;
		    for (j in 1:length(beta.names)){
    		    if(lhs[[1]][[i]] == beta.names[[j]]){
    		        isStateVariables[[i]] = FALSE;
    		        break;
    		    }
		    }
		}
		
		# Start outputing the functions in converted_function.h
		#ret_head = "#include <iostream>\n#include <armadillo>\nusing namespace std;\nusing namespace arma;\nvoid function_arma_hello_world(void) {\n\t\tarma::mat a(2,2);\n\ta(0, 0) = 1;\n\ta(1, 1) = 2;\n\ta(0, 1) = -3;\n\ta(1, 0) = -2;\n\tprintf(\"hello world!\\n\");\n\ta.print();\n}\n\n"
		
		#ret_head = "#include <iostream>\n#include <armadillo>\nusing namespace std;\nusing namespace arma;\n\n\n"
		ret = "#// [[Rcpp::depends(RcppArmadillo)]]\n\n#include <RcppArmadillo.h>\nusing namespace std;\nusing namespace arma;\nusing namespace Rcpp;\n\n"
		
		# structure prototype
		ret= paste0(ret, "\n\nstruct C_INFDS{\n//public:\n\tint Nx, NxState, Ny, Nbeta, Ntheta, Nb, Nmu, NLambda, Nbpar, Nu, Npar0, totalT, Nsubj, Neta, Nbetax, N, MAXGIB;\n\tarma::mat Z, G,  tspan, H, U1, allDelta, thetatild, sytild, EStild, EItild, Iytild, bAdaptParams, OMEGAb, bacc, Sigmae, dSigmaede, dSigmaede2, Sigmab, dSigmabdb, dSigmabdb2, start, startpars, b, sy, ES, EI, Iy, lowBound, upBound, y0, SigmaEta, mu, dmudparMu, dmudparMu2, Lambda, dLambdparLamb, dLambdparLamb2, P0, par, trueb, Tfilter, tidx, lens, ICb;\n\tdouble omega, maxT, delt, gainpara, gainparb, errtrol, errtrol1, gainpara1, gainparb1, setAccept, scaleb;\n\tint alp, maxIterStage1, MAXITER, KKO, IT, isInfo;\n\tarma::field<arma::mat> fulldt, timeDiscrete, Deltat, meanY, dXstarAll, dXstarAll2, dXtildthetafAll, dXtildthetafAll2, tobs, Y;\n\tarma::cube Xtild;\n\tarma::Mat<double> allT;\n};\n\n")
		
		# function C extern
		ret= paste0(ret, "\n\n#ifdef __cplusplus\nextern \"C\" {\n#endif\n\narma::mat hello_world();\n\narma::mat dynfunICM(const int isPar, const arma::mat &xin, arma::vec &i, const int t, const int isStart, struct C_INFDS &InfDS);\narma::cube dfdxFreeICM(const int isPar, arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);\narma::cube dfdparFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);\narma::cube dfdx2FreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);\narma::cube dfdxdpFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);\narma::cube dfdpdxFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);\narma::cube dfdpar2FreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);\nvoid VDPMeas(arma::vec x, int Ny, int Nx, int NxState, arma::mat InfDS_Lambda, arma::mat mu, arma::mat *yPred, arma::mat *Jy);\nvoid setParsFreeICwb(C_INFDS &InfDS);\n\narma::mat rowProjection(arma::mat data, arma::mat index);\narma::vec span_vec(int start, int end, int step);\n\n#ifdef __cplusplus\n}\n#endif\n")
		
		
		#hello_world
		ret = paste0(ret, "\narma::mat hello_world() {\n\tarma::mat m1 = arma::eye<arma::mat>(3, 3);\n\tRprintf(\"Hello!!\\n\");\n\treturn m1;\n}\n\n")
		
		#span_vec
		ret = paste0(ret, "\narma::vec span_vec(int start, int end, int step){\n\tarma::vec x;\n\tint i;\n\n\tx.set_size((end - start)/step + 1);\n\tx[0] = start;\n\tfor(i = 1; i < (end - start)/step + 1; i++)\n\t\n\tx[i]= x[i-1] + step;\n\treturn x;\n}\n")
		
		#row_projection
		ret = paste0(ret, "\narma::mat rowProjection(arma::mat data, arma::mat index){\n\t//index is a row vector\n\tint i;\n\tarma::mat output;\n\toutput.reset();\n\tfor (i = 0; i < (int)index.n_cols; i++){\n\t\toutput.insert_rows(output.n_rows, data.row(index(i)-1));\n\t}\n\treturn output;\n}\n\n")
		
		#calculateTheta
		ret = paste0(ret, "\narma::mat calculateTheta(const int isPar, const arma::mat &y, arma::vec &i, struct C_INFDS &InfDS){\n\t\n\tif(isPar == 1)\n\t\tprintf(\"Should NOT happen! isPar = 1\");\n\t\n\tarma::mat b, betai, temp, Hi, thetaf, par;\n\t\n\tthetaf = arma::zeros<arma::mat>(InfDS.Ntheta, y.n_cols);\n\tpar = InfDS.par;\n\n\tif (isPar==0){\n\t\ttemp = span_vec(1, InfDS.Nbeta, 1).t();//Nbeta is the number of fixed effects parameters\n\t\tpar = rowProjection(par, temp);\n\t}\n\t\n\tb = arma::zeros<arma::mat>(thetaf.n_rows, thetaf.n_cols) ;\n\tbetai = arma::zeros<arma::mat>(InfDS.Ntheta, 1);// beta\n\tfor (int ii = 0; ii < int(i.n_elem); ii++){\n\t\tint current_i = int(i(ii));   \n\t\t\n\t\tHi = InfDS.H.rows(1+(current_i-1)*InfDS.Ntheta-1, current_i*InfDS.Ntheta-1);\n\t\t\n\t\tif (isPar == 1)\n\t\t\tbetai = y.col(ii).rows(InfDS.NxState,InfDS.NxState + Hi.n_cols-1);\n\t\telse\n\t\t\tbetai = reshape(par,InfDS.Nbeta, 1);\n\t\t\n\t\tif (InfDS.Nb > 0)\n\t\tb.col(ii) = reshape(InfDS.b.row(current_i-1),InfDS.Nb,1);\n\n\t\tthetaf.col(ii) = Hi * betai + InfDS.Z * b.col(ii);\n\t}\n\n\treturn thetaf;\n}\n\n")
		
		
		
		#----------------------------------------------------------------------------------------------
		# output code for function dynfun
		ret = paste0(ret, "arma::mat dynfunICM(const int isPar, const arma::mat &xin, arma::vec &i, const int t, const int isStart, struct C_INFDS &InfDS){\n\n\t//local parameters\n\tarma::mat y, r, thetaf;\n\t\n\t// if i is empty, traverse all vectors\n\tif(i.is_empty()){\n\t\ti = span_vec(1, InfDS.Nsubj, 1);\n\t}\t\n\n\t//input parameters\n\ty = xin;\n\tr.set_size(InfDS.Nx, int(i.n_elem));\n\tr.clear();\n\tr = y;\n")
		

		# Judge whether we needs calculateTheta
		c_i <- lapply(rhs, function(x){grep(paste0("thetaf(0,s)"),x, fixed = TRUE)})
		if(length(c_i[[1]]) > 0)	
			ret=paste0(ret,"\n\tthetaf=calculateTheta(isPar, y, i,InfDS);\n\n")
		
		# [todo] replace the initial condition by information from prep.initial
		# ret = paste0(ret,"\tif (isStart==1){\n\t\tr.zeros();\n\t\tint row, s;\n\t\tfor (s = 0; s < int(i.n_elem); s++){\n\t\t\tfor (row = 0; row < InfDS.NxState; row++)\n\t\t\t\tr(row, s) = thetaf(row +1, s);\n\n\t\t\tif (isPar == 1){\n\t\t\t\tfor (row = InfDS.NxState; row < InfDS.NxState + InfDS.Nbeta; row++){\n \t\t\t\t\tr(row, s) = y(row, s);\n \t\t\t\t}\n \t\t\t}\n\t\t}\n\n\t}\n")
		
		ret = paste0(ret, "\t\n\tr.zeros();\n\tint s;\n\tfor (s = 0; s < int(i.n_elem); s++){")
        for (i in 1:n){
			for (j in 1:length(lhs[[1]])){
			    # gsub (a, b, c) : in c replace a with b
				if(isStateVariables[[i]] == TRUE)
			    rhs[[1]][[i]]=gsub(paste0("\\<",lhs[[1]][[j]],"\\>"),paste0("y(",j-1,", s)"),rhs[[1]][[i]])
			}
            if(isStateVariables[[i]] == TRUE)
                ret=paste(ret,paste0("\t\t\tr(",i-1,", s)= ",rhs[[1]][[i]],";"),sep="\n")

        }

	    #ret=paste0(ret,"\n\t\t\tif (isPar == 1){\n\t\t\t\tfor (row = InfDS.NxState; row < InfDS.NxState + InfDS.Nbeta; row++){\n\t\t\t\t\tr(row, s) = 0;\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\n\treturn r;\n}\n")
		ret=paste0(ret,"\n\t}\n\treturn r;\n}\n")
	    ret=paste0(ret,"\n//----------------\n")

	    #----------------------------------------------------------------------------------------------
		
		# Output dfdxFreeICM
	    ret=paste0(ret, "arma::cube dfdxFreeICM(const int isPar, arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){\n\tarma::mat Hi, b, bb, thetaf, rowNum, y;\n\tarma::cube r;\n\n\t// if i is empty, traverse all vectors\n\tif(i.is_empty()){\n\t\ti = span_vec(1, InfDS.Nsubj, 1);\n\t}\n\n\ty = xin;\n\tr = arma::zeros<arma::cube>(InfDS.Nx, InfDS.Nx, y.n_cols) ;\n\tthetaf = arma::zeros<arma::mat>(InfDS.Ntheta, y.n_cols) ;\n")
		
		# Judge whether we needs calculateTheta
		c_i <- lapply(rhsj, function(x){grep(paste0("thetaf(0,s)"),x, fixed = TRUE)})
		if(length(c_i[[1]]) > 0)	
			ret=paste0(ret,"\n\tthetaf=calculateTheta(isPar, y, i,InfDS);\n\n")

		# repalce the state variable with corresponding InfDS variable
        for (i in 1:length(jacob[[1]])){
            for (j in 1:length(lhs[[1]])){
				rhsj[[1]][[i]]=gsub(paste0("\\<",lhs[[1]][[j]],"\\>"),paste0("y(",j-1,",s)"),rhsj[[1]][[i]])
            }
        }

		# [todo] replace the initial condition by information from prep.initial
		#ret=paste(ret, "\tif (isStart==1){\n\t\tif (isPar == 0)\n\t\t\t; //Undefined case in dfdxParIC\n\t\telse{\n\t\t\tint s;\n\t\t\tfor (s = 0; s < int(y.n_cols); s++){\n\t\t\t\tr.slice(s)(2, 2)=1 ;\n\t\t\t\tr.slice(s)(3, 3)=1 ;\n\t\t\t\tr.slice(s)(4, 4)=1 ;\n\t\t\t\tr.slice(s)(5, 5)=1 ;\n\t\t\t\tr.slice(s)(6, 6)=1 ;\n\t\t\t} \n\t\t}\n\t}\n")
		
	    ret=paste(ret,"\t\n\tint s;\n\n\tfor (s = 0; s < int(y.n_cols); s++){")
	    for (i in 1:length(jacob[[1]])){
	        col <- (i-1)%/%n  
	        row <- (i-1)%%n 
	        if(all(isStateVariables[col+1]==TRUE,isStateVariables[row+1]==TRUE)){
	            if(nchar(rhsj[[1]][[i]]) > 1 || !grepl("0",rhsj[[1]][[i]], useBytes = 1))
	                ret=paste(ret,paste0("\t\tr.slice(s)(",row,",",col,") = ",rhsj[[1]][[i]],";"),sep="\n")
	        }
	    }

	    # output(isPar)
	    # ret = paste(ret, "\n\t\t\tif(isPar == 1){")
	    # for (i in 1:length(jacob[[1]])){
	        # col <- (i-1)%/%n  
	        # row <- (i-1)%%n 
	        # if(any(isStateVariables[col+1]==FALSE, isStateVariables[row+1]==FALSE)){
	            # if(nchar(rhsj[[1]][[i]]) > 1 || !grepl("0",rhsj[[1]][[i]], useBytes = 1))
                    # ret=paste(ret,paste0("\t\t\t\tr.slice(s)(",row,",",col,") = ",rhsj[[1]][[i]],";"),sep="\n")
	        # }
	    # }
	    # ret = paste(ret, "\n\t\t\t\tr.slice(s) = r.slice(s).t();\n\t\t\t}\n\t\t}\n\treturn r;\n}\n")
		ret = paste(ret, "\n\t}\n\treturn r;\n}\n")
 	    ret=paste0(ret,"\n//----------------\n")
		
		#----------------------------------------------------------------------------------------------
		#Output dfdparFreeIC
		ret=paste0(ret, "arma::cube dfdparFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){\n\tarma::mat y ;\n\tarma::cube r;\n\n\t// if i is empty, traverse all vectors\n\tif(i.is_empty()){\n\t\ti = span_vec(1, InfDS.Nsubj, 1);\n\t}\n\ty = xin;\n\tr = arma::zeros<arma::cube>(InfDS.Ntheta, InfDS.NxState, y.n_cols);\n")
		
		# Judge whether we needs calculateTheta
		c_i <- lapply(rhsp, function(x){grep(paste0("thetaf(0,s)"),x, fixed = TRUE)})
		if(length(c_i[[1]]) > 0)	
			ret=paste0(ret,"\n\tthetaf=calculateTheta(0, y, i,InfDS);\n\n")
			
		for (i in 1:length(rhsp[[1]])){
            for (j in 1:length(lhs[[1]])){
				rhsp[[1]][[i]]=gsub(paste0("\\<",lhs[[1]][[j]],"\\>"),paste0("y(",j-1,",s)"),rhsp[[1]][[i]])
            }
        }
		ret = paste(ret, "\tint s;\n\tfor (s = 0; s < int(y.n_cols); s++){\t\t")
		for (i in 1:length(rhsp[[1]])){
	        col <- (i-1)%/%n.theta 
	        row <- (i-1)%%n.theta
	        
			#following line for debugging
			#ret = paste0(ret, paste0(i," ", row, " ", col, "\n"))
	        if(nchar(rhsp[[1]][[i]]) > 1 || !grepl("0",rhsp[[1]][[i]], useBytes = 1))
                ret=paste(ret,paste0("\t\tr.slice(s)(",row,",",col,") = ",rhsp[[1]][[i]],";"),sep="\n")
	        
	    }
		ret = paste(ret, "\n\t}\n\treturn r;\n}\n")
		ret=paste0(ret,"\n//----------------\n")
		
		#----------------------------------------------------------------------------------------------
		# Output dfdx2
		ret=paste0(ret, "arma::cube dfdx2FreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){\n\t// i and t are dummy variables\n\tarma::mat thetaf, y;\n\tarma::cube r;\n\n\t// if i is empty, traverse all vectors\n\tif(i.is_empty()){\n\t\ti = span_vec(1, InfDS.Nsubj, 1);\n\t}\n\n\ty = xin ;  \n\tr = arma::zeros<arma::cube>(InfDS.Nx * InfDS.Nx, InfDS.Nx, y.n_cols) ;\n\tthetaf = arma::zeros<arma::mat>(InfDS.Ntheta, y.n_cols) ;\n")
		
		# Judge whether we needs calculateTheta
		c_i <- lapply(rhsx2, function(x){grep(paste0("thetaf(0,s)"),x, fixed = TRUE)})
		if(length(c_i[[1]]) > 0)	
			ret=paste0(ret,"\n\tthetaf=calculateTheta(0, y, i,InfDS);\n\n")
			
		for (i in 1:length(rhsx2[[1]])){
            for (j in 1:length(lhs[[1]])){
				rhsx2[[1]][[i]]=gsub(paste0("\\<",lhs[[1]][[j]],"\\>"),paste0("y(",j-1,",s)"),rhsx2[[1]][[i]])
            }
        }
		
		ret = paste(ret, "\tint s;\n\tfor (s = 0; s < int(y.n_cols); s++){\t\t")
		for (i in 1:length(rhsx2[[1]])){
	        row <- (i-1)%/%n.state
	        col <- (i-1)%%n.state
	        
			#following line for debugging
			#ret = paste0(ret, paste0(i," ", row, " ", col, "\n"))
	        if(nchar(rhsx2[[1]][[i]]) > 1 || !grepl("0",rhsx2[[1]][[i]], useBytes = 1))
				ret=paste(ret,paste0("\t\tr.slice(s)(",row,",",col,") = ",rhsx2[[1]][[i]],";"),sep="\n")
	        
	    }
		ret = paste(ret, "\n\t}\n\treturn r;\n}\n")
		ret=paste0(ret,"\n//----------------\n")
		
		#----------------------------------------------------------------------------------------------
		# Output dfdxdp
		ret=paste0(ret, "arma::cube dfdxdpFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){\n\t// i and t are dummy variables\n\tarma::mat y;\n\tarma::cube r;\n\n\t// if i is empty, traverse all vectors\n\tif(i.is_empty()){\n\t\ti = span_vec(1, InfDS.Nsubj, 1);\n\t}\n\n\ty = xin ;  \n\tr = arma::zeros<arma::cube>(InfDS.Nx * InfDS.Nx, InfDS.Ntheta, y.n_cols);\n")
		
		#browser()
		# Judge whether we needs calculateTheta
		c_i <- lapply(rhsxp, function(x){grep(paste0("thetaf(0,s)"),x, fixed = TRUE)})
		if(length(c_i[[1]]) > 0)	
			ret=paste0(ret,"\n\tthetaf=calculateTheta(0, y, i,InfDS);\n\n")
			
		for (i in 1:length(rhsxp[[1]])){
            for (j in 1:length(lhs[[1]])){
				rhsxp[[1]][[i]]=gsub(paste0("\\<",lhs[[1]][[j]],"\\>"),paste0("y(",j-1,",s)"),rhsxp[[1]][[i]])
            }
        }
		
		ret = paste(ret, "\tint s;\n\tfor (s = 0; s < int(y.n_cols); s++){\t\t")
		for (i in 1:length(rhsxp[[1]])){
	        row <- (i-1)%/%n.theta
	        col <- (i-1)%%n.theta
	        
			#following line for debugging
			#ret = paste0(ret, paste0(i," ", row, " ", col, "\n"))
	        if(nchar(rhsxp[[1]][[i]]) > 1 || !grepl("0",rhsxp[[1]][[i]], useBytes = 1))
				ret=paste(ret,paste0("\t\tr.slice(s)(",row,",",col,") = ",rhsxp[[1]][[i]],";"),sep="\n")
	        
	    }
		ret = paste(ret, "\n\t}\n\treturn r;\n}\n")
		ret=paste0(ret,"\n//----------------\n")
		
		
		#----------------------------------------------------------------------------------------------
		# Output dfdpdx
		ret=paste0(ret, "arma::cube dfdpdxFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){\n\t// i and t are dummy variables\n\tarma::mat  y;\n\tarma::cube r;\n\n\t// if i is empty, traverse all vectors\n\tif(i.is_empty()){\n\t\ti = span_vec(1, InfDS.Nsubj, 1);\n\t}\n\n\ty = xin ;  \n\tr = arma::zeros<arma::cube>(InfDS.Nx*InfDS.Ntheta, InfDS.Nx, y.n_cols) ;\n")
		
		# Judge whether we needs calculateTheta
		c_i <- lapply(rhspx, function(x){grep(paste0("thetaf(0,s)"),x, fixed = TRUE)})
		if(length(c_i[[1]]) > 0)	
			ret=paste0(ret,"\n\tthetaf=calculateTheta(0, y, i,InfDS);\n\n")
			
		for (i in 1:length(rhspx[[1]])){
            for (j in 1:length(lhs[[1]])){
				rhspx[[1]][[i]]=gsub(paste0("\\<",lhs[[1]][[j]],"\\>"),paste0("y(",j-1,",s)"),rhspx[[1]][[i]])
            }
        }
		
		ret = paste(ret, "\tint s;\n\tfor (s = 0; s < int(y.n_cols); s++){\t\t")
		for (i in 1:length(rhspx[[1]])){
	        row <- (i-1)%/%n.state
	        col <- (i-1)%%n.state
	        
			#following line for debugging
			#ret = paste0(ret, paste0(i," ", row, " ", col, "\n"))
	        if(nchar(rhspx[[1]][[i]]) > 1 || !grepl("0",rhspx[[1]][[i]], useBytes = 1))
				ret=paste(ret,paste0("\t\tr.slice(s)(",row,",",col,") = ",rhspx[[1]][[i]],";"),sep="\n")
	        
	    }
		ret = paste(ret, "\n\t}\n\treturn r;\n}\n")
		ret=paste0(ret,"\n//----------------\n")

		
		#----------------------------------------------------------------------------------------------
		#Output dfdp2
		ret=paste0(ret, "arma::cube dfdpar2FreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS){\n\t// i and t are dummy variables\n\tarma::mat y ;\n\tarma::cube r;\n\n\t// if i is empty, traverse all vectors\n\tif(i.is_empty()){\n\t\ti = span_vec(1, InfDS.Nsubj, 1);\n\t}\n\n\ty = xin ;  \n\tr = arma::zeros<arma::cube>(InfDS.NxState*InfDS.Ntheta, InfDS.Ntheta, y.n_cols) ;\n\n")
		
		# Judge whether we needs calculateTheta
		c_i <- lapply(rhsp2, function(x){grep(paste0("thetaf(0,s)"),x, fixed = TRUE)})
		if(length(c_i[[1]]) > 0)	
			ret=paste0(ret,"\n\tthetaf=calculateTheta(0, y, i,InfDS);\n\n")
			
		for (i in 1:length(rhsp2[[1]])){
            for (j in 1:length(lhs[[1]])){
				rhsp2[[1]][[i]]=gsub(paste0("\\<",lhs[[1]][[j]],"\\>"),paste0("y(",j-1,",s)"),rhsp2[[1]][[i]])
            }
        }
		
		ret = paste(ret, "\tint s;\n\tfor (s = 0; s < int(y.n_cols); s++){\n\t\t;")
		for (i in 1:length(rhsp2[[1]])){
	        row <- (i-1)%/%n.theta
	        col <- (i-1)%%n.theta
	        
			#following line for debugging
			#ret = paste0(ret, paste0(i," ", row, " ", col, "\n"))
	        if(nchar(rhsp2[[1]][[i]]) > 1 || !grepl("0",rhsp2[[1]][[i]], useBytes = 1))
				ret=paste(ret,paste0("\t\tr.slice(s)(",row,",",col,") = ",rhsp2[[1]][[i]],";"),sep="\n")
	        
	    }
		ret = paste(ret, "\n\t}\n\treturn r;\n}\n")
		
		#----------------------------------------------------------------------------------------------
		#Output setParsFreeICwb
		
		#browser()
		
		
		# ret = paste(ret, "\n\tstartB = InfDS.Nbeta + InfDS.Nmu + InfDS.NLambda + InfDS.Ny + 1; \n")
		# Sigmab = object@random.params.inicov
		# sigmab.names <- unique(c(object$random.params.inicov))
		# ret = paste0(ret, "\tInfDS.Sigmab = arma::zeros<arma::mat>(", nrow(Sigmab), ", ", ncol(Sigmab), ");\n")
		# startB = min((1:length(param.names))[ param.names %in% object@random.params.inicov])+1
		# for(i in 1: nrow(Sigmab)){
			# for(j in 1: ncol(Sigmab)){
				# if(Sigmab[i, j] > 0){
					# #print(param.names[Lambda[i, j]])
					# for(k in (1:length(param.names))[ param.names %in% object@random.params.inicov] ){
						# if(k - startB < 0)
							# expr = gsub(param.names[k], paste0("exp(InfDS.par(startB - ", -(k-startB),"))"), expr)
						# else
							# expr = gsub(param.names[k], paste0("exp(InfDS.par(startB + ", (k-startB),"))"), expr)
					# }
					# ret = paste0(ret, "\tInfDS.Sigmab(", i-1, ", ", j-1, ") = ", expr , ";\n")
				# }
			# }
			
		# }
		
		ret=paste0(ret,"\n//----------------\n")
		ret=paste0(ret, "void setParsFreeICwb(C_INFDS &InfDS){\n\tint startM, startL, startE, startB;\n\tarma::mat D, L;\n\tInfDS.par = real(InfDS.par);\n\n\n")

		ret = paste(ret, "\n\tstartB = InfDS.Nbeta + InfDS.Nmu + InfDS.NLambda + InfDS.Ny + 1; \n")
		
		
		ret = paste0(ret, "\tInfDS.Sigmab = arma::zeros<arma::mat>(", nrow(Sigmab), ", ", ncol(Sigmab), ");\n")
		for(i in 1: nrow(Sigmab)){
			for(j in 1: ncol(Sigmab)){
				expr.char = deparse(Sigmab[i,j][[1]])
				if(length(expr.char) > 0 && expr.char[1] != "0"){
					#print(param.names[Lambda[i, j]])
					for(k in (1:length(known.vars))){
					    if(k == 1)
							expr.char = gsub(known.vars[k], paste0("InfDS.par(startB-", -(k-2),")"), expr.char )
						else
							expr.char = gsub(known.vars[k], paste0("InfDS.par(startB+", k-2,")"), expr.char )
					}
					ret = paste0(ret, "\tInfDS.Sigmab(", i-1, ", ", j-1, ") = ", expr.char  , ";\n")
				}
			}
			
		}
		
		ret = paste0(ret, "\tInfDS.dSigmabdb = arma::zeros<arma::mat>(", nrow(dSigmabdb), ", ", ncol(dSigmabdb), ");\n")
		for(i in 1: nrow(dSigmabdb)){
			for(j in 1: ncol(dSigmabdb)){
				expr.char = deparse(dSigmabdb[i,j][[1]])
				if(length(expr.char) > 0 && expr.char[1] != "0"){
					#print(param.names[Lambda[i, j]])
					for(k in (1:length(known.vars))){
						if(k == 1)
							expr.char = gsub(known.vars[k], paste0("InfDS.par(startB-", -(k-2),")"), expr.char )
						else
							expr.char = gsub(known.vars[k], paste0("InfDS.par(startB+", k-2,")"), expr.char )
					}
					ret = paste0(ret, "\tInfDS.dSigmabdb(", i-1, ", ", j-1, ") = ", expr.char  , ";\n")
				}
			}
			
		}
		
		ret = paste0(ret, "\tInfDS.dSigmabdb2 = arma::zeros<arma::mat>(", nrow(dSigmabdb2), ", ", ncol(dSigmabdb2), ");\n")
		for(i in 1: nrow(dSigmabdb2)){
			for(j in 1: ncol(dSigmabdb2)){
				expr.char = deparse(dSigmabdb2[i,j][[1]])
				if(length(expr.char) > 0 && expr.char[1] != "0"){
					#print(param.names[Lambda[i, j]])
					for(k in (1:length(known.vars))){
						if(k == 1)
							expr.char = gsub(known.vars[k], paste0("InfDS.par(startB-", -(k-2),")"), expr.char )
						else
							expr.char = gsub(known.vars[k], paste0("InfDS.par(startB+", k-2,")"), expr.char )
					}
					ret = paste0(ret, "\tInfDS.dSigmabdb2(", i-1, ", ", j-1, ") = ", expr.char  , ";\n")
				}
			}
			
		}
		
		# ret = paste0(ret, "\tInfDS.dSigmabdb2 = arma::zeros<arma::mat>(", nrow(dSigmabdb2), ", ", ncol(dSigmabdb2), ");\n")
		# for(i in 1: nrow(dSigmabdb2)){
			# for(j in 1: ncol(dSigmabdb2)){
				# if(dSigmabdb2[i, j] > 0){
					# #print(param.names[Lambda[i, j]])
					# for(k in (1:length(param.names))[ param.names %in% object@random.params.inicov] ){
						# if(k - startB < 0)
							# expr = gsub(param.names[k], paste0("InfDS.par(startB - ", -(k-startB),")"), expr)
						# else
							# expr = gsub(param.names[k], paste0("InfDS.par(startB + ", (k-startB),")"), expr)
					# }
					# ret = paste0(ret, "\tInfDS.dSigmabdb2(", i-1, ", ", j-1, ") = ", expr , ";\n")
				# }
			# }
			
		#}
		
		
		#ret contains all the transferred code, currently uses ret_head for testing
        #object@c.string <- ret
		object@c.string <- ret
        return(object)
	}
)

#------------------------------------------------------------------------------
# Define method for dynrNoise class

# TODO This is probably the easiest armadillo code to write
#  It will be the most similar to the corresponding writeCcode method.

setMethod("writeArmadilloCode", "dynrNoise",
	function(object, covariates, param.names, dSigmaede, dSigmaede2){
		#browser()
		ret <- ""
		ret = paste(ret, "\n\tstartE = InfDS.Nbeta + InfDS.Nmu + InfDS.NLambda + 1;\n")
		#Sigmae = object$params.observed[[1]]
		if (length(object$params.observed) > 0){ 
		  Sigmae <- object$params.observed[[1]]
		} else { Sigmae <- matrix(0L, nrow=0, ncol=0) }
		
		startE = min(Sigmae[Sigmae>0]) + 1
		ret = paste0(ret, "\tInfDS.Sigmae = arma::zeros<arma::mat>(", nrow(Sigmae), ", ", ncol(Sigmae), ");\n")
		for(i in 1: nrow(Sigmae)){
			for(j in 1: ncol(Sigmae)){
				if(Sigmae[i, j] > 0){
					#print(param.names[Lambda[i, j]])
					if(Sigmae[i, j] - startE < 0)
						ret = paste0(ret, "\tInfDS.Sigmae(", i-1, ", ", j-1, ") = exp(InfDS.par(startE - ", -(Sigmae[i, j] - startE), '));\n')
					else
						ret = paste0(ret, "\tInfDS.Sigmae(", i-1, ", ", j-1, ") = exp(InfDS.par(startE + ", Sigmae[i, j] - startE, '));\n')
				}
			}
			
		}
		
		ret = paste0(ret, "\tInfDS.dSigmaede = arma::zeros<arma::mat>(", nrow(dSigmaede), ", ", ncol(dSigmaede), ");\n")
		for(i in 1: nrow(dSigmaede)){
			for(j in 1: ncol(dSigmaede)){
				expr = as.character(dSigmaede[i,j])
				if(expr != "0"){
					for(k in min(Sigmae[Sigmae>0]):max(Sigmae[Sigmae>0])){
						if(k - startE < 0)
							expr = gsub(param.names[k], paste0("InfDS.par(startE - ", -(k-startE),")"), expr)
						else
							expr = gsub(param.names[k], paste0("InfDS.par(startE + ", (k-startE),")"), expr)
					}
					ret = paste0(ret, "\tInfDS.dSigmaede(", i-1, ", ", j-1, ") = ", expr , ";\n")
				}
			}
			
		}
		
		ret = paste0(ret, "\tInfDS.dSigmaede2 = arma::zeros<arma::mat>(", nrow(dSigmaede2), ", ", ncol(dSigmaede2), ");\n")
		for(i in 1: nrow(dSigmaede2)){
			for(j in 1: ncol(dSigmaede2)){
				expr = as.character(dSigmaede2[i,j])
				if(expr != "0"){
					for(k in min(Sigmae[Sigmae>0]):max(Sigmae[Sigmae>0])){
						if(k - startE < 0)
							expr = gsub(param.names[k], paste0("InfDS.par(startE - ", -(k-startE),")"), expr)
						else
							expr = gsub(param.names[k], paste0("InfDS.par(startE + ", (k-startE),")"), expr)
					}
					ret = paste0(ret, "\tInfDS.dSigmaede2(", i-1, ", ", j-1, ") = ", expr , ";\n")
				}
			}
			
		}
		
		ret=paste0(ret,"}\n//----------------\n")
		
		object@c.string <- ret
		return(object)
	}
)

#------------------------------------------------------------------------------
# Define method for dynrInitial class

setMethod("writeArmadilloCode", "dynrInitial",
	function(object, covariates){
		ret <- ""
		object@c.string <- ret
		return(object)
	}
)

#------------------------------------------------------------------------------
# Define method for dynrRegime class

setMethod("writeArmadilloCode", "dynrRegimes",
	function(object, covariates){
		ret <- ""
		object@c.string <- ret
		return(object)
	}
)

#------------------------------------------------------------------------------
# Define method for dynrRegime class

setMethod("writeArmadilloCode", "dynrRandom",
	function(object, covariates){
		ret <- ""
		object@c.string <- ret
		return(object)
	}
)


#------------------------------------------------------------------------------
# End
