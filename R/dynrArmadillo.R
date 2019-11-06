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
	function(object, covariates, show=TRUE) { 
		return(standardGeneric("writeArmadilloCode")) 
	})


#------------------------------------------------------------------------------
# Define method for dynrMeasurement class

setMethod("writeArmadilloCode", "dynrMeasurement",
	function(object, covariates){
		ret <- ""
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
	function(object, covariates){
		#browser()
		formula <- object$formula
		formula2 <- object$formula2
		jacob <- object$jacobianOriginal
		dfdtheta<- object$dfdtheta
		dfdx2<- object$dfdx2
		dfdxdtheta<- object$dfdxdtheta
		dfdthetadx<- object$dfdthetadx
		dfdtheta2<- object$dfdtheta2
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
		
		#Parse theta.formula
		fmlt = processFormula(theta.formula)
		lhst = fmlt[[1]][1]
		rhst = fmlt[[1]][2]
		
		
		# Replace theta_i in formula with thetaf
		# - thetaf is calculated by calculateTheta()
		# - the variable name of thetaf is from the LHS of theta.formula
		# - in jacobian (dfdx), LHS of theta.formula is already replaced by RHS of theta.formula in rhsj (to get correct differentiation), thus, we replace the RHS of theta formula by thetaf
		rhs <- lapply(rhs, function(x){gsub(paste0(lhst),paste0("thetaf(0,s)"),x, fixed = TRUE)})
		rhsj <- lapply(rhsj, function(x){gsub(paste0(lhst),paste0("thetaf(0,s)"),x, fixed = TRUE)})
		rhsp <- lapply(rhsp, function(x){gsub(paste0(lhst),paste0("thetaf(0,s)"),x, fixed = TRUE)})
		rhsx2 <- lapply(rhsx2, function(x){gsub(paste0(lhst),paste0("thetaf(0,s)"),x, fixed = TRUE)})
		rhsxp <- lapply(rhsxp, function(x){gsub(paste0(lhst),paste0("thetaf(0,s)"),x, fixed = TRUE)})
		rhspx <- lapply(rhspx, function(x){gsub(paste0(lhst),paste0("thetaf(0,s)"),x, fixed = TRUE)})
		rhsp2 <- lapply(rhsp2, function(x){gsub(paste0(lhst),paste0("thetaf(0,s)"),x, fixed = TRUE)})
		
		
		
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
		
		ret_head = "#include <iostream>\n#include <armadillo>\nusing namespace std;\nusing namespace arma;\n\n\n"
		#----------------------------------------------------------------------------------------------
		# output code for function dynfun
		ret = "arma::mat dynfunICM(const int isPar, const arma::mat &xin, arma::vec &i, const int t, const int isStart, struct C_INFDS &InfDS){\n\n\t//local parameters\n\tarma::mat y, r, thetaf;\n\t\n\t// if i is empty, traverse all vectors\n\tif(i.is_empty()){\n\t\ti = span_vec(1, InfDS.Nsubj, 1);\n\t}\t\n\n\t//input parameters\n\ty = xin;\n\tr.set_size(InfDS.Nx, int(i.n_elem));\n\tr.clear();\n\tr = y;\n"
		
		ret = paste0(ret_head, ret)

		# Judge whether we needs calculateTheta
		c_i <- lapply(rhs, function(x){grep(paste0("thetaf(0,s)"),x, fixed = TRUE)})
		if(length(c_i[[1]]) > 0)	
			ret=paste0(ret,"\n\tthetaf=calculateTheta(isPar, y, i,InfDS);\n\n")
		
		# [todo] replace the initial condition by information from prep.initial
		ret = paste0(ret,"\tif (isStart==1){\n\t\tr.zeros();\n\t\tint row, s;\n\t\tfor (s = 0; s < int(i.n_elem); s++){\n\t\t\tfor (row = 0; row < InfDS.NxState; row++)\n\t\t\t\tr(row, s) = thetaf(row +1, s);\n\n\t\t\tif (isPar == 1){\n\t\t\t\tfor (row = InfDS.NxState; row < InfDS.NxState + InfDS.Nbeta; row++){\n \t\t\t\t\tr(row, s) = y(row, s);\n \t\t\t\t}\n \t\t\t}\n\t\t}\n\n\t}\n")
		
		ret = paste0(ret, "\telse{\n\t\tr.zeros();\n\t\tint row, s;\n\t\tfor (s = 0; s < int(i.n_elem); s++){")
        for (i in 1:n){
			for (j in 1:length(lhs[[1]])){
			    # gsub (a, b, c) : in c replace a with b
				if(isStateVariables[[i]] == TRUE)
			    rhs[[1]][[i]]=gsub(paste0("\\<",lhs[[1]][[j]],"\\>"),paste0("y(",j-1,", s)"),rhs[[1]][[i]])
			}
            if(isStateVariables[[i]] == TRUE)
                ret=paste(ret,paste0("\t\t\tr(",i-1,", s)= ",rhs[[1]][[i]],";"),sep="\n")

        }

	    ret=paste0(ret,"\n\t\t\tif (isPar == 1){\n\t\t\t\tfor (row = InfDS.NxState; row < InfDS.NxState + InfDS.Nbeta; row++){\n\t\t\t\t\tr(row, s) = 0;\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t}\n\treturn r;\n}\n")
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
		ret=paste(ret, "\tif (isStart==1){\n\t\tif (isPar == 0)\n\t\t\t; //Undefined case in dfdxParIC\n\t\telse{\n\t\t\tint s;\n\t\t\tfor (s = 0; s < int(y.n_cols); s++){\n\t\t\t\tr.slice(s)(2, 2)=1 ;\n\t\t\t\tr.slice(s)(3, 3)=1 ;\n\t\t\t\tr.slice(s)(4, 4)=1 ;\n\t\t\t\tr.slice(s)(5, 5)=1 ;\n\t\t\t\tr.slice(s)(6, 6)=1 ;\n\t\t\t} \n\t\t}\n\t}\n")
		
	    ret=paste(ret,"\telse{\n\t\tint s;\n\n\t\tfor (s = 0; s < int(y.n_cols); s++){")
	    for (i in 1:length(jacob[[1]])){
	        col <- (i-1)%/%n  
	        row <- (i-1)%%n 
	        if(all(isStateVariables[col+1]==TRUE,isStateVariables[row+1]==TRUE)){
	            if(nchar(rhsj[[1]][[i]]) > 1 || !grepl("0",rhsj[[1]][[i]], useBytes = 1))
	                ret=paste(ret,paste0("\t\t\tr.slice(s)(",row,",",col,") = ",rhsj[[1]][[i]],";"),sep="\n")
	        }
	    }

	    # output(isPar)
	    ret = paste(ret, "\n\t\t\tif(isPar == 1){")
	    for (i in 1:length(jacob[[1]])){
	        col <- (i-1)%/%n  
	        row <- (i-1)%%n 
	        if(any(isStateVariables[col+1]==FALSE, isStateVariables[row+1]==FALSE)){
	            if(nchar(rhsj[[1]][[i]]) > 1 || !grepl("0",rhsj[[1]][[i]], useBytes = 1))
                    ret=paste(ret,paste0("\t\t\t\tr.slice(s)(",row,",",col,") = ",rhsj[[1]][[i]],";"),sep="\n")
	        }
	    }
	    ret = paste(ret, "\n\t\t\t\tr.slice(s) = r.slice(s).t();\n\t\t\t}\n\t\t}\n\t}\n\treturn r;\n}\n")
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
	function(object, covariates){
		ret <- ""
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
