#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2017-09-14
# Filename: dynrTaste.R
# Purpose: Compute shocks and chi-squared diagnostics following
#  Chow, Hamaker, and Allaire (2009).  Using Innovative Outliers to
#    Detect Discrete Shifts in Dynamics in Group-Based State-Space Models.
#    Multivariate Behavioral research.  44:465–496.
#    DOI: 10.1080/00273170903103324
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
require(dynr)
data(NonlinearDFAsim)
EMGdata <- dynr.data(NonlinearDFAsim, id='id', time='time', 
                     observed=c('y1','y2','y3','y4','y5','y6'))

#---- (3) Specify recipes for all model pieces ----

#---- (3a) Measurement ----
recMeas <- prep.measurement(
  values.load = matrix(c(1, 0,
                         0.5, 0,
                         0.5, 0,
                         0, 1,
                         0, 0.5,
                         0, 0.5), nrow=6, byrow=TRUE),
  params.load=matrix(c(0, 0,
                       'lambda_P2', 0,
                       'lambda_P3', 0,
                       0, 0,
                       0, 'lambda_N5',
                       0, 'lambda_N6'), nrow=6, byrow=TRUE),
  #values.int = ,
  #params.int = ,
  obs.names = c('y1','y2','y3','y4','y5','y6'),
  state.names = c('P','N'))

#---- (3b) Dynamic and measurement noise cov structures----

recNoise <- prep.noise(
  values.latent=matrix(c(0.5, 0,
                         0, 0.5), nrow=2, byrow=TRUE),
  params.latent=matrix(c('psi_P', 0,
                         0, 'Psi_N'), nrow=2, byrow=TRUE),
  values.observed=diag(rep(0.5, 6)),
  params.observed=diag(paste0('theta_', 1:6)) )

#---- (3d) Initial condition specification ----

recIni <- prep.initial(
  values.inistate=c(0.5, 0.5),
  params.inistate=c(0, 0),
  values.inicov=matrix(c(0.5, 0,
                         0, 0.5), nrow=2, byrow=TRUE),
  params.inicov=matrix(c(0, 0,
                         0, 0), nrow=2) )


#---- (3e) Dynamic model ----

recDyn <- prep.matrixDynamics(
  values.dyn=matrix(c(0.5, 0,
                      0, 0.5), nrow=2, byrow=TRUE),
  params.dyn=matrix(c('beta_P', 0,
                      0, 'beta_N'), nrow=2, byrow=TRUE),
  isContinuousTime = FALSE)

#---- (4a) Create model  ----

rsmod <- dynr.model(
  measurement = recMeas,
  dynamics = recDyn,
  noise = recNoise,
  initial = recIni,
  data = EMGdata,
  outfile = "RSLinearDiscrete.c")

plotFormula(dynrModel = rsmod, ParameterAs = rsmod$param.names,
            printDyn = TRUE, printMeas = TRUE) +
  ggtitle("(A)")+
  theme(plot.title = element_text(hjust = 0.5, vjust=0.01, size=16)) 


#---- (4c) Create model and cook it all up  ----

yum <- dynr.cook(rsmod, debug_flag=TRUE)
#yum2 <- dynr.cook(rsmod, debug_flag=FALSE)
F <- yum@residual_cov # array [p, p, t] ?
v <- yum@innov_vec # matrix [p, t] ?
P_1 <- yum@error_cov_predicted # array, [q, q, t]
eta <- yum@eta_filtered # matrix, [q, t]
eta_1 <- yum@eta_predicted # matrix, [q, t]
Lambda <- matrix(1, 1, 1)
#---- (5) Serve it! ----

##############
x <- yum
model <- rsmod

dynr.taste <- function(cookDebug, dynrModel, alpha=0.95) {
  #TODO check for non-regime switching
	x <- cookDebug
	model <- dynrModel
	
	latentDim <- dim(x$error_cov_predicted)[1]
	observedDim <- dim(x$residual_cov)[1]
	timeDim <- dim(x$residual_cov)[3]
	personID <- unique(model$data$id)
	
	# fitted parameters
	coefx <- coef(x)
	
	# Measurement matrix
	# Lambda in Chow, Hamaker, and Allaire
	# number matches the parameters in coefx
	# e.g., 0 0
	#       3 0
	#       4 0
	# 3, 4th elements in coefx match
	Lambda_p <- model$measurement$params.load[[1]]
	Lambda <- model$measurement$values.load[[1]]
	Lambda_idx <- which(Lambda_p !=0, arr.ind=TRUE)
	Lambda[Lambda_idx] <- coefx[ Lambda_p[Lambda_idx] ]
	
	# Dynamics matrix
	# B in Chow, Hamaker, and Allaire
	B_p <- model$dynamics$params.dyn[[1]]
	B <- model$dynamics$values.dyn[[1]]
	B_idx <- which(B_p !=0, arr.ind=TRUE)
	B[B_idx] <- coefx[ B_p[B_idx] ]
	
	# Array of inverse covariance matrices (i.e. information matrices) for the observed variables
	# F^-1  in Chow, Hamaker, and Allaire
	F_inv <- array(apply(x$residual_cov, 3, solve), 
	               c(observedDim, observedDim, timeDim))
	
	# Compute Kalman gain for every person/time combination
	P_pred <- x$error_cov_predicted
	t_Lambda <- t(Lambda)
	v <- x$innov_vec
	K <- array(NA, c(observedDim, latentDim, timeDim))
	obsChi <- numeric(timeDim)
	
	for(i in 1:timeDim){
		F_inv_i <- F_inv[,,i]
		K[,,i] <- P_pred[,,i] %*% t_Lambda %*% F_inv_i
		v_i <- v[,i] # matrix(v[,i], ncol=1)
		obsChi[i] <- t(v_i) %*% F_inv_i %*% v_i
	}
	
	#qqplot(qchisq(ppoints(timeDim), df=observedDim), obsChi)
	#qqline(obsChi, distribution = function(p) qchisq(p, df = observedDim), prob = c(0.1, 0.6))
	
	# For each person, loop backward from their final time to their first time
	# computing r and N
	r <- matrix(NA, latentDim, timeDim)
	N <- array(NA, c(latentDim, latentDim, timeDim))
	model_tstart <- model$data$tstart
	latChi <- numeric(timeDim)
	for(j in 1:length(personID)){
		beginTime <- model_tstart[j] + 1
		endTime <- model_tstart[j+1]
		# set endTime r and N to 0
		r[,endTime] <- 0
		N[,,endTime] <- 0
		for(i in (endTime-1):beginTime){
			# Probably could/should compute shocks in this same loop!!!
			ri <- matrix(r[,i+1], latentDim, 1)
			Ni <- matrix(N[,,i+1], latentDim, latentDim)
			Ki <- matrix(K[,,i+1], latentDim, observedDim)
			Li <- B - Ki %*% Lambda
			obsInfI <- matrix(F_inv[,,i+1], observedDim, observedDim) #TODO check for off by one index
			vi <- matrix(v[,i+1], nrow=observedDim, ncol=1) #TODO check for off by one index
			ui <- obsInfI %*% vi - t(Ki) %*% ri # TODO save ui? It is used in s_ij Eq 26.
			# Could also just compute shocks in this loop
			rnew <- t_Lambda %*% ui + t(B) %*% ri
			r[,i] <- rnew
			Nnew <- t_Lambda %*% obsInfI %*% Lambda + t(Li) %*% Ni %*% Li
			N[,,i] <- Nnew
			Ninv <- try(solve(Nnew), silent=TRUE)
			if(class(Ninv) == "try-error"){Ninv <- MASS::ginv(Nnew)}
			latChi[i] <- t(rnew) %*% Ninv %*% rnew
		}
	}
	#qqplot(qchisq(ppoints(timeDim), df=latentDim), latChi)
	#qqline(latChi, distribution = function(p) qchisq(p, df = latentDim), prob = c(0.1, 0.6))
	
	#qqplot(qchisq(ppoints(timeDim), df=observedDim+latentDim), obsChi + latChi)
	#qqline(obsChi + latChi, distribution = function(p) qchisq(p, df = observedDim+latentDim), prob = c(0.1, 0.6))
	
	
	# Return list of r array and N array
	# perhaps automatically plot
	# perhaps add plot=TRUE argument
	
	# t-test
	t_values <- matrix(NA, latentDim, timeDim)
	for (t in 1:timeDim) {
	  t_values[,t] <- r[,t] / sqrt(diag(N[,,t]))  
	}

	# one-sided test?
	t_start <- (model_tstart + 1)[-length(model_tstart)]
	t_end <- model_tstart[-1]
	t_df <- t_end - model_tstart[-length(model_tstart)]
	
	t_values_split <- lapply(1:length(t_end), function(i) {
		t_values[, t_start[i]:t_end[i]]
	})
	
	t_test_split <- lapply(1:length(t_end), function(i) {
		critical_value <- qt(alpha, df=(t_df[i]-latentDim))
		# shock point index: latent * time
		which(t_values_split[[i]] > critical_value, arr.ind=FALSE)  
	})
	
	## test observed chi-square. index time points
	chi_test_obs <- which(obsChi > qchisq(alpha, observedDim))
	## test latent chi-square. index time points
	chi_test_lat <- which(latChi > qchisq(alpha, latentDim))
	## test combined chi-square
	chi_test_comb <- which( (obsChi + latChi) > qchisq(alpha, (observedDim + latentDim)))
	
	list
		r=r, N=N, 
		chi_obs=obsChi, chi_lat=latChi,
		t_test=t_test_split,
		chi_test_obs=chi_test_obs, chi_test_lat=chi_test_lat,
		chi_test_comb=chi_test_comb)
}











