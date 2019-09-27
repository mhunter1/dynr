##' A function that takes the data and some Bsplines specification parameters, and return a matrix of values that contain the
##' generalized cross-validation values at the specified range of lambda values
##' @param theTimes The time points at which derivative estimation are requested
##' @param norder #Order of Bsplines - usually 2 higher than roughPenaltyMax
##' @param roughPenaltyMax  Penalization order. Usually set to 2 higher than the highest-order derivatives desired
##' @param  dataMatrix Data of size total number of time points x  total number of subjects
##' @param LowLambda Lower limit of lambda values to be tested. Here, lambda is a positive smoothing parameter, larger means more smoothing)
##' @param upLambda Upper limit of lambda
##' @param lambdaInt = The interval of lambda values to be tested.
##' @param isPlot A binary flag on whether to plot the gcv values (0 = no, 1 = yes)
##' 
##' @return a data frame containing:
##' 1. lambda values; 2. edf (effective degrees of freedom); 3. GCV (average generalized cross-validation values)
##' 
##' @references
##' Chow, S-M. (2019). Practical Tools and Guidelines for Exploring and Fitting Linear 
##' and Nonlinear Dynamical Systems Models. Multivariate Behavioral Research. https://www.nihms.nih.gov/pmc/articlerender.fcgi?artid=1520409
##' 
##' Chow, S-M., *Bendezu, J. J., Cole, P. M., & Ram, N. (2016). A Comparison of Two-
##' Stage Approaches for Fitting Nonlinear Ordinary Differential Equation (ODE) 
##' Models with Mixed Effects. Multivariate Behavioral Research, 51, 154-184. Doi: 10.1080/00273171.2015.1123138.
##' @examples
##' matt = plotGCV(theTimes,norder,roughPenaltyMax,out2,lambdaLow, 
##' lambdaHi,lambdaBy,isPlot)
plotGCV = function(theTimes,norder,roughPenaltyMax,dataMatrix,lowLambda,upLambda,lambdaInt,isPlot){
  rng =  range(theTimes)
  nt   =  length(theTimes)
  nbasis =  nt + norder - 2
  oscbasis = create.bspline.basis(rng, nbasis,norder, theTimes) #create a system of spline basis functions by inputting range, number of basis functions, order, knot #points
  OscfdPar = fdPar(oscbasis, roughPenaltyMax, 10^-5) #Use an arbitrary lambda value here to set up the fda object
  oscfd1 = smooth.basis(theTimes,dataMatrix,OscfdPar)
  #lamout2 = plotGCVRMSE.fd(lowLambda,upLambda, lambdaInt, theTimes, dataMatrix, OscfdPar,isPlot)
  lamvec = seq(lamlow, lamhi, lamdel)
  lamout = matrix(0,length(lamvec),3)
  m = 0
  for (lambda in lamvec) {
    m = m + 1
    lamout[m,1] = lambda
    fdParobj$lambda = lambda
    #Setting up the basis functions and other properties of the fda object
    smoothlist = smooth.basis(theTimes, y, fdParobj)
    #Extract the smoothed or approximated curve
    xfd = smoothlist$fd 
    # Get the degrees of freedom of the smoothing curve
    lamout[m,2] = smoothlist$df
    # Here I extract the mean of the N gcv values (e.g., across participants)
    lamout[m,3] =mean(smoothlist$gcv)  
  }
  if (isPlot==1){
    cat("lambda, deg. freedom, gcv\n")
    for (i in 1:m) {
      cat(format(round(lamout[i,],2)))
      cat("\n")
    }
    #par(mfrow=c(2,1))
    plot(lamvec, lamout[,2], type="b",xlab="lambda",ylab="df")
    title("Effective degrees of freedom")
    plot(lamvec, lamout[,3], type="b",xlab="lambda",ylab="Mean gcv")
    title("Mean GCV")
    par(mfrow=c(1,1))
  }
  colnames(lamout) = c("lambda","edf","GCV")
  return(data.frame(lamout2))
}



##' A function that takes the data and some Bsplines specification parameters and returns
##' smoothed estimated derivatives at the order specified by derivOrder
##' @param theTimes The time points at which derivative estimation are requested
##' @param norder #Order of Bsplines - usually 2 higher than roughPenaltyMax
##' @param roughPenaltyMax  Penalization order. Usually set to 2 higher than the highest-order derivatives desired
##' @param lambda = a positive smoothing parameter: larger --> more smoothing
##' @param dataMatrix Data of size total number of time points x  total number of subjects
##' @param derivOrder The order of the desired derivative estimates
##' 
##' @return a list containing:
##' 1. out (a matrix of derivative estimates); 2. basisCoef (basis coefficients); 3. basis 2 (basis functions)
##' 
##' @references
##' Chow, S-M. (2019). Practical Tools and Guidelines for Exploring and Fitting Linear 
##' and Nonlinear Dynamical Systems Models. Multivariate Behavioral Research. https://www.nihms.nih.gov/pmc/articlerender.fcgi?artid=1520409
##' 
##' Chow, S-M., *Bendezu, J. J., Cole, P. M., & Ram, N. (2016). A Comparison of Two-
##' Stage Approaches for Fitting Nonlinear Ordinary Differential Equation (ODE) 
##' Models with Mixed Effects. Multivariate Behavioral Research, 51, 154-184. Doi: 10.1080/00273171.2015.1123138.
##' @examples
##' x = getdx(theTimes,norder,roughPenaltyMax,sp,out2,0)[[1]] #Smoothed level
##' dx = getdx(theTimes,norder,roughPenaltyMax,sp,out2,1)[[1]] #Smoothed 1st derivs
##' d2x = getdx(theTimes,norder,roughPenaltyMax,sp,out2,2)[[1]] #Smoothed 2nd derivs
getdx <- function(theTimes,norder,roughPenaltyMax,lambda,dataMatrix,derivOrder){
  rng =  range(theTimes)
  nt   =  length(theTimes)
  nbasis =  nt + norder - 2
  oscbasis = create.bspline.basis(rng, nbasis,norder, theTimes) #create a system of #spline basis functions by inputting range, number of basis functions, order, knot #points
  OscfdPar = fdPar(oscbasis, roughPenaltyMax, lambda)
  oscfd1 = smooth.basis(theTimes,dataMatrix,OscfdPar)
  oscfd1.fd = oscfd1$fd
  basisCoef=coef(oscfd1.fd) #Get basis coefficients, c, for all b and i
  basis2 = coef(OscfdPar) #Get basis function matrix, Phi
  out = eval.fd(theTimes,oscfd1.fd,derivOrder)
  return(list(out,basisCoef,basis2))
}

