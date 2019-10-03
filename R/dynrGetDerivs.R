##' A function to evaluate the generalized cross-validation (GCV) values 
##' associated with derivative estimates via Bsplines at a range of specified 
##' smoothing parameter (lambda) values
##' @param theTimes The time points at which derivative estimation are requested
##' @param norder Order of Bsplines - usually 2 higher than roughPenaltyMax
##' @param roughPenaltyMax  Penalization order. Usually set to 2 higher than the highest-order derivatives desired
##' @param dataMatrix Data of size total number of time points x  total number of subjects
##' @param LowLambda Lower limit of lambda values to be tested. Here, lambda is a 
##' positive smoothing parameter, with larger values resulting in greater smoothing)
##' @param upLambda Upper limit of lambda
##' @param lambdaInt The interval of lambda values to be tested.
##' @param isPlot A binary flag on whether to plot the gcv values (0 = no, 1 = yes)
##' 
##' @return A data frame containing:
##' 1. lambda values; 2. edf (effective degrees of freedom); 
##' 3. GCV (Generalized cross-validation value as averaged across units (e.g., subjects))
##' 
##' @references
##' Chow, S-M. (2019). Practical Tools and Guidelines for Exploring and Fitting Linear 
##' and Nonlinear Dynamical Systems Models. Multivariate Behavioral Research. https://www.nihms.nih.gov/pmc/articlerender.fcgi?artid=1520409
##' 
##' Chow, S-M., *Bendezu, J. J., Cole, P. M., & Ram, N. (2016). A Comparison of Two-
##' Stage Approaches for Fitting Nonlinear Ordinary Differential Equation (ODE) 
##' Models with Mixed Effects. Multivariate Behavioral Research, 51, 154-184. Doi: 10.1080/00273171.2015.1123138.
##' @examples
##' #outMatrix = plotGCV(theTimes,norder,roughPenaltyMax,out2,lambdaLow, 
##' #lambdaHi,lambdaBy,isPlot)
plotGCV = function(theTimes,norder,roughPenaltyMax,dataMatrix,lowLambda,upLambda,lambdaInt,isPlot){
  rng =  range(theTimes)
  nt   =  length(theTimes)
  nbasis =  nt + norder - 2
  oscbasis = create.bspline.basis(rng, nbasis,norder, theTimes) #create a system of spline basis functions by inputting range, number of basis functions, order, knot #points
  fdParobj = fdPar(oscbasis, roughPenaltyMax, 10^-5) #Use an arbitrary lambda value here to set up the fda object
  #oscfd1 = smooth.basis(theTimes,dataMatrix,fdParobj)
  #lamout2 = plotGCVRMSE.fd(lowLambda,upLambda, lambdaInt, theTimes, dataMatrix, fdParobj,isPlot)
  lamvec = seq(lowLambda, upLambda, lambdaInt)
  lamout = matrix(0,length(lamvec),3)
  m = 0
  for (lambda in lamvec) {
    m = m + 1
    lamout[m,1] = lambda
    fdParobj$lambda = lambda
    #Setting up the basis functions and other properties of the fda object
    smoothlist = smooth.basis(theTimes, dataMatrix, fdParobj)
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
  return(data.frame(lamout))
}



##' A wrapper function to call functions in the fda package to obtain
##' smoothed estimated derivatives at a specified order
##' @param theTimes The time points at which derivative estimation are requested
##' @param norder Order of Bsplines - usually 2 higher than roughPenaltyMax
##' @param roughPenaltyMax  Penalization order. Usually set to 2 higher than the highest-order derivatives desired
##' @param lambda A positive smoothing parameter: larger --> more smoothing
##' @param dataMatrix Data of size total number of time points x  total number of subjects
##' @param derivOrder The order of the desired derivative estimates
##' 
##' @return A list containing:
##' 1. out (a matrix containing the derivative estimates
##' at the specified order that matches the dimension of dataMatrix); 
##' 2. basisCoef (estimated basis coefficients); 3. basis2 (basis functions)
##' 
##' @references
##' Chow, S-M. (2019). Practical Tools and Guidelines for Exploring and Fitting Linear 
##' and Nonlinear Dynamical Systems Models. Multivariate Behavioral Research. https://www.nihms.nih.gov/pmc/articlerender.fcgi?artid=1520409
##' 
##' Chow, S-M., *Bendezu, J. J., Cole, P. M., & Ram, N. (2016). A Comparison of Two-
##' Stage Approaches for Fitting Nonlinear Ordinary Differential Equation (ODE) 
##' Models with Mixed Effects. Multivariate Behavioral Research, 51, 154-184. Doi: 10.1080/00273171.2015.1123138.
##' @examples
##' #x = getdx(theTimes,norder,roughPenaltyMax,sp,out2,0)[[1]] #Smoothed level
##' #dx = getdx(theTimes,norder,roughPenaltyMax,sp,out2,1)[[1]] #Smoothed 1st derivs
##' #d2x = getdx(theTimes,norder,roughPenaltyMax,sp,out2,2)[[1]] #Smoothed 2nd derivs
getdx <- function(theTimes,norder,roughPenaltyMax,lambda,dataMatrix,derivOrder){
  rng =  range(theTimes)
  nt   =  length(theTimes)
  nbasis =  nt + norder - 2
  oscbasis = create.bspline.basis(rng, nbasis,norder, theTimes) #create a system of #spline basis functions by inputting range, number of basis functions, order, knot #points
  fdParobj = fdPar(oscbasis, roughPenaltyMax, lambda)
  oscfd1 = smooth.basis(theTimes,dataMatrix,fdParobj)
  oscfd1.fd = oscfd1$fd
  basisCoef=coef(oscfd1.fd) #Get basis coefficients, c, for all b and i
  basis2 = coef(fdParobj) #Get basis function matrix, Phi
  out = eval.fd(theTimes,oscfd1.fd,derivOrder)
  return(list(out,basisCoef,basis2))
}


#---- Functions by Marco Bachl to perform Johnson-Neyman and plots of simple slopes/region of significance ----
#at https://rpubs.com/bachl/jn-plot
#Function to implement the Johnson-Neyman technique
jnt <- function(.lm, predictor, moderator, alpha=.05) {
  b1 = coef(.lm)[predictor]
  b3 = coef(.lm)[stringi::stri_startswith_fixed(names(coef(.lm)), paste0(predictor,":")) | stringi::stri_endswith_fixed(names(coef(.lm)), paste0(":",predictor))]
  se_b1 = coef(summary(.lm))[predictor, 2]
  se_b3 = coef(summary(.lm))[stringi::stri_startswith_fixed(names(coef(.lm)), paste0(predictor,":")) | stringi::stri_endswith_fixed(names(coef(.lm)), paste0(":",predictor)), 2]
  print(paste0("b1 = ", b1))
  print(paste0("\nb3 = ", b3))
  print(paste0("\nseb1 = ", se_b1))
  print(paste0("seb3 = ", se_b3))
  COV_b1b3 = vcov(.lm)[predictor, stringi::stri_startswith_fixed(names(coef(.lm)), paste0(predictor,":")) | stringi::stri_endswith_fixed(names(coef(.lm)), paste0(":",predictor))]
  t_crit = qt(1-alpha/2, .lm$df.residual)
  # see Bauer & Curran, 2005
  a = t_crit^2 * se_b3^2 - b3^2
  b = 2 * (t_crit^2 * COV_b1b3 - b1 * b3)
  c = t_crit^2 * se_b1^2 - b1^2
  jn = c(
    (-b - sqrt(b^2 - 4 * a * c)) / (2 * a),
    (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)
  )
  JN = sort(unname(jn))
  JN = JN[JN>=min(.lm$model[,moderator]) & JN<=max(.lm$model[,moderator])]
  JN
}

##' A Function to plot simple slopes and region of significance.
##' @param .lm A regression object from running a linear model of the form: 
##' lm(y~ x1+x2+x1:x2), yielding: y = b0 + b1*x1 + b2*x2 + b3*x1*x2 + residual. 
##' In this case, one may rewrite the lm as y = b0 + (b1+b3*x2)*x1 + b2*x2 + residual,
##' where (b1+b3*x2) is referred to as the simple slope of x1, x1 is the predictor,
##' and x2 is the moderator whose values yield different simple slope values for x1. 
##' @param predictor The independent variable for which simple slope is requested
##' @param moderator The moderator whose values affect the simple slopes of the predictor. 
##' Appears on the horizontal axis.
##' @param alpha The designated alpha level for the Johnson-Neyman technique
##' @param jn A binary flag requesting the Johnson-Neyman test (T or F)
##' @param title0 Title for the plot
##' @param predictorLab Label for the predictor
##' @param moderatorLab Label for the moderator
##' 
##' @return A region of significance plot with simple slopes of the predictor on
##' the vertical axis, and values of the moderator on the horizontal axis.
##' 
##' @references 
##' Adapted from functions written by Marco Bachl to perform the Johnson-Neyman test 
##' and produce a plot of simple slopes and region of significance available at:
##' https://rpubs.com/bachl/jn-plot
##' 
##' @examples
##' # g = lm(y~x1:x2,data=data) 
##' # theta_plot(g, predictor = "x1", moderator = "x2", 
##' #           alpha = .05, jn = T, title0=" ",
##' #           predictorLab = "x1", moderatorLab = "x2")
theta_plot <- function(.lm, predictor, moderator, 
                       alpha=.05, jn=F,title0,
                       predictorLab, moderatorLab) {
  theme_set(theme_minimal())
  .data = tibble::tibble(b1 = coef(.lm)[predictor],
                     b3 = coef(.lm)[stringi::stri_startswith_fixed(names(coef(.lm)), paste0(predictor,":")) | stringi::stri_endswith_fixed(names(coef(.lm)), paste0(":",predictor))],
                     Z = quantile(.lm$model[,moderator], seq(0,1,.01)),
                     theta = b1 + Z * b3,
                     se_b1 = coef(summary(.lm))[predictor, 2],
                     COV_b1b3 = vcov(.lm)[predictor, stringi::stri_startswith_fixed(names(coef(.lm)), paste0(predictor,":")) | stringi::stri_endswith_fixed(names(coef(.lm)), paste0(":",predictor))],
                     se_b3 = coef(summary(.lm))[stringi::stri_startswith_fixed(names(coef(.lm)), paste0(predictor,":")) | stringi::stri_endswith_fixed(names(coef(.lm)), paste0(":",predictor)), 2],
                     se_theta = sqrt(se_b1^2 + 2 * Z * COV_b1b3 + Z^2 * se_b3^2),
                     ci.lo_theta = theta+qt(alpha/2, .lm$df.residual)*se_theta,
                     ci.hi_theta = theta+qt(1-alpha/2, .lm$df.residual)*se_theta)
  if (jn) {
    JN = jnt(.lm=.lm, predictor=predictor, moderator=moderator, alpha=alpha)
    JN_lines = geom_vline(xintercept=JN, linetype=2)
    JN_regions = ifelse(length(JN) == 0, 
                        "no significance regions or \nentire range of the moderator is in the significance region", 
                        paste(round(JN,2), collapse = "; "))
    Xlab = paste0(moderatorLab, " (JN Significance Regions: ", JN_regions,")")
  }
  else {
    Xlab = moderator
    JN_lines = NULL
  }
 # .data magrittr::%>%
    ggplot(.data,aes(Z, theta, ymin=ci.lo_theta, ymax=ci.hi_theta)) + 
    geom_ribbon(alpha = .2) + geom_line() + 
    ggtitle(paste(title0, "Simple slope of", predictorLab, 
                  "as function of", moderatorLab)) + 
    geom_hline(yintercept=0, linetype=2) + 
    labs(x = Xlab, y= "Simple slope") + 
    JN_lines +
    theme(axis.text=element_text(size=13),
          axis.title=element_text(size=14,face="bold"))
  
}

