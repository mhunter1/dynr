##' A function to evaluate the generalized cross-validation (GCV) values 
##' associated with derivative estimates via Bsplines at a range of specified 
##' smoothing parameter (lambda) values
##' @param theTimes The time points at which derivative estimation are requested
##' @param norder Order of Bsplines - usually 2 higher than roughPenaltyMax
##' @param roughPenaltyMax  Penalization order. Usually set to 2 higher than the highest-order derivatives desired
##' @param dataMatrix Data of size total number of time points x  total number of subjects
##' @param lowLambda Lower limit of lambda values to be tested. Here, lambda is a 
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
plotGCV = function(theTimes,norder,roughPenaltyMax,dataMatrix,lowLambda,upLambda,lambdaInt,isPlot){
  # Pseudo-example from vignette
  # outMatrix = plotGCV(theTimes,norder,roughPenaltyMax,out2,lambdaLow, 
  # lambdaHi,lambdaBy,isPlot)
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
##' 
##' @examples
##' data("LinearOsc")
##' # Number of subjects is 10
##' numP <- length(unique(LinearOsc$ID))
##' # Number of time points is 100
##' numT <- max(table(LinearOsc$ID))
##' out2 <- matrix(LinearOsc$x, ncol=numP, byrow=FALSE)
##' theTimes <- LinearOsc$theTimes[1:numT]
##' # Order of Bsplines - usually 2 higher than roughPenaltyMax
##' norder <- 6
##' # Penalization order
##' roughPenaltyMax <- 4 
##' # Pick lambda value that gives the low GCV
##' # Could/should use plotGCV instead
##' sp <- 1/2
##' # Smoothed level
##' x <- getdx(theTimes, norder, roughPenaltyMax, sp, out2, 0)[[1]]
##' # Smoothed 1st derivs
##' dx <- getdx(theTimes, norder, roughPenaltyMax, sp, out2, 1)[[1]]
##' # Smoothed 2nd derivs
##' d2x = getdx(theTimes, norder, roughPenaltyMax, sp, out2, 2)[[1]]
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

##' A function to plot simple slopes and region of significance.
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
theta_plot <- function(.lm, predictor, moderator, 
                       alpha=.05, jn=F,title0,
                       predictorLab, moderatorLab) {
  # Here is one pseudo-example taken from the vignette
  # Not to be confused with a real example ##' @examples
  # g = lm(y~x1:x2,data=data) 
  # theta_plot(g, predictor = "x1", moderator = "x2", 
  #           alpha = .05, jn = T, title0=" ",
  #           predictorLab = "x1", moderatorLab = "x2")
  theme_set(theme_minimal())
  b1 = NULL; b3 = NULL; Z = NULL; theta = NULL; se_b1 = NULL; COV_b1b3 = NULL
  se_b3 = NULL; se_theta = NULL; ci.lo_theta=NULL; ci.hi_theta = NULL
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

##' A Function to plot the flow or velocity field for a one or two dimensional 
##' autonomous ODE system from the phaseR package written by Michael J. Grayling. 
##' @note The phaseR package was taken off cran as off 10/1/2019 so we are 
##' exporting some selected functions from phaseR_2.0 published on 8/20/2018.
##' For details of these functions please see original documentations on the
##' phaseR package.
##' 
##' @param deriv A function computing the derivative at a point for the ODE system 
##' to be analysed. For examples see the phaseR package guide.
##' @param xlim A vector of length two setting the lower and upper limits of the variable 
##' to be plotted on the horizontal axis (usually the first variable returned by the function deriv)
##' @param ylim A vector of length two setting the lower and upper limits of the variable 
##' to be plotted on the vertical axis (usually the second variable returned by the function deriv)
##' @param parameters Parameters of the ODE system, to be passed to deriv. 
##' Supplied as a vector; the order of the parameters can be found from the 
##' deriv file. Defaults to NULL.
##' @param system Set to either "one.dim" or "two.dim" to indicate the type of system 
##' being analysed. Defaults to "two.dim".
##' @param points Sets the density of the line segments to be plotted. Defaults to 11.
##' @param col Sets the color of the plotted line segments. Defaults to "gray". 
##' Should be a vector of length one. Will be reset accordingly if it is a vector of the wrong length. 
##' @param arrow.type Sets the type of line segments plotted. Options include:
##' "proportional" = the length of the line segments reflects the magnitude of the derivative. 
##' "equal" the line segments take equal lengths, simply reflecting the gradient of 
##' the derivative(s). Defaults to "equal".
##' @param arrow.head Sets the length of the arrow heads. Passed to arrows. Defaults to 0.05.
##' @param frac Sets the fraction of the theoretical maximum length line segments can take 
##' without overlapping, that they can actually attain. In practice, frac can be set to greater 
##' than 1 without line segments overlapping. 
##' @param add Logical. Defaults to TRUE.
##' TRUE = the flow field is added to an existing plot; FALSE = a new plot is created. 
##' @param xlab Label for the x-axis of the resulting plot. Defaults to "x".
##' @param ylab	Label for the y-axis of the resulting plot. Defaults to "y".
##' @param state.names State names for ode functions that do not use positional states
##' @param ... Additional arguments to be passed to either plot or arrows.
##' 
##' @return Returns a list with the following components:
##' add, arrow.head, arrow.type, col, deriv, dx, dy, frac, parameters, points,
##' system, x, xlab, xlim, y, ylab, ylim. Most of these components correspond simply
##' to their original input values. 
##' 
##' The only new elements are:
##' 
##' dx = A matrix. In the case of a two dimensional system, the values of the 
##' derivative of the first dependent derivative at all evaluated points.
##' 
##' dy = A matrix. In the case of a two dimensional system, the values of the derivative of the second dependent variable at all evaluated points. In the case of a one dimensional system, the values of the derivative of the dependent variable at all evaluated points.
##'
##' x	= A vector. In the case of a two dimensional system, the values of the first dependent variable at which the derivatives were computed. In the case of a one dimensional system, the values of the independent variable at which the derivatives were computed.
##'
##' y	= A vector. In the case of a two dimensional system, the values of the second dependent variable at which the derivatives were computed. In the case of a one dimensional system, the values of the dependent variable at which the derivatives were computed.
##'
##'  
##' @references 
##' Grayling, Michael J. (2014). phaseR: An R Package for Phase Plane Analysis of Autonomous
##' ODE Systems. The R Journal, 6(2), 43-51. DOI: 10.32614/RJ-2014-023. Available at
##' https://doi.org/10.32614/RJ-2014-023
dynr.flowField <- function (deriv, xlim, ylim, parameters = NULL, system = "two.dim", 
          points = 21, col = "gray", arrow.type = "equal", arrow.head = 0.05, 
          frac = 1, add = TRUE, xlab = "x", ylab = "y", state.names = c("x", 
                                                                        "y"), ...) {
  # Pseudo-example from vignette
  #Osc <- function(t, y, parameters) {
  #  dy <- numeric(2)
  #  dy[1] <- y[2]
  #  dy[2] <- parameters[1]*y[1]+parameters[2]*dy[1]   
  #  return(list(dy))
  #}
  #
  #param <- coef(g)
  #dynr.flowField(Osc, xlim = c(-3, 3), 
  #                  ylim = c(-3, 3),
  #                  xlab="x", ylab="dx/dt",
  #                  main=paste0("Oscillator model"),
  #                  cex.main=2,
  #                  parameters = param, 
  #                  points = 15, add = FALSE,
  # col="blue",
  # arrow.type="proportional",
  # arrow.head=.05)  
  if ((!is.vector(xlim)) | (length(xlim) != 2)) {
    stop("xlim is not a vector of length 2 as required")
  }
  if (xlim[2] <= xlim[1]) {
    stop("xlim[2] is less than or equal to xlim[1]")
  }
  if ((!is.vector(ylim)) | (length(ylim) != 2)) {
    stop("ylim is not a vector of length 2 as required")
  }
  if (ylim[2] <= ylim[1]) {
    stop("ylim[2] is less than or equal to ylim[1]")
  }
  if (points <= 0) {
    stop("points is less than or equal to zero")
  }
  if (!(system %in% c("one.dim", "two.dim"))) {
    stop("system must either be set to \"one.dim\" or \"two.dim\"")
  }
  if (is.vector(col) == FALSE) {
    stop("col is not a vector as required")
  }
  if (length(col) > 1) {
    col <- col[1]
    print("Note: col has been reset as required")
  }
  if (!(arrow.type %in% c("proportional", "equal"))) {
    stop("arrow.type must either be set to \"proportional\" or \"equal\"")
  }
  if (arrow.head <= 0) {
    stop("arrow.head is less than or equal to zero")
  }
  if (frac <= 0) {
    stop("frac is less than or equal to zero")
  }
  if (!is.logical(add)) {
    stop("add must be logical")
  }
  x <- seq(from = xlim[1], to = xlim[2], length = points)
  y <- seq(from = ylim[1], to = ylim[2], length = points)
  dx <- matrix(0, ncol = points, nrow = points)
  dy <- matrix(0, ncol = points, nrow = points)
  xmax.length <- x[2] - x[1]
  ymax.length <- y[2] - y[1]
  if (add == FALSE) {
    plot(1, xlim = c(xlim[1] - xmax.length, xlim[2] + xmax.length), 
         ylim = c(ylim[1] - ymax.length, ylim[2] + ymax.length), 
         type = "n", xlab = xlab, ylab = ylab, ...)
  }
  if (system == "one.dim") {
    for (i in 1:points) {
      dy[1, i] <- deriv(0, setNames(c(y[i]), state.names[1]), 
                        parameters)[[1]]
    }
    for (i in 2:points) {
      dy[i, ] <- dy[1, ]
    }
    abs.dy <- abs(dy)
    abs.dy.non <- abs.dy[which(abs.dy != 0)]
    max.abs.dy <- max(abs(dy))
    coefficient <- frac * min(xmax.length, ymax.length)/(2 * 
                                                           sqrt(2) * max(sqrt(2 * abs.dy.non/(abs.dy.non + (1/abs.dy.non))), 
                                                                         sqrt(2 * (1/abs.dy.non)/(abs.dy.non + (1/abs.dy.non)))))
    for (i in 1:points) {
      for (j in 1:points) {
        if (dy[i, j] != 0) {
          factor <- sqrt(2/(abs.dy[i, j] + (1/abs.dy[i, 
                                                     j])))
          y.shift <- coefficient * factor * sqrt(abs.dy[i, 
                                                        j])
          x.shift <- coefficient * factor/sqrt(abs.dy[i, 
                                                      j])
          if (dy[i, j] < 0) {
            y.shift <- -y.shift
          }
        }
        if (dy[i, j] == 0) {
          y.shift <- 0
          x.shift <- coefficient * sqrt(2)
        }
        if (arrow.type == "proportional") {
          if (dy[i, j] != 0) {
            prop <- abs.dy[i, j]/max.abs.dy
            y.shift <- y.shift * prop
            x.shift <- x.shift * prop
          }
          if (dy[i, j] == 0) {
            x.shift <- y.shift * mean(abs.dy)/max.abs.dy
          }
        }
        arrows(x[i] - x.shift, y[j] - y.shift, x[i] + 
                 x.shift, y[j] + y.shift, length = arrow.head, 
               col = col, ...)
      }
    }
    return(list(add = add, arrow.head = arrow.head, arrow.type = arrow.type, 
                col = col, deriv = deriv, dy = dy, frac = frac, parameters = parameters, 
                points = points, system = system, x = x, xlab = xlab, 
                xlim = xlim, y = y, ylab = ylab, ylim = ylim))
  }
  else {
    for (i in 1:length(x)) {
      for (j in 1:length(y)) {
        df <- deriv(0, setNames(c(x[i], y[j]), state.names), 
                    parameters)
        dx[i, j] <- df[[1]][1]
        dy[i, j] <- df[[1]][2]
      }
    }
    abs.dx <- abs(dx)
    abs.dy <- abs(dy)
    abs.dx.non <- abs.dx[which((abs.dx != 0) & (abs.dy != 
                                                  0))]
    abs.dy.non <- abs.dy[which((abs.dx != 0) & (abs.dy != 
                                                  0))]
    max.length <- max(sqrt(dx^2 + dy^2))
    coefficient <- frac * min(xmax.length, ymax.length)/(2 * 
                                                           sqrt(2) * max(sqrt(2 * (abs.dy.non/abs.dx.non)/((abs.dy.non/abs.dx.non) + 
                                                                                                             (abs.dx.non/abs.dy.non))), sqrt(2 * (abs.dx.non/abs.dy.non)/((abs.dy.non/abs.dx.non) + 
                                                                                                                                                                            (abs.dx.non/abs.dy.non)))))
    for (i in 1:points) {
      for (j in 1:points) {
        if ((dx[i, j] != 0) | (dy[i, j] != 0)) {
          if ((dx[i, j] != 0) & (dy[i, j] != 0)) {
            factor <- sqrt(2/((abs.dy[i, j]/abs.dx[i, 
                                                   j]) + (abs.dx[i, j]/abs.dy[i, j])))
            y.shift <- coefficient * factor * sqrt(abs.dy[i, 
                                                          j]/abs.dx[i, j])
            x.shift <- coefficient * factor/sqrt(abs.dy[i, 
                                                        j]/abs.dx[i, j])
            if (dy[i, j] < 0) {
              y.shift <- -abs(y.shift)
            }
            if (dx[i, j] < 0) {
              x.shift <- -abs(x.shift)
            }
          }
          if ((dx[i, j] == 0) & (dy[i, j] != 0)) {
            y.shift <- coefficient * sqrt(2)
            x.shift <- 0
            if (dy[i, j] < 0) {
              y.shift <- -abs(y.shift)
            }
          }
          if ((dx[i, j] != 0) & (dy[i, j] == 0)) {
            y.shift <- 0
            x.shift <- coefficient * sqrt(2)
            if (dx[i, j] < 0) {
              x.shift <- -abs(x.shift)
            }
          }
          if (arrow.type == "proportional") {
            prop <- sqrt((abs.dx[i, j]^2 + abs.dy[i, 
                                                  j]^2))/max.length
            y.shift <- y.shift * prop
            x.shift <- x.shift * prop
          }
          arrows(x[i] - x.shift, y[j] - y.shift, x[i] + 
                   x.shift, y[j] + y.shift, length = arrow.head, 
                 col = col, ...)
        }
      }
    }
  }
  return(list(add = add, arrow.head = arrow.head, arrow.type = arrow.type, 
              col = col, deriv = deriv, dx = dx, dy = dy, frac = frac, 
              parameters = parameters, points = points, system = system, 
              x = x, xlab = xlab, xlim = xlim, y = y, ylab = ylab, 
              ylim = ylim))
}


##' A Function to perform numerical integration of the chosen ODE system, 
##' for a user-specified set of initial conditions. Plots the resulting solution(s) 
##' in the phase plane. This function from the phaseR package written by 
##' Michael J. Grayling. 
##' @note The phaseR package was taken off cran as off 10/1/2019 so we are 
##' exporting some selected functions from phaseR_2.0 published on 8/20/2018.
##' For details of these functions please see original documentations on the
##' phaseR package.
##' 
##' @param deriv A function computing the derivative at a point for the specified 
##' ODE system. See the phaseR package guide for more examples.
##' @param y0	The initial condition(s) (ICs). In one-dimensional system, this can 
##' either be a single number indicating a single IC or a vector indicating multiple 
##' ICs. In two-dimensional system, this can either be a vector of length two 
##' reflecting the location of the two dependent variables initially, 
##' or it can be matrix where each row reflects a different set of ICs. 
##' Alternatively this can be left blank and the user can use locator to specify initial condition(s) on a plot. In this case, for one dimensional systems, all initial conditions are taken at tlim[1], even if not selected so on the graph. Defaults to NULL.
##' @param n If y0 is left NULL so initial conditions can be specified using 
##' locator, n sets the number of initial conditions to be chosen. Defaults to NULL.
##' @param tlim Sets the limits of the independent variable for which the solution 
##' should be plotted. Should be a vector of length two. If tlim[2] > tlim[1], 
##' then tstep should be negative to indicate a backwards trajectory.
##' @param tstep The step length of the independent variable, used 
##' in numerical integration. Defaults to 0.01.
##' @param parameters	Parameters of the ODE system, to be passed to deriv. 
##' Supplied as a vector; the order of the parameters can be found from the deriv file. 
##' Defaults to NULL.
##' @param system	Set to either "one.dim" or "two.dim" to indicate the type of 
##' system being analysed. Defaults to "two.dim".
##' @param col The color(s) to plot the trajectories in. Will be reset accordingly if it is a vector not of the length of the number of initial conditions. Defaults to "black".
##' @param add Logical.  Defaults to TRUE.
##' TRUE = the trajectories added to an existing plot; FALSE = a new plot is created.
##' @param state.names State names for the ODE functions that do not use positional states
##' @param ... Additional arguments to be passed to either plot or arrows.
##' 
##' @return Returns a list with the following components:
##' add, col, deriv, n, parameters, system, tlim, tstep, t, x, y, ylab, y0. 
##' Most of these components correspond simply to their original input values. 
##' 
##' The only new elements are:
##' t = A vector containing the values of the independent variable at each integration step.
##'
##' x	= In the two dimensional system case, a matrix whose columns are the 
##' numerically computed values of the first dependent variable for each set of ICs.
##'
##' y	= In the two dimensional system case, a matrix whose columns are the numerically computed values of the second dependent variable for each initial condition. In the one dimensional system case, a matrix whose columns are the numerically computed values of the dependent variable for each initial condition.
##'
##' y0	= As per input, but converted to a matrix if supplied as a vector initially.
##' 
##' 
##' @references 
##' Grayling, Michael J. (2014). phaseR: An R Package for Phase Plane Analysis of Autonomous
##' ODE Systems. The R Journal, 6(2), 43-51. DOI: 10.32614/RJ-2014-023. Available at
##' https://doi.org/10.32614/RJ-2014-023
dynr.trajectory <- function (deriv, y0 = NULL, n = NULL, tlim, tstep = 0.01, parameters = NULL, 
          system = "two.dim", col = "black", add = TRUE, state.names = c("x", 
                                                                         "y"), ...) {
  # Pseudo-example from vignette
  #Osc <- function(t, y, parameters) {
  #  dy <- numeric(2)
  #  dy[1] <- y[2]
  #  dy[2] <- parameters[1]*y[1]+parameters[2]*dy[1]   
  #  return(list(dy))
  #}
  #
  #param <- coef(g)
  #dynr.flowField(Osc, xlim = c(-3, 3), 
  #                  ylim = c(-3, 3),
  #                  xlab="x", ylab="dx/dt",
  #                  main=paste0("Oscillator model"),
  #                  cex.main=2,
  #                  parameters = param,
  #                  points = 15, add = FALSE,
  # col="blue",
  # arrow.type="proportional",
  # arrow.head=.05)  
  #IC <- matrix(c(-2, -2), ncol = 2, byrow = TRUE)  #Initial conditions
  # phaseR::trajectory(Osc, y0 = IC, parameters = param,tlim=c(0,10))
  if (tstep == 0) {
    stop("tstep is equal to 0")
  }
  if (tlim[1] == tlim[2]) {
    stop("tlim[1] is equal to tlim[2]")
  }
  if ((tlim[1] > tlim[2]) & (tstep > 0)) {
    stop("tstep must be negative if tlim[1] > tlim[2]")
  }
  if ((tlim[1] < tlim[2]) & (tstep < 0)) {
    stop("tstep must be positive if tlim[1] < tlim[2]")
  }
  if (!(system %in% c("one.dim", "two.dim"))) {
    stop("system must either be set to one.dim or two.dim")
  }
  if (!is.vector(col)) {
    stop("col is not a vector as required")
  }
  if (!is.logical(add)) {
    stop("add must be logical")
  }
  if (is.null(y0) & is.null(n)) {
    stop(paste("Both y0 and n cannot be NULL"))
  }
  if (!is.null(y0) & !is.null(n)) {
    warning("n is non-NULL whilst y0 has also been specified")
  }
  if (is.null(y0) & (add == FALSE)) {
    stop(paste("y0 cannot be null and add set to FALSE"))
  }
  if (is.null(y0)) {
    y0 <- locator(n = n)
    if (system == "two.dim") {
      re.set <- matrix(0, ncol = 2, nrow = n)
      for (i in 1:n) {
        re.set[i, ] <- c(y0$x[i], y0$y[i])
      }
      y0 <- re.set
    }
    if (system == "one.dim") {
      re.set <- numeric(n)
      for (i in 1:n) {
        re.set[i] <- y0$y[i]
      }
      y0 <- re.set
    }
  }
  if ((!is.vector(y0)) & (!is.matrix(y0))) {
    stop("y0 is neither a number, vector or matrix as required")
  }
  if (is.vector(y0)) {
    y0 <- as.matrix(y0)
  }
  if ((system == "one.dim") & (all(dim(y0) > 1))) {
    stop("For system equal to \"one.dim\" y0 must contain either a vector or a matrix where either nrow(y0) or ncol(y0) is one")
  }
  if ((system == "two.dim") & (!any(dim(y0) == 2))) {
    stop("For system equal to \"two.dim\" y0 must contain either a vector of length two or a matrix where either nrow(y0) or ncol(y0) is two")
  }
  if (system == "one.dim") {
    if (ncol(y0) > nrow(y0)) {
      y0 <- t(y0)
    }
    state.names <- state.names[1]
  }
  else {
    if ((nrow(y0) == 2) & (ncol(y0) != 2)) {
      y0 <- t(y0)
    }
  }
  if (nrow(y0) > length(col)) {
    col <- rep(col, nrow(y0))
    message("Note: col has been reset as required")
  }
  else if (nrow(y0) < length(col)) {
    col <- col[1:nrow(y0)]
    message("Note: col has been reset as required")
  }
  t <- seq(from = tlim[1], to = tlim[2], by = tstep)
  x <- matrix(0, nrow = length(t), ncol = nrow(y0))
  if (system == "two.dim") {
    y <- matrix(0, nrow = length(t), ncol = nrow(y0))
  }
  method <- ifelse(tstep > 0, "ode45", "lsoda")
  for (i in 1:nrow(y0)) {
    phase.trajectory <- deSolve::ode(times = t, y = setNames(c(y0[i, 
                                                         ]), state.names), func = deriv, parms = parameters, 
                            method = method)
    x[, i] <- phase.trajectory[, 2]
    if (system == "two.dim") {
      y[, i] <- phase.trajectory[, 3]
    }
    if ((add == FALSE) & (i == 1)) {
      if (system == "one.dim") {
        plot(t, x[, i], col = col[i], type = "l", ...)
      }
      else {
        plot(x[, i], y[, i], col = col[i], type = "l", 
             ...)
      }
    }
    else {
      if (system == "one.dim") {
        lines(t, x[, i], col = col[i], type = "l", ...)
      }
      else {
        lines(x[, i], y[, i], col = col[i], type = "l", 
              ...)
      }
    }
  }
  if (system == "one.dim") {
    points(rep(tlim[1], nrow(y0)), y0, col = col, ...)
    return(list(add = add, col = col, deriv = deriv, n = n, 
                parameters = parameters, system = system, t = t, 
                tlim = tlim, tstep = tstep, y = x, y0 = y0))
  }
  else {
    points(y0[, 1], y0[, 2], col = col, ...)
    return(list(add = add, col = col, deriv = deriv, n = n, 
                parameters = parameters, system = system, t = t, 
                tlim = tlim, tstep = tstep, x = x, y = y, y0 = y0))
  }
}
