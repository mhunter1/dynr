#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Illustration I: Linear oscillator example
# Author: Sy-Miin Chow
# Last modified: 9/7/2018
# The simulation model features:
# dx1(t)/dt = x2(t)
# dx2(t)/dt = eta*x1(t) + zeta*x2(t)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
require('dynr') 
options(scipen=999)

#Variables in this data set are: 
#ID = ID of participants
#x = true score 
#theTimes = time index

#Loading simulated data generated using the linear oscillator model
data("LinearOsc")
load('LinearOsc.rda')


#Structure the one indicator for n individuals into a matrix with n columns
#If different individuals have different number of rows, consider doing
#this step separately for each individual

n = 50 #Number of subjects
T = 100 #Number of time points
out2 = matrix(out1$x,ncol=n,byrow=FALSE)
theTimes = out1$theTimes[1:T]
norder = 6 #Order of Bsplines - usually 2 higher than roughPenaltyMax
roughPenaltyMax = 4 #penalization order 
#  #lambdaLow, lambdaHi, lambdaBy = specify an interval of lambda (a positive smoothing parameter,
#                                 larger means more smoothing) to be tested from lambdaLow to
#                                 lambdaHi, as separated by lambdaBy
#isPlot = a binary flag on whether to plot the gcv values (0 = no, 1 = yes)
lambdaLow = 1e-10
lambdaHi = 2
lambdaBy = .1
isPlot = 1 #Whether to plot the average GCV values as a function of lambda values
matt = plotGCV(theTimes,norder,roughPenaltyMax,out2,lambdaLow, lambdaHi,lambdaBy,isPlot)

#Extract the lambda value that gives the minimum GCV
sp = matt[matt[,"GCV"] == min(matt[,"GCV"]),"lambda"]

#Extract smoothed level, first and second derivative estimates at the lambda value selected above
x = getdx(theTimes,norder,roughPenaltyMax,sp,out2,0)[[1]] #Smoothed level
dx = getdx(theTimes,norder,roughPenaltyMax,sp,out2,1)[[1]] #Smoothed 1st derivs
d2x = getdx(theTimes,norder,roughPenaltyMax,sp,out2,2)[[1]] #Smoothed 2nd derivs

#Put level and derivative estimates into a data frame
dxall = data.frame(time = rep(theTimes,n),
                   x = matrix(x,ncol=1,byrow=FALSE), 
                   dx = matrix(dx,ncol=1,byrow=FALSE),
                   d2x = matrix(d2x,ncol=1,byrow=FALSE))


g = lm(d2x~x+dx-1,data=dxall)


#Component + residuals plot to show the association between smoothed d2x and smoothed x
#after partialling out the effect of smoothed dx
par(mgp=c(2.5,0.5,0))
car::crPlots(g,terms=~x,
        main=c(""),layout=c(1,1),
        cex.lab=1.3,cex.axis=1.2,
        xlab=expression(hat(eta)[i](t)),
        ylab=expression(paste("Component+Residuals ", "  ",d^2,hat(eta)[i](t)/dt^2))
)

#Component + residuals plot to show the association between smoothed d2x and smoothed dx
#after partialling out the effect of smoothed x
par(mgp=c(2.5,0.5,0))
car::crPlots(g,terms=~dx,
        main=c(""),layout=c(1,1),
        cex.lab=1.3,cex.axis=1.2,
        xlab=expression(paste(d, hat(eta)[i](t)/dt)),
        ylab=expression(paste("Component+Residuals ", "  ",d^2,hat(eta)[i](t)/dt^2))
)

