#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2017-09-27
# Filename: fern.R
# Purpose: Create data from a Barnsley fern, add noise, and fit a dynr model
#  to them.
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------

fern <- function(x){
	p <- c(.02, .15, .13, .70)
	A1 <- matrix(c(0, 0, 0, .27), 2, 2)
	B1 <- matrix(c(.5, 0), 2, 1)
	A2 <- matrix(c(-.139, .246, .263, .224), 2, 2)
	B2 <- matrix(c(.57, -.036), 2, 1)
	A3 <- matrix(c(.17, .222, -.215, .176), 2, 2)
	B3 <- matrix(c(.408, .0893), 2, 1)
	A4 <- matrix(c(.781, -.032, .034, .739), 2, 2)
	B4 <- matrix(c(.1075, .27), 2, 1)
	u <- matrix(1, 1, 1)
	AL <- list(A1, A2, A3, A4)
	BL <- list(B1, B2, B3, B4)
	r <- sample(size=1, x=1:4, prob=p)
	AL[[r]]%*%x + BL[[r]]%*%u
}

fern2 <- function(x){
	p <- c(.02, .15, .13, .70)
	#p <- c(.01, .07, .07, .85)
	A1 <- matrix(c(0, 0, 0, .16), 2, 2)
	B1 <- matrix(c(0, 0), 2, 1)
	A2 <- matrix(c(.20, .23, -.26, .22), 2, 2)
	B2 <- matrix(c(0, 1.6), 2, 1)
	A3 <- matrix(c(-.15, .26, .28, .24), 2, 2)
	B3 <- matrix(c(0, .44), 2, 1)
	A4 <- matrix(c(.85, -.04, .04, .85), 2, 2)
	B4 <- matrix(c(0, 1.6), 2, 1)
	u <- matrix(1, 1, 1)
	AL <- list(A1, A2, A3, A4)
	BL <- list(B1, B2, B3, B4)
	r <- sample(size=1, x=1:4, prob=p)
	AL[[r]]%*%x + BL[[r]]%*%u
}


tree <- function(x){
	r <- runif(1)
	y <- c(0, 0)
	if(r < .10){
		y <- c(.05*x[1], .60*x[2])
	} else if(r < .20){
		y <- c(.05*x[1],
			-.5*x[2]+1.0)
	} else if(r < .40){
		y <- c(.46*x[1]-.15*x[2],
			.39*x[1]+.38*x[2]+.6)
	} else if(r < .60){
		y <- c(.47*x[1]-.15*x[2],
			.17*x[1]+.42*x[2]+1.1)
	} else if(r < .80){
		y <- c(.43*x[1]+.28*x[2],
			-.25*x[1]+.45*x[2]+1.0)
	} else {
		y <- c(.42*x[1]+.26*x[2],
			-.35*x[1]+.31*x[2]+.7)
	}
	return(y)
}

#http://mathworld.wolfram.com/BarnsleysTree.html
btree <- function(z){
	c <- complex(real=0.6, imaginary=1.1)
	c*(z - sign(Re(z)))
}


nIter <- 10000 #30000
x <- matrix(0, nIter, 2)


x[1,] <-  c(0.5, 0) #runif(2)

for(i in 1:(nIter-1)){
	x[i+1,] <- fern2(x[i,]) # or fern2
}

plot(x, pch='.', xlim=range(x), ylim=range(x))


# 25% leaves a blob
# Add 5% Gaussian White noise
maxRange <- max(apply(apply(x, 2, range), 2, diff))
pctNoise <- matrix(rnorm(nIter*2, mean=0, sd=0.05*maxRange/3), nIter, 2)
# Plus/minus 3 SDs is about the range of the noise
xn <- x + pctNoise
plot(xn, pch='.', xlim=range(xn), ylim=range(xn))

colnames(xn) <- c('x1', 'x2')
xn <- cbind(id=1, time=1:nrow(xn), xn)



#------------------------------------------------------------------------------
# Build dynr model

require(dynr)

da <- dynr.data(xn, observed=c('x1', 'x2'))

lo <- prep.loadings(list(eta1='x1', eta2='x2'), params=c())
#no <- togo.noise('zero', 2, '', 'diagonal', 2, '')
no <- prep.noise(matrix(0, 2, 2), matrix(0, 2, 2), diag(.2, 2), diag(c('n1', 'n2'), 2))
nocon <- prep.noise(matrix(0, 2, 2), matrix(0, 2, 2), diag(.2, 2), diag('mnoise', 2))
dlab <- outer(paste0('m', 1:2), 1:2, paste0)
dlab <- lapply(1:4, function(x) matrix(paste0(dlab, '_', x), 2, 2))
dval <- lapply(1:4, function(x) diag(x/10, nrow=2))
ilab <- lapply(1:4, function(x) matrix(paste0(paste0('i', 1:2), '_', x), 2, 1))
ival <- lapply(1:4, function(x) matrix(x/10, nrow=2, ncol=1))
ilabcon <- lapply(1:4, function(x) matrix(c('fixed', paste0(paste0('i', 2), '_', x)), 2, 1))
ivalcon <- lapply(1:4, function(x) matrix(c(0, x/10), nrow=2, ncol=1))
dy <- prep.matrixDynamics(
	params.dyn=dlab,
	values.dyn=dval,
	params.int=ilab,
	values.int=ival,
	isContinuousTime=FALSE)
dycon <- prep.matrixDynamics(
	params.dyn=dlab,
	values.dyn=dval,
	params.int=ilabcon,
	values.int=ivalcon,
	isContinuousTime=FALSE)
re <- prep.regimes(
	values=matrix(c(0, 1, 1, 1,
					0, 1, 1, 1,
					0, 1, 1, 1,
					0, 1, 1, 1), 4, 4, byrow=TRUE),
	params=matrix(c('fixed', 'g12', 'g13', 'g14',
					'fixed', 'g22', 'g23', 'g24',
					'fixed', 'g32', 'g33', 'g34',
					'fixed', 'g42', 'g43', 'g44'), 4, 4, byrow=TRUE))
recon <- prep.regimes(
	values=matrix(c(0, 1, 1, 1,
					0, 1, 1, 1,
					0, 1, 1, 1,
					0, 1, 1, 1), 4, 4, byrow=TRUE),
	params=matrix(c('fixed', 'p2', 'p3', 'p4',
					'fixed', 'p2', 'p3', 'p4',
					'fixed', 'p2', 'p3', 'p4',
					'fixed', 'p2', 'p3', 'p4'), 4, 4, byrow=TRUE))
ic <- prep.initial(values.inistate=c(.5, 0), params.inistate=c('fixed', 'fixed'), values.inicov=diag(1, 2), params.inicov=diag('fixed', 2), values.regimep=c(0, 1, 1, 1), params.regimep=c('fixed', 'p2', 'p3', 'p4'))

m <- dynr.model(dynamics=dy, measurement=lo, noise=no, regimes=re, initial=ic, data=da, outfile='fern.c')

nr <- 4
m$lb <- c(rep(-2, 4*nr), rep(-10, 2*nr), rep(NA, 2), rep(NA, nr*nr-nr + nr-1))
m$ub <- c(rep(2, 4*nr), rep(10, 2*nr), rep(NA, 2), rep(NA, nr*nr-nr + nr-1))

#mc <- dynr.cook(m)



mcon <- dynr.model(dynamics=dycon, measurement=lo, noise=nocon, regimes=recon, initial=ic, data=da, outfile='fernCon.c')

nr <- 4
mcon$lb <- c(rep(-2, 4*nr), rep(-5, 1*nr), rep(NA, 1), rep(-6, nr-1))
mcon$ub <- c(rep(2, 4*nr), rep(5, 1*nr), rep(NA, 1), rep(15, nr-1))


mconc <- dynr.cook(mcon)




#------------------------------------------------------------------------------
# Serve the results

plot(mconc$eta_smooth_final[1,], mconc$eta_smooth_final[2,], pch='.', xlim=range(mconc$eta_smooth_final), ylim=range(mconc$eta_smooth_final))



pdf('plotFern.pdf', height=12*2.5, width=9*2.5)
plot(xn[,c('x1', 'x2')], pch='.', xlim=range(xn[,c('x1')]), ylim=range(xn[,c('x2')]), col='grey', xlab='', ylab='', xaxt='n', yaxt='n', bty='n', cex=.1)
points(x[,1], x[,2], col=89, pch='.', cex=.2)
points(mconc$eta_smooth_final[1,], mconc$eta_smooth_final[2,], pch='.', col=51, cex=.2)
dev.off()

printex(mcon,ParameterAs=coef(mconc), printInit=TRUE, printRS=TRUE,
        outFile="fernTree2.tex")
tools::texi2pdf("fernTree2.tex")
system(paste(getOption("pdfviewer"), "fernTree2.pdf"))


#------------------------------------------------------------------------------
# TODO add coefficient checking against true parameters in fern2


#------------------------------------------------------------------------------
# End
