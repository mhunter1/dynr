options(digits = 12)
library(ipoptr)
#------use .C()-------
dyn.load(file.path("~/Dropbox/Projects/SUSP15_IPOPT/CDAEKFR",paste("ipoptR",.Platform$dynlib.ext,sep="")))
dynamic<-function(n=10000,thin=500,seed=trunc(runif(1)*1e6))
{
  tmp=.C("gibbs",as.integer(n),as.integer(thin),
         as.integer(seed),x=as.double(1:n),
         y=as.double(1:n))
  mat=cbind(1:n,tmp$x,tmp$y) 
  colnames(mat)=c("Iter","x","y")
  mat
}

eval_f <- function(x) {
  tmp=.C("brekfis_obj_R",negloglike=as.double(0),params=as.double(x))
  return(tmp$negloglike)
}
## Gradient of Rosenbrock Banana function
eval_grad_f <- function(x) {
  return( c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
             200 * (x[2] - x[1] * x[1]) ) )
}
# initial values
x0 <- c(-0.1,.1,.1,1,1,2.9,-1,-1,-1,-.7)
#library(pryr)
#sexp_type(x0)
eval_f(x0)

res <- ipoptr( x0=x0,p
               eval_f=eval_f,
               eval_grad_f=eval_grad_f )
print( res )



#------use .Call()-------
dyn.load(file.path("~/Dropbox/Projects/SUSP15_IPOPT/CDAEKFR",paste("ipoptRcall",.Platform$dynlib.ext,sep="")))
eval_f <- function(x) {
  tmp=.Call("brekfis_obj_RcallC",x)
  return(tmp)
}
eval_f(x0)
