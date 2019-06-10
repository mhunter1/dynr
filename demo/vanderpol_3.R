#------------------------------------------------------------------------------
# Author: Hui-Ju Hung
# Date: 2019-05-23
# Filename: vanderpol_3.R
# Purpose: Moel script for dynr/SAEM gateway function
# Note: Workable for developer dynr on arma branch
#------------------------------------------------------------------------------

#---------
#function for parsing theta.formula (remove intercept terms and random names
prep.thetaFormula <- function(formula, intercept.names, random.names){
    
    
    
    #fun
    fml=lapply(formula, as.character)
    lhs=lapply(fml,function(x){x[[2]]})
    rhs=lapply(fml,function(x){x[[3]]})
    
    
    for(i in 1:length(formula)){
        formula[[i]]=as.character(formula[[i]])
        for (j in 1:length(intercept.names)){
            # gsub (a, b, c) : in c replace a with b
            rhs[[i]]=gsub(paste0(intercept.names[j]),paste0("0"),rhs[[i]], fixed = TRUE)
        }
        
        for (j in 1:length(random.names)){
            # gsub (a, b, c) : in c replace a with b
            rhs[[i]]=gsub(paste0(random.names[j]),paste0("0"),rhs[[i]], fixed = TRUE)
        }
        
        formula[[i]]=as.formula(paste0(lhs[[i]], ' ~ ', rhs[[i]]))
        #formula[[i]]=paste0(lhs[[i]], ' ~ ', rhs[[i]])
    }
    
    
    return(formula)
}
#print(prep.thetaFormula(model@dynamics@theta.formula, intercept.names = c("mu1", "mu2"),random.names = c('b_zeta', 'b_x1', 'b_x2')))

#---
#line 9- line 121: model specification
library('dynr')
state.names = c('x1', 'x2')
beta.names = c('zeta0', 'zeta1', 'zeta2', 'mu1', 'mu2')
covariate.names = c('u1', 'u2')
theta.names = c('zeta_i', 'zeta_i_2', 'zeta_i_3')
#[todo] Z*b's b
random.names = c('b_zeta', 'b_x1', 'b_x2')
intercept.names = c('mu1', 'mu2', 'mu3')


NxState = length(state.names)
Nbeta = length(beta.names)
Nx = NxState + Nbeta
Ntheta = length(theta.names)



#random data
N = 100
T = 300
vdpData <- data.frame(id=rep(1:N,each=T), time=rep(seq(0.005,1.5,by=0.005),N),
                      y=rnorm(100*300),
                      u1 = rnorm(100*300), u2 = rnorm(100*300))
colnames(vdpData) <- c("id","time","y","u1", "u2") # try
data <- dynr.data(vdpData, id="id", time="time",
                  observed=c('y'),
                  covariates=c('u1','u2'))


meas <- prep.measurement(
    values.load=matrix(c(1,0), 1, 2),
    params.load=matrix(c('fixed'), 1, 2),
    obs.names = c('y'),
    state.names=state.names)


initial <- prep.initial(
    values.inistate=c(3, 1),
    params.inistate=c("mu1", "mu2"),
    values.inicov=matrix(c(.5,.2,
                           .2,.6),ncol=2,byrow=T), 
    params.inicov=matrix(c('v10','c120',
                           'c120','v20'),ncol=2,byrow=T)
)

mdcov <- prep.noise(
    values.latent=diag(0, 2),
    params.latent=diag(c("fixed","fixed"), 2),
    values.observed=diag(rep(0.3,2)),
    params.observed=diag(c("var_1","var_2"),2)
)


formula=
    list(x1 ~ x2,
         x2 ~ -61.68503 * x1 + zeta_i * (1 - x1^2) * x2,
         zeta0 ~0,
         zeta1 ~0,
         zeta2 ~0,
         mu1 ~0,
         mu2 ~0
    )
# theta.formula  = list (zeta_i ~ 1 * zeta0  + u1 * zeta1 + u2 * zeta2 + 1*0,
# zeta_i_2 ~ 0,
# zeta_i_3 ~ 0)

theta.formula = list( zeta_i ~ 1 * zeta0 + u1 * zeta1 + u2 * zeta2 + 1 * b_zeta,
                      zeta_i_2 ~ 1 * mu1 + 1 * b_x1,
                      zeta_i_3 ~ 1 * mu2 + 1 * b_x2)

theta.formula2 = prep.thetaFormula(theta.formula, intercept.names, random.names)





dynm<-prep.formulaDynamics(formula=formula,
                           startval=c(zeta0=-1,
                                           zeta1=.5,
                                           zeta2=.2),
                                isContinuousTime=FALSE,
								state.names=state.names,
								theta.formula=theta.formula2,
								theta.names=theta.names,
								beta.names=beta.names,
								saem=TRUE)



model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=data,armadillo=TRUE,
                    outfile="VanDerPol.c")

# to do consist the formula in Line 71 and here
#model@dynamics@theta.formula = list( zeta_i ~ 1*zeta_0 + u1*zeta_1 + u2*zeta_2 + 1*b_zeta,
#                     zeta_i_2 ~ 1*mu1 + 1*b_x1,
#                      zeta_i_3 ~ 1*mu2 + 1*b_x2)

fitted_model <- dynr.cook(model, saem=TRUE, optimization_flag = TRUE, hessian_flag = TRUE, verbose=TRUE, debug_flag=TRUE)
print(fitted_model)

#---------
#following: parsing model to get the H and Z matrices

# # TODO interpret strict as character with match.arg.  strict=c('allowOne' let's people fly with x without writing 1*x, 'convent' let's no one get away with anything, 'hippie' is very loose
# formula2matrix <- function(formula, variables, strict=TRUE, col.match=TRUE, process=TRUE){
#     if(process){
#         pf <- dynr:::processFormula(list(formula))[[1]]
#         lhs <- pf[1]
#         rhs <- pf[2]
#         lhs2 <- sapply(lhs, strsplit, split=' + ', fixed=TRUE)[[1]]
#         rhs2 <- sapply(rhs, strsplit, split=' + ', fixed=TRUE)[[1]]
#         rhs3 <- strsplit(rhs2, split=" * ", fixed=TRUE)
#         if(strict == FALSE){
#             rhs4 <- lapply(rhs3, function(x){if(length(x) == 1) c("1", x) else x})
#         } else {rhs4 <- rhs3}
#     } else {
#         lhs2 <- formula[[1]]
#         rhs4 <- formula[[2]]
#     }
#     lens <- sapply(rhs4, length)
#     if(any(lens > 2)){
#         stop(paste0("I spy with my little eye terms in the right hand side of your formula with more than two parts.\n",
#                     "This function can only handle linear formulas.\nOffending term has: ",
#                     paste(rhs4[lens > 2], collapse=", ")))
#     }
#     if(any(lens < 2)){
#         stop(paste0("I spy with my little eye terms in the right hand side of your formula with less than two parts.\n",
#                     "Perhaps you forgot to multiply by 1.\nOffending term has: ",
#                     paste(rhs4[lens < 2], collapse=", ")))
#     }
#     
#     eleNames <- sapply(rhs4, function(x) {x[!(x %in% variables)]} )
#     ncol <- ifelse(col.match, length(variables), length(eleNames))
#     cnam <- if(col.match) variables else eleNames
#     
#     rmat <- matrix("0", nrow=length(lhs2), ncol=ncol, dimnames=list(lhs2, cnam))
#     for(aterm in 1:length(rhs4)){
#         x <- rhs4[[aterm]]
#         if(sum(x %in% colnames(rmat)) == 2){stop("Both parts of term are in the 'variables'")}
#         # TODO warning/error if strict = TRUE and trying to overwrite nonzero element
#         # if strict = FALSE try to write next nonzero element that matches
#         if(x[2] %in% colnames(rmat)){
#             if(strict==TRUE & !all(rmat[, x[2]] %in% "0")){
#                 warning(paste0("Overwriting element in column ", x[2], ". It was ", rmat[1, x[2]], " and now will be ", x[1], "."))
#             }
#             rmat[, colnames(rmat) %in% x[2]] <- x[1]
#         } else if(x[1] %in% colnames(rmat)){
#             if(strict==TRUE & !all(rmat[, x[1]] %in% "0")){
#                 warning(paste0("Overwriting element in column ", x[1], ". It was ", rmat[1, x[1]], " and now will be ", x[2], "."))
#             }
#             rmat[, colnames(rmat) %in% x[1]] <- x[2]
#         }
#     }
#     return(rmat)
# }
# 
# 
# # Example useage
# #formula2matrix(theta1 ~ x1 + x2 + mylabel*x3, variables=c('x1', 'x2', 'x3'), strict=FALSE)
# 
# #formula2matrix(theta1 ~ x1 + x2 + mylabel*x3, variables=c('1', 'mylabel'), strict=FALSE)
# 
# #formula2matrix(theta1 + theta2 ~ x1 + x2 + mylabel*x3, variables=c('x1', 'x2', 'x3'), strict=FALSE)
# 
# ## Factor model
# #formula2matrix(x1 + x2 + x3 + x4 + x5 ~ F, variables="F", strict=FALSE)
# #t(formula2matrix(F ~ x1 + x2 + x3 + x4 + x5, variables=paste0('x', 1:5), strict=FALSE))
# #t(formula2matrix(F ~ x1 + load2*x2 + load3*x3 + load4*x4 + load5*x5, variables=paste0('x', 1:5), strict=FALSE))
# 
# ## Regression among latent variables
# #formula2matrix(F1 ~ F2 + F3, variables=c("F2", "F3"), strict=FALSE)
# 
# 
# formula2design <- function(..., covariates, random.names){
#     dots <- list(...)
#     
#     pf <- list(dynr:::processFormula(dots))
#     lhs <- unlist(lapply(pf,function(x){lapply(x,"[[",1)})[[1]])
#     rhs <- lapply(pf,function(x){lapply(x,"[[",2)})[[1]]
#     rhs2 <- sapply(rhs, strsplit, split=' + ', fixed=TRUE)
#     rhs3 <- lapply(rhs2, strsplit, split=" * ", fixed=TRUE)
#     #rhs4 <- lapply(rhs3, function(x){if(length(x) == 1) c("1", x) else x})
#     eleNames <- unique(unlist(sapply(rhs3, function(rlist){ sapply(rlist, function(x) {x[!(x %in% covariates) & !(x %in% random.names)]} ) } )))
#     
#     dots <- cbind(lhs, rhs3)
#     
#     fmat <- t(apply(dots, 1, formula2matrix, variables=eleNames, strict=TRUE, process=FALSE))
#     rmat <- t(apply(dots, 1, formula2matrix, variables=random.names, strict=TRUE, process=FALSE))
#     dimnames(fmat) <- list(lhs, eleNames)
#     dimnames(rmat) <- list(lhs, random.names)
#     return(list(fixed=fmat, random=rmat))
# }
# 
# 
# # Example useage
# #formula2design(
# #	theta1 ~ 1*zeta_0 + u1*zeta_1 + u2*zeta_2 + 1*b_zeta,
# #	theta2 ~ mu_x1*1 + 1*b_x1,
# #	theta3 ~ mu_x2*1 + 1*b_x2,
# #	covariates=c('1', 'u1', 'u2'),
# #	random.names=c('b_zeta', 'b_x1', 'b_x2'))
# 
# 
# 
# multiformula2matrix <- function(..., variables){
#     dots <- list(...)
#     
#     pf <- list(dynr:::processFormula(dots))
#     lhs <- unlist(lapply(pf,function(x){lapply(x,"[[",1)})[[1]])
#     rhs <- lapply(pf,function(x){lapply(x,"[[",2)})[[1]]
#     rhs2 <- sapply(rhs, strsplit, split=' + ', fixed=TRUE)
#     rhs3 <- lapply(rhs2, strsplit, split=" * ", fixed=TRUE)
#     #rhs4 <- lapply(rhs3, function(x){if(length(x) == 1) c("1", x) else x})
#     
#     dots <- cbind(lhs, rhs3)
#     
#     fmat <- t(apply(dots, 1, formula2matrix, variables=variables, strict=TRUE, process=FALSE))
#     dimnames(fmat) <- list(lhs, variables)
#     return(fmat)
# }
# 
# # Example useage
# # Regression among latent variables
# #multiformula2matrix(F1 ~ a*F2 + b*F3, F2 ~ 0*F2, F3 ~ c*F2, variables=c("F1", "F2", "F3"))
# 
# matrix2ArmadilloCode <- function(variable.name, matrix.input){
#     str = paste0(variable.name, " = \"")
#     dim(matrix.input)
#     for(i in c(1:dim(matrix.input)[1])){
#         for(j in c(1:dim(matrix.input)[2])){
#             str = paste0(str, " ", as.character(matrix.input[i,j]))
#         }
#         str = paste0(str, ";")
#     }
#     str = paste0(str, "\";")
#     return(str)
# }
# 
# 
# r =formula2design( 
# 	theta1 ~ 1*zeta_0 + u1*zeta_1 + u2*zeta_2 + 1*b_zeta,
# 	theta2 ~ 1*mu1 + 1*b_x1,
# 	theta3 ~ 1*mu2 + 1*b_x2,
# 	covariates=c('1', 'u1', 'u2'),
# 	random.names=c('b_zeta', 'b_x1', 'b_x2'))
# 
# #covariate.names=covariate.names,beta.names=beta.names,
# r$fixed
# r$random
# 
# #todo: 
# #1. InfDS.U1(InfDS.Nsubj * InfDS.Ncovariate), [u1, u2] for each subject
# #2. read [u1, u2] for each subject, fill in into u1 u2 in r$fixed and create InfDS.H 
# #3. IndDS.Z = r$random
# 
# Z= apply(r$random, 1, as.numeric)
# 
# H = matrix(nrow=0, ncol=0)
# #U1 <- read.csv("U1.csv", header=FALSE)
# #for (line in c(1: dim(U1)[1])){
# for (line in c(1: as.integer(dim(vdpData)[1]/T))){
# #for (line in c(1: 2)){ 
#     #u1 <- U1[line,1]    
#     #u2 <- U1[line,2]
#     
#     U = c(1:length(covariate.names))
#     for (u in c(1:length(covariate.names))){
#         U[u] = mean(vdpData[covariate.names[u]][((line-1)*T+1):(line*T),])
#     }
#     
#     temp <- r$fixed
#     for (i in c(1:dim(r$fixed)[1])){
#         for (j in c(1:dim(r$fixed)[2])){
#             for (u in c(1:length(covariate.names)))
#                 if(r$fixed[i,j] %in% covariate.names[u]){
#                     #r$fixed[i,j] <- U1[line,u]
#                     r$fixed[i,j] <- U[u]
#                 }
#             # if(r$fixed[i,j] %in% "u1"){
#             #     r$fixed[i,j] <- u1
#             # }
#             # else if (r$fixed[i,j] %in% "u2"){
#             #     r$fixed[i,j] <- u2
#             # }
#             # else{
#             #     r$fixed[i,j] = as.numeric(r$fixed[i,j])
#             # }
#         }
#     }
#     
#     
#     r$fixed
#     
#     if(dim(H)[1] == 0 && dim(H)[2] == 0){
#         H = t(apply(r$fixed, 1, as.numeric))
#     } else{
#         H=rbind(H,t(apply(r$fixed, 1, as.numeric)))
#     }
#     r$fixed<- temp
# }
# mean(vdpData['u1'][1:T,])
# mean(vdpData['u2'][1:T,])
# H[1:3,] 
# mean(vdpData['u1'][601:900,])
# mean(vdpData['u2'][601:900,])
# H[7:9,]  
# 
# matrix2ArmadilloCode("InfDS.Z", Z)
# matrix2ArmadilloCode("InfDS.H", H)
# 
# #ask/discuss: 
# #1. How to obtain InfDS.U1 from the model?
# #2. How to obtain the information in main.cpp from model?
# write(paste0(matrix2ArmadilloCode("InfDS.Z", Z), '\n', 
#             matrix2ArmadilloCode("InfDS.H", H)),  
#       file( "temp.txt", "w" ))
