dynr.mi=function(model,m=5,aux.variable=c("x1","x2"),imp.obs=0,imp.exo=0,imp.aux=0){    #multiple lag; #factor  #get variable names

  # m=2
  # aux.variable=c("x1","x2")
  # imp.obs=TRUE
  # imp.exo=TRUE
  # 
  
  
data=model@data$original.data
k=length(model@param.names) #- length(model@initial@paramnames)             #number of parameters estimated
nolag=1


ynames=model@data$observed.names
xnames=model@data$covariate.names
y=data[,colnames(data)==ynames]
x=data[,colnames(data)==xnames]
ID=model@data$id
time=model@data$time

au=data[,colnames(data)==aux.variable]  

dataforlag=cbind(ID,y,x)

#dataforlag=ts(dataforlag)

#head(stats::lag(as.xts(dataforlag),2))
#head(lag(dataforlag,2))

datalag <- 
  dataforlag %>%
  group_by(ID) %>%
  mutate_all(lag)

dataformice=cbind(dataforlag[,-1],datalag[,-1],au)
dataformice=data.frame(dataformice)
# dataformice[,3]=as.factor(dataformice[,3])    #identify factors
# dataformice[,7]=as.factor(dataformice[,7])

colnames(dataformice)=c()
m=m
imp=mice(dataformice,m=m)


pmcarqhat=matrix(NA, nrow=m,ncol=k) #parameter estimates from each imputation
pmcaru=array(NA,dim=c(k,k,m)) #vcov of par estimates from each imputation

for (j in 1:m){
  
  completedata=complete(imp,action=j)
  colnames(completedata)=c(ynames,xnames,
                          paste("lag",ynames,sep=''),
                          paste("lag",xnames,sep=''),
                          colnames(au))
  
  if (imp.obs==TRUE){
  imp.data.obs= completedata[,colnames(completedata)==ynames]  
  }else{
    imp.data.obs=y
  }
 
  if (imp.exo==TRUE){
    imp.data.exo= completedata[,colnames(completedata)==xnames]  
  }else{
    imp.data.exo=x
  }
  
  newdata=cbind(ID,time,imp.data.obs,imp.data.exo)
  
  save(newdata,file="test.rdata")
  
  colnames(newdata)=c("ID","Time","wp","hp","ca","cn")
  
  data <- dynr.data(newdata, id="ID", time="Time", 
                    observed=c("wp","hp"),covariates=c("ca","cn"))
  
 modelnew=model
 modelnew@data=data
 #   model <- dynr.model(dynamics=model@dynamics, measurement=model@measurement,
 #                      noise=model@noise, initial=model@initial, data=data, transform=model@transform,
 #                     outfile=paste("trial4",i,".c",sep=""))
  
   
   # model <- dynr.model(dynamics=dynm, measurement=meas,
   #                     noise=mdcov, initial=initial, data=data,#transform=trans,
   #                     outfile=paste("trial2",i,".c",sep=""))
   
   
  trial <- dynr.cook(modelnew)  #names(trial) get names of the params
  summary(trial)
  
  #getting parameter estimates
  pmcarqhat[j,]=coef(trial)[1:k]
  pmcaru[, ,j]= vcov(trial)[c(1:k),c(1:k)]
}

pqbarmcarimpute <- apply(pmcarqhat, 2, mean) 
pubarmcarimpute <- apply(pmcaru, 1:2, mean)
#ubar <- apply(u, c(2, 3), mean)
pe.mcarimpute <- pmcarqhat - matrix(pqbarmcarimpute, nrow = m, ncol = k, byrow = TRUE)
pb.mcarimpute <- (t(pe.mcarimpute) %*% pe.mcarimpute)/(m - 1)
pvcovmcarimpute <- pubarmcarimpute + (1 + 1/m) * pb.mcarimpute #vcov for estimates
psemcarimpute=sqrt(diag(pvcovmcarimpute))

pqbarmcarimpute+2*psemcarimpute
pqbarmcarimpute-2*psemcarimpute

return(cbind(pqbarmcarimpute,psemcarimpute))
}

mi.cook(rawdata,m=5,original.data=ddata,aux.variable=c("x1","x2"),imp.y=TRUE,imp.x=TRUE,k,nolag=1)




############################################################################################################
mcartemphp=mcarmisshp[1,]
mcartempwp=mcarmisswp[1,]
mcartempca=mcarmissca[1,]
mcartempcn=mcarmisscn[1,]
mcartempx1=x1full[1,]
mcartempx2=x2full[1,]

mcarmisscnlag1=cbind(rep(NA,n),mcarmisscn[,c(1:(nt-1))])
mcartempcnlag1=mcarmisscnlag1[1,]
mcarmisscalag1=cbind(rep(NA,n),mcarmissca[,c(1:(nt-1))])
mcartempcalag1=mcarmisscalag1[1,]
mcarmisshplag1=cbind(rep(NA,n),mcarmisshp[,c(1:(nt-1))])
mcartemphplag1=mcarmisshplag1[1,]
mcarmisswplag1=cbind(rep(NA,n),mcarmisswp[,c(1:(nt-1))])
mcartempwplag1=mcarmisswplag1[1,]

for (i in 2:n){
  mcartemphp=c(mcartemphp,mcarmisshp[i,])
  mcartempwp=c(mcartempwp,mcarmisswp[i,])
  mcartempca=c(mcartempca,mcarmissca[i,])
  mcartempcn=c(mcartempcn,mcarmisscn[i,])
  mcartempx1=c(mcartempx1,x1full[i,])
  mcartempx2=c(mcartempx2,x2full[i,])
  mcartempcnlag1=c(mcartempcnlag1,mcarmisscnlag1[i,])
  mcartempcalag1=c(mcartempcalag1,mcarmisscalag1[i,])
  mcartemphplag1=c(mcartemphplag1,mcarmisshplag1[i,])
  mcartempwplag1=c(mcartempwplag1,mcarmisswplag1[i,])
}

mcarformice=cbind(mcartempwp,mcartemphp,mcartempca,mcartempcn,mcartempwplag1,mcartemphplag1,mcartempcalag1,mcartempx1,mcartempx2)#,mcartempcnlag1,mcartempx3,mcartempx4)
mcarformicef=data.frame(mcarformice)
mcarformicef[,3]=as.factor(mcarformice[,3])
mcarformicef[,7]=as.factor(mcarformice[,7])
levels(mcarformicef[,3])
mcarimp=mice(mcarformicef,m=m)

#initial list to store outputs 
k=11     #number of parameters estimated  
pmcarqhat=matrix(NA, nrow=m,ncol=k) #parameter estimates from each imputation
pmcaru=array(NA,dim=c(k,k,m)) #vcov of par estimates from each imputation

for (j in 1:m){

mcarnew=complete(mcarimp,action=j)
mcarcnforfit=matrix(mcarnew[,4],nrow=n,byrow=TRUE)
mcarcaforfit=matrix(mcarnew[,3],nrow=n,byrow=TRUE)

#create dataset for dynr
ddata=matrix(,ncol=6,nrow=0)

#MCAR data
for (i in 1:n){
  temp=cbind(rep(i,nt),seq(1:nt),mcarcaforfit[i,],mcarcnforfit[i,],mcarmisswp[i,],mcarmisshp[i,])
  ddata=rbind(ddata,temp)
}

colnames(ddata)=c("ID","Time","ca","cn","wp","hp")

data <- dynr.data(ddata, id="ID", time="Time", 
                  observed=c("wp","hp"),covariates=c("ca","cn"))

meas <- prep.measurement(
  values.load=matrix(c(1,0,
                       0,1),ncol=2,byrow=T), # starting values and fixed values
  params.load=matrix(rep("fixed",4),ncol=2),
  state.names=c("wp","hp"),
  obs.names=c("wp","hp")
)

 initial <- prep.initial(
   values.inistate=c(-.5,-.9),
   params.inistate=c('mu_wp', 'mu_hp'),
  values.inicov=matrix(c(1,-0.3,
                         -0.3,1),byrow=T,ncol=2),
  params.inicov=matrix(c("v_11","c_21",
                         "c_21","v_22"),byrow=T,ncol=2))

# initial <- prep.initial(
#   values.inistate=rep(0,4),
#   params.inistate=rep('fixed',4), 
#   values.inicov=diag(rep(.1,4)), 
#   params.inicov=diag(rep('fixed',4))) 


 mdcov <- prep.noise(
   values.latent=matrix(c(0.5,-0.05,
                          -0.05,0.5),byrow=T,ncol=2), 
   params.latent=matrix(c("v_wp","c_hw",
                          "c_hw","v_hp"),byrow=T,ncol=2), 
   values.observed=diag(rep(0,2)), 
   params.observed=diag(c('fixed','fixed'),2)) 

 #wp(t)=a*wp(t-1)+b*hp(t-1)+c*ca +d*cn+we(t)
 #hp(t)=a1*hp(t-1)+b1*wp(t-1)+c1*ca +d1*cn+he(t) 
 
 formula =list(
   list(wp ~ a*wp + b*hp + c*ca + d*cn,
        hp ~ a1*hp + b1*wp +c1*ca + d1*cn
   ))
 # 
 # a=0.4 ; b=-0.3 ; 
 # b1=-0.2 ; a1=.3
 # 
 # c=.3 ; c1=.3
 # d=-.5 ; d1=-.4
 dynm  <- prep.formulaDynamics(formula=formula,
                               startval=c(a = .4, b = -.3, b1=-.2, a1=.3, 
                                          c = .3, c1=.3, d=-.5, d1=-.4
                               ), isContinuousTime=FALSE) #,jacobian=jacob
 

 model <- dynr.model(dynamics=dynm, measurement=meas,
                           noise=mdcov, initial=initial, data=data,#transform=trans,
                           outfile=paste("trial2",i,".c",sep=""))

# You can uncomment these lines and it will now work with dynr version v0.1.8-36-gf96ade0
# model$ub <- c(a=1.5, b=1.5, b1=1.5, a1=1.5, c=5, c1=5, d=5, d1=5, v_wp=2, c_hw=NA, v_hp=2, mu_wp=5, mu_hp=5, v_11=5, c_21=NA, v_22=5)

# model$lb <- c(a=-2, b=-2, b1=-2, a1=-2, c=-5, c1=-5, d=-5, d1=-5, v_wp=1e-4, c_hw=NA, v_hp=1e-4, mu_wp=-5, mu_hp=-5, v_11=1e-4, c_21=NA, v_22=1e-4)


 
 trial <- dynr.cook(model)
 summary(trial)
 
 #getting parameter estimates
 pmcarqhat[j,]=coef(trial)[1:11]
 pmcaru[, ,j]= vcov(trial)[c(1:11),c(1:11)]
}
 
 pqbarmcarimpute <- apply(pmcarqhat, 2, mean) 
 pubarmcarimpute <- apply(pmcaru, 1:2, mean)
 #ubar <- apply(u, c(2, 3), mean)
 pe.mcarimpute <- pmcarqhat - matrix(pqbarmcarimpute, nrow = m, ncol = k, byrow = TRUE)
 pb.mcarimpute <- (t(pe.mcarimpute) %*% pe.mcarimpute)/(m - 1)
 pvcovmcarimpute <- pubarmcarimpute + (1 + 1/m) * pb.mcarimpute #vcov for estimates
 psemcarimpute=sqrt(diag(pvcovmcarimpute))
 
 pqbarmcarimpute+2*psemcarimpute
 pqbarmcarimpute-2*psemcarimpute
 
 return(cbind(pqbarmcarimpute,psemcarimpute))
 }

      
