#---------------------------------------------------------------------
# Author: Dongjun You
# Date: 2018-05-15
# Filename: ShockDemo1.R
# Purpose: An illustrative example of using dynr.taste to detect
#   outliers from a state space model
#---------------------------------------------------------------------
library(dynr)
library(data.table)

beta2 <- matrix(c(0.8, -0.2,
                 -0.2,  0.7), ncol=2, byrow=TRUE)
lambda2 <- matrix(c(1.0, 0.0,
                    0.9, 0.0,
                    0.8, 0.0, 
                    0.0, 1.0, 
                    0.0, 0.9, 
                    0.0, 0.8), ncol=2, byrow=TRUE)
psi2 <- matrix(c(0.3, -0.1,  
                -0.1,  0.3), ncol=2, byrow=TRUE)
theta2 <- diag(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2), ncol=6, nrow=6)

true_params <- c(.8, -.2, -.2, .8,# beta
                 .9, .8, .9, .8,# lambda
                 .3, -.1, .3,# psi
                 .2, .2, .2, .2, .2, .2)# theta

# yy2 <- gen_taste(time=60, nsubj=100, nlat=2, nobs=6,
#                  lambda=lambda2, theta=theta2, beta=beta2, psi=psi2, 
#                  nshockLat=3, nshockObs=6, shockMag=2.5, gen_method="random")
# write.table(yy2$y, file="ShockData.txt", row.names=FALSE)
# write.table(yy2$shocks$shockL, file="ShockLat.txt", row.names=FALSE)
# write.table(yy2$shocks$shockO, file="ShockObs.txt", row.names=FALSE)

data(ShockData1)
data(ShockLat1)
data(ShockObs1)
data_shk <- dynr.data(ShockData, id='id', time='time',
                      observed=c('V1','V2','V3','V4','V5','V6'))
# data_shk <- dynr.data(yy2$y, id='id', time='time',
#                       observed=c('V1','V2','V3','V4','V5','V6'))

meas_shk <- prep.measurement(
  values.load=matrix(c(1.0, 0.0,
                       0.9, 0.0,
                       0.8, 0.0, 
                       0.0, 1.0, 
                       0.0, 0.9, 
                       0.0, 0.8), ncol=2, byrow=TRUE),
  params.load=matrix(c('fixed','fixed',
                       'l_21','fixed',
                       'l_31','fixed',
                       'fixed','fixed',
                       'fixed','l_52',
                       'fixed','l_62'), ncol=2, byrow=TRUE),
  state.names=c("eta_1",'eta_2'),
  obs.names=c('V1','V2','V3','V4','V5','V6') )

nois_shk <- prep.noise(
  values.latent=matrix(c(0.3, -0.1,  
                         -0.1, 0.3), ncol=2, byrow=TRUE),
  params.latent=matrix(c('psi_11','psi_12',
                         'psi_12','psi_22'), ncol=2, byrow=TRUE),
  values.observed=diag(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2), ncol=6, nrow=6),
  params.observed=diag( paste0('e_', 1:6), 6) )

init_shk <- prep.initial(
  values.inistate=c(0,0),
  params.inistate=c('mu_1','mu_2'),
  values.inicov=matrix(c(0.3, -0.1,  
                         -0.1,  0.3), ncol=2, byrow=TRUE),
  params.inicov=matrix(c('c_11','c_12',
                         'c_12','c_22'), ncol=2, byrow=TRUE) )

dynm_shk <- prep.matrixDynamics(
  values.dyn=matrix(c(0.8, -0.2,
                      -0.2,  0.7), ncol=2, byrow=TRUE),
  params.dyn=matrix(c('b_11','b_12',
                      'b_21','b_22'), ncol=2, byrow=TRUE),
  isContinuousTime=FALSE)

model_shk <- dynr.model(dynamics=dynm_shk, measurement=meas_shk,
                        noise=nois_shk, initial=init_shk,
                        data=data_shk, outfile=paste0("model_shk.c"))
cook_shk <- dynr.cook(model_shk, debug_flag=TRUE)

# detect outliers
taste_shk <- dynr.taste(model_shk, cook_shk, conf.level=.99)

# apply the detected outliers, and re-fit the model
taste2_shk <- dynr.taste2(model_shk, cook_shk, taste_shk,
                          outlierTest="t",
                          newOutfile="taste2_shk.c")

# (optional) detect outliers from the re-fitted model
taste_taste2_shk <- dynr.taste(taste2_shk$new_dynrModel,
                               taste2_shk$new_dynrCook, 
                               conf.level=.99)

true_params <- c(.8, -.2, -.2, .8,# beta
                 .9, .8, .9, .8,# lambda
                 .3, -.1, .3,# psi
                 .2, .2, .2, .2, .2, .2)# theta

# compare true parameter and estimated ones
#  cook: parameter estimates when the data include true outliers
#  taste2: parameter estimates after applying the estimated outliers by the function 'dynr.taste2'
# check the how closely recover the true parameter from 'dynr.taste2', i.e., after applying outliers
data.frame(true=true_params, cook=coef(cook_shk)[1:17],
           taste2=coef(taste2_shk$new_dynrCook)[1:17])

# helper function to combine detected innovative outliers (delta) using t-test
delta_t_L <- function(tasteOut) {
  shocks <- lapply(tasteOut, function(tasteout_i) {
    data.table::data.table(id=tasteout_i$taste[["id"]],
                           time=tasteout_i$taste[["time"]],
                           tasteout_i$delta_t$delta.L)
  } )
  data.table::rbindlist(shocks)
}

# helper function to combine detected additive outliers (delta) using t-test
delta_t_O <- function(tasteOut) {
  shocks <- lapply(tasteOut, function(tasteout_i) {
    data.table::data.table(id=tasteout_i$taste[["id"]],
                           time=tasteout_i$taste[["time"]],
                           tasteout_i$delta_t$delta.O)
  } )
  data.table::rbindlist(shocks)
}

# the number of outliers from the null model
# the number of innovative outliers 
chk_t_L <- delta_t_L(taste_shk)
eta1_t <- chk_t_L[, .(time, id, shock=d_eta_1!=0)][shock==TRUE,]
eta1_t[, .N]# 163
table(eta1_t[, time])
eta2_t <- chk_t_L[, .(time, id, shock=d_eta_2!=0)][shock==TRUE,]
eta2_t[, .N]# 137
table(eta2_t[, time])
# the number of additive outliers
chk_t_O <- delta_t_O(taste_shk)
v1_t <- chk_t_O[, .(time, id, shock=d_V1!=0)][shock==TRUE,]
v1_t[, .N]# 115
table(v1_t[, time])
v2_t <- chk_t_O[, .(time, id, shock=d_V2!=0)][shock==TRUE,]
v2_t[, .N]# 119
table(v2_t[, time])
v3_t <- chk_t_O[, .(time, id, shock=d_V3!=0)][shock==TRUE,]
v3_t[, .N]# 114
table(v3_t[, time])
v4_t <- chk_t_O[, .(time, id, shock=d_V4!=0)][shock==TRUE,]
v4_t[, .N]# 122
table(v4_t[, time])
v5_t <- chk_t_O[, .(time, id, shock=d_V5!=0)][shock==TRUE,]
v5_t[, .N]# 118
table(v5_t[, time])
v6_t <- chk_t_O[, .(time, id, shock=d_V6!=0)][shock==TRUE,]
v6_t[, .N]# 110
table(v6_t[, time])

# the number of outliers from the updated model applied with the estimated outliers (compare with the results from the above null model, and check the reduced number of outliers )
# the number of innovative outliers 
chk2_t_L <- delta_t_L(taste_taste2_shk)
eta1_2_t <- chk2_t_L[, .(time, id, shock=d_eta_1!=0)][shock==TRUE,]
eta1_2_t[, .N]# 50
table(eta1_2_t[, time])
eta2_2_t <- chk2_t_L[, .(time, id, shock=d_eta_2!=0)][shock==TRUE,]
eta2_2_t[, .N]# 68
table(eta2_2_t[, time])
# the number of additive outliers
chk2_t_O <- delta_t_O(taste_taste2_shk)
v1_2_t <- chk2_t_O[, .(time, id, shock=d_V1!=0)][shock==TRUE,]
v1_2_t[, .N]# 49
table(v1_2_t[, time])
v2_2_t <- chk2_t_O[, .(time, id, shock=d_V2!=0)][shock==TRUE,]
v2_2_t[, .N]# 44
table(v2_2_t[, time])
v3_2_t <- chk2_t_O[, .(time, id, shock=d_V3!=0)][shock==TRUE,]
v3_2_t[, .N]# 44
table(v3_2_t[, time])
v4_2_t <- chk2_t_O[, .(time, id, shock=d_V4!=0)][shock==TRUE,]
v4_2_t[, .N]# 45
table(v4_2_t[, time])
v5_2_t <- chk2_t_O[, .(time, id, shock=d_V5!=0)][shock==TRUE,]
v5_2_t[, .N]# 50
table(v5_2_t[, time])
v6_2_t <- chk2_t_O[, .(time, id, shock=d_V6!=0)][shock==TRUE,]
v6_2_t[, .N]# 49
table(v6_2_t[, time])

shockLat <- yy2$shocks$shockL
shockObs <- yy2$shocks$shockO

# calculate true and false detection ratio for each individual 
l_id <- lapply(1:100, function(idd) {
  rr <- taste_shk[[idd]][["taste"]]
  chi_L <- which(rr[["chi.L.shk"]])
  chi_O <- which(rr[["chi.O.shk"]])
  t_L_1 <- which(rr[["t.L.shk.eta_1"]])
  t_L_2 <- which(rr[["t.L.shk.eta_2"]])
  
  t_O_1 <- which(rr[["t.O.shk.V1"]])
  t_O_2 <- which(rr[["t.O.shk.V2"]])
  t_O_3 <- which(rr[["t.O.shk.V3"]])
  t_O_4 <- which(rr[["t.O.shk.V4"]])
  t_O_5 <- which(rr[["t.O.shk.V5"]])
  t_O_6 <- which(rr[["t.O.shk.V6"]])
  
  sk_L_1 <- shockLat[id==idd, .(time_L, lat)][lat==1, time_L]
  sk_L_2 <- shockLat[id==idd, .(time_L, lat)][lat==2, time_L]
  
  sk_OO <- shockObs[id==idd, .(time_O, obs)]
  sk_O_1 <- sk_OO[obs==1, time_O]; sk_O_2 <- sk_OO[obs==2, time_O]
  sk_O_3 <- sk_OO[obs==3, time_O]; sk_O_4 <- sk_OO[obs==4, time_O]
  sk_O_5 <- sk_OO[obs==5, time_O]; sk_O_6 <- sk_OO[obs==6, time_O]
  
  detect_L_1 <- if (lsk_L_1 <- length(sk_L_1)) {t_L_1 %in% sk_L_1} else {NA}
  detect_L_2 <- if (lsk_L_2 <- length(sk_L_2)) {t_L_2 %in% sk_L_2} else {NA}
  
  detect_O_1 <- if (lsk_O_1 <- length(sk_O_1)) {t_O_1 %in% sk_O_1} else {NA}
  detect_O_2 <- if (lsk_O_2 <- length(sk_O_2)) {t_O_2 %in% sk_O_2} else {NA}
  detect_O_3 <- if (lsk_O_3 <- length(sk_O_3)) {t_O_3 %in% sk_O_3} else {NA}
  detect_O_4 <- if (lsk_O_4 <- length(sk_O_4)) {t_O_4 %in% sk_O_4} else {NA}
  detect_O_5 <- if (lsk_O_5 <- length(sk_O_5)) {t_O_5 %in% sk_O_5} else {NA}
  detect_O_6 <- if (lsk_O_6 <- length(sk_O_6)) {t_O_6 %in% sk_O_6} else {NA}
  
  T_detect_L_1 = sum(detect_L_1)/lsk_L_1; T_detect_L_2 = sum(detect_L_2)/lsk_L_2
  T_detect_O_1 = sum(detect_O_1)/lsk_O_1; T_detect_O_2 = sum(detect_O_2)/lsk_O_2
  T_detect_O_3 = sum(detect_O_3)/lsk_O_3; T_detect_O_4 = sum(detect_O_4)/lsk_O_4
  T_detect_O_5 = sum(detect_O_5)/lsk_O_5; T_detect_O_6 = sum(detect_O_6)/lsk_O_6
  
  data.table(
    T_detect_L_1 = T_detect_L_1, F_detect_L_1 = sum(!detect_L_1)/(60 - lsk_L_1),
    T_detect_L_2 = T_detect_L_2, F_detect_L_2 = sum(!detect_L_2)/(60 - lsk_L_2),
    T_detect_O_1 = T_detect_O_1, F_detect_O_1 = sum(!detect_O_1)/(60 - lsk_O_1),
    T_detect_O_2 = T_detect_O_2, F_detect_O_2 = sum(!detect_O_2)/(60 - lsk_O_2),
    T_detect_O_3 = T_detect_O_3, F_detect_O_3 = sum(!detect_O_3)/(60 - lsk_O_3),
    T_detect_O_4 = T_detect_O_4, F_detect_O_4 = sum(!detect_O_4)/(60 - lsk_O_4),
    T_detect_O_5 = T_detect_O_5, F_detect_O_5 = sum(!detect_O_5)/(60 - lsk_O_5),
    T_detect_O_6 = T_detect_O_6, F_detect_O_6 = sum(!detect_O_6)/(60 - lsk_O_6)
  )
})

# average outlier detection rates
# T_ : true detection rates, F_ : false detection rates
rbindlist(l_id)[, lapply(.SD, mean, na.rm=TRUE)]

