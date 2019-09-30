#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2016-10-31
# Filename: RSLinearDiscreteYang.R
# Purpose: Show regime-switching 
# measurement and dynamics in linear discrete time model
# with comparison with results reported in Yang & Chow  (2010).
#------------------------------------------------------------------------------


#---- (1) Load packages ----
require(dynr)

#---- (2) Read in data ----
# create dynr data object
data(EMG)
EMGdata <- dynr.data(EMG, id='id', time='time', 
                     observed='iEMG', covariates='SelfReport')

#---- (3) Specify recipes for all model pieces ----

#---- (3a) Measurement ----
recMeas <- prep.measurement(
  values.load = rep(list(matrix(1, 1, 1)), 2),
  values.int = list(matrix(3, 1, 1), matrix(4, 1, 1)),
  params.int = list(matrix('mu_1', 1, 1), matrix('mu_2', 1, 1)),
  values.exo = list(matrix(0, 1, 1), matrix(1, 1, 1)),
  params.exo = list(matrix('fixed', 1, 1), matrix('beta_2', 1, 1)),
  obs.names = c('iEMG'),
  state.names = c('eta'),
  exo.names = c("SelfReport"))

#---- (3b) Dynamic and measurement noise cov structures----

recNoise <- prep.noise(
  values.latent=matrix(.5, 1, 1),
  params.latent=matrix('dynNoise', 1, 1),
  values.observed=matrix(0, 1, 1),
  params.observed=matrix('fixed', 1, 1))

# ---- (3c) Regimes-switching model ----

recReg <- prep.regimes(
  values = matrix(c(.7, -1, 0, 0), 2, 2),
  params = matrix(c('c11', 'c21', 'fixed', 'fixed'), 2, 2))

#---- (3d) Initial condition specification ----

recIni <- prep.initial(
  values.inistate=matrix(0, 1, 1),
  params.inistate=matrix('fixed', 1, 1),
  values.inicov=matrix(1, 1, 1),
  params.inicov=matrix('fixed', 1, 1),
  values.regimep=c(10, 0),
  params.regimep=c('fixed', 'fixed'))


#---- (3e) Dynamic model ----

recDyn <- prep.matrixDynamics(
  values.dyn = list(matrix(.5, 1, 1), matrix(.1, 1, 1)),
  params.dyn = list(matrix('phi_1', 1, 1), matrix('phi_2', 1, 1)),
  isContinuousTime = FALSE)

#---- (4a) Create model  ----

rsmod <- dynr.model(
  dynamics = recDyn,
  measurement = recMeas,
  noise = recNoise,
  initial = recIni,
  regimes = recReg,
  data = EMGdata,
  outfile = "RSLinearDiscrete.c")

rsmod$ub[c('phi_1', 'phi_2')] <- 1.1

#---- (4b) Check model specification  ----

printex(rsmod,
        ParameterAs = rsmod$param.names,
        printInit = TRUE, printRS = TRUE,
        outFile = "RSLinearDiscreteYang.tex")
#tools::texi2pdf("RSLinearDiscreteYang.tex")
#system(paste(getOption("pdfviewer"), "RSLinearDiscreteYang.pdf"))

#---- (4c) Create model and cook it all up  ----

yum <- dynr.cook(rsmod, verbose = FALSE)

#---- (5) Serve it! ----

summary(yum)


# Overview multicomponent plot quickly showing
#  (1) latent states with regimes
#  (2) regime histogram
#  (3) typeset model specification
plot(yum, dynrModel = rsmod, style = 1, textsize = 5)

# Figure 1 (B) from R Journal paper
# Plot showing the latent state estimates along with
#  which regime is most likely.
dynr.ggplot(yum, dynrModel = rsmod, style = 1,
            names.regime = c("Deactivated", "Activated"),
            title = "(B) Results from RS-AR model", numSubjDemo = 1,
            shape.values = c(1),
            text = element_text(size=16),
            is.bw = TRUE)

# Figure 2 (A) from R Journal paper
# Plot of model equations in terms of parameter names
plotFormula(dynrModel = rsmod, ParameterAs = names(rsmod),
            printDyn = TRUE, printMeas = TRUE) +
  ggtitle("(A)")+
  theme(plot.title = element_text(hjust = 0.5, vjust=0.01, size=16)) 

# Figure 2 (B) from R Journal paper
# Plot of model equations in terms of parameter values
plotFormula(dynrModel = rsmod, ParameterAs = coef(yum),
            printDyn = TRUE, printMeas = TRUE) +
  ggtitle("(B)")+
  theme(plot.title = element_text(hjust = 0.5, vjust=0.01, size=16)) 

