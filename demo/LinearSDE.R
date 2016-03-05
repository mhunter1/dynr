#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2016-02-25
# Filename: LinearSDE.R
# Purpose: Translate LinearSDE.R into the user-spec for the same model
#  in dynr.
#------------------------------------------------------------------------------



require(dynr)
require(KFAS)

meas <- dynr.matrixLoadings(values=matrix(c(1,0), 1, 2), params=matrix(0, 1, 2))

ecov <- dynr.matrixErrorCov(diag(c(0.00001,1)), diag(c(0, 3)), diag(1.5,1), diag(4, 1))
ecov$c.string
ecov$startval

dP_dt <- dynr.dP_dt



