#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2016-02-25
# Filename: LinearSDE.R
# Purpose: Translate LinearSDE.R into the user-spec for the same model
#  in dynr.
#------------------------------------------------------------------------------



require(dynr)

meas <- dynr.matrixLoadings(values=matrix(c(1,0), 1, 2), params=matrix(0, 1, 2))

ecov <- dynr.matrixErrorCov(diag(c(-10, log(1))), diag(c(0, 3)), diag(log(1.5), 1), diag(4, 1))

dP_dt <- dynr.dP_dt



