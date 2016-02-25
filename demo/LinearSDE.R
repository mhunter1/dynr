#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2016-02-25
# Filename: LinearSDE.R
# Purpose: Translate LinearSDE.R into the user-spec for the same model
#  in dynr.
#------------------------------------------------------------------------------



require(dynr)

meas <- dynr.matrixLoadings(values=matrix(c(1,0), 1, 2), params=matrix(0, 1, 2))



dP_dt <- dynr.dP_dt



