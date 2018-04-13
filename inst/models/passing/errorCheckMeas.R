#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2018-04-12 15:08:22
# Filename: errorCheckMeas.R
# Purpose: Check that errors in measurement are caught and reported properly.
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Load packages

require(dynr)


#------------------------------------------------------------------------------
# Check that the dimensions of load, exo, and int match with what is in 
#  obs.names, state.names, and exo.names

# Error
testthat::expect_error(
	me1 <- prep.measurement(values.load=matrix(c(1, 0), 1, 1),
		obs.names=c('y'),
		state.names=c('x1', 'x2')),
	regexp="Matrix values.load has 1 columns and 2 state.names (i.e., column names).\nHint: they should match.\nNot even the King of Soul could set these column names.",
	fixed=TRUE)

# Error
testthat::expect_error(
	me2 <- prep.measurement(values.load=matrix(c(1, 0), 1, 1),
		obs.names=c('y1', 'y2'),
		state.names=c('x1')),
	regexp="Matrix values.load has 1 rows and 2 obs.names (i.e., row names).\nHint: they should match.\nNot even the King of Pop could set these row names.",
	fixed=TRUE)

# Error
testthat::expect_error(
	me3 <- prep.measurement(values.load=matrix(c(1, 0), 1, 1),
		values.exo=matrix(c(1, 0), 1, 1),
		obs.names=c('y'),
		state.names=c('x1'),
		exo.names=c('u1', 'u2')),
	regexp="Matrix values.exo has 1 columns and 2 exo.names (i.e., column names).\nHint: they should match.\nNot even the King of Soul could set these column names.",
	fixed=TRUE)

# Error
testthat::expect_error(
	me4 <- prep.measurement(values.load=matrix(c(1, 0), 1, 2),
		values.int=matrix(c(1, 1), 2, 1),
		obs.names=c('y'),
		state.names=c('x1', 'x2')),
	regexp="Matrix values.int has 2 rows and 1 obs.names (i.e., row names).\nHint: they should match.\nNot even the King of Pop could set these row names.",
	fixed=TRUE)

