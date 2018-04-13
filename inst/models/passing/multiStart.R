#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2017-11-29
# Filename: multiStart.R
# Purpose: Check that multiple starting values are handled well.
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------


require(dynr)


# Sy-Miin found that this gave an incorrect error
#  No error should be given.
initial <- prep.initial(
	values.inistate=c(-1,.5),
	params.inistate=c('mu_readLevel',
					  'mu_readSlope'),
	values.inicov=matrix(c(.2,.02,
						   .02,.1),byrow=TRUE,ncol=2),
	params.inicov=matrix(c("v_11","c_12",
						   "c_12","v_22"),byrow=TRUE,ncol=2))


#------------------------------------------------------------------------------


