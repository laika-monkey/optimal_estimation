#!/usr/bin/env python
#############################################################################_80
# file:		optimal_estimation.py
# author:	Walter R. Sessions, 2016.01.12
# purpose:	do post processing on cases	
################################################################################
################################################################################

from libtools import dbg
from subprocess import call
import numpy as np
from numpy import matrix
from numpy.linalg import inv, svd
from scipy.linalg import sqrtm
from scipy.io import netcdf
import os
import math
import time
import sys
from optimal_estimation import add_ssps
from oe_real_hs3 import namelist, check_missing, write_retrieved_qsub
from oe_real_hs3 import plot_jeff_qsub, plot_real_qsub

#_select flight segements VVVVVVVV change that to shift
from flight_namelists import experiments

#_pull out desired experiment
experiment = 'HS3'
segments   = experiments[experiment] 

DIR_LOG  = os.environ['LOG']
DIR_PROD = os.environ['PRODUCTS']
DIR_LBL  = 'LBL-RTM_hs3' 
DIR_TMP  = '/data/wsessions/TMP'
DEBUG    = 1

out_labels =	[
##	'NVACLIMO1_tau-ref_M2_DV27_CLD_NOCTP_ATP-NAAPS_3KM',		#_done
##	'NVACLIMO1_tau-ref_M2_DV27_CLD_CTP-SHIS_ATP-NAAPS_3KM',
##	'NVACLIMO1_tau-ref_M2_DV27_CLD_NOCTP_NOATP_3KM',
##	'NVACLIMO1_tau-ref_M2_DV27_CLD_CTP-SHIS_NOATP_3KM',
#	'NVACLIMO1_tau-ref_M2_DV27_NOCLD_NOCTP_NOATP_3KM',
	'NVACLIMO1_tau-ref_M2_DV27_NOCLD_NOCTP_ATP-NAAPS_3KM',
	]


if __name__ == '__main__':
	from hs3_utils import Flight_segment
	from qsubmissions import lblrtm_hs3_lim, lbldis_hs3_lim
	from hs3_utils import find_shiscplgdas
	import re

	############################################################################
	#_PLOTTING ONLY
	if 0:	#_WITH RETRIEVALS
	  ppath = namelist.get('dir_plot')
	  namelist.update({' dir_plot' : os.path.join(ppath,out_label,'DASHBOARD')})
	  for dtg, values in segments.iteritems():
		#_plot comparisons between OD from SHIS and CPL 
		# SEPARATE! CANNOT RUN WITH real_case_qsub()
		dbg(dtg)
		final = 'lbldis_output.{0}.final.cdf'.format(out_label)
		clear = 'lbldis_output.{0}.clear.cdf'.format(out_label)
		namelist.update({	'lbldis_output'	: final, 
							'lbldis_clear'	: clear, })
		namelist.update(values)
		flight = Flight_segment(dtg)
		flight.plot(**namelist)			#_replot everything

	if 0:	#_GENERATE SUMMARY
	 for out_label in out_labels:	
	  for dtg, values in segments.iteritems():
		namelist['dtg'] = dtg
		final = 'lbldis_output.{0}.final.cdf'.format(out_label)
		clear = 'lbldis_output.{0}.clear.cdf'.format(out_label)
		namelist.update({	'lbldis_output'	: final, 
							'lbldis_clear'	: clear,
							'out_label'		: out_label })
		namelist.update(values)

		dbg(''); dbg(dtg)
	#	try:
	#	write_retrieved_qsub(**namelist)
		plot_jeff_qsub(**namelist)
		plot_real_qsub(**namelist)
	#	except:
	#		dbg((dtg, 'FAILED TO PLOT'))


	#_FULL CAMPAIGN SUMMARIES
	if 0:
		#_send to method that reads in all retrieval data
		# with indiv methods to handle plotting	
		plot_real_summary(**namelist)


	if 1: #_SEE WHAT'S COMPLETE
	 for out_label in out_labels:	
	  dbg(out_label)
	  for dtg, values in segments.iteritems():
		final = 'lbldis_output.{0}.final.cdf'.format(out_label)
		clear = 'lbldis_output.{0}.clear.cdf'.format(out_label)
		namelist.update({	'lbldis_output'	: final, 
							'lbldis_clear'	: clear,
							'out_label'		: out_label })
		namelist.update(values)
		check_missing(dtg, **namelist)		
			

	if 0:
	  for dtg, values in segments.iteritems():
		#_plot comparisons between OD from SHIS and CPL 
		# SEPARATE! CANNOT RUN WITH real_case_qsub()
		plot_real(**namelist)
