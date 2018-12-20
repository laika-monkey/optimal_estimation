#!/usr/bin/env python
#############################################################################_80
# file:		optimal_estimation.py
# author:	Walter R. Sessions, 2014.10.04
# purpose:	Attempts to use optimal estimation techniques taken from Clive
#			Rodgers's Inverse Methods for Atmospheric Sounding, Vol II to
#			retrieve profile information from hyperspectral infrared imagers
#			such as S-HIS and AIRS. Code also includes cases to run tests on
#			pseudo-obs.
#
#			Much of this code is intended to be split out into separate 
#			files, such as the namelist and hs3_2013 dictionaries. All
#			in here just for testing.
################################################################################
################################################################################


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
from flight_namelists import experiments
segments = experiments['HS3']
from libtools import dbg


if __name__ == '__main__':
	from hs3_utils import Flight_segment
	from qsubmissions import lblrtm_hs3_lim, lbldis_hs3_lim
	from hs3_utils import find_shiscplgdas
	import re
	from oe_real_hs3 import write_retrieved_qsub
	from oe_real_hs3 import plot_jeff_qsub 
	from oe_real_hs3 import plot_jeffs_qsub 
	from oe_real_hs3 import plot_jeffz_qsub 
	from oe_real_hs3 import plot_real_qsub
	from oe_real_hs3 import namelist
	
	#_yank out label
##	namelist = {}
	out_labels = [

		#_LAST TWO
	#	'NVACLIMO1_tau-ref_M2_DV26_CLD_CTP-NAAPS_3KM',
	##	'NVACLIMO1_tau-ref_M2_DV27_NOCLD_CTP-NAAPS_3KM',
	#	'NVACLIMO1_tau-ref_M2_DV27_CLD_NOCTP_NOATP_3KM',
	#	'NVACLIMO1_tau-ref_M2_DV27_CLD_NOCTP_ATP-NAAPS_3KM',
		'NVACLIMO1_tau-ref_M2_DV27_CLD_CTP-SHIS_ATP-NAAPS_3KM',
	###	'NVACLIMO1_tau-ref_M2_DV27_CLD_CTP-SHIS_NOATP_3KM',


##		'NVACLIMO1_tau-ref_M2_DV27_CLD_NOCTP_3KM',
	##	'NVACLIMO1_tau-ref_M2_DV26_NOCLD_CTP-NAAPS_3KM',

	##	'NVACLIMO1_tau-ref_M2_DV27_NOCLD_CTP-NAAPS_3KM',
	#	'NVACLIMO1_tau-ref_M2_DV26_NOCLD_NOCTP_3KM',
	##	'NVACLIMO1_tau-ref_M2_DV27_NOCLD_NOCTP_3KM',

##	'NVACLIMO1_tau-ref_M2_DV27_CLD_CTP-SHIS_ATP-NAAPS_3KM'
##	'NVACLIMO1_tau-ref_M2_DV27p_CLD_CTP-SHIS_ATP-NAAPS_3KM'
		]

	############################################################################
	#_real_#####################################################################
	############################################################################
	#_LBL-RTM, REAL FLIGHT LOOP (must be run before retrieval)
	if 0:	#_GENERATE SUMMARY
	 for out_label in out_labels:
	  for dtg, values in segments.iteritems():
		namelist['dtg'] = dtg
		final = 'lbldis_output.{0}.final.cdf'.format(out_label)
		clear = 'lbldis_output.{0}.clear.cdf'.format(out_label)
		namelist.update({	'lbldis_output'	: final, 
							'lbldis_clear'	: clear,
							'out_label'		: out_label,
							'smooth'		: False })
		namelist.update(values)
	
		dbg(''); dbg(dtg)
		write_retrieved_qsub(**namelist)
	#	plot_jeff_qsub(**namelist)
	#	plot_real_qsub(**namelist)

	#_PLOT IN ONE
	if 1:
	  for dtg, values in segments.iteritems():
		namelist['dtg'] = dtg
		namelist.update({	'out_labels'	: out_labels,
							'smooth'		: False,
							'bt_thresh'		: 7.,  })
		namelist.update(values)
	
		dbg(''); dbg(dtg)
	#	plot_jeffs_qsub(**namelist)
		plot_jeffz_qsub(**namelist)

	if 0:	#_JUST WANT DASHBOARDS OF ONE SEGMENT FOR ALL EXPS
	  ppath = namelist.get('dir_plot')
	  exps = []
	  exps = ['NVACLIMO1_M2_DV{0}_{1}_{2}_3KM'.format(d,c,p)
				for d in ['26','27'] 
				for c in ['CLD','NOCLD'] 
				for p in ['CTP-NAAPS','NOCTP']]
	  for e in exps:
	    namelist.update({'dir_plot' : os.path.join(ppath,e,'DASHBOARD')})
	    for dtg, values in segments.iteritems():
			if dtg != '20130824170000': continue
			#_plot comparisons between OD from SHIS and CPL 
			# SEPARATE! CANNOT RUN WITH real_case_qsub()
			dbg(dtg)
			final = 'lbldis_output.{0}.final.cdf'.format(e)
			clear = 'lbldis_output.{0}.clear.cdf'.format(e)
			namelist.update({	'lbldis_output'	: final, 
								'lbldis_clear'	: clear,
								'out_label' : e })
			namelist.update(values)
			flight = Flight_segment(dtg)
			try:
				flight.plot(**namelist)			#_replot everything
			except:
				print 'cannot plot flight for', out_label

	#_PLOTTING ONLY
	if 0:	#_WITH RETRIEVALS
	  ppath = namelist.get('dir_plot')
	  namelist.update({'dir_plot' : os.path.join(ppath,out_label,'DASHBOARD')})
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
	if 0:	#_NO RETRIEVALS
	  for dtg, values in segments.iteritems():
		#_plot comparisons between OD from SHIS and CPL 
		# SEPARATE! CANNOT RUN WITH real_case_qsub()
		dbg(dtg)
		namelist.update(values)
		namelist['out_label'] = 'EVENT-ONLY'	#_change to whatever
		flight = Flight_segment(dtg)
		flight.plot_flight(**namelist)	#_flight only
