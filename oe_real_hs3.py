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

#_select flight segements VVVVVVVV change that to shift
from flight_namelists import experiments

if 'DISPLAY' not in os.environ:
	import matplotlib
	matplotlib.use('Agg')

#_pull out desired experiment
experiment = 'HS3'
segments = experiments[experiment]
'''
To turn off uncertainties
	'model'			set 'uncertainty' to 0
	'instrument'	set 'uncertainty' to 0
'''

DIR_LOG		= os.environ['LOG']
DIR_PROD	= os.environ['PRODUCTS']
DIR_LBL		= 'LBL-RTM_hs3' 
DIR_TMP		= os.environ['WORK'] + '/kwargs'
DIR_WORK	= os.environ['WORK']
DIR_SBATCH	= os.environ['WORK'] + '/batch_scripts'
DEBUG		= 1

'''
NVACLIMO1_tau-ref_M2_DV26_CLD_CTP-SHIS_ATP-NAAPS_3KM

TESTS:
NVACLIMO1_tau-ref_M2_DV26_CLD_CTP-NAAPS_3KM			done, 01.12	
NVACLIMO1_tau-ref_M2_DV27_CLD_CTP-NAAPS_3KM			done, 01.12

NVACLIMO1_tau-ref_M2_DV26_CLD_NOCTP_3KM				done, 01.12
NVACLIMO1_tau-ref_M2_DV27_CLD_NOCTP_3KM				done, 01.12
	
NVACLIMO1_tau-ref_M2_DV26_NOCLD_CTP-NAAPS_3KM		resubbed 12.31, w/ all
NVACLIMO1_tau-ref_M2_DV27_NOCLD_CTP-NAAPS_3KM		resubbed 12.31, w/ all	

NVACLIMO1_tau-ref_M2_DV26_NOCLD_NOCTP_3KM			resubbed 12.31, w/ all
NVACLIMO1_tau-ref_M2_DV27_NOCLD_NOCTP_3KM			resubbed 12.31, w/ all

NVACLIMO1_tau-ref_M3_DV27_CLD_CTP-NAAPS_3KM

NVACLIMO1_tau-ref_M2_DV27_CLD_CTP-NAAPS_1KM			resubbed 12.18, w/ all
'''

#_used instead of config_file
namelist = {
	#_where to get temperary files
	'hostname'		: 'iris.ssec.wisc.edu',

	#_if directory exists, skip
	'rerun'			: False,  
  
	#_smooth summary lines
	'smooth'		: False, #True,

	#_correct biases
	'bias_correct'	: False, 

	#_ok, see if this breaks everything
	'experiment'	: experiment,	#_used for LBL I/O naming 
##	'out_label'		: 'NVACLIMO1misr',	
	'out_label'		: 'NVACLIMO1',	
#	'out_label'		: 'KAOLINITE',	
#	'out_label'		: 'GYPSUM',	
#	'out_label'		: 'QUARTZ',	
#	'out_label'		: 'DUST',	
								# the loop designates it to be used
								# only as a kind of top level label
								# when looping over flight segments.

	#_these are separate so that uncertainty
	# can be looped via jv without attempts to retrieve
	# and covar matrix generated
	'state_vars'	: ['tau','ref'],	#_vars to retrieve
	'uncert_vars'	: ['z','thickness','surf_emis','surf_temp','rel_hum'],

	#_location of naaps data, only used if set in dynamic CTP
##	'dir_naaps'		: os.path.join(DIR_WORK, 'NRL', 'NVA_CLIMO1'),
	'dir_naaps'		: os.path.join(DIR_PROD, 'NRL', 'NVA_CLIMO1'),
##	'dir_naaps'		: os.path.join(DIR_PROD, 'NRL', 'NVA_CLIMO1misr'),

	#_location of input radiances and output directory for fov plots
	'dir_shis'		: os.path.join(DIR_PROD, experiment.lower()), 
	'dir_out'		: os.path.join(DIR_WORK, experiment.lower() + '_out'), 
	'dir_plot'		: os.path.join(DIR_PROD, 'plots'), 

	#_input sources (incomplete)
	'surf_temp'			: 'GDAS',	#_GDAS || LBL-ATM || LBL-ATM_old
	'profile_source'	: 'GDAS',	#_GDAS || SHIS || ECMWF || (only gdas works)
									# used for rel_hum jacobian
	'surf_temp_var'		: 'skt',	#_sst || skt, only used for ECMWF
##	'surf_temp_src'		: 'ECMWF',	#_GDAS || LBL-ATM || ECMWF || SHIS (n/a)
									# NVM, use apriori['surf_temp'] in same way,
									# which is built to permit arbitrary sfc_tmp
									# values to be passed.

	#_use DR cloud top pressure for layer
	'dynamic_CTP'		: False, #['SHIS'],  # ['SHIS'], False,
	'dynamic_ATP'		: ['NAAPS'], # ['NAAPS'], False,

	#_use previous field of view as apriori if uncertainty under limit
#	'persistence'		: { 'tau' : 0.2, 'ref' : 0.2 },	
	'persistence'		: False, 

	#_format of lblrtm directory, [experiment]_[dtg0.dtg1]_fov[###]
	'dir_lblrtm_fmt':os.path.join(DIR_PROD,DIR_LBL,'{4}/{3}_{0}.{1}_fov{2:06d}'),

	# 2. VARIABLES
#	'dv'			: -26,			#_microwindow (window, no O3+ 1100)
	'dv'			: -27,			#_microwindow (window, no O3)

	'obs_angle'		: 0,			# 180 (zenith); 0 (nadir) 

	#_variables that are in JV but not SV are sources of uncertainty
	# that we are not attempting to back out

##	#_REAL
##	############################################################################
##	############################################################################
##	'dtg'			: '20130821021538', 
##	'fidx'			: range(425, 426),
##	#_END_REAL
##	############################################################################
##	############################################################################

	'apriori' : [
					{	'type'		: 'aerosol',
						'tau'		: 0.01, 
						'ref'		: 0.5,		#_in the middle of af dust 
						'z'			: 3.0,
						'z_top'		: 3.01,
						'thickness' : 0.5,		#_if thickness retrieved, z_top
												# is ignored
						'ssp_db'	: 'ssp_db.mie_gypsum.lognormal_sigma_0p699',
						'ref_wn'	: -1,		#_-1 : geo limit, or waven 900.
						},
					{	'type'		: 'aerosol',
						'tau'		: 0.01, 
						'ref'		: 0.5,	
						'z'			: 3.0,
						'z_top'		: 3.01,
						'thickness' : 0.5,	
						'ssp_db' : 'ssp_db.mie_kaolinite.lognormal_sigma_0p699',
						'ref_wn'	: -1,	
						},
			#		{	'type'		: 'aerosol',
			#			'tau'		: 0.01, 
			#			'ref'		: 0.5,	
			#			'z'			: 2.0,
			#			'z_top'		: 2.01,
			#			'thickness' : 0.5,	
			#			'ssp_db'	: 'ssp_db.mie_quartz.lognormal_sigma_0p699',
			#			'ref_wn'	: -1,	
			#			},
			##		{	'type'		: 'aerosol',
			##			'tau'		: 0.01, 
			##			'ref'		: 0.5,	
			##			'z'			: 3.0,
			##			'z_top'		: 3.01,
			##			'thickness' : 0.5,	
			##			'ssp_db'	: 'ssp_db.shettle_dust.gamma_sigma_0p100',
			##			'ref_wn'	: -1,	
			##			},
		#			{	'type'		: 'cloud',
		#				'tau'		: 0.01, #1.0, 
		#				'ref'		: 1.0,
		#				'z'			: 1.00,
		#				'z_top'		: 1.01,
		#				'thickness' : 0.5,	
		#				'ssp_db'	: 'ssp_db.mie_wat.gamma_sigma_0p100',
		#				'ref_wn'	: -1,
		#				},
				],


	#_by how much are we perturbing jacobian
	'sensitivity'	: {
						'tau'		: 0.01,	#_optical depth
						'ref'		: 0.10, # 0.1,	#_effective radius
						'z'			: 0.51,	#_cloud height
						'thickness' : 1,	#_number of layers (this gets weird)
						'surf_temp'	: 0.5,	#_surface temperature
						'surf_emis'	: 0.1,	#_surface emissivity
						'rel_hum'	: 5.0,	#_percent relative humidity
						'ssp_dist'	: 0.1,	#_portion to increment portion of OD
						},

	#_basically stating that we're dumb ONLY USE FOR STATE VARS (2.19.15)
	# for uncert_vars, these are multiplied and can be used as switches
	#_for jacobian	: 1|0
	#_for state		: unit values
	'uncertainty'	: {				#_uncertainty in a priori values
						'tau'		: 10.,
						'ref'		: 1,
						'z'			: 1,
						'thickness' : 1,	#_layers
						'surf_temp'	: 1,	#_K
						'surf_emis'	: 1,	#_fractional 
						'rel_hum'	: 1,	#_percent relative humidity
						'model'		: 0,	#_flat factor to apply, leave off 
						'instrument': 1,	#_0 for off, 1 for on
						'ssp_dist'	: 0.5,	
						},

	#_ability to add arbitrary layers to run (NOT FULLY IMPLEMENTED)
	'static'		: True,						#_Try adding static layer?
	'static_layer'	: [{'type'		: 'water_cloud',
						'tau'		: 1.0, 
						'ref'		: 1.0,		#_in the middle of af dust 
						'z'			: 1.00,
						'z_top'		: 1.01,
						'thickness' : 0.5,		#_if thickness retrieved, z_top
												# is ignored
						'dbnum'		: 5,
						'ref_wn'	: -1,}],	#_-1 : geo limit, or waven 900.

	'sensor' : 'SHIS',		# define the sensor used for the retrieval.
							# Possible options are SHIS and AIRS
							# NOTE: AIRS not implemented yet.

	#_how many processors to take over
	'nproc'	: 10,

	#_remove input/output files in lbldis
	'clean'	: True,	

	'NETCDF_LIB' : '/opt/netcdf4/4.1.3-intel-14.0-2/lib/',	# netcdf lib 
															# for lbldis
	'ssp_db_files' : [],
	}
				
#_add these for total column/surf values
namelist['apriori'][0].update({ 
			'surf_emis'	: [[100, 3000], [0.985, 0.985]],
			'rel_hum'	: 1.0, })

#_add in spectral property file options
add_ssps(**namelist)

#_add levels at sensitivity of z
namelist.update({ 'dummy_levs' : namelist['sensitivity']['z'] })

#_if retrieving thickness, overwrite z_top
if 'thickness' in namelist['uncert_vars']:
	for layer in namelist['apriori']:
		layer['z_top'] = layer['z'] + layer['thickness']


#############################################################################_80
#_real/simulated_cases_#########################################################
################################################################################


def real_case_sbatch(fidx=None, dir_plot='.', bias_correct=False,
	hostname=False, **kwargs):
	'''
	Given a list of field of view indices and a format for
	the output from previously run LBL-RTM, submits jobs to
	the queue to run LBL-DIS. A ton of options are hidden in 
	kwargs.

	dtg				str(14),	dtg within range of flight
	dir_lblrtm_fmt	str,		format to be used for input lblrtm dir
	out_label		str,		used for output imagery	
	experiment		str,		unused
	'''

	from hs3_utils import Flight_segment as F
	from libtools import mkdir_p, capture_stdout
	import pickle
	import re

	#_get process number
	pid = os.getpid()
	rerun = kwargs.get('rerun')

	#_pull out some vars
	out_label = kwargs.get('out_label')
	dir_lblrtm_fmt = kwargs.get('dir_lblrtm_fmt')
	experiment = kwargs.get('experiment')

	#_create directory for plots
	ppath = os.path.join(dir_plot, out_label)
	mkdir_p(ppath)
	kwargs.update({'dir_plot' : ppath})

	#_do entire flight if not specified
	flight = F(**kwargs)
	fidx = xrange(flight.size) if fidx is None else fidx

	#_check if bias correction desired
	if bias_correct:
		bias = _read_bias(**kwargs)
		kwargs.update({'bias_correction' : bias})

	#_build sbatch script
	scpt = os.path.expanduser('~/lib/run_oe_hs3.py')

	#_last job ID number
	last_jid = False

	#_Loop over each fov (in this case they're just variable TAU/REF/Z opts)
	for ii, i in enumerate(fidx):

		#_build output directory name. Args really should be passed
		# if lblrtm_fmt is a thing. A half effort otherwise
		dir_lblrtm = dir_lblrtm_fmt.format(	flight.dtg0, flight.dtg1, i,
											experiment, flight.dtg0[2:])
		kwargs.update({'dir_lblrtm' : dir_lblrtm})

		fov_label = '{0}.{1}'.format(dir_lblrtm.split('/')[-1], out_label)
		out = \
'''#!/bin/bash
#SBATCH --job-name={2}
#SBATCH --partition=all
#SBATCH --share
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/odyssey/scratch/%u/logs/{0}.{2}-%A.txt
#SBATCH --nice=1000
module purge
module load license_intel
module load impi
module load intel/15.0-2
module load hdf5/1.8.14
module load netcdf4/4.3.3
module load anaconda27/base_2015a_sharppy13
echo ${{HOSTNAME}}
export TMPDIR={1}/${{SLURM_JOB_NAME}}.${{SLURM_JOB_ID}}
mkdir -p $TMPDIR
'''.format(out_label, '{0}/tmp/'.format(os.environ['WORK']), fov_label)
		out += '\n'.join(['export {0}={1}'.format(var, os.environ[var]) \
			for var in ['PYTHONPATH', 'PRODUCTS', 'WORK', 'PATH', 'LOG']])

		#_generate kwarg file that scrpt will load to get settings
		fmt = 'kwgs.{4}.{2}.{3}.{0:04d}.{1:07d}'
		fmt_args = (i, pid, flight.dtg0, flight.dtg1, out_label)
		file_kwarg = os.path.join(DIR_TMP, fmt.format(*fmt_args))
		pickle.dump(kwargs, open(file_kwarg, 'wb')) 

		#_add leading hostname to kwargs
		if hostname:
			file_kwarg = '{0}:{1}'.format(hostname, file_kwarg)

		#_build actual script call string
		out += '\nsource activate my_root'
		out += '\necho `which python`'
		out += '\n' + ' '.join((scpt, str(i), file_kwarg))
		out += '\nrm -rf $TMPDIR'

		#_KLUDGE~DELETE. Skips retrievals already performed.
	#	outfile = os.path.join(dir_lblrtm, kwargs.get('lbldis_output'))
		outfile = os.path.join(dir_lblrtm, kwargs.get('lbldis_output'))
		outfile = os.path.join(dir_lblrtm, kwargs.get('lbldis_input'))
		if os.path.exists(outfile) and not rerun:
			dbg((outfile, 'exists... skipping'))
			continue

		#_write batch script
		sname = '{1}/sbatch.{0}'.format(fov_label, DIR_SBATCH)
		with open(sname, 'w') as f:
			print out
			f.write(out)

		#_if using persistence, only run after previous one has exited
		#_Every tenth restarts process.
		dependency = ''
		if type(kwargs['persistence']) == dict and last_jid and ii % 10:
			#_first time through, will not pass because last_jid iti == False
			dependency = '--dependency=afterany:{0:d}'.format(last_jid)

		#_put it all together and submit
		cmd = ' '.join(('sbatch', dependency, sname))
		#_GET BATCH ID AND PASS TO NEXT AS OK IF PERSISTENCE SET
		dbg(cmd)
	#	os.system(cmd)
		stdout = capture_stdout(cmd)

		#_get ID
		last_jid = int(re.search('Submitted batch job (\d+)', stdout[0]).group(1))


def _read_bias(fname='/data/wsessions/csbias.dat', **kwargs):
	''' read in bias file and put in to dict keyed by mw number '''
	with open(fname, 'r') as f:
		for i, line in enumerate(f.readlines()):

			#_skip comment
			if not i:
				bias = {}
				continue

			tmp = [t.strip() for t in line.split(',')]
			if len(tmp) == 1:
				key = float(tmp[0].strip())
				bias.update({ key : {} })
			elif len(tmp) == 2:
				bias[key].update({ float(tmp[0]) : float(tmp[1]) })

	return bias


def write_retrieved_cod_qsub(dtg=None, out_label='HS3',
	state_vars=['tau'], experiment='HS3',**kwargs):
	'''
	Produce .retrieved.cod file based upon flight object for previously
	run retrieval. This method was necessitated when moving from real_case()
	to real_case_qsub() because the latter does not generate one
	at runtime.

	dtg			str{14},			Full date-time-group for flight segement	
	flight		hs3.Flight_segment,	Object with all the goodies.	
	'''

	def sum_tau(lbldis):
		''' sum all lbldis optical depth outputs. MOVE THIS METHOD TO CLASS '''
		#_keywords in aerosol ssp discrimination
		ssp_cloud = ['wat', 'ice']

		#_initialize tau
		tau = 0

		#_pull out ssp labels
		ssp_labels = [ssp.split('.')[1] for ssp in lbldis.ssp_db]
		ssp_labels = [ssp.split('_')[1] for ssp in ssp_labels]

		#_loop over layers
		for i, layer in enumerate(lbldis.layers):

			#_check if aerosol ssp
			if ssp_labels[layer['dbnum']] in ssp_cloud:
				tau += layer['tau'][0]	
		
		#_return total AOD
		return tau

	from hs3_utils import Flight_segment
	from lblrtm_utils import Lbldis_input
	from time import sleep

	fidx = kwargs.get('fidx', None)
	dir_lblrtm_fmt = kwargs.get('dir_lblrtm_fmt', 'lbldis_out_{0}')
	
	#_DELETE FOR NOW ONLY ALLOW TAU BECAUSE
	# 2016.06.09 Revisiting this, I do not know why this was restricted
	#			it may be because we're only writing one habits tau out...
	#			which is dumb... 
	state_vars = ['tau']

	#_read in flight segment file
	flight = Flight_segment(dtg=dtg, **kwargs)

	#_if no indices are specified, attempt all
	fidx = range(flight.size) if fidx is None else fidx
	nfov = len(fidx)

	#_kludge to speed this up
	missing_no, total = check_missing(dtg, **kwargs)
	if missing_no == total:
		dbg('No retrievals for {0}'.format(dtg))
		return False

	else:
		pass

	#_name of output file
	file_out = flight.FLIGHT_file.replace('.nc',
		'.{0}.retrieved'.format(out_label))

	#_get length of state dimension
	nv = len(state_vars)

	#_initialize header labels and their format 
	header = ['fov']
	fmt = ['{0:10d}']			#_FOV NUMBER
	fmt_header = ['{0:>10s}']	#_fov label
	for j, v in enumerate(state_vars):
		fmt.append('{{{0:d}:10.5f}}'.format(j+1))	#_ fov, tau, ref, tau, ref
		fmt_header.append('{{{0:d}:>10s}}'.format(j+1))
		header.append(v)

	#_CPL header, for length of tau_type, tau_532
	for j in range(20): #state_vars):	#_uncertainty format (?)
		fmt.append('{{{0:d}:10.5f}}'.format(nv+j+1))

	#_put them together
	fmt_header = ','.join(fmt_header) + '\n' 
	fmt = ','.join(fmt) + '\n'

	dbg(file_out)

	#_init output list
	lines_out = []

	with open(file_out,'w') as f:

		#_write header
		f.write(fmt_header.format(*header))

		#_loop over indices and write output
		for i in fidx:

			#_build path names
			dtgfov = (flight.dtg0, flight.dtg1, i, experiment, flight.dtg0[2:])
			dir_lblrtm = dir_lblrtm_fmt.format(*dtgfov)
			lbldis_out = 'lbldis_input.{0}.final'.format(out_label)
			lbldis_out = os.path.join(dir_lblrtm, lbldis_out)

			#_initialize output line tuple with fov number
			tup = [i]

			#_read in final lbldis output file
			try:
				lbldis_out = Lbldis_input(lbldis_out)
		
				#_sum up cloud taus
				tau = sum_tau(lbldis_out)

			except IOError:
				dbg('Missing {0}'.format(lbldis_out), 3)
				tau = -9999.
		#		continue

			tup.extend([tau])	#_this is where a generator with ref would go
			tup.extend(flight[i].CPL_tau_532)
			tup.extend(flight[i].CPL_tau_type)

			#_write line to file
			lines_out.append(tup) 
		#	f.write(fmt.format(*tup))

		#_write output
		[f.write(fmt.format(*tup)) for tup in lines_out]
		nret = len(lines_out)
	
		#_give a status
		dbg('Retrievals vs. Total FOVS {0:>6d} / {1:>6d}'.format(nret, nfov))
		sleep(1)

	#_if no retrievals, delete file
	if len(lines_out) == 0:
		dbg('No retrievals for {0}'.format(dtg))
		return False

	#_return output filename
	return file_out

		

def write_retrieved_qsub(dtg=None, out_label='HS3',
	state_vars=['tau'], experiment='HS3',**kwargs):
	'''
	Produce .retrieved file based upon flight object for previously
	run retrieval. This method was necessitated when moving from real_case()
	to real_case_qsub() because the latter does not generate one
	at runtime.

	dtg			str{14},			Full date-time-group for flight segement	
	flight		hs3.Flight_segment,	Object with all the goodies.	
	'''
	def sum_tau(lbldis):
		''' sum all lbldis optical depth outputs. MOVE THIS METHOD TO CLASS '''
		#_keywords in aerosol ssp discrimination
		ssp_aero = ['dust','sand','gypsum','kaolinite','quartz']

		#_initialize tau
		tau = 0

		#_pull out ssp labels
		ssp_labels = [ssp.split('.')[1] for ssp in lbldis.ssp_db]
		ssp_labels = [ssp.split('_')[1] for ssp in ssp_labels]

		#_loop over layers
		for i, layer in enumerate(lbldis.layers):

			#_check if aerosol ssp
			if ssp_labels[layer['dbnum']] in ssp_aero:
			##	dbg(('AEROSOL', ssp_labels[layer['dbnum']]))
				tau += layer['tau'][0]	
		
		#_return total AOD
		return tau

	from hs3_utils import Flight_segment
	from lblrtm_utils import Lbldis_input
	from time import sleep

	fidx = kwargs.get('fidx', None)
	dir_lblrtm_fmt = kwargs.get('dir_lblrtm_fmt', 'lbldis_out_{0}')
	
	#_DELETE FOR NOW ONLY ALLOW TAU BECAUSE
	state_vars = ['tau']

	#_read in flight segment file
	flight = Flight_segment(dtg=dtg, **kwargs)

	#_if no indices are specified, attempt all
	fidx = range(flight.size) if fidx is None else fidx
	nfov = len(fidx)

	#_kludge to speed this up
	missing_no, total = check_missing(dtg, **kwargs)
	if missing_no == total:
		dbg('No retrievals for {0}'.format(dtg))
		return False
	else:
		pass

	#_name of output file
	file_out = flight.FLIGHT_file.replace('.nc',
		'.{0}.retrieved'.format(out_label))

	#_get length of state dimension
	nv = len(state_vars)

	#_initialize header labels and their format 
	header = ['fov']
	fmt = ['{0:10d}']			#_FOV NUMBER
	fmt_header = ['{0:>10s}']	#_fov label
	for j, v in enumerate(state_vars):
		fmt.append('{{{0:d}:10.5f}}'.format(j+1))	#_ fov, tau, ref, tau, ref
		fmt_header.append('{{{0:d}:>10s}}'.format(j+1))
		header.append(v)

	#_CPL header, for length of tau_type, tau_532
	for j in range(20): #state_vars):	#_uncertainty format (?)
		fmt.append('{{{0:d}:10.5f}}'.format(nv+j+1))

	#_put them together
	fmt_header = ','.join(fmt_header) + '\n' 
	fmt = ','.join(fmt) + '\n'

	dbg(file_out)

	#_init output list
	lines_out = []

	#_write header
	with open(file_out,'w') as f:

		#_write header
		f.write(fmt_header.format(*header))

		#_loop over indices and write output
		for i in fidx:

			#_build path names
			dtgfov = (flight.dtg0, flight.dtg1, i, experiment, flight.dtg0[2:])
			dir_lblrtm = dir_lblrtm_fmt.format(*dtgfov)
			lbldis_out = 'lbldis_input.{0}.final'.format(out_label)
			lbldis_out = os.path.join(dir_lblrtm, lbldis_out)

			#_initialize output line tuple with fov number
			tup = [i]

			#_read in final lbldis output file
			try:
				lbldis_out = Lbldis_input(lbldis_out)
		
				#_sum up aerosol taus
				tau = sum_tau(lbldis_out)

			except IOError:
				dbg('Missing {0}'.format(lbldis_out), 3)
				tau = -9999.
		#		continue

			tup.extend([tau])	#_this is where a generator with ref would go
			tup.extend(flight[i].CPL_tau_532)
			tup.extend(flight[i].CPL_tau_type)

			#_write line to file
			lines_out.append(tup) 
		#	f.write(fmt.format(*tup))

		#_write output
		[f.write(fmt.format(*tup)) for tup in lines_out]
		nret = len(lines_out)
	
		#_give a status
		dbg('Retrievals vs. Total FOVS {0:>6d} / {1:>6d}'.format(nret, nfov))
		sleep(1)

	#_if no retrievals, delete file
	if len(lines_out) == 0:
		dbg('No retrievals for {0}'.format(dtg))
		return False

	#_return output filename
	return file_out


def real_case(dtg=None, fidx=None, out_label='', desc='', experiment='hs3', 
	dir_lblrtm_fmt='lblrtm_{0}.{1}', **kwargs):
	'''
	dtg				str(14),	dtg within range of flight
	dir_lblrtm_fmt	str,		format to be used for input lblrtm dir
	out_label		str,		used for output imagery	
	'''
	from hs3_utils import Flight_segment as F
	from numpy import diag, matrix
	from lblrtm_utils import microwindow_average, microwindow_stdev
	from time import sleep
	import matplotlib.pyplot as plt
	from libtools import shrink_ticks, setup_groups, mkdir_p
	from os import fork, waitpid
	from multiprocessing import Process, Pipe

	raise RuntimeError, 'stop using this one'

	#_simulated input file name
	flight = F(dtg=dtg)

	file_tru = flight.FLIGHT_file.replace('.nc', '.txt')	#_true value file
	file_out = flight.FLIGHT_file.replace('.nc',			#_retr output
		'.{0}.retrieved'.format(out_label))

	#_get values being tested
	state_vars = kwargs.get('state_vars')

	#_create directory for plots
	ppath = 'PLOTS_{0}'.format(out_label)
	mkdir_p(ppath)

	#_print header
	header = ['fov']

	#_output format for results
	fmt = ['{0:10d}']			#_FOV NUMBER
	fmt_header = ['{0:>10s}']	#_fov label
	for j, v in enumerate(state_vars):
		fmt.append('{{{0:d}:10.5f}}'.format(j+1))
		fmt_header.append('{{{0:d}:>10s}}'.format(j+1))
		header.append(v)

	#_for length of tau_type, tau_532
	for j in range(20): #state_vars):	#_uncertainty format (?)
		nv = len(state_vars)
		fmt.append('{{{0:d}:10.5f}}'.format(nv+j+1))

	fmt = ','.join(fmt)
	fmt += '\n'
	fmt_header = ','.join(fmt_header) + '\n' 

	#_write header
	with open(file_out,'w') as f:
		f.write(fmt_header.format(*header))

	#_do entire flight if not specified
	fidx = xrange(flight.size) if fidx is None else fidx

	#_split fov objects into forkable groups
	groups = setup_groups(fidx, **kwargs)

	#_Loop over each fov (in this case they're just variable TAU/REF/Z opts)
	for fovs in groups:
		nfov = len(fovs)
		pipe = [None]*nfov
		pope = [None]*nfov
		thread = [None]*nfov
		line_out = [None]*nfov
	
		for j, i in enumerate(fovs): #_list index = j, flight index = i
			dbg(('TESTING FOV NUMBER', i)) 
			fov = flight[i]
			wav = flight.SHIS_wavenumber

			#_update path to lblrtm run
			dtgfov = (flight.dtg0, flight.dtg1, i, experiment)
			kwargs.update({'dir_lblrtm' : dir_lblrtm_fmt.format(*dtgfov)}) 
	
			#_init pipe
			pipe[j], pope[j] = Pipe()

			#_initialize plotting, for first ten instances only
			fig, ax = plt.subplots(3,2)
			fig.suptitle(desc, size='xx-small')
			pname = '{1}_{0:06d}.png'.format(i, out_label) 
			pname = os.path.join(ppath, pname)
			kwargs.update({	'fig' : fig, 'ax' : ax, 'pname' : pname })
						
			if os.path.exists(pname):
				os.unlink(pname)

			#_calculate instrument covariance matrix
			dv = kwargs.get('dv')

			#_when using microwindows, average observation data
			if dv < 0:
				std, d		= microwindow_stdev(fov.SHIS_radiances, wav, dv) 
				nesr, d		= microwindow_average(flight.SHIS_HBB_NESR, wav, dv,
								error=True)
				y, wave_obs = microwindow_average(fov.SHIS_radiances, wav, dv)
	
				#_std deviation within microwindows applied to ref
				# blackbody 
				cov	= matrix(diag(nesr) + diag(std))	
			else:
				cov = matrix(diag(flight.SHIS_HBB_NESR.copy()))

			#_figure out how to best organize this
			# Returns converged solution and all.
			# Maybe switch to just all and pop
			args = (fov, wav, cov)
			kwgs = kwargs.copy()
			kwgs['pipe'] = pope[j]
			thread[j] = Process(target=optimal_estimation, args=args, 
								kwargs=kwgs)
			thread[j].start()

		#_collect output
		for j, i in enumerate(fovs):

			#_get output from optimal estimation
			x, all_x, x_std = pipe[j].recv() 

			#_init with fov number
			tup_retr = [i]

			#_format retrieved values for output
			tup_retr.extend([x[v] for v in state_vars])

			dbg(('FOV', j, i))
			#_retrieved into out line
			tup_retr.extend(flight[i].CPL_tau_532)
			tup_retr.extend(flight[i].CPL_tau_type)
			line_out[j] = fmt.format(*tup_retr)

			#_close thread
			thread[j].join()

		#_write output as we go	
		with open(file_out,'a') as f:
			[f.write(l) for l in line_out] 

	dbg(file_out)
	return file_out

################################################################################
#_plotting_#####################################################################
################################################################################

	
def get_all_residuals(out_label, fractional=False, threshold=0.15, **kwargs):
	'''
	read in all differences between shis-oe retrieval and cpl

	Residuals are SHIS - CPL, CPL optical depths include aerosol AND pbl types.
	fractional	bool,	Report residuals as fraction of CPL aod
	'''
	from numpy import append, vstack, empty, linspace 
	from hs3_utils import Flight_segment as F
	from libtools import epoch2iso, shrink_ticks

	residuals = empty(0)

	#_read in all retrievals and get residuals
	for dtg, values in segments.iteritems():
		dbg('')

		#_generate input file name
		flight = F(dtg=dtg, **kwargs) 
		fname = flight.FLIGHT_file.replace('.nc',
			'.{0}.retrieved'.format(out_label))

		#_read in data
		data = read_real_retrieval(fname, **kwargs)
		if data.size == 0:
			dbg('No data for {0}'.format(dtg))
			continue

		#_put CPL data into arrays	
		cpl_typ = vstack(data.CPL_type)	
		cpl_tau = vstack(data.CPL_tau)	

		#_trim down to number of layers possibly defined in type
		cpl_tau = cpl_tau[:,:5]
		cpl_typ = cpl_typ[:,::2]

		#_find fovs with any aerosol
		idx_pbl = (cpl_typ == 1).any(1)	#_index of PBL type
		idx_ele = (cpl_typ == 2).any(1) #_index of AEROSOL type
		idx_aer = append(idx_pbl[:,None], idx_ele[:,None], axis=1).any(1)#_aero
		idx_cld = (cpl_typ == 3).any(1) * idx_aer	#_aerosol w/ cloud
		idx_clr = (cpl_typ != 3).all(1) * idx_aer	#_aerosol w/o cloud

		#_find where not missing
		idx_here = cpl_tau >= 0

		#_get column totals of elevated and pbl aersol when not missing
		col_pbl = ((cpl_typ == 1) * idx_here * (cpl_tau)).sum(1) #max(1) #sum(1)
		col_ele = ((cpl_typ == 2) * idx_here * (cpl_tau)).sum(1) #max(1) #sum(1)
		cpl_aod = col_ele + col_pbl
		
		shis_aod = data.tau[cpl_aod > 1e-5] 
		cpl_aod = cpl_aod[cpl_aod > 1e-5] 

		print data.tau.size, shis_aod.size, cpl_aod.size, 'TESTING'

		#_get indices wherea bove threshold
		idx = cpl_aod > threshold

		#_get residual difference between CPL and S-HIS retrievals
		if not fractional:
			residuals = append(residuals, shis_aod[idx] - cpl_aod[idx])
		else:
			residuals = append(residuals, (data.tau - cpl_aod) / cpl_aod)
	
	return residuals


def plot_jeffz_qsub(dtg=None, fidx=None, out_label='', dz_cpl=25, smooth=True,
	dir_plot='.', out_labels=None, shis_thresh=1.25, cpl_thresh=0.05,
	bt_thresh=5, 
	plot_cpl=True,
	plot_aod=True, plot_modis=True, 
	plot_cod=False,
	plot_naaps=True, 
	plot_cmask=True,
	gradient_kibosh=False,	#_leave this off.  removes arb spikes.
	**kwargs):

	'''
	fname	str,	path to output file of optimal estimation
		
	comparison plot between CPL AOD retrieval and SHIS AOD retrieval
	'''
	import matplotlib.pyplot as plt
	from numpy import append, array, arange, vstack, linspace, mean, max, round
	from numpy import zeros, diff
	from hs3_utils import Flight_segment as F
	from hs3_utils import read_cod
	from matplotlib.cm import spectral
	from libtools import epoch2iso, shrink_ticks, mkdir_p
	from numpy.ma import masked_where
	from libgeo import p2z
	from scipy.stats import spearmanr as R
	from libnva import get_naapsaod_track, get_naapsconc_xsect
	import mpl_toolkits.axisartist as AA
	from mpl_toolkits.axes_grid1 import host_subplot
	from netCDF4 import Dataset
	from numpy.ma import masked_where
	from matplotlib.cm import bwr
	import ccplot_cmap as ccc

	if out_labels is None:
		dbg('mooooop')
		return

	#_initialize
	ncol = 1
	nrow = sum([plot_cpl, plot_aod, plot_cod, plot_naaps, plot_cmask]) 
	fig = plt.figure()

	row = 1 #_pretty easy to clean this up by using fig.axes[index] later
	if plot_aod:
		ax0 = fig.add_subplot(nrow, ncol, row)	#_shis
		row += 1
	if plot_cpl:
		ax1 = fig.add_subplot(nrow, ncol, row)	#_cpl
		row += 1
	if plot_naaps:
		ax2 = fig.add_subplot(nrow, ncol, row)	#_naaps_conc
		row += 1
	if plot_cmask:
		ax3 = fig.add_subplot(nrow, ncol, row, aspect=26)	#_cloud_mask
		row += 1
	if plot_cod:
		ax4 = fig.add_subplot(nrow, ncol, row)	#_cloud_mask
		row += 1	
	
	#_read in flight data
	flight = F(dtg=dtg, **kwargs) 
	en, yn, xn = flight.CPL_epoch, flight.CPL_latitude, flight.CPL_longitude 
	naaps = get_naapsaod_track(en, yn, xn, **kwargs)
	
	#_get start and end of area in cpl
	if fidx is None:
		f0 = 0
		f1 = -1 
		fidx = range(flight.size)
	else:
		f0 = min(fidx)
		f1 = max(fidx)

	#_set top header with dates
	st = '{0} - {1}'.format(epoch2iso(flight[fidx[f0]].CPL_epoch),
							epoch2iso(flight[fidx[f1]].CPL_epoch))
	fig.suptitle(st)
	
	letters = 'abcdefghijklmnopqrstuvwxyz'
	labs = []

	fname = flight.FLIGHT_file.replace('.nc', '.retrieved')

	color = ['red','blue','green','violet','orange']
	if plot_cpl:
		#_CPL nonsense
		top = 900; nz=6
		top = 350; nz=6
		top = 225; nz=6
		cpl_plt, cpl_x, cpl_y = flight.map_cpl(**kwargs) 
		cpl_max = cpl_plt.shape[1]
		cpl_nx, cpl_ny = cpl_plt.shape

		CTP = masked_where(flight.SHIS_CTP <= 0, flight.SHIS_CTP)
		CTZ = p2z(CTP)
		CTZ[CTZ <= 0] = -9999

		#_testing
		cmfile = '/home/wsessions/lib/cmaps/calipso-backscatter.cmap'
		cmap, norm, ticks = ccc.loadcolormap(cmfile, 'caliop')
	##	cb = ax1.pcolormesh(cpl_x[f0:f1], cpl_y[:top], cpl_plt.T[:top,f0:f1],
	##		vmin=0,vmax=1e-4, cmap=cmap, zorder=0, norm=norm)
		cb = ax1.pcolormesh(cpl_x[f0:f1], cpl_y[:top], cpl_plt.T[:top,f0:f1],
			vmin=0,vmax=1e-4, cmap=spectral, zorder=0)

		xlim = ax1.xaxis.get_data_interval()
		ax1.set_xlim((cpl_x[f0], cpl_x[f1])) #xlim)
		ax1.set_ylim(0, cpl_y[:top][-1])
		ax1.set_yticks(linspace(0, cpl_y[:top][-1], nz))
		ax1.set_yticklabels(['{0:4.1f}'.format(vv) for vv in 
			linspace(0, 6, nz)])
		##	linspace(0, 5, nz)]) #_jeff reviews 2016

		ax1.set_xticks(arange(cpl_x[f0], cpl_x[f1], len(cpl_x[f0:f1])/5))
		ax1.set_xticklabels([epoch2iso(ttt)[-8:] for
			ttt in flight.CPL_epoch[f0:f1:len(cpl_x[f0:f1])/5]])

		ax1.set_ylabel('Height (km)', size='xx-small')
		ax1.set_xlabel('532 nm Backscatter', size='xx-small')
		ax1.plot([f0, f0], [0, cpl_y[-1]], 'r-')
		ax1.plot([f1, f1], [0, cpl_y[-1]], 'y-')
				
	if plot_naaps:
		#_NAAPS nonsense
		xsect, naaps_z, valid_dt, valid_km = get_naapsconc_xsect(en, yn, xn,
			heights=True, valid_time=True, **kwargs)
		print cpl_x.shape, naaps_z.shape, xsect.shape
		cs = ax2.pcolormesh(cpl_x, naaps_z/1e3, xsect)
		ax2.set_xticks(arange(cpl_x[f0], cpl_x[f1], len(cpl_x[f0:f1])/5))
		ax2.set_xticklabels([epoch2iso(ttt)[-8:] for
			ttt in flight.CPL_epoch[f0:f1:len(cpl_x[f0:f1])/5]])
		ax2.set_xlim(cpl_x[f0], cpl_x[f1]) 

		dt_hour = valid_dt / 3600.
		dt_km	= valid_km # / 1e3

		ax2.get_yaxis().set_visible(False)

		#_setup distance
		if 0:
			#_split y axis on right side
			ax2_dt = ax2.twinx()
			ax2_dx = ax2.twinx()

			r0, = ax2_dt.plot(dt_hour, 'w-', linewidth=1)
			r1, = ax2_dx.plot(dt_km, 'k-', linewidth=1)
			ax2.set_ylabel('Distance to Grid (km)', size='xx-small')
			ax2_dt.set_ylabel('Time to Forecast (hr)', size='xx-small')
			ax2_dt.set_yticks(linspace(-3, 3, 7))
		
			ax2_dx.get_yaxis().set_visible(False)
			ax2_dx.set_ylim(0, 80)
			ax2.set_yticklabels([round(w) for w in linspace(0, 80, 9)])
		else:
			pass	

		naaps_dummy_height = 5
		ax2.set_yticks(linspace(0, naaps_dummy_height, 9))	#_faking this a bit	
		ax2.set_ylim(0, naaps_dummy_height)

	#_all cod should be the same? noo
	if plot_cmask: 
		if len(out_labels) > 1: raise RuntimeError, 'ONLY PLOTS ONE'
		cod = read_cod(flight, dtg, fidx=fidx, out_label=out_labels[0],**kwargs)

		#_cloud mask nonsense
		mask_cpl = cod.COD_CPL > cpl_thresh
		mask_shi = cod.COD_SHIS > shis_thresh
		mask_btd = find_BT_diff(flight, fidx=fidx, bt_thresh=bt_thresh,**kwargs)
		nt = len(mask_cpl)
		cpl = arange(nt)[mask_cpl]
		shi = arange(nt)[mask_shi]
		btd = arange(nt)[mask_btd]
		ncpl = len(cpl); nshi = len(shi); nbtd = len(btd) 

		#_label vlines
		xloc = len(fidx) / 30
		bbox = { 
			'boxstyle'	: 'round',
			'facecolor' : 'white', 
			'alpha'		: .75, 
			'edgecolor' : None
			}
		props = {
			'bbox'			: bbox, 
			'size'			: 'xx-small', 
			'verticalalignment' : 'center',
			}

		#_plot masks
		if True:
			cloud_mask = zeros((2, nt))
		##	cloud_mask[0] = mask_cpl
			cloud_mask[0] = mask_shi
			cloud_mask[1] = mask_btd
		##	ax3.text(xloc, 0.5,'CPL > {0:4.2f}'.format(cpl_thresh), **props)  
			ax3.text(xloc, 0.5,'SHIS > {0:4.2f}'.format(shis_thresh), **props) 
			ax3.text(xloc, 1.5,'${{\delta}}B_{{\\nu}}T$ > {0:3.1f}K'.format(
				bt_thresh), **props)  
		
		#_test of variaous BvT settings
		if False:
			#_how man to test
			dB = arange(3,9)
			ndB = dB.size
			cloud_mask = zeros((ndB, nt))
			tmp = arange(nt)
			
			for iii, dBvT in enumerate(dB):
				
				cloud_mask[iii] = find_BT_diff(flight, fidx=fidx, 
									bt_thresh=dBvT, **kwargs)
				ax3.text(xloc, 0.5 + iii,'{0:3.1f}'.format(dBvT), **props)  


		#_another one off
		if False:
			cloud_mask = zeros((3, nt))
			mask_all = append(mask_shi[:,None], mask_btd[:,None], 1).any(1)
			cloud_mask[0] = mask_shi
			cloud_mask[1] = mask_btd
			cloud_mask[0] = mask_shi
		#	cloud_mask[0] = mask_all
		#	cloud_mask[1] = mask_all
		#	cloud_mask[2] = mask_all
			ax3.text(xloc, 1.5,'Red = Masked'.format(shis_thresh), **props) 
			
		#_one off crap, ignore
		if False:
			cloud_mask = zeros((8, nt))
			cloud_mask[0] = cod.COD_SHIS > 0.25 
			cloud_mask[1] = cod.COD_SHIS > 0.50
			cloud_mask[2] = cod.COD_SHIS > 0.75
			cloud_mask[3] = cod.COD_SHIS > 0.75
			cloud_mask[4] = cod.COD_SHIS > 1.25
			cloud_mask[5] = cod.COD_SHIS > 1.50
			cloud_mask[6] = cod.COD_SHIS > 1.75
			cloud_mask[7] = cod.COD_SHIS > 2.00

		#_plot mask
		ax3.pcolor(cloud_mask, zorder=0, cmap=bwr)
		ax3.set_xlim(0, len(fidx))

		ax3.set_xticks(arange(cpl_x[f0], cpl_x[f1], len(cpl_x[f0:f1])/5))
		ax3.set_xticklabels([epoch2iso(ttt)[-8:] for
			ttt in flight.CPL_epoch[f0:f1:len(cpl_x[f0:f1])/5]])

		ax3.yaxis.set_visible(False)
		ax3.set_ylabel('Cloud Masks', size='xx-small')

	if plot_aod:
		#_DELETE
		jeff_dict = {}

		#_plot for each run type
		for n, out_label in enumerate(out_labels):

			#_generate input file name
			fname = flight.FLIGHT_file.replace('.nc',
				'.{0}.retrieved'.format(out_label))
	
			#_read in data
			data = read_real_retrieval(fname, **kwargs)

			#_give different colors based upon if cloud present, number oflayers
			x = data.tau
			y = arange(x.size)
	
			#_smooth out data
			if smooth:
				x = array([ mean(x[max([n-2,0]):n+3]) for n in range(x.size) ])

			#_do some neighbor screening
			if gradient_kibosh:
				#_only look at cases above 0.5
				idx_mag = x > 0.1 #25 #_ orig 0.5
				
				#_generate an array removing outliers
				idx_dif = abs(diff(x)) > 0.25 #_ orig 2.5
				idx_dif = [d0*d1 for d0, d1 in zip(idx_dif[:-1], idx_dif[1:])]
				idx_dif = append(append([False], idx_dif), [False])

				print 'P('
				print idx_dif*idx_mag
				print (idx_dif*idx_mag).sum(), idx_dif.size	
				print ')P'

				#_mask out these spikes
				x = masked_where(idx_dif*idx_mag, x)
						
			#_if producing cloud mask, use it
			if plot_cmask:
				mask_x = append(mask_shi[:,None], mask_btd[:,None], 1).any(1)
			##	mask_x = mask_shi
				x = masked_where(mask_x, x)
			##	x = x[mask_x == False]
			##	y = y[mask_x == False]
			##	x = x[mask_shi*mask_btd == False]
			##	y = y[mask_shi*mask_btd == False]

		##	ax0.scatter(y, x, color=color[n], label=out_label, zorder=0,
			ax0.scatter(y, x, color=color[n], label='UW-SSEC/S-HIS', zorder=3,
				marker='x', s=2)

			#_DELETE
			jeff_dict['UW'] = x

		##	ax0.plot(x, color=color[n], linestyle='-', linewidth='0.35',
		##		label=out_label, zorder=0, marker='x', markevery=50, ms=3)

		### CPL_stuff
		cpl_typ = vstack(data.CPL_type)	
		cpl_tau = vstack(data.CPL_tau)	

		#_trim down to number of layers possibly defined in type
		cpl_tau = cpl_tau[:,:5]
		cpl_typ = cpl_typ[:,::2]

		#_find where not missing
		idx_here = cpl_tau >= 0

		#_get column totals of elevated and pbl aersol when not missing
		col_pbl = ((cpl_typ == 1) * idx_here * (cpl_tau)).sum(1) #max(1) #sum(1)
		col_ele = ((cpl_typ == 2) * idx_here * (cpl_tau)).sum(1) #max(1) #sum(1)
		cpl_aod = col_ele + col_pbl
	
		y = cpl_aod
		if smooth:
			y = array([ mean(y[max([n-2,0]):n+3]) for n in range(y.size) ])

		#_DELETE
		jeff_dict['CPL'] = y
		jeff_dict['NAAPS'] = naaps
		jeff_dict['time'] = [epoch2iso(t2) for t2 in flight.CPL_epoch[:]]

		ax0.plot(y, 'k-', marker='x', linewidth=0.6, ms=3, label='CPL',
			markevery=50, zorder=1)
		ax0.plot(naaps, 'k-', marker='o', linewidth=0.6, label='NAAPS', ms=3,
			markevery=50, zorder=1)
		ax0.set_xticklabels(data.fov)
		ax0.grid(True)
	#	ax0.get_xaxis().set_visible(False)
		
		#_DELETE 2016.06.23 plot manual aod values
		ext_prof = vstack(flight.CPL_ext_532)
		import scipy.integrate
		ax0.plot(scipy.integrate.trapz(ext_prof)*30., linewidth=0.6, label='CPL?',
			zorder=4, color='magenta')


		ax0.set_xticks(arange(cpl_x[f0], cpl_x[f1], len(cpl_x[f0:f1])/5))
		ax0.set_xticklabels([epoch2iso(ttt)[-8:] for
			ttt in flight.CPL_epoch[f0:f1:len(cpl_x[f0:f1])/5]])

		ax0.set_xlim([0,x.size])
		ax0.set_ylim([0, 2.0])

		#_plot MODIS AOD if available
		if plot_modis:
	
			#_read in flight data
			fmod = flight.FLIGHT_file.replace('CPL.GDAS','MODIS')
			with Dataset(fmod, 'r') as cdf_mod:
	
				#_mask missing
				aod_mod = cdf_mod.variables['aod'][:]
				dis_mod = cdf_mod.variables['distance'][:]
				aod_mod = masked_where(aod_mod <= 0, aod_mod)
				aod_mod = masked_where(dis_mod > 3e4, aod_mod)

				x_mod = range(aod_mod.size)
				ax0.scatter(x_mod, aod_mod, s=4, marker='d', color='cyan',
					zorder=2, label='MODIS (T&A)')

				#_DELETE
				jeff_dict['MODIS'] = aod_mod

		props = dict(boxstyle='round', facecolor='white', alpha=.5)
		leg = ax0.legend(fontsize=8, loc=2, frameon=False)
		leg.get_frame().set_alpha(0.5)
				
		ax = ax0.twinx()
		ax.yaxis.set_visible(False)

		#_build label box
		props = dict(boxstyle='round', facecolor='white', alpha=.9)

		ax0.set_ylabel('AOD (532 nm)', size='xx-small')		

		#_DELETE dump file for jeff
		jn = fname.replace('.retrieved','aod_dump.csv').split('/')[-1]
		if False:
		  with open(jn, 'w') as f:
			order = ['time', 'UW', 'CPL', 'MODIS', 'NAAPS']
			f.write(','.join(order) + '\n')

			for iii, ttt in enumerate(jeff_dict['time']):
				tmp0 = [jeff_dict['time'][iii]]
				tmp1 = []
				for l in order[1:]:
					try:
						tmp1.append('{0:7.5f}'.format(jeff_dict[l][iii]))
					except ValueError:
						print jeff_dict[l][iii]
						tmp1.append('-9999')
				tmp0.extend(tmp1)# 1 = ['{0:7.5f}'.format(n) for n in tmp]
				f.write(','.join(tmp0) + '\n')
	

	if plot_cod:
		if len(out_labels) > 1: raise RuntimeError, 'ONLY PLOTS ONE'
		cod = read_cod(flight, dtg, fidx=fidx, out_label=out_labels[0], **kwargs)
		#_plot COD for shis and CPL
		ax4.plot(cod.COD_SHIS, color='red', linewidth=.4, zorder=2,label='SHIS')
		ax4.plot(cod.COD_CPL, color='black', linewidth=.4, zorder=3,label='CPL')
		ax4.set_ylim(0, 5)
		ax4.set_xlim(0, len(fidx))
		ax4.set_ylabel('COD', size='xx-small')

		#_add legend
		props = dict(boxstyle='round', facecolor='white', alpha=.5)
		ax4.legend(fontsize=6, loc=2)

	# swap with fig.axes loop
##	[shrink_ticks(axis) for axis in [ax0, ax1, ax2, ax3, ax2_dt, ax2_dx]]
	[shrink_ticks(axis) for axis in fig.axes]

	#_which
	www = ''.join([plot_aod*'.aod', plot_cpl*'.cpl', plot_cod*'.cod',
			plot_naaps*'.naaps', plot_cmask*'.cmask'])	

	#_name file and save
	dbg(("TESTING", out_labels, len(out_labels), out_labels[0]))
#	pname = fname.replace('.retrieved','jeffz2.{0}.png'.format(www)).replace(
#		out_label, '').split('/')[-1]
	pname = fname.replace('.retrieved','.jeffz2{0}.png'.format(www)).split('/')[-1]

	if gradient_kibosh:
		pname = pname.replace('.png', '.kibosh.png')		

	dir_out = dir_plot
	mkdir_p(dir_out)
	pname = os.path.join(dir_out, pname)
	plt.savefig(pname)
	plt.close()
	dbg(pname)


def plot_jeffs_qsub(dtg=None, fidx=None, out_label='', dz_cpl=25, smooth=True,
	dir_plot='.', out_labels=None, thresh=0.5, thresh_bt=5, **kwargs):
	'''
	fname	str,	path to output file of optimal estimation
		
	comparison plot between CPL AOD retrieval and SHIS AOD retrieval
	'''
	import matplotlib.pyplot as plt
	from numpy import append, array, arange, vstack, linspace, mean, max
	from hs3_utils import Flight_segment as F
	from hs3_utils import read_cod
	from matplotlib.cm import spectral
	from libtools import epoch2iso, shrink_ticks, mkdir_p
	from numpy.ma import masked_where
	from libgeo import p2z
	from scipy.stats import spearmanr as R
	from hs3_utils import read_cod

	if out_labels is None:
		dbg('mooooop')
		return

	#_initialize
	ncol = 1
	nrow = len(out_labels) + 1 #_cpl 
	fig = plt.figure()
	ax1 = fig.add_subplot(nrow, 1, nrow)	#_cpl

	flight = F(dtg=dtg, **kwargs) 
	#_get start and end of area in cpl
	if fidx is None:
		f0 = 0
		f1 = -1 #9999
	else:
		f0 = min(fidx)
		f1 = max(fidx)
	
	letters = 'abcdefghijklmnopqrstuvwxyz'
	labs = []
		
	#_get mask for BTD
	mask_btd = find_BT_diff(flight, fidx=fidx, **kwargs)

	for n, out_label in enumerate(out_labels):
	
		cod = read_cod(flight, dtg, fidx=fidx, out_label=out_label, **kwargs)
		mask_shi = cod.COD_SHIS > 1.25
		mask_idx = append(mask_shi[:,None], mask_btd[:,None], 1).any(1)

		#_generate input file name
		fname = flight.FLIGHT_file.replace('.nc',
			'.{0}.retrieved'.format(out_label))
	
		#_read in data
		data = read_real_retrieval(fname, **kwargs)
	
		cpl_typ = vstack(data.CPL_type)	
		cpl_tau = vstack(data.CPL_tau)	

		#_trim down to number of layers possibly defined in type
		cpl_tau = cpl_tau[:,:5]
		cpl_typ = cpl_typ[:,::2]

		#_find fovs with any aerosol
		idx_pbl = (cpl_typ == 1).any(1)	#_index of PBL type
		idx_ele = (cpl_typ == 2).any(1) #_index of AEROSOL type
		idx_aer = append(idx_pbl[:,None], idx_ele[:,None], axis=1).any(1)
		idx_cld = (cpl_typ == 3).any(1) * idx_aer	#_aerosol w/ cloud
		idx_clr = (cpl_typ != 3).all(1) * idx_aer	#_aerosol w/o cloud

		#_find where not missing
		idx_here = cpl_tau >= 0

		#_get column totals of elevated and pbl aersol when not missing
		col_pbl = ((cpl_typ == 1) * idx_here * (cpl_tau)).sum(1) #max(1) #sum(1)
		col_ele = ((cpl_typ == 2) * idx_here * (cpl_tau)).sum(1) #max(1) #sum(1)
		cpl_aod = col_ele + col_pbl
		ax0 = fig.add_subplot(nrow, 1, n+1)	#_line

		#_plot MODIS AOD if available
		if True:  
			from netCDF4 import Dataset
	
			#_read in flight data
			fmod = flight.FLIGHT_file.replace('CPL.GDAS','MODIS')
			with Dataset(fmod, 'r') as cdf_mod:
	
				#_mask missing
				aod_mod = cdf_mod.variables['aod'][:]
				dis_mod = cdf_mod.variables['distance'][:]
				aod_mod = masked_where(aod_mod <= 0, aod_mod)
				aod_mod = masked_where(dis_mod > 3e4, aod_mod)

				x_mod = range(aod_mod.size)
				ax0.scatter(x_mod, aod_mod, s=4, marker='d', color='cyan',
					zorder=1, label='MODIS (T&A)')

		#_give different colors based upon if cloud present, number of layers
		x = data.tau
		y = cpl_aod
	
		#_smooth out data
		if smooth:
			y = array([ mean(y[max([n-2,0]):n+3]) for n in range(y.size) ])
			x = array([ mean(x[max([n-2,0]):n+3]) for n in range(x.size) ])

		print "TESTING", x.size, mask_idx.size
		ax0.plot(y, 'k-', label='CPL', zorder=0, linewidth=0.6)
		xx = arange(x.size)[mask_idx == False]
		xy = x[mask_idx == False]
		ax0.scatter(xx, xy, color='red', marker='o', s=2, label='S-HIS',
			zorder=2)

		if n == 0:
			leg = ax0.legend(fancybox=True, loc='upper left', fontsize=6)
			leg.get_frame().set_alpha(0.9)
		ax0.set_xticklabels(data.fov)
		ax0.grid(True)
		ax0.get_xaxis().set_visible(False)
		ax0.set_xlim([0,x.size])
		ax0.set_ylim([0, 2.0])

		ax = ax0.twinx()
	#	labs.append('{0}={1}'.format(letters[n], '_'.join(out_label.split('_')[2:])))
		lab='{0}'.format('_'.join(out_label.split('_')[2:]))
	##	lab='r={0:4.2f}, {1}'.format(R(x,y)[0], '_'.join(out_label.split('_')[2:]))
		ax.set_ylabel(letters[n], size='xx-small', rotation=270)
		ax.yaxis.set_visible(False)

		#_build label box
		props = dict(boxstyle=None, facecolor='white', alpha=.9)
		ax0.text(x.size*0.2, 1.75, lab, color='k', verticalalignment='top',
	##		size='xx-small')
			size='xx-small', bbox=props)

		shrink_ticks(ax0)
		
	fig.text(0.04, 0.6, 'Optical Thickness (532 nm)', va='center', size='x-small',
		rotation='vertical')

	#_CPL nonsense
	top = 225; nz=6
	cpl_plt, cpl_x, cpl_y = flight.map_cpl(**kwargs) 
	cpl_max = cpl_plt.shape[1]
	cpl_nx, cpl_ny = cpl_plt.shape
	CTP = masked_where(flight.SHIS_CTP <= 0, flight.SHIS_CTP)
	CTZ = p2z(CTP)
	CTZ[CTZ <= 0] = -9999
	cb = ax1.pcolormesh(cpl_x[f0:f1], cpl_y[:], cpl_plt.T[:,f0:f1],
		vmin=0,vmax=1e-4, cmap=spectral, zorder=0)
#	cb = ax1.pcolormesh(cpl_x[f0:f1], cpl_y[:top], cpl_plt.T[:top,f0:f1],
#		vmin=0,vmax=1e-4, cmap=spectral, zorder=0)
	xlim = ax1.xaxis.get_data_interval()
	ax1.set_xlim((cpl_x[f0], cpl_x[f1])) #xlim)
	ax1.set_ylim(0, cpl_y[:][-1])
	ax1.set_yticks(linspace(0, cpl_y[:][-1], nz))
#	ax1.set_ylim(0, cpl_y[:top][-1])
#	ax1.set_yticks(linspace(0, cpl_y[:top][-1], nz))
	ax1.set_yticklabels(['{0:4.1f}'.format(vv) for vv in 
		linspace(0, 18, nz)])
	ax1.set_xticks(arange(cpl_x[f0], cpl_x[f1], len(cpl_x[f0:f1])/5))
	ax1.set_xticklabels([epoch2iso(ttt)[-8:] for
		ttt in flight.CPL_epoch[f0:f1:len(cpl_x[f0:f1])/5]])
	shrink_ticks(ax1)
	ax1.set_ylabel('Height (km)', size='xx-small')
	ax1.set_xlabel('532 nm Backscatter', size='xx-small')
	ax1.plot([f0, f0], [0, cpl_y[-1]], 'r-')
	ax1.plot([f1, f1], [0, cpl_y[-1]], 'y-')
			
	fig.suptitle('Retrieval Comparison', size='xx-small')
#	fig.suptitle(','.join(labs), size='xx-small')

	#_move everything up for legend
##	leg = ax.

	#_plot mask data 
	cod = read_cod(flight, dtg, fidx=fidx, out_label=out_label, **kwargs)
	mask_cpl = cod.COD_CPL > thresh
	mask_shi = cod.COD_SHIS > thresh
	mask_btd = find_BT_diff(flight, fidx=fidx, thresh_bt=thresh_bt, **kwargs)
	nt = len(mask_cpl)
	cpl = arange(nt)[mask_cpl]
	shi = arange(nt)[mask_shi]
	btd = arange(nt)[mask_btd]
	ncpl = len(cpl); nshi = len(shi); nbtd = len(btd) 

	#_name file and save
	pname = fname.replace('.retrieved','jeffz.png').replace(out_label,
			 '').split('/')[-1]
	dir_out = dir_plot
	mkdir_p(dir_out)
	pname = os.path.join(dir_out, pname)
	plt.savefig(pname)
	plt.close()
	dbg(pname)


def find_BT_diff(flight, fidx=None, bt_thresh=5, surf_temp=-1, **kwargs):
	'''
	flag where in fidx there is a difference of thresh between the
	surface temperature and the average observed radiance
	'''
	from lblrtm_utils import get_surf_temp 
	from libgeo import planck_inv
	from lblrtm_utils import microwindow_average, microwindow_stdev
	from numpy import array

	if fidx is None:
		fidx = range(len(flight))

	#_get wavenumbers for observations
	wave = flight.SHIS_wavenumber

	bt_diff = [None] * len(fidx)
	for i in fidx:

		#_pull out desired field of view
		fov = flight[i]	
		surf_temp = get_surf_temp(fov, surf_temp_src='GDAS', **kwargs)
	
		#_get average brightness temperature for microwindows
		r, w = microwindow_average(fov.SHIS_radiances, wave, -26)
		bt = planck_inv(r, w*100, domain='wavenumber')[-3:].mean()	
	##	dbg(('USED FOR dBvT', w[-3:]))
	
		#_set flag based on thresh
		bt_diff[i] = True if surf_temp - bt > bt_thresh else False

	return array(bt_diff)


def plot_jeff_qsub(dtg=None, fidx=None, out_label='', dz_cpl=25, smooth=True,
	dir_plot='.', **kwargs):
	'''
	fname	str,	path to output file of optimal estimation
		
	comparison plot between CPL AOD retrieval and SHIS AOD retrieval
	'''
	import matplotlib.pyplot as plt
	from numpy import append, array, arange, vstack, linspace, mean, max
	from hs3_utils import Flight_segment as F
	from matplotlib.cm import spectral
	from libtools import epoch2iso, shrink_ticks, mkdir_p
	from numpy.ma import masked_where
	from libgeo import p2z

	#_generate input file name
	flight = F(dtg=dtg, **kwargs) 
	fname = flight.FLIGHT_file.replace('.nc',
		'.{0}.retrieved'.format(out_label))

	#_get start and end of area in cpl
	if fidx is None:
		f0 = 0
		f1 = -1
	else:
		f0 = min(fidx)
		f1 = max(fidx)
	
	#_read in data
	data = read_real_retrieval(fname, **kwargs)
	
	cpl_typ = vstack(data.CPL_type)	
	cpl_tau = vstack(data.CPL_tau)	

	#_trim down to number of layers possibly defined in type
	cpl_tau = cpl_tau[:,:5]
	cpl_typ = cpl_typ[:,::2]

	#_find fovs with any aerosol
	idx_pbl = (cpl_typ == 1).any(1)	#_index of PBL type
	idx_ele = (cpl_typ == 2).any(1) #_index of AEROSOL type
	idx_aer = append(idx_pbl[:,None], idx_ele[:,None], axis=1).any(1) #_all aero
	idx_cld = (cpl_typ == 3).any(1) * idx_aer	#_aerosol w/ cloud
	idx_clr = (cpl_typ != 3).all(1) * idx_aer	#_aerosol w/o cloud

	#_find where not missing
	idx_here = cpl_tau >= 0

	#_get column totals of elevated and pbl aersol when not missing
	col_pbl = ((cpl_typ == 1) * idx_here * (cpl_tau)).sum(1) #max(1) #sum(1)
	col_ele = ((cpl_typ == 2) * idx_here * (cpl_tau)).sum(1) #max(1) #sum(1)
	cpl_aod = col_ele + col_pbl
##	cpl_aod = append(col_pbl[:,None], col_ele[:,None], axis=1).max(1)

	#_initialize
	fig = plt.figure()
	ax0 = fig.add_subplot(211)	#_line
	ax1 = fig.add_subplot(212)	#_cpl
	
	#_give different colors based upon if cloud present, number of layers
	x = data.tau
	y = cpl_aod
	
	#_smooth out data
	if smooth:
		y = array([ mean(y[max([n-2,0]):n+3]) for n in range(y.size) ])
		x = array([ mean(x[max([n-2,0]):n+3]) for n in range(x.size) ])

	ax0.plot(y, 'r-', label='CPL')
	ax0.plot(x, 'k-', label='S-HIS')
	leg = ax0.legend(fancybox=True, loc='upper left')
	leg.get_frame().set_alpha(0.5)
	ax0.set_xticklabels(data.fov)
	ax0.grid(True)
	ax0.get_xaxis().set_visible(False)
	ax0.set_ylabel('Optical Thickness (532 nm)', size='xx-small')
	ax0.set_xlim([0,x.size])
	ax0.set_ylim([0, 5.0])
	shrink_ticks(ax0)

	#_CPL nonsense
	cpl_plt, cpl_x, cpl_y = flight.map_cpl(**kwargs) 
	cpl_max = cpl_plt.shape[1]
	cpl_nx, cpl_ny = cpl_plt.shape
	CTP = masked_where(flight.SHIS_CTP <= 0, flight.SHIS_CTP)
	CTZ = p2z(CTP)
	CTZ[CTZ <= 0] = -9999
	cb = ax1.pcolormesh(cpl_x[f0:f1], cpl_y[:450], cpl_plt.T[:450,f0:f1], vmin=0,vmax=1e-4,
			cmap=spectral, zorder=0)
	ax1.scatter(cpl_x, CTZ, marker='x', linewidth=0.5, s=4, color='yellow',
			zorder=1)
	xlim = ax1.xaxis.get_data_interval()
	ax1.set_xlim((cpl_x[f0], cpl_x[f1])) #xlim)
	ax1.set_ylim(0, cpl_y[:450][-1])
	ax1.set_yticks(linspace(0, cpl_y[:450][-1], 10))
	ax1.set_yticklabels(['{0:4.1f}'.format(vv) for vv in 
		linspace(0, 9, 10)])
	ax1.set_xticks(arange(cpl_x[f0], cpl_x[f1], len(cpl_x[f0:f1])/5))
	ax1.set_xticklabels([epoch2iso(ttt)[-8:] for
		ttt in flight.CPL_epoch[f0:f1:len(cpl_x[f0:f1])/5]])
	shrink_ticks(ax1)
	ax1.set_ylabel('Height (km)', size='xx-small')
	ax1.set_xlabel('532 nm Backscatter', size='xx-small')
	ax1.plot([f0, f0], [0, cpl_y[-1]], 'r-')
	ax1.plot([f1, f1], [0, cpl_y[-1]], 'y-')
			
	fig.suptitle('Retrieval Comparison', size='xx-small')

	#_name file and save
	pname = fname.replace('.retrieved','.jeff.png').split('/')[-1]
	dir_out = os.path.join(dir_plot, out_label)
	mkdir_p(dir_out)
	pname = os.path.join(dir_out, pname)
	plt.savefig(pname)
	plt.close()
	dbg(pname)



def plot_real_qsub(dtg=None, fidx=None, out_label='', dz_cpl=25, smooth=True,
	dir_plot='.', **kwargs):
	'''
	fname	str,	path to output file of optimal estimation
		
	comparison plot between CPL AOD retrieval and SHIS AOD retrieval
	'''
	import matplotlib.pyplot as plt
	from numpy import append, array, arange, vstack, linspace, mean, max
	from hs3_utils import Flight_segment as F
	from matplotlib.cm import spectral
	from libtools import epoch2iso, shrink_ticks
	from numpy.ma import masked_where
	from libgeo import p2z
	from libtools import mkdir_p

	#_generate input file name
	flight = F(dtg=dtg, **kwargs) 
	fname = flight.FLIGHT_file.replace('.nc','.{0}.retrieved'.format(out_label))

	#_get start and end of area in cpl
	if fidx is None:
		f0 = 0
		f1 = 9999
	else:
		f0 = min(fidx)
		f1 = max(fidx)
	
	#_read in data
	data = read_real_retrieval(fname, **kwargs)
	
	cpl_typ = vstack(data.CPL_type)	
	cpl_tau = vstack(data.CPL_tau)	

	#_trim down to number of layers possibly defined in type
	cpl_tau = cpl_tau[:,:5]
	cpl_typ = cpl_typ[:,::2]

	#_find fovs with any aerosol
	idx_pbl = (cpl_typ == 1).any(1)	#_index of PBL type
	idx_ele = (cpl_typ == 2).any(1) #_index of AEROSOL type
	idx_aer = append(idx_pbl[:,None], idx_ele[:,None], axis=1).any(1) #_all aero
	idx_cld = (cpl_typ == 3).any(1) * idx_aer	#_aerosol w/ cloud
	idx_clr = (cpl_typ != 3).all(1) * idx_aer	#_aerosol w/o cloud

	#_find where not missing
	idx_here = cpl_tau >= 0

	#_get column totals of elevated and pbl aersol when not missing
	col_pbl = ((cpl_typ == 1) * idx_here * (cpl_tau)).sum(1) #max(1) #sum(1)
	col_ele = ((cpl_typ == 2) * idx_here * (cpl_tau)).sum(1) #max(1) #sum(1)
	cpl_aod = col_ele + col_pbl
##	cpl_aod = append(col_pbl[:,None], col_ele[:,None], axis=1).max(1)

	#_initialize
	fig = plt.figure()
	axe = fig.add_subplot(221)	#_scatter
	ax0 = fig.add_subplot(222)	#_line
	ax1 = fig.add_subplot(223)	#_cpl
	
	#_give different colors based upon if cloud present, number of layers
	x = data.tau
	y = cpl_aod

	#_plot cloudy aerosol cases
	axe.scatter(x[idx_cld], y[idx_cld], marker='o', color='r', label='cloud', s=2)
	
	#_plot clear aerosol cases
	axe.scatter(x[idx_clr], y[idx_clr], marker='x', color='b', label='clear', s=2)
	axe.legend()

	#_label
	axe.set_xlabel('S-HIS retrieval (red=start, yellow=stop)', size='xx-small')
	axe.set_ylabel('CPL retrieval', size='xx-small')
	axe.grid(True)

	#_set title
	axe.set_title('{0}'.format(smooth*'SMOOTHED'))	

	#_smooth out data
	if smooth:
		y = array([ mean(y[max([n-2,0]):n+3]) for n in range(y.size) ])
		x = array([ mean(x[max([n-2,0]):n+3]) for n in range(x.size) ])

	ax0.plot(y, 'r-', label='cpl')
	ax0.plot(x, 'k-', label='shis')
	ax0.legend()
	ax0.set_xticklabels(data.fov)
#	ax1.plot(x, 'k-')
#	ax1.plot(x[idx_cld], 'b-')
#	ax1.plot(x[idx_clr], 'r-')

	axe.set_xlim(0, 1) #4)
	axe.set_ylim(0, 1) #4)

	#_CPL nonsense
	cpl_plt, cpl_x, cpl_y = flight.map_cpl(**kwargs) 
	cpl_max = cpl_plt.shape[1]
	cpl_nx, cpl_ny = cpl_plt.shape
	CTP = masked_where(flight.SHIS_CTP <= 0, flight.SHIS_CTP)
	CTZ = p2z(CTP)
	CTZ[CTZ <= 0] = -9999
	cb = ax1.pcolormesh(cpl_x, cpl_y, cpl_plt.T, vmin=0, vmax=1e-4,
			cmap=spectral, zorder=0)
	ax1.scatter(cpl_x, CTZ, marker='x', linewidth=0.5, s=4, color='yellow',
			zorder=1)
	xlim = ax1.xaxis.get_data_interval()
	ax1.set_xlim(xlim)
	ax1.set_ylim(0, cpl_y[-1])
	ax1.set_yticks(linspace(0, cpl_y[-1], 11))
	ax1.set_yticklabels(['{0:4.1f}'.format(vv) for vv in
		linspace(0, cpl_y[-1], 11)*dz_cpl/1e3])
	ax1.set_xticks(arange(0, cpl_nx, cpl_nx/5))
	ax1.set_xticklabels([epoch2iso(ttt)[-8:] for
		ttt in flight.CPL_epoch[::cpl_nx/5]])
	shrink_ticks(ax1)
	ax1.set_ylabel('just assume 0-18km', size='xx-small')
##	ax1.set_ylabel('HEIGHTS TOO LOW ~1KM', size='xx-small')
	ax1.plot([f0, f0], [0, cpl_y[-1]], 'r-')
	ax1.plot([f1, f1], [0, cpl_y[-1]], 'y-')
##	ax1.plot([f0, f0], [0, cpl_ny], 'r-')
##	ax1.plot([f1, f1], [0, cpl_ny], 'y-')

	#_name file and save
	pname = fname.replace('.retrieved','.png').split('/')[-1]
	dir_out = os.path.join(dir_plot, out_label)
	mkdir_p(dir_out)
	pname = os.path.join(dir_out, pname)
	plt.savefig(pname)
	plt.close()
	dbg(pname)


def plot_real(fname=None, dtg=None, out_label='', **kwargs):
	'''
	fname	str,	path to output file of optimal estimation
		
	comparison plot between CPL AOD retrieval and SHIS AOD retrieval
	'''
	import matplotlib.pyplot as plt
	from numpy import append, vstack
	from hs3_utils import Flight_segment as F
	
	if fname is None:
		flight = F(dtg=dtg) 
		fname = flight.FLIGHT_file.replace('.nc',
				'.{0}.retrieved'.format(out_label))

	#_read in data
	data = read_real_retrieval(fname, **kwargs)

	cpl_typ = vstack(data.CPL_type)	
	cpl_tau = vstack(data.CPL_tau)	

	#_trim down to number of layers possibly defined in type
	cpl_tau = cpl_tau[:,:5]
	cpl_typ = cpl_typ[:,::2]

	#_find fovs with any aerosol
	idx_pbl = (cpl_typ == 1).any(1)
	idx_ele = (cpl_typ == 2).any(1)
	idx_aer = append(idx_pbl[:,None], idx_ele[:,None], axis=1).any(1)
	idx_cld = (cpl_typ == 3).any(1) * idx_aer
	idx_clr = (cpl_typ != 3).all(1) * idx_aer

	#_find where not missing
	idx_here = cpl_tau >= 0

	#_get column totals of elevated and pbl aersol when not missing
	col_pbl = ((cpl_typ == 1) * idx_here * (cpl_tau)).max(1) #sum(1)
	col_ele = ((cpl_typ == 2) * idx_here * (cpl_tau)).max(1) #sum(1)
##	cpl_aod = col_ele + col_pbl
	cpl_aod = append(col_pbl[:,None], col_ele[:,None], axis=1).max(1)

	#_initialize
	fig = plt.figure()
	axe = fig.add_subplot(221)
	ax0 = fig.add_subplot(222)
	ax1 = fig.add_subplot(223)
	
	#_give different colors based upon if cloud present, number of layers
	x = data.tau
	y = cpl_aod

	#_plot cloudy aerosol cases
	axe.scatter(x[idx_cld], y[idx_cld], marker='o', color='b')
	
	#_plot clear aerosol cases
	axe.scatter(x[idx_clr], y[idx_clr], marker='x', color='r')

	#_label
	axe.set_xlabel('S-HIS retrieval', size='xx-small')
	axe.set_ylabel('CPL retrieval', size='xx-small')
	axe.grid(True)

	ax0.plot(y, 'r-')
	ax1.plot(x, 'k-')
#	ax1.plot(x[idx_cld], 'b-')
#	ax1.plot(x[idx_clr], 'r-')

	axe.set_xlim(0, 4)
	axe.set_ylim(0, 4)

	#_name file and save
	pname = fname.replace('.retrieved','.png')
	plt.savefig(pname)
	dbg(pname)


def plot_close(pname=None, fig=None, **kwargs):
	''' call to plot figure in case not plotted elsewhere '''
	import matplotlib.pyplot as plt
	import os, time
	if not os.path.exists(pname):
		dbg(('SAVING', pname))
		fig.tight_layout()
		plt.savefig(pname)
	time.sleep(5)


def plot_DELETE(r_obs, n_obs, r_lbl, n_lbl, fig=None, clddef={}, surf_temp=0, dv=0,
	pname='default.png', label='', last=False, r_clr=None, n_clr=None, sup=None,
	ax=None, **kwargs):
	'''
	create plot of r values in brightness temperature space
	To have multiple on a single plot, pass axes through ax
	'''
	import os
	import matplotlib.pyplot as plt
	from libgeo import planck_inv
	from numpy import array, tile
	from libtools import shrink_ticks
	import os

	#_check if clear sky passed
	clear = False if r_clr is None else True
	first = True if len(ax.lines) == 0 else False

	#_kill matrices
	r_obs = r_obs if r_obs.ndim == 1 else array(r_obs).squeeze()
	r_lbl = r_lbl if r_lbl.ndim == 1 else array(r_lbl).squeeze()

	try:
		tau = clddef[0]['tau'][0]
	except:
		tau = clddef[0]['tau']

	ref = clddef[0]['ref']
	zkm = clddef[0]['z']
	ztp = clddef[0]['z_top']
	thk = ztp - zkm 
	tit = label + ' '
	tit += '{0:4.2f} tau, {1:4.2f} ref, {2:4.2f} z, {3:4.2f} thk '.format(tau,
															 ref, zkm, thk)
	tit += '(red = obs, green = lbl)'

	#_check if in brightness temperature space
	t		= tile(surf_temp, n_obs.size)
	t_obs	= planck_inv(r_obs/1e3, n_obs*100, domain='wavenumber')	
	t_lbl	= planck_inv(r_lbl/1e3, n_lbl*100, domain='wavenumber')	

	#_line settings
	arg = { 'linewidth' : 0.3 }
	if dv < 0:
		arg.update({'marker' : 'x' })

	#_clear and obs should be plotted first, then skipped for alpha
	if first: 
		#_plot obs in red
		ax.plot(n_obs, t_obs, color='r', **arg)

		#_if clear sky radiances are passed
		if clear:
			#_if last, setup axes and the like
			print 'CLEARALPHATEST', last, clear
			#_output from lbldis is many flavored if tau vs taus given
			if r_clr.ndim != 1:
				r_clr = array(r_clr).squeeze()

			#_plot clear sky radiances in blue
			t_clr = planck_inv(r_clr/1e3, n_clr*100, domain='wavenumber')
			ax.plot(n_clr, t_clr, color='#CCCCFF', **arg)

	#_plot this iteration's radiances
	ax.plot(n_lbl, t_lbl, color='g', **arg)

	#_if last, setup axes and the like
	print 'ALPHATEST', last, clear
	if last:

		#_set supertitle if not set
		if fig is not None and sup is not None:
			fig.suptitle(sup, size='xx-small')

		#_loop over lines, skip clear sky and obs, set alphas
		iskip = 2 if clear else 1
		alphas = linspace(0.2, 1, len(ax.lines)-iskip)
		print 'PREALPHAS', alphas[i], iskip
		for i, line in enumerate(ax.lines):
			print 'ALPHAS', alphas[i], iskip
			if i < iskip:
				continue
			#_fade old lines
			plt.setp(line, 'alpha', alphas[i])

		#_shrink axis tick marks
		shrink_ticks(ax) 
		ax.grid(True)
		ax.set_xlim(n_lbl[0]-20, n_lbl[-1]+20)
		ax.set_title(tit, size='xx-small')

		#_do something where the list of values tried is somewhere


################################################################################
#_i/o_##########################################################################
################################################################################


def read_real_retrieval(fname, **kwargs):
	'''
	Read in output file produced by optimal estimation test
	on real SHIS data. File should be comma deliminated and 
	
	FOV#, TAU, (other state vars), 10*tau_532, 10*layertype 

	TYPES,
	0	Missing	
	1	PBL (aerosol)
	2	Elevated (aerosol)
	3	Cloud
	4	Indeterminate	

	fname	str,	full path to file

	'''
	from numpy import array, recarray, ndarray
	
	dbg(('opening', fname))
	with open(fname, 'r') as f:
		lines = [l.strip() for l in f.readlines()]

		#_pull out state variables
		state_vars = [l.strip() for l in lines[0].split(',')[1:]]
		nvar = len(state_vars)
		del lines[0]

		#_initialize output recarray
		nfov = len(lines)
		dtype = [('fov','i4'), ('CPL_type',ndarray),('CPL_tau',ndarray)]
		[dtype.append((sv, 'f4')) for sv in state_vars]
		data = recarray((nfov,), dtype)

		#_loop over lines and fill arrays
		for k, line in enumerate(lines):
			cols = line.split(',')
			fov = cols[0]
			cols = cols[1:]
		
			#_add field of view number
			data[k].fov = fov

			#_loop over state variables
			for i, variable in enumerate(state_vars):
				setattr(data[k], variable, cols[i])

			#_add tau and type
			data[k].CPL_tau = array([float(c) for c in cols[nvar:nvar+10]]) 
			data[k].CPL_type = array([float(c) for c in cols[nvar+10:]])
	
	return data 


################################################################################
#_tools_########################################################################
################################################################################


def dbg(msg, l=1, err=False):
	''' 
	if global debug is set to true, be more verbose 
	msg : str, Message to be printed
	l   : int, Debug level of message.  Set higher for lower level 
	        messages.  As debug increases, noisiness should also.
	'''
	import inspect
	import re

	if hasattr(msg, '__iter__'):
		msg = ' '.join([str(m) for m in msg])
	elif type(msg) != str:
		msg = str(msg)

	if DEBUG >= l:
		curf    = inspect.currentframe()
		calf    = inspect.getouterframes(curf,2)
		file, line, method = calf[1][1:4]
		file     = '.'.join(file.split('/')[-1].split('.')[:-1])
		scream  = '[%s.%s.%i] %s' % (file,method,line,msg) \
			if not re.search('^\s*$', msg) else ' ' 

		if 0: #len(scream) > 80: shit makes it hard to cp
			ii = 80
			ss = ''
			nl = 70
			sx = '\n        + '
			ss = scream[:ii] + sx 
			for i in range(len(scream) / nl -1 ):
				ss += scream[ii:ii+nl] + sx 
				ii += nl 
			else:
				ss += scream[ii:]
			print ss	
			scream = ss
			exit()

		if not err:
			print scream
		else:
			raise RuntimeError, scream


def gen_outlabel(dv=None, apriori=None, dynamic_CTP=None, out_label=None,
	dynamic_ATP=None, cld=False, bias_correct=False, state_vars=None,
	persistence=False, **kwargs):
	'''
	generate an out label based on settings to prevent forgetting to mod later
	out_label	str,	The arbitrary part that this overwrites
	'''
	n_habit = []
	zkm = [] 
	for i, layer in enumerate(apriori):

		#_check if there is a cloud
		if layer['type'] == 'cloud':
			cld = True

		else:
			#_increment habit type
			n_habit.append(layer['ssp_db'])
			zkm.append(layer['z_top'])

	n_habit = len(set(n_habit))
	zkm = '-'.join(['{0}KM'.format(int(z)) for z in set(zkm)])

	#_is there a cloud layer?
	cld = '{0}CLD'.format('NO'*(not cld))

	#_what are we using for cloud top pressure
	try:
		if cld:
			ctp = 'CTP-{0}'.format('-'.join(dynamic_CTP))
		else:
			ctp = 'NOCTP'
	except:
		ctp = 'NOCTP' 
	try:
		atp = 'ATP-{0}'.format('-'.join(dynamic_ATP))
	except:
		atp = 'NOATP' 

	#_bias correction?
	b = bias_correct * 'b'
	
	#_use persistence?
	p = 'p' * (type(persistence) == dict)

	#_state vars being retrieved
	sv = '-'.join(state_vars)

	#_update out_label
	arg = (out_label, n_habit, abs(dv), cld, ctp, zkm, b, sv, atp, p)
	return '{0}_{7}_M{1}_DV{2}{6}{9}_{3}_{4}_{8}_{5}'.format(*arg)


def check_missing(dtg, fidx=None, dir_lblrtm_fmt=None, lbldis_output='',
	lbldis_input='', lbldis_clear='', dir_plot='.', **kwargs):
	from hs3_utils import find_shiscplgdas as find_scg
	import re
	from hs3_utils import Flight_segment as F

	#_get 
	fname = find_scg(dtg, **kwargs)
	redtg = re.search('SHIS.CPL.GDAS.COLLOC.(\d+).(\d+).nc', fname)

	#_pull out start and end dtg
	dtg0 = '20{0}'.format(redtg.group(1))
	dtg1 = '20{0}'.format(redtg.group(2))

	total = 0

	if fidx is None:
		flight = F(file_seg=fname, **kwargs)		
		fidx = range(len(flight))

	#_loop over fidx to get
	missing_count = 0
	missing_count_in = 0
	missing_count_clr = 0
	missing_count_pn = 0
	for i in fidx:
		#_build file path
		dir_lblrtm = dir_lblrtm_fmt.format(dtg0, dtg1, i, experiment, dtg0[2:])
		lbldis = os.path.join(dir_lblrtm, lbldis_output)
		lbldis_in = os.path.join(dir_lblrtm, lbldis_input)
		lbldis_clr = os.path.join(dir_lblrtm, lbldis_clear)
	#	pn = os.path.join(dir_plot, pname.format(dtg0, i))

		#_check if missing
		if not os.path.exists(lbldis):
			missing_count += 1	
		if not os.path.exists(lbldis_in):
			missing_count_in += 1
		if not os.path.exists(lbldis_clr):
			missing_count_clr += 1
	#	if not os.path.exists(pn):
	#		missing_count_pn += 1
		#_wait what	
		total += 1

	#_report
	dbg('{0} missing out {1:>6d} / {2:>6d}'.format(dtg, missing_count, len(fidx)))
	dbg('{0} missing in  {1:>6d} / {2:>6d}'.format(dtg, missing_count_in, len(fidx)))
	dbg('{0} missing clr {1:>6d} / {2:>6d}'.format(dtg, missing_count_clr, len(fidx)))
#	dbg('{0} missing pn  {1:>6d} / {2:>6d}'.format(dtg, missing_count_pn, len(fidx)))
	return missing_count, len(fidx)


def plot_sbatch(fidx=[], **kwargs):
	''' launch plotting of sbatch crap '''
	from libtools import dummy_sbatch
	from pickle import dump
	from os import unlink

	if len(fidx) == 0:

		from hs3_utils import Flight_segment	
		fidx = xrange(Flight_segment(**kwargs).size)
	
	dtg = kwargs.get('dtg')
	for i in fidx:

		kw_fn = '/home/wsessions/kwargs/plot_sbatch_{0}_{1:05d}'.format(dtg, i)	
		sb_fn = dummy_sbatch('/home/wsessions/lib/plot_sbatch.py {0}'.format(kw_fn),
				unique_id='{0:05d}'.format(i))	

		kwargs.update({'fidx' : [i]})
		dump(kwargs, open(kw_fn, 'wb'))

		cmd = 'sbatch {0}'.format(sb_fn)
		os.system(cmd)
		unlink(sb_fn)


################################################################################
#_SHELL_########################################################################
################################################################################


if __name__ == '__main__':
	from hs3_utils import Flight_segment
	from qsubmissions import lblrtm_hs3_lim, lbldis_hs3_lim
	from hs3_utils import find_shiscplgdas
	import re

	#_yank out label
	out_label = gen_outlabel(**namelist)
	
	namelist.update({'out_label' : out_label})
	dbg(out_label)

	############################################################################
	#_real_#####################################################################
	############################################################################
	#_LBL-RTM, REAL FLIGHT LOOP (must be run before retrieval)
	if 0:
	  #_DELETE, LIMIT TO THIS CASE
	  for dtg, values in segments.iteritems():
		#_lblrtm_hs3_lim currently requires that there
		# be a kwargs pickle available to get default opts.
		# This is dumb.
		dbg(dtg)
		namelist['dtg'] = dtg
		namelist.update(values)
		file_seg = find_shiscplgdas(**namelist)
		lblrtm_hs3_lim(file_flight=file_seg, **namelist)		


	#_OPTIMAL ESTIMATION, REAL FLIGHT LOOP
	if 1:
	  if namelist['surf_temp'] == -1:
		#_for some reason, -1 breaks in the simulated retrievals
		raise RuntimeError, 'PUT IT BACK TO GDAS'
	  for dtg, values in segments.iteritems():
		namelist['dtg'] = dtg
		final = 'lbldis_output.{0}.final.cdf'.format(out_label)
		input = 'lbldis_input.{0}.final'.format(out_label)
		clear = 'lbldis_output.{0}.clear.cdf'.format(out_label)
		postd = 'posterior.{0}.pk'.format(out_label)
		namelist.update({	'lbldis_output'	: final, 
							'lbldis_clear'	: clear,
							'lbldis_input'	: input,
							'posterior_dump' : postd, })
		namelist.update(values)
	
		#_run retrieval over real flight segment	
		file_out = real_case_sbatch(**namelist)


	if 0:	#_GENERATE SUMMARY
	  for dtg, values in segments.iteritems():
		namelist['dtg'] = dtg
		final = 'lbldis_output.{0}.final.cdf'.format(out_label)
		clear = 'lbldis_output.{0}.clear.cdf'.format(out_label)
		namelist.update({	'lbldis_output'	: final, 
							'lbldis_clear'	: clear, })
		namelist.update(values)
	
		dbg(''); dbg(dtg)
	#	try:
		write_retrieved_qsub(**namelist)
		#	plot_jeffz_qsub(**namelist)
		#	plot_jeffz_qsub(**namelist)
		#	plot_real_qsub(**namelist)
	#	except:
	#		dbg((dtg, 'FAILED TO PLOT'))


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
	if 0:	#_WITH RETRIEVALS SBATCH
	  ppath = namelist.get('dir_plot')
	  namelist.update({'dir_plot' : os.path.join(ppath,out_label,'DASHBOARD')})
	  for dtg, values in segments.iteritems():
		#_plot comparisons between OD from SHIS and CPL 
		# SEPARATE! CANNOT RUN WITH real_case_qsub()
		dbg(dtg)
		final = 'lbldis_output.{0}.final.cdf'.format(out_label)
		clear = 'lbldis_output.{0}.clear.cdf'.format(out_label)
		namelist.update({	'lbldis_output'	: final, 
							'lbldis_clear'	: clear,
							'dtg'			: dtg })
		namelist.update(values)
		plot_sbatch(**namelist)

	if 0:	#_NO RETRIEVALS
	  for dtg, values in segments.iteritems():
		#_plot comparisons between OD from SHIS and CPL 
		# SEPARATE! CANNOT RUN WITH real_case_qsub()
		dbg(dtg)
		namelist.update(values)
		namelist['out_label'] = 'EVENT-ONLY'	#_change to whatever
		flight = Flight_segment(dtg)
		flight.plot_flight(**namelist)	#_flight only


	#_FULL CAMPAIGN SUMMARIES
	if 0:
		#_send to method that reads in all retrieval data
		# with indiv methods to handle plotting	
		plot_real_summary(**namelist)


	if 0: #_SEE WHAT'S COMPLETE
	  dbg(out_label)
	  for dtg, values in segments.iteritems():
		final = 'lbldis_output.{0}.final.cdf'.format(out_label)
		clear = 'lbldis_output.{0}.clear.cdf'.format(out_label)
		input = 'lbldis_input.{0}.final'.format(out_label)
		pname = '{0}/DASHBOARD/dashboard_{0}_{{0}}_f{{1:04d}}.png'.format(out_label)
		namelist.update({	'lbldis_output'	: final, 
							'lbldis_clear'	: clear,
							'lbldis_input'	: input,
							'pname'			: pname,  })
		namelist.update(values)
		check_missing(dtg, **namelist)		
			

	if 0:
	  for dtg, values in segments.iteritems():
		#_plot comparisons between OD from SHIS and CPL 
		# SEPARATE! CANNOT RUN WITH real_case_qsub()
		plot_real(**namelist)
