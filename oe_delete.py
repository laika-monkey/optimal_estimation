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

#_select flight segements VVVVVVVV change that to shift
from flight_namelists import experiments

#_pull out desired experiment
experiment = 'HS3'
segments   = experiments[experiment] 

print 'WITHT HE NEW RDR, NO SHIS PROFILE DATA, MAKE SURE THAT IT WORKS'
print 'why is THICKNESS in .jacobians as a method twice'
print 'MAKE SURE SURF_EMIS IS GETTING TO WRITE_LBLDIS'
print 'LAYER LOOPING OF APRIORI WILL ONLY WORK IF ALL STATE_VARS ARE IN LAYERS'
print 'SO DONT RETRIEVE SURF_TEMP'

'''
To turn off uncertainties
	'model'			set 'uncertainty' to 0
	'instrument'	set 'uncertainty' to 0
'''

# DIR_LOG  = os.path.expanduser('~/qsub_logs/')
DIR_LOG  = os.environ['LOG']
DIR_PROD = os.environ['PRODUCTS']
DIR_LBL  = 'LBL-RTM_hs3' #_generally want to use just LBL-RTM
DIR_TMP  = '/data/wsessions/TMP'
DEBUG    = 1


#_used instead of config_file
namelist = {
	#_if directory exists, skip
	'rerun'			: False,  

	#_ok, see if this breaks everything
	'experiment'	: experiment,	#_used for LBL I/O naming 
##	'out_label'		: 'ERR_total',		#_for sim, name what is being retreived?
	'out_label'		: 'GDAS',	#_use for retr output and plots
##	'out_label'		: 'GDASmisr',	#_use for retr output and plots
#_fucked up and deleted some of the NVACLIMO when it wasn't working,
# now the naming is back to GDAS/GDASmisr for no good reason, but they mean
# the same
##	'out_label'		: 'NVACLIMO1misr',	#_use for retr output and plots
##	'out_label'		: 'NVACLIMO1',	#_use for retr output and plots
									# the loop designates it to be used
									# only as a kind of top level label
									# when looping over flight segments.

	#_these are separate so that uncertainty
	# can be looped via jv without attempts to retrieve
	# and covar matrix generated
##_Don't retrieve thickness without updating levels
## First, give more level resolution.  Second, update sensitivity
## based upon TAPE7 at1mosphere.
	'state_vars'	: ['tau'],	#_vars to retrieve
	'uncert_vars'	: ['ref','z','thickness','surf_emis','surf_temp','rel_hum'],

	#_location of naaps data, only used if set in dynamic_CTP
	'dir_naaps'		: os.path.join(DIR_PROD, 'NRL', 'NVA_CLIMO1'),
##	'dir_naaps'		: os.path.join(DIR_PROD, 'NRL', 'NVA_CLIMO1misr'),

	#_location of input radiances and output directory for fov plots
	'dir_shis'		: os.path.join(DIR_PROD, experiment.lower()), 
	'dir_out'		: os.path.join(DIR_PROD, experiment.lower()), 

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
	'dynamic_CTP'		: False, #['NAAPS'], # ['NAAPS','SHIS'], False,

	#_format of lblrtm directory, [experiment]_[dtg0.dtg1]_fov[###]
	'dir_lblrtm_fmt': os.path.join(DIR_PROD, DIR_LBL,'{3}_{0}.{1}_fov{2:06d}'),

	# 2. VARIABLES
	'dv'			: -26,			#_microwindow (window, no O3)
##	'dv'			: -4,			#_this was used for simtests
##	'dv'			: -7,			#_microwindow (window)
##	'dv'			: 0.5,			#_resultion of spectrum (cm-1)
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

	#_SIMULATED
	############################################################################
	############################################################################
	#_label used for experiment (image output, flight file input for sim case)
	# fsim overwrites sim_experiment if both passed. If only sim_exp, fsim is
	# generated in methods as needed.
##	'sim_experiment'	: 'tau_ref',
	'fsim'	: 'SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.nc',
##	'fsim'	: 'SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.NOISE.0004485.nc', RUN
##	'fsim'	: 'SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.NOISE.0021783.nc', NOT
##	'fsim'	: 'SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.NOISE.0030252.nc', NOT
	#_noise [0003934, 0010976, 0014129, 0018558, 0022304, 0024691, 0032751


	#_directory of simulated run static case
	'dir_lblrtm'		: '/data/wsessions/LBL-RTM/test_fov',
##	'out_label'			: ['TAU_BD-BASELINE','TAU_BD-Z', 'TAU_BD-THICKNESS',
##							'TAU_BD-EMISSIVITY', 'TAU_BD-TEMPERATURE',
##							'TAU_BD-REL_HUM'],
##	#_END_SIMULATED
##	############################################################################
##	############################################################################

##  ssps = ['ssp_db.shettle_dust.gamma_sigma_0p100',
##          'ssp_db.shettle_sand.gamma_sigma_0p100',
##			'ssp_db.mie_gypsum.lognormal_sigma_0p699',
##			'ssp_db.mie_kaolinite.lognormal_sigma_0p699',
##			'ssp_db.mie_quartz.lognormal_sigma_0p699',	]
##			'ssp_db.mie_wat.gamma_sigma_0p100'	]

	#_actual state in test (0.5, 0.5, 3, 1.0)
	#_a priori values for state vars
	'apriori_simulated' : [
					{	'type'		: 'aerosol', #_only use with simulated test
						'tau'		: 0.01, 
					#	'tau'		: 0.01, 
						'ref'		: 0.5,	
						'z'			: 3.0,
						'z_top'		: 3.01,
						'thickness' : 0.5,	
						'rel_hum'	: 1.0,		#_1 == apriori is GDAS
						'ssp_db' : 'ssp_db.shettle_dust.gamma_sigma_0p100',
						'ref_wn'	: 900,	#_to match dummy_profiles.py	
						'surf_emis'	: [[100, 3000], [0.985, 0.985]],
						},
				],
	'apriori' : [
					{	'type'		: 'aerosol',
						'tau'		: 0.01, 
						'ref'		: 0.5,		#_in the middle of af dust 
						'z'			: 3.0,
						'z_top'		: 3.01,
						'thickness' : 0.5,		#_if thickness retrieved, z_top
												# is ignored
						'rel_hum'	: 1.0,		#_1 == apriori is GDAS
						'ssp_db'	: 'ssp_db.mie_gypsum.lognormal_sigma_0p699',
						'ref_wn'	: -1,		#_-1 : geo limit, or waven 900.
						'surf_emis'	: [[100, 3000], [0.985, 0.985]],
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
				#	{	'type'		: 'aerosol',
				#		'tau'		: 0.01, 
				#		'ref'		: 0.5,	
				#		'z'			: 3.0,
				#		'z_top'		: 3.01,
				#		'thickness' : 0.5,	
				#		'ssp_db'	: 'ssp_db.mie_quartz.lognormal_sigma_0p699',
				#		'ref_wn'	: -1,	
				#		},
					{	'type'		: 'cloud',
						'tau'		: 0.01, #1.0, 
						'ref'		: 1.0,
						'z'			: 1.00,
						'z_top'		: 1.01,
						'thickness' : 0.5,	
						'ssp_db'	: 'ssp_db.mie_wat.gamma_sigma_0p100',
						'ref_wn'	: -1,
						},
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
	'nproc'	: 20,

	#_remove input/output files in lbldis
	'clean'	: True,	

	'NETCDF_LIB' : '/opt/netcdf4/4.1.3-intel-14.0-2/lib/',	# netcdf lib 
															# for lbldis
	'ssp_db_files' : [],
	}

#_add in spectral property file options
#_single scattering property file list
dir_ssp = os.path.join(os.environ['PRODUCTS'], 'lbldis_inputs')
ssps = [layer['ssp_db'] for layer in namelist['apriori']]
[namelist['ssp_db_files'].append(os.path.join(dir_ssp, f)) for f in ssps]
[layer.update({'dbnum' : ssps.index(layer['ssp_db'])}) \
	for layer in namelist['apriori']]

#_add levels at sensitivity of z
namelist.update({ 'dummy_levs' : namelist['sensitivity']['z'] })

#_if retrieving thickness, overwrite z_top
if 'thickness' in namelist['uncert_vars']:
	for layer in namelist['apriori']:
		layer['z_top'] = layer['z'] + layer['thickness']


#############################################################################_80
#_real/simulated_cases_#########################################################
################################################################################


def real_case_qsub(fidx=None, rerun=True, **kwargs):
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
	from libtools import mkdir_p
	import pickle

	#_get process number
	pid = os.getpid()

	#_pull out some vars
	out_label = kwargs.get('out_label')
	dir_lblrtm_fmt = kwargs.get('dir_lblrtm_fmt')
	experiment = kwargs.get('experiment')

	#_create directory for plots
	ppath = 'PLOTS_{0}'.format(out_label)
	mkdir_p(ppath)

	#_do entire flight if not specified
	flight = F(**kwargs)
	fidx = xrange(flight.size) if fidx is None else fidx

	#_qsubmission prefix		
	env = ','.join(['='.join((var, os.environ[var])) \
		for var in ['PYTHONPATH', 'PRODUCTS', 'PATH', 'LOG']])
	qsub = ' '.join((   'qsub -v', env, '-o', DIR_LOG, '-e', DIR_LOG,
		'-cwd -S /opt/ShellB3/bin/python'))
	scpt = os.path.expanduser('~/lib/run_oe.py')

	#_Loop over each fov (in this case they're just variable TAU/REF/Z opts)
	for i in fidx:

		#_build output directory name. Args really should be passed
		# if lblrtm_fmt is a thing. A half effort otherwise
		dir_lblrtm = dir_lblrtm_fmt.format(	flight.dtg0, flight.dtg1, i,
											experiment)
		kwargs.update({'dir_lblrtm' : dir_lblrtm})

		#_KLUDGE~DELETE. Skips retrievals already performed.
		outfile = os.path.join(dir_lblrtm, kwargs.get('lbldis_output'))
		if os.path.exists(outfile) and not rerun:
			continue

		#_generate kwarg file that scrpt will load to get settings
		fmt = 'kwgs.{4}.{2}.{3}.{0:04d}.{1:07d}'
		fmt_args = (i, pid, flight.dtg0, flight.dtg1, out_label)
		file_kwarg = os.path.join(DIR_TMP, fmt.format(*fmt_args))
		pickle.dump(kwargs, open(file_kwarg, 'wb')) 

		#_put it all together and submit
		cmd = ' '.join((qsub, scpt, str(i), file_kwarg))
		dbg(cmd)
		os.system(cmd)


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
		dbg(('testing', missing_no, total, nfov))

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
			dtgfov = (flight.dtg0, flight.dtg1, i, experiment)
			dir_lblrtm = dir_lblrtm_fmt.format(*dtgfov)
			lbldis_out = 'lbldis_input.{0}.final'.format(out_label)
			lbldis_out = os.path.join(dir_lblrtm, lbldis_out)

			#_read in final lbldis output file
			try:
				lbldis_out = Lbldis_input(lbldis_out)
			except IOError:
				dbg('Missing {0}'.format(lbldis_out), 3)
				continue
		
			#_sum up aerosol taus
			tau = sum_tau(lbldis_out)

			#_initialize output line tuple with fov number
			tup = [i]
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


def simulated_test_qsub(fsim=None, rerun=True, **kwargs):
	'''
	test against simulated data

	Simulated runs were produced using dummy_profiles.py script
	The values for tau/ref can be found in the hs3 SIMULATED file with
	the .txt extension

	experiment	str,	used for selecting the input SIMULATED data file
	out_label	str,	used for naming of output files, plots, retrievals

	OUTPUT
	Produces a SIMULATED.retrieved text file to be compared with SIMULATED.txt

	If the simulated scenario label has any variables that are
	not in the state_var list, then the apriori is updated to be 
	the exact value.

	THINGS REQUIRED FOR THIS TEST:
		DUMMY S-HIS DATA FROM dummy_profiles.py, output in $PRODUCTS/hs3
		APPROPRIATELY SPACED OUT LBL-RTM LEVELS, in $PRODUCTS/LBL-RTM
	'''
	import matplotlib.pyplot as plt
	from hs3_utils import Flight_segment as F
	from numpy import diag, matrix
	from lblrtm_utils import microwindow_average, microwindow_stdev
	from time import sleep
	from libtools import shrink_ticks, setup_groups, mkdir_p
	from os import fork, waitpid
	from multiprocessing import Process, Pipe
	from pickle import dump

	#_announce
	dbg('RUNNING SIMULATED TEST')
	dbg('uses one particular atmospheric profile from test_fov')
	pid = os.getpid()
	out_label = kwargs.get('out_label')
	sim_experiment = kwargs.get('sim_experiment')
	dir_shis = kwargs.get('dir_shis')
	state_vars = kwargs.get('state_vars')

	#_simulated input file name
	# Either have it be named here from sim_exp, or 
	# be explicit with filename and it will yank out the label
	if fsim is None:
		fsim = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0:s}.nc'.format(sim_experiment)
	else:
		import re
		res = re.search('SIMULATED.(.*?).nc', fsim)
		sim_experiment = res.group(1)		
		kwargs.update({'sim_experiment' : sim_experiment})
	
	#_file_seg==simulated SHIS, file_tru==state vars for file_seg,
	# file_out==where to dump retrieval values
	file_seg = os.path.join(dir_shis, fsim)
	kwargs.update({'file_seg' : file_seg, 'plot' : False})
	flight = F(file_seg=file_seg)

	#_pull out true values
	true_values = _get_true_values(flight.notes) 
	meh = flight.notes.split('\n')

	#_qsubmission prefix		
	env = ','.join(['='.join((var, os.environ[var])) \
		for var in ['PYTHONPATH', 'PRODUCTS', 'PATH', 'LOG']])
	qsub = ' '.join((   'qsub -v', env, '-o', DIR_LOG, '-e', DIR_LOG,
		'-cwd -S /opt/ShellB3/bin/python'))
	scpt = os.path.expanduser('~/lib/run_oe.py')

	#_create directory for plots
	ppath = 'PLOTS_{0}'.format(out_label)
	mkdir_p(ppath)

	#_Loop over each fov (in this case they're just variable TAU/REF/Z opts)
	for i, fov in enumerate(flight):
		
		#_initialize plotting, for first ten instances only
		kwargs.update({	'sup'	: meh[i],
						'comment' : meh[i]	})

		#_name of final output file
		lbldis_output = 'lbldis_output.SIMULATED.' + \
			'{0:04d}.{1}.{2}.final.cdf'.format(i, sim_experiment, out_label)
		lbldis_clear = 'lbldis_output.SIMULATED.' + \
			'{0:04d}.{1}.{2}.clear.cdf'.format(i, sim_experiment, out_label)
		posterior_pickle = 'posterior.SIMULATED.' + \
			'{0:04d}.{1}.{2}.pk'.format(i, sim_experiment, out_label)

		#_KLUDGE~DELETE. Skips retrievals already performed.
		outfile = os.path.join(kwargs.get('dir_lblrtm'), lbldis_output) 
		if os.path.exists(outfile) and not rerun:
			continue

		#_figure out how to best organize this
		# Returns converged solution and all.
		# Maybe switch to just all and pop
		kwgs = kwargs.copy()
		kwgs.update({	'lbldis_output'  : lbldis_output,
						'lbldis_clear'   : lbldis_clear,
						'posterior_dump' : posterior_pickle })

		#_update the apriori with true values if not in the state
		#_ONE LEVEL LIMITATION IN PLACE WRS
		for variable, value in true_values.iteritems():
			if variable in state_vars:
				continue

			#_else put in true value
			kwgs['apriori'][0][variable] = value[i]
	
			if variable == 'z':
				kwgs['apriori'][0]['z_top'] = value[i]+.01
	
		#_write kwarg pickle file
		fkwgs = 'kwargs.f{0:04d}.{1}.{2}.{3}.pk'.format(i, sim_experiment,
			out_label, pid)
		fkwgs = os.path.join(DIR_TMP, fkwgs)
		dump(kwgs, open(fkwgs, 'wb'))

		#_submit job
		dbg(('TESTING FOV NUMBER', i))
		cmd = ' '.join((qsub, scpt, str(i), fkwgs))
		dbg(cmd)
		os.system(cmd)


def _get_true_values(string):
	import re
	o = {}
##	opt = [a.split(',') for a in string.split('\n')]
	for a in string.split('\n'):
		for b in a.split(','):
			res = re.search('(.*?)-([\d.]+)', b)
			varname = '{0}'.format(res.group(1))
			value = float('{0}'.format(res.group(2)))
			try:
				o[varname].append(value)
			except KeyError:
				o.update({varname : [value]})
	return o


def add_ssps(apriori=None, **kwargs): 
	''' pull out ssp files form apriori and put them in kwargs '''
	#_setup path to ssp databases
	dir_ssp = os.path.join(os.environ['PRODUCTS'], 'lbldis_inputs')

	#_pull out all ssp files from layers
	ssps = [layer['ssp_db'] for layer in apriori]

	#_add ssp file to list of all (add a check later to prevent dupes)
	[kwargs['ssp_db_files'].append(os.path.join(dir_ssp, f)) for f in ssps]

	#_add database number associated with ssp
	[layer.update({'dbnum' : ssps.index(layer['ssp_db'])}) for layer in apriori]


def write_simulated_qsub(sim_experiment='tau', out_label=None, dir_lblrtm=None,
	fsim=None, **kwargs):
	''' go back through output from simulated_test_qsub and make output '''
	import matplotlib.pyplot as plt
	from hs3_utils import Flight_segment as F
	from numpy import diag, matrix
	from lblrtm_utils import microwindow_average, microwindow_stdev
	from lblrtm_utils import Lbldis_input
	from time import sleep
	from libtools import shrink_ticks, setup_groups, mkdir_p
	from os import fork, waitpid
	from multiprocessing import Process, Pipe
	from pickle import load
	import re

	#_announce
	dbg('WRITE SIMULATED RESULTS')
	dir_shis = kwargs.get('dir_shis')

	#_simulated input file name
	if fsim is None:
		fsim = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0:s}.nc'.format(sim_experiment)
	else:
		res = re.search('SIMULATED.(.*?).nc', fsim)
		sim_experiment = res.group(1)		
		kwargs.update({'sim_experiment' : sim_experiment})

	#_file_seg==simulated SHIS, file_tru==state vars for file_seg,
	# file_out==where to dump retrieval values
	file_seg = os.path.join(dir_shis, fsim)
	file_out = file_seg.replace('.nc', '.{0}.retrieved'.format(out_label))	
	flight = F(file_seg=file_seg)

	#_get values being tested
	state_vars = kwargs.get('state_vars')

	#_get true values from flight attributes
	true_vals = flight.notes.split('\n')

	#_read in all values used for simulated data generation
##	[re.search('(.*?)-[\d.]+', t).group(1) for t in true_vals[0]]
	header = [re.search('(.*?)-[\d.]+', t).group(1) \
		for t in true_vals[0].split(',')]
	header.extend(state_vars)	#_get known state, then add retrieved

	#_this is supposed to get the real values up front of the lines...
	lines = []
	for truth in true_vals:
		values = [re.search('.*?-([\d.]+)', t).group(1) \
			 for t in truth.split(',')]
		lines.append([float('{0}'.format(v)) for v in values])


	#_output format for results
	fmt = []
	nv = len(state_vars)
	buf = len(header)
	for j, v in enumerate(header):
		fmt.append('{{{0:d}:10.5f}}'.format(j))
	for j, v in enumerate(state_vars):	#_uncertainty format (?)
		fmt.append('{{{0:d}:10.5f}}'.format(buf+j))
	##	fmt.append('{{{0:d}:10.5f}}'.format(nv+j+1))
	fmt = ','.join(fmt) 
	fmt += '\n'

	#_collect output
	line_out = [None] * flight.size
	for i, fov in enumerate(flight):

		#_generate output filename
		lbldis_input = 'lbldis_input.SIMULATED.' + \
			'{0:04d}.{1}.{2}.final'.format(i, sim_experiment, out_label)
		posterior = 'posterior.SIMULATED.' + \
			'{0:04d}.{1}.{2}.pk'.format(i, sim_experiment, out_label)
		fname = os.path.join(dir_lblrtm, lbldis_input)
		pname = os.path.join(dir_lblrtm, posterior)

		try:
			input = Lbldis_input(fname)

			#_pull out retrieved layer
			x = input.layers[0]		#_pull out first layer (dict = {'tau' 'ref'}
			x['tau'] = x['tau'][0]
			x_std = load(open(pname, 'rb'))
	
			#_format retrieved values for output
			tup_retr = [x[v] for v in state_vars]
			tup_retr.extend([x_std[v] for v in state_vars])
			lines[i].extend(tup_retr)
	
		except IOError:
			dbg('WARNING: missing {0}'.format(fname))
			tup_retr = [-9999 for v in state_vars]
			tup_retr.extend([-9999 for v in state_vars])
			lines[i].extend(tup_retr)
			
		if len(input.layers) != 1:
			raise RuntimeError, 'not yet implemented for more than one layer'

		line_out[i] = fmt.format(*lines[i])

	#_write output as we go
	dbg(file_out)

	#_initialize output (fmt produces header with state variable names)
	with open(file_out, 'w') as f:
		fmt =','.join(['{{{0:d}:>10s}}'.format(i) for i,v in enumerate(header)])
		f.write(fmt.format(*header) + '\n')
		[f.write(l) for l in line_out] 


def simulated_test(sim_experiment=None, out_label=None, fsim=None, **kwargs):
	'''
	test against simulated data

	Simulated runs were produced using dummy_profiles.py script
	The values for tau/ref can be found in the hs3 SIMULATED file with
	the .txt extension

	experiment	str,	used for selecting the input SIMULATED data file
	out_label	str,	used for naming of output files, plots, retrievals

	OUTPUT
	Produces a SIMULATED.retrieved text file to be compared with SIMULATED.txt

	THINGS REQUIRED FOR THIS TEST:
		DUMMY S-HIS DATA FROM dummy_profiles.py, output in $PRODUCTS/hs3
		APPROPRIATELY SPACED OUT LBL-RTM LEVELS, in $PRODUCTS/LBL-RTM
	'''
	import matplotlib.pyplot as plt
	from hs3_utils import Flight_segment as F
	from numpy import diag, matrix
	from lblrtm_utils import microwindow_average, microwindow_stdev
	from time import sleep
	from libtools import shrink_ticks, setup_groups, mkdir_p
	from os import fork, waitpid
	from multiprocessing import Process, Pipe

	#_announce
	dbg('RUNNING SIMULATED TEST')
	dbg('uses one particular atmospheric profile from test_fov')
	dir_test = '/data/wsessions/hs3'

	#_simulated input file name
	if fsim is None:
		fsim = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0:s}.nc'.format(sim_experiment)
	else:
		import re
		res = re.search('SIMULATED.(.*?).nc', fsim)
		sim_experiment = res.group(1)		

	#_file_seg==simulated SHIS, file_tru==state vars for file_seg,
	# file_out==where to dump retrieval values
	file_seg = os.path.join(dir_test, fsim)
	file_tru = file_seg.replace('.nc', '.txt')			#_true value file
	file_out = file_seg.replace('.nc', '.retrieved')	#_retr output
	file_out = file_out.replace(sim_experiment, '{0}.{1}'.format(
		sim_experiment, out_label))	#_retr output
	flight = F(file_seg=file_seg)

	#_get values being tested
	state_vars = kwargs.get('state_vars')

	#_read in all values used for simulated data generation
	lines = [l.replace('\n','').replace(',',' ') for l in 
			open(file_tru, 'r').readlines()]
	header = lines[0]; del lines[0]
	header = header.split()
	header.extend(state_vars)

	lines = [[float(ll) for ll in l.split()] for l in lines]

	#_create directory for plots
	ppath = 'PLOTS_{0}'.format(out_label)
	mkdir_p(ppath)

	#_initialize output (fmt produces header with state variable names)
	f = open(file_out, 'w')
	fmt = ','.join(['{{{0:d}:10s}}'.format(i) for i,v in enumerate(header)])
	f.write(fmt.format(*header) + '\n')
	f.close() #_to keep it updating, reopen...

	#_output format for results
	fmt = []
	for j, v in enumerate(header):
		fmt.append('{{{0:d}:10.5f}}'.format(j))

	for j, v in enumerate(state_vars):	#_uncertainty format (?)
		nv = len(state_vars)
		fmt.append('{{{0:d}:10.5f}}'.format(nv+j+1))

	fmt = ','.join(fmt) 
	fmt += '\n'

	#_split fov objects into forkable groups
	groups = setup_groups(range(flight.size), **kwargs)

	#_Loop over each fov (in this case they're just variable TAU/REF/Z opts)
	for fovs in groups:
		nfov = len(fovs)
		pipe = [None]*nfov
		pope = [None]*nfov
		thread = [None]*nfov
		line_out = [None]*nfov
	
		for j, i in enumerate(fovs): #_list index = j, flight index = i
			dbg(('TESTING FOV NUMBER', i, lines[i]))
			fov = flight[i]
			wav = flight.SHIS_wavenumber

			#_init pipe
			pipe[j], pope[j] = Pipe()

			#_initialize plotting, for first ten instances only
			fig, ax = plt.subplots(3,2)
			fig.suptitle(lines[i], size='xx-small')
			pname = '{2}_{1}_{0:06d}.png'.format(i, out_label, sim_experiment) 
			pname = os.path.join(ppath, pname)
			kwargs.update({	'fig' : fig, 'ax' : ax, 'pname' : pname })
						
			if os.path.exists(pname):
				os.unlink(pname)

			#_calculate instrument covariance matrix
			dv = kwargs.get('dv')

			#_when using microwindows, average observation data
			if dv < 0:
				std, d		= microwindow_stdev(fov.SHIS_radiances, wav, dv) 
				nesr, d		= microwindow_average(flight.SHIS_HBB_NESR, wav, dv)
				y, wave_obs = microwindow_average(fov.SHIS_radiances, wav, dv)
	
				#_std deviation within microwindows applied to ref
				# blackbody 
				cov	= matrix(diag(nesr) + diag(std))	
	
			else:
				cov = matrix(diag(flight.SHIS_HBB_NESR.copy()))

##			#_throw in an option to turn off error from obs platform
##			uncert = kwargs.get('uncertainty', {'model':False})
##			if not uncert['model']:
##				cov[:] = 0.0	

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

			#_format retrieved values for output
			tup_retr = [x[v] for v in state_vars]
			tup_retr.extend([x_std[v] for v in state_vars])
			lines[i].extend(tup_retr)
			line_out[j] = fmt.format(*lines[i])

			#_close thread
			thread[j].join()

		#_write output as we go	
		dbg(('TESTING', file_out))
		with open(file_out,'a') as f:
			[f.write(l) for l in line_out] 


################################################################################
#_retrieval_####################################################################################################################################################


def optimal_estimation(fov, wave_obs, Sr, surf_temp=-1, dir_lblrtm='.', 
	ref_wn=900, L_curve_index=4, L_curve_size=151, form='rodgers',
	max_iter=25, out_type=dict, pipe=None, apriori={}, dynamic_CTP=False,
	posterior_dump=None, **kwargs):
	'''
	phys_info
	fov		hs3_utils.Flight_segment record
			Contains collocated GDAS, SHIS, and CPL data

	form    str,    carissimo or rodgers, defines the form of OE
                    as following either:

	                CARISSIMO, 2005
	                The physical retrieval methodology for IASI: the d-IASI code

	                RODGERS, 2000
	                Inverse Methods for Atmospheric Sounding
		
					TURNER, 2014
					Information Content and Uncertainties in Thermodynamic
					Profiles and Liquid Cloud Properties Retrieved from the
					Ground-Based Atmospheric Emitted Radiance Interferometer
					(AERI)
	  
                    Both including the usage of Gamma found in Turner, 2014
	L_curve_index	int,	range over which to do L-curve calculations

	Sr		matrix,	instrument error covariance matrix

	ret_vars	list,	contains strings for var to actually do a retrieval upon

	out_type	type, dict will output a dictionary by state var,
						defaults to matrix
	surf_temp	float,	Surface temperature value to use. Set to -1 to 
						use LBL-ATM surface temperature, set to -2 to use
						GDAS (not implemented yet), or >0 to be arbitrary
	dynamic_CTP		bool||list,	if false, use static apriori cloud top pressure
								otherwise use what's available based on
								priority defined by list. If nothing available,
								will default to static.
	posterior_dump	str,		for qsub write. dumps the dictionary with 
								posterior uncertainties to a pickle. Arg is
								name of file	
	apriori		list(dicts),	definition of initial layers
	'''
	from shis import convolve2shis
	from lblrtm_utils import microwindow_average, get_surf_temp
	from lblrtm_utils import check_reff_v_ssp as check_ssp
	from numpy import append, array, arange, c_, diag, eye, matrix, power, sqrt
	from numpy import zeros, diff, tile, trace, linspace, concatenate
	from scipy.linalg import sqrtm
	from libtools import strictly_increasing
	from copy import deepcopy
	from hs3_utils import get_naaps_ctp
	from pickle import dump

	#_for testing
	import matplotlib.pyplot as plt
	from libgeo import planck, p2z
	from tape7 import tape7

	############################################################################
	#_LINEAR_LEAST_SQUARES_#####################################################
	############################################################################

	#_pull out sensitivity dictionary
	sensitivity = kwargs.get('sensitivity')
	uncertainty = kwargs.get('uncertainty') 
	state_vars = kwargs.get('state_vars', ['tau'])
	uncert_vars = kwargs.get('uncert_vars', ['tau'])
	out_label = kwargs.get('out_label', '')

	nstate = len(state_vars)
	nlayer = len(apriori)

	#_get height of observer from LBL-RTM to limit convergence
	t7 = tape7(os.path.join(dir_lblrtm, 'TAPE7'))
	max_z = t7.H1 - sensitivity['z']*2

	#_set surface temperature, if -1, use from LBL-ATM prof
	surf_temp = get_surf_temp(fov, surf_temp_src=surf_temp,
		dir_lblrtm=dir_lblrtm, **kwargs)
	
	#_define an a priori cloud layer dictionary. 
	cld = deepcopy(apriori) #_first guess

	#_use static cloud top height or use DR SHIS guess
	#_convert pressure to AGL... sorta...
	dbg(dynamic_CTP)
	for layer in cld:
		''' MAKE THIS A METHOD ELSEWHERE '''
		if dynamic_CTP:
			for ctp_source in dynamic_CTP:
				if ctp_source == 'NAAPS':
					arg = (fov.SHIS_epoch,fov.SHIS_latitude,fov.SHIS_longitude) 
					CTP, CBP = get_naaps_ctp(*arg, **kwargs)
					CTZ = p2z(CTP) / 1e3 
					CBZ = p2z(CBP) / 1e3 

					#_if apprioriate layer and value NAAPS, use
					if layer['type'] == 'aerosol' and	(type(CTZ) == list and max(CTZ) > 0) or \
														(type(CTZ) == float and CTZ > 0): 
				#	if layer['type'] == 'aerosol' and CTZ > 0:
						layer['z_top'] = max(CTZ) #CTZ
						layer['z'] = max(CBZ) #CTZ-0.1
						dbg(('using just the top layer', max(CTZ), max(CBZ)))
					else:
						continue

				elif ctp_source == 'SHIS':
					CTZ = p2z(fov.SHIS_CTP) / 1e3

					#_if apprioriate layer and value present, use SHIS DR 
					if layer['type'] == 'aerosol' and CTZ > 0:
						layer['z_top'] = CTZ
						layer['z'] = CTZ-0.1
					else:
						continue
					
				#_set, move on	
				break			

	#_update cloud into kwargs and add analysis sfc temperature
	kwargs.update({'clddef' : cld, 'surf_temp' : surf_temp})

	#_If using full obs spectrum, convolve lbldis to retrieval
	# If using microwindows, average retrieval to channels
	obs_rads = fov.SHIS_radiances
	y_std    = None	#_not currently an option.  Use microwindows 
	dv       = kwargs.get('dv')

	#_when using microwindows, average observation data
	if dv < 0:
		y, wave_obs = microwindow_average(obs_rads, wave_obs, dv)
	else:
		y = obs_rads[:]

	#_MAKE THE NEXT THREE OPERATIONS NOT NECESSARY
	#_and convert to matrix
	y = matrix(y)

	#_make sure matrix is [m,1]
	if y.shape[1] != 1:
		y = y.T

	#_maybe I don't know
	y *= 1e3

	#_generate a clear sky run based on surf temp guess (for conv plots)
	kwargs['r_clr'], kwargs['n_clr'] = clear_sky_radiance(dir_lblrtm, **kwargs)

	#_get values for uncertainty jacobians
	dbg('launching uncertainty calculations')
	Ki, Fx, wave_lbl = jacobian_forker(dir_lblrtm,
		jacobian_vars=uncert_vars, **kwargs)
	dbg('finishing uncertainty calculations')

	#_a priori and definition of the vector of the first guess
	# v and vg in Carismo
	xi = matrix([[a[v] for v in state_vars] for a in apriori]).flatten().T
	xa = matrix([[a[v] for v in state_vars] for a in apriori]).flatten().T

	#_get sizes of measurement and state space
	m = len(y)
	n = len(xi)

	'''
	#_NOT GAMMA DIST, but a regularization factor. (see pveg about range) 
	# See Tikhonov-Miller / Phillips-Twomey / Constrained linear invers. method
	#_CARISSIMO, 2.2.3
	gamma = 10.**linspace(L_curve_index, L_curve_index, L_curve_size) 
	gamma = tile(eye(n), (L_curve_size, 1, 1)) * gamma[:,None,None]

	#_See Carissimo 2005, 2.2 Inversion Methodology
	N			= len(y)			#_size of DOF
	chi2alpha	= N + 2*sqrt(2.*N)	# condition to be <= in eq(3) 
									#_assumes 2sig tolerance interval
	chi2n		= array([1e4, 1e3])	# used in conditional escape below 
	chi2steps	= 1000. 
	'''

	#_(v-vg).T*L*(v-vg), L = B.T = apriori covar
	#_(xi-xa).T*L*(xi-xa)
	############################################################################
	#_optimal_estimation_#######################################################
	############################################################################

	#_initialize collection of all iterations
	all_xi	= xi.copy()		#_keep a running list of state iterations
	all_J	= []			#_list of cost function values
	all_di2	= []			#_DELETE WRS all convergence criteria

	#_arbitrary gamma values for testing
	gamma_arb = [1000, 300, 100, 30, 10, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1]
	
	#_i => OD, Reff (2xMeas)
	#_define covariance matrices
	Sm = matrix(diag(Fx.A.squeeze() * uncertainty['model']))
	Sr = Sr * uncertainty['instrument']
	Sa = array([[uncertainty[unc]	for unc in state_vars] \
									for i in range(nlayer)])
	Sa = matrix(diag(Sa.flatten()))
##	Sa = matrix(diag([[uncertainty[unc]	\
##		for unc in state_vars]			\
##		for i in range(nlayer)]))
##	Sa = matrix(diag(array([uncertainty[unc] for unc in state_vars])))
									#_uncertainty covar of apriori
##	St = sensitivity['surf_temp']	#_include these even though not used?
##	St = uncertainty['surf_temp']	#_include these even though not used?
##	Sh = sensitivity or sensitivity ['z']
	#_s-1 precision/determination of information in prior
	#_total covariance  (dF/dstate * d_unc_state => variance in F by state)
	# uncertainty of SHIS comes from sensor
	# uncertainty of model appears... arbitrary? Do I just trust docs?
	#_assumes uncorrelated error.  Used to normalize ellipse?
	#_Se should include all error that would prevent the Forward model
	# from preproducing y from x.
	Se = Sr + Sm

	#_hrm..... noooooo
	for lidx in range(nlayer):
 	  for jv in uncert_vars:

		#_uncertainty of state vars in Sa, so skip here
		if jv in state_vars or uncertainty[jv] == 0 or	\
			(lidx != 0 and (jv[:4] == 'surf' or jv == 'rel_hum')):	#_skip surface vars
			continue

		#_initialize jacobian
		#_DELETE
		Ktmp = matrix(Ki[lidx][jv]) * uncertainty[jv]

		#_covar for
		Stmp = Ktmp * sensitivity[jv] * Ktmp.T
		Stmp = matrix(diag(diag(Stmp)))

		Se = Se + Stmp

	''' look at eigen decomposition of covariance matrices '''
	'''
	LOOK BELOW FOR OTHER QUOTE and WHY IS SENSITIVITY TO TAU/REF NOT HERE

	CHECK DIM V2 TO SEE IF HE'S INCLUDES THE DESIRED ONES IN Sm and THE NON
	IN Se
	'''
	for oe in xrange(max_iter):
		#_update cloud for jacobians to current guess (x_n), updates at end 
		[[layer.update({sv : xi[i+j*nstate,0]})	\
			for i, sv in enumerate(state_vars)]	\
			for j, layer in enumerate(kwargs['clddef'])]

		#_Fx == Forward Model Obs at Current Guess
		#_calculate jacobian around first guess, get K(jac), F(x), dim of F(x)  
		K, Fx, wave_lbl = jacobian_forker(dir_lblrtm,
			jacobian_vars=state_vars, **kwargs)

		#_plot current iteration DELETE::WRS
		plot(y.copy(), wave_obs, Fx.copy(), wave_lbl, label='OE', **kwargs)

		#_jacobians (build matrix from jacobian dictionary)
		for lidx in range(nlayer):
			
			for i, sv in enumerate(state_vars):
				tmp = K[lidx][sv] if K[lidx][sv].ndim == 2 \
								else K[lidx][sv][:,None]

				#_if first state variable and first layer, initialize
				if not lidx and not i:
					Ki = tmp
				else:	
					Ki = append(Ki, tmp, 1)

		Ki = matrix(Ki)
##		Ki = array([[K[lidx][v].flatten() \
##			for v in state_vars] \
##			for lidx in range(nlayer)])
##		Ki = matrix([[K[lidx][v] \
##			for v in state_vars] \
##			for lidx in range(nlayer)]).T
#_#		#_s-1 precision/determination of information in prior
#_#		#_total covariance  (dF/dstate * d_unc_state => variance in F by state)
#_#		# uncertainty of SHIS comes from sensor
#_#		# uncertainty of model appears... arbitrary? Do I just trust docs?
#_#		#_assumes uncorrelated error.  Used to normalize ellipse?
#_#		#_Se should include all error that would prevent the Forward model
#_#		# from preproducing y from x.
#_#		Se = Sr
#_#		for jv in jacobian_vars:
#_#			#_uncertainty of state vars in Sa
#_#			if jv in state_vars or uncertainty[jv] == 0:
#_#				continue
#_#
#_#			#_covar for
#_#			Stmp = K[jv] * sensitivity[jv] * K[jv].T
#_#			Se = Se + Stmp 
##		Se = Sr + Sm # + Kt.T*St*Kt	# + Kh*Sh*Kh.T
##		Se = Sr + Sm # + Kt*St*Kt.T	# + Kh*Sh*Kh.T	#_orig?
##		sum([K[sv] * uncertainty[sv] * K[sv].T for sv in state_vars]) 

		#_choose between formulation of Carissimo (2005) or Rodgers (2000)
		if form == 'carissimo':
			pass
			'''
			#_change of variables to make the problem dimless, easier to
			# manipulate and more efficient.
			#_CARISSIMO ET AL. (2005) - The physical retrieval methodology for
			#							IASI: the delta-IASI code
			Yn	= (y - Fx) + Ki*(xi - xa)#_current x minus first guess minus xa?
	
			#_the inverse of the variance of all errors * dF/dstate
			Ji	= matrix(inv(sqrtm(Se)) * Ki)	#_CARISSIMO eq.18~
			G	= Ji*sqrtm(Sa)			#_Gn = Jn*sqrt(B), B=a priori covar mtrx
			U,D2,V = svd(G.T*G)
			V	= U						#_matrix being factorized is symmetric
			z	= inv(sqrtm(Se))*Yn		#_same reference

			#_new radiance values
			# (SHIS minus LBL-DIS rads + dF/dx(first guess minux a priori)
	
			#_HANSEN (1992) - 	Analysis of discrete ill-posed problems by means
			#					of the L-curve.
			curv = []
			for g in range(0): #range(L_curve_size):	#_can this be indexed
	
				#_CARISSIMO eq.20: Solve for u == inv(Sa)*x_n+1
				# u_n+1 = V inv(gI+D2)V.T G.T z
				u = V * inv(gamma[g] + D2) * V.T * G.T * z
	
				#_first and second derivative of u
				u1 = -V * inv(gamma[g]*eye(n) + D2) * V.T*u
				u2 = -2. * inv(gamma[g]*eye(n) + D2) * V.T*u1

##				Should all be dimensionless
##				a(gamma)= u.Tu
##				a(g)'	= u.T'u + u.Tu'
##				a(g)"	= u.T"u + 2u.T'u' + u.Tu" 	
##	
##				b(gamma)= (Gu-z).T(Gu-z)
##				b(g)'	= (Gu').T(Gu-z) + (Gu-z).T(Gu')
##				b(g)"	= (Gu").T(Gu-z) + 2(Gu').T(Gu') + (Gu-z).T(Gu")

				a1 = u1.T*u + u.T*u1				#_first der of param a
				a2 = u2.T*u + 2.*u1.T*u1 + u.T*u2	#_first der of param a
				b1 = (G*u1).T*(G*u-z) + (G*u-z).T*(G*u1)
				b2 = (G*u2).T*(G*u-z) + 2.*(G*u1).T*(G*u1) + (G*u-z).T*(G*u2)
	
				#_Curvature of the L-curve may be then expressed in terms of gamma
				# CARISSIMO eq.23
				cgamma	= abs(a1*b2 - a2*b1) / power((a1**2 + b1**2), 1.5)
				curv	= append(curv, cgamma)
##				this is currently a mess.  dimensions wrong? 
	
			#_search for the best gamar regularization param for the current
			# iter  Gamma is a ratio of prior information to obs, so g > 1
			#
			# should an array of gamma 1000, 300, 100, 30... show up, it's from 
			# TURNER 2014
			gammamax = 1 # curv.max()

			#_build giant matrices, calculate all at once.
			#_optimal estimation rearranged according to Carissimo 2005
			# eq.21
			u = V * inv(gammamax*eye(n) + D2) * V.T * G.T * z
	
			#_??? state vector for current iteration, w/ dims
			# Sa*u == x_n+1  (is this square root correct and is it just a step)
			#_carissimo above eq.20, u_n+1 definition
			xoe = sqrtm(Sa) * u + xa 

			#_approximate the current state vector
			xoe[0] = round(xoe[0].real, 3)
			xoe[1] = round(xoe[1].real, 2)

			#_chisquared test stuff
			chi2current	= (y-Fx).T * inv(Se) * (y-Fx)
			chi2n		= append(chi2n, chi2current)
			chi2steps	= (chi2n[-2] - chi2n[-1]) / chi2n[-2]
	
			print chi2current.shape, chi2n[-1], N, n, gam
			#_set exit conditions
			#_CHARISSIMO (3), 2.2.1 Iterization
			# (R - Fv).T*S^-1*(R-Fv) <= chi2alpha
			if chi2n[-1] < chi2alpha and gam == 1:		#_and di2 < 0.01:
	##		if chi2current < chi2alpha and gam == 1:	#_and di2 < 0.01:
				out_status = 1	#_good
				print 'GAMMA STATIC'
				break
			elif chi2steps <= 0.1:
				out_status = -1	#_bad, oscillation or slow convergeance
				break
			elif chi2n[-1] >= chi2n[-2]:
				out_status = -3	#_bad, residuals increasing
				break
			else:
				#_still converging... this is superfluous.
				xi = xoe.copy()
				continue
			'''

		#_RODGER_###############################################################
		elif 'rodgers':
			try:
				gam = gamma_arb[oe]	# curv.max()
			except IndexError:
				gam = 1.

			#_RODGERS eq.5.9
			#_n-form
			dbg((xa.shape, Sa.shape, Ki.shape, Se.shape, y.shape, Fx.shape))
			dbg('DOCTOR')
			xoe = xa + inv(gam*inv(Sa) + Ki.T*inv(Se)*Ki) * Ki.T \
					* inv(Se) * (y-Fx + Ki*(xi-xa)) 
			#_m-form
			#	= xa + Sa * Ki.T * inv(Ki*Sa*Ki.T + Se) * (y-Fx + Ki*(xi-xa))

			#_ensure output is within bounds
			xoe = check_state_boundaries(xoe, **kwargs)

			#_add to previously computed retrievals
			all_xi = concatenate((all_xi, xoe), axis=1)

			#_TURNER 2014 eq.3-5, calculating the posterior error covar
			B	= gam*inv(Sa) + Ki.T*inv(Se)*Ki
			A	= inv(B) * Ki.T*inv(Se)*Ki		#_averaging kernel	
			Sxi = inv(B) * (gam**2 * inv(Sa) + Ki.T*inv(Se)*Ki) * inv(B)
			degf= trace(A)

			#_RODGERS eq.5.29
			# di2 = (x_i - x_i+1)' * S^-1 * (x_i - x_i+1) << n
			di2 = (all_xi[:,oe]-all_xi[:,oe+1]).T * inv(Sxi) \
				* (all_xi[:,oe]-all_xi[:,oe+1])
			all_di2.append(di2)	#_DELETE WRS for testing

			#_update initial point
			xi = xoe.copy()
			dbg(('CURRENT STATE','\n',xi))
	
			#_MAP / Minimum cost function.  Find CF
			J = (xi-xa).T*inv(Sa)*(xi-xa) + (Fx-y).T*inv(Se)*(Fx-y)
			all_J.append(J)

##			crit_val = (xi - xoe).T*inv(Sa)*(xi - xoe) 
##			if crit_val < n/100:
			dbg(('TESTING_CONV', oe, di2, J, di2 < 0.01, J.A[0] < 0.01))

##			if di2 < 0.001 and gam == 1:	
			if J.A[0] < 0.01 and gam == 1:	
				#_convergeance
				break

		elif 'jacobs':
			#_SEE SECTION 7 FOR JACOBS TREATMENT USING ADJOINT
			pass

	#_update final iteration and rerun with no cloud and converged cloud
	[[layer.update({sv : xi[i+j*nstate,0]})			\
		for i, sv in enumerate(state_vars)]	\
		for j, layer in enumerate(kwargs['clddef'])]
	lbldis_final(dir_lblrtm, **kwargs)

	#_make sure plot is complete DELETE WRS
	plot(y.copy(), wave_obs, Fx.copy(), wave_lbl,
		label='OE_LAST', last=True, **kwargs)
	dbg(os.getcwd())
	try:
		plot_close(**kwargs)
	except IOError:
		dbg('closing image failed')

	#_prepare output type	
	if out_type == dict:

		#_initialize output dictionary
		out = {}
		xi = array(xi)
	
		#_create dictionary for retrieval uncertainty
		posterior_dict = {}
		posterior_std = sqrt(diag(Sxi))
	
		#_build output jacobian dictionary keyed with state_var order
		fname_unc = 'uncertainties_LOOP-FOR-LAYERS_{0}.txt'.format(out_label)
		with open(os.path.join(dir_lblrtm, fname_unc), 'w') as f:
			for i, v in enumerate(state_vars):
				out.update({v : xi[i,0]})
				posterior_dict.update({v : posterior_std[i]})

				f.write('{0:10s} {1:10.5f}\n'.format(v, posterior_std[i]))

		#_dump posterior
		if posterior_dump is not None:
			pkname = os.path.join(dir_lblrtm, posterior_dump)
			dump(posterior_dict, open(pkname, 'wb'))
			
		if pipe is None:
			return out, all_xi, posterior_dict

		else:
			pipe.send((out, all_xi, posterior_dict))
			pipe.close()

	#_don't use this... it's functionally broken
	else:
		return xi	


################################################################################
################################################################################
################################################################################


def del_init_layers(dbnum=0, ssp_db_files='./dbfile', apriori={}, **kwargs):
	'''
	given the namelist, initialize layers

	This has been added as an initial step toward 
	bimodal jacobian juggling. Probably going to work 
	to merge this with the boundary checks after each
	iteration.

	The potentially more sensible option is to define 
	each layer individually, not as a fraction of a single
	layer. Or do mulitple retrievals over individual ssps.

	Regardless, this isn't a great idea.

	2015.04.20

	It was a bad idea, deprecated

	2015.05.10
	'''
	ssp_dist = apriori['ssp_dist']

	#_loop over db files and define
	if type(dbnum) == list:
		cld = []
		for i, db in enumerate(dbnum):
			cld.append({ 'dbnum' : db, 'ref_wn' : kwargs.get('ref_wn', -1) })
			cld[i].update(apriori)
			cld[i]['tau'] *= ssp_dist

			#_locked into bimodal
			ssp_dist = 1 - ssp_dist
	
			if i > 1:
				raise RuntimeError, 'Only currently supporting bimodal dists.'
	
	else:
		#_fract is 1, so no multiply
		cld = [{ 'dbnum' : dbnum, 'ref_wn' : kwargs.get('ref_wn', -1) }]
		cld[0].update(apriori)	#_first guess

	return cld


def jacobian_array(K, state_vars=['tau'], **kwargs):
	''' put jacobian dictionary into a matrix '''
	#_loop over jacobian_vars and put into a 2D array
	from numpy import matrix
	for i, sv in enumerate(state_vars):
		if not i:
			tmp = K[sv][:,None]
		else:
			tmp = append(tmp, K[sv][:,None], axis=1)

	#_update return value
	return matrix(tmp)


def jacobians(dir_lblrtm, sensitivity={}, jacobian_vars=['tau', 'ref'],
	clddef=[], output_array=True, pipe=None, apriori_rads=None, nsend=1,
	clean=False, lidx=0, **kwargs):
	'''
	Builds an MxN matrix for the jacobian describing the relationship 
	of the forward model and the state variables linearized around
	certain values.  

	fov				hs3_utils.Flight_segment record containing single FOV for
							collocated SHIS, CPL, and GDAS data

	sensitivity		dict,	dictionary of the sensitivities of various
							state variables
	jacobian_vars	list,	strings of keys to sensitivity.dict above
							for the variables that will make up the jacobian

	clddef			list,	list for cloud layers, see write_lbldis_params
							gets updated here to that method.
	output_array	bool,	loosey goosey way of doing this.  When true,
							output Jacobian is in a single array with the
							column order following jacobian_vars.  Otherwise,
							output as 1D arrays of dictionary labeled by 
							jacobian_vars
	lidx			int,	layer index within clddef to perturb

	Each run is performed depending on the presence of the label in the
	jacobian_vars list.  Do not try to loop over that and be slick about it.
	There are too many variable specific things that will come into play
	down the road to do that.

	The Jacobian will be linearized around the cloud values passed in
	clddef, then shifted according to the corresponding change values
	found in sensitivity
	'''
	from lblrtm_utils import write_lbldis_parameters as lbldis_write
	from numpy import append, matrix, zeros
	from copy import deepcopy
	import time
	from shis import convolve2shis
	from tape7 import tape7

	dbg(dir_lblrtm)
	dbg(jacobian_vars)
	
	#_get processing ID
	pid = os.getpid()

	#_lbldis gets persnickety about this
	if 'LB_LIBRARY_PATH' not in os.environ:
		os.environ['LD_LIBRARY_PATH'] = '/opt/netcdf-4.1.2/lib/'

	#_create a shortcut
	path = os.path.join

	#_move to working directory
	dir_pwd = os.getcwd()
	os.chdir(dir_lblrtm)

	#_check to make sure tau values are in list
	for cloud in clddef:
		if type(cloud['tau']) != list:
			cloud['tau'] = [cloud['tau']]

	#_initialize dictionary
	jacobian = {}

	#############################
	#_optical_depth_sensitivity_#
	#############################
	#_always do optical depth, but keep tabbed for symmetry
	if 'tau' in jacobian_vars:	

		#_make a copy to futz with here
		cloud = deepcopy(clddef) 

		#_add second instance of lbldis run to tau field in clddef
	##  2015.04.21 WRS removed. Let write_lbldis split this?
	##	also with the inclusion of lidx, I don't want every tau arb inc'd 
	##	dOD = sensitivity['tau'] / len(cloud)	#_split inc among layers(1)
	##	for cld in cloud:
	##		cld['tau'].append(cld['tau'][0] + dOD) #_make sure this ok 
		dOD = sensitivity['tau']
		for l, cld in enumerate(cloud):
			cld['tau'].append(cld['tau'][0] + dOD * (lidx == l))

		#_create a name for this run based upon state variables tested
		label = 'tau_{0:5.3f}_{1:07d}'.format(cloud[lidx]['tau'][-1], pid)
		fname = path(dir_lblrtm, 'lbldis.input.' + label)
		oname = 'lbldis.out.{0}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kwargs)

		if call(['lbldis', fname, '0', oname]):
			dbg('lbldis failed to run {0}'.format(os.getcwd))
			raise RuntimeError, 'lbldis failed to run\n{0}\n{1}'.format(fname, oname)

		ncdf = netcdf.netcdf_file(oname+'.cdf','r')
		rads = ncdf.variables['radiance'][:] #_[nchan, ninstances]
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave_tmp = wave.copy()
			wave, tmp0 = convolve2shis(wave_tmp, rads[:,0])
			wave, tmp1 = convolve2shis(wave_tmp, rads[:,1])
	
			n_rads = zeros((wave.size, 2))
			n_rads[:,0] = tmp0.copy().squeeze()
			n_rads[:,1] = tmp1.copy().squeeze()
			rads = n_rads.copy()

		#_find first column of Jacobian (dMeasurement / dState)
		try:
			jacobian['tau'] = ((rads[:,1] - rads[:,0]) / dOD)
		except IndexError:
			jacobian['tau'] = ((apriori_rads - rads[:,0]) / dOD)

		#_keep a copy of the apriori rads
		apriori_rads = rads[:,0].copy()

		if clean:
			[os.unlink(ff) for ff in [fname, '{0}.cdf'.format(oname)]]

		if pipe is not None:
		##	for i in range(nsend):
		##		pipe.put((jacobian, matrix(apriori_rads).T, wave))
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			sys.exit(0)	#_don't exit? Need apriori below?
	
	################################
	#_effective_radius_sensitivity_#
	################################
	elif 'ref' in jacobian_vars: 
		
		#_make a copy to futz with here
		cloud = deepcopy(clddef)
	##	for cld in cloud:
	##		cld['ref'] = cld['ref'] + sensitivity['ref']
		cloud[lidx]['ref'] = cloud[lidx]['ref'] + sensitivity['ref']

		label = 'ref_{0:5.3f}_{1:07d}'.format(cloud[lidx]['ref'], pid)
		fname = path(dir_lblrtm, 'lbldis.input.{0}'.format(label))
		oname = 'lbldis.out.{0}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kwargs)
 
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run', fname
		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_this version assumes that rads is from one instance
		tmp = (rads - apriori_rads) / sensitivity['ref']
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T
		jacobian['ref'] = tmp 

		if clean:
			[os.unlink(ff) for ff in [fname, '{0}.cdf'.format(oname)]]

		if pipe is not None:
			t= [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)

	##############################
	#_cloud_position_sensitivity_#
	##############################
	elif 'z' in jacobian_vars: 

		#_make a copy to futz with here
		cloud = deepcopy(clddef)
		cloud[lidx]['z'] = cloud[lidx]['z'] + sensitivity['z']
		cloud[lidx]['z_top'] = cloud[lidx]['z_top'] + sensitivity['z']

		label = 'z_{0:5.3f}_{1:07d}'.format(cloud[lidx]['z'], pid)
		fname = path(dir_lblrtm, 'lbldis.input.{0}'.format(label))
		oname = 'lbldis.out.{0:s}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kwargs)
 
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run'
		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_find third column of jacobians
		#_this version assumes that rads is from on instance
		tmp = (rads - apriori_rads) / sensitivity['z']
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T

		jacobian['z'] = tmp 

		if clean:
			[os.unlink(ff) for ff in [fname, '{0}.cdf'.format(oname)]]

		if pipe is not None:
			t = [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)
	
	###########################
	#_cloud_depth_sensitivity_#
	###########################
	elif 'thickness' in jacobian_vars: 

		#_make a copy to futz with here
		cloud = deepcopy(clddef)
		thickness = cloud[lidx]['thickness'] + sensitivity['thickness']
		cloud[lidx]['z_top'] = cloud[lidx]['z'] + thickness 

		label = 'thickness_{0:5.3f}_{1:07d}'.format(thickness, pid)
		fname = path(dir_lblrtm, 'lbldis.input.{0}'.format(label))
		oname = 'lbldis.out.{0:s}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kwargs)
 
		#call(['lbldis', fname, '0', oname])
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run'
		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_this version assumes that rads is from on instance
		tmp = (rads - apriori_rads) / sensitivity['thickness']
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T

		jacobian['thickness'] = tmp 

		if clean:
			[os.unlink(ff) for ff in [fname, '{0}.cdf'.format(oname)]]

		if pipe is not None:
			t = [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)

	###################################
	#_surface_temperature_sensitivity_#
	###################################
	elif 'surf_temp'in jacobian_vars:
		if lidx:
			sys.exit(0)

		#_make a copy to futz with here
		kw = kwargs.copy()
##		kw['surf_temp'] = clddef[0]['surf_temp'] + sensitivity['surf_temp']
		kw['surf_temp'] = kwargs['surf_temp'] + sensitivity['surf_temp']
		cloud = deepcopy(clddef)
	##	for cld in cloud:
	##		cld['surf_temp'] = cld['surf_temp'] + sensitivity['surf_temp']

		label = 'surf_temp_{0:5.3f}_{1:07d}'.format(kw['surf_temp'], pid)
		fname = path(dir_lblrtm, 'lbldis.input.' + label)
		oname = 'lbldis.out.{0:s}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kw)
 
		#call(['lbldis', fname, '0', oname])
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run'
		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_find third column of jacobians
		#_this version assumes that rads is from on instance
		tmp = (rads - apriori_rads) / sensitivity['surf_temp']
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T

		jacobian['surf_temp'] = tmp 

		if clean:
			[os.unlink(ff) for ff in [fname, '{0:s}.cdf'.format(oname)]]

		if pipe is not None:
			t = [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)
		pass

	##################################
	#_surface_emissivity_sensitivity_# 
	##################################
	elif 'surf_emis' in jacobian_vars:
		if lidx:
			sys.exit(0)
		#_make a copy to futz with here
		kw = kwargs.copy()
		kw['surf_emis'] = clddef[0]['surf_emis'][:]
		for k, e in enumerate(kw['surf_emis'][1]):
			kw['surf_emis'][1][k] = min((e + sensitivity['surf_emis'], 1))
		cloud = deepcopy(clddef)

		label = 'surf_emis_{0:5.3f}_{1:07d}'.format(kw['surf_emis'][1][0], pid)
		fname = path(dir_lblrtm, 'lbldis.input.' + label)
		oname = 'lbldis.out.{0:s}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kw)
 
		#call(['lbldis', fname, '0', oname])
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run'
		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_find third column of jacobians
		#_this version assumes that rads is from on instance
		tmp = (rads - apriori_rads) / sensitivity['surf_emis']
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T

		jacobian['surf_emis'] = tmp 

		if clean:
			[os.unlink(ff) for ff in [fname, '{0:s}.cdf'.format(oname)]]

		if pipe is not None:
			t = [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)

	###############################
	#_cloud_thickness_sensitivity_#
	###############################
	#_WRS WHY IS THIS IN HERE TWICE
	elif 'thickness' in jacobian_vars: 
		dbg("MAKE SURE THERE IS A WAY FOR THICKNESS TO BE SHRUNK")

		#_make a copy to futz with here
		cloud = deepcopy(clddef)
		for cld in cloud:
			#_need to pass current layer base and top, get returned new base/top
			dz_orig = cld['z_top'] - cld['z']
			cld['z'], cld['z_top'], dz = increment_thickness(dir_lblrtm, **cld)

		#_change in thickness
		dz -= dz_orig

		label = 'thickness_{0:5.3f}_{1:07d}'.format(dz, pid)
		fname = path(dir_lblrtm, 'lbldis.input.' + label)
		oname = 'lbldis.out.{0:s}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kwargs)
 
		#call(['lbldis', fname, '0', oname])
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run'
		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_find third column of jacobians
		#_this version assumes that rads is from on instance
		tmp = (rads - apriori_rads) / dz 
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T

		jacobian['thickness'] = tmp 

		if pipe is not None:
			t = [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)

	################################
	#_effective_radius_sensitivity_#
	################################
	elif 'ssp_dist' in jacobian_vars: 
		if lidx:
			sys.exit(0)
		
		#_make a copy to futz with here
		cloud = deepcopy(clddef)
		for cld in cloud:
			cld['ref'] = cld['ref'] + sensitivity['ref']

		label = 'ref_{0:5.3f}_{1:07d}'.format(cloud[0]['ref'], pid)
		fname = path(dir_lblrtm, 'lbldis.input.{0}'.format(label))
		oname = 'lbldis.out.{0}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kwargs)
 
		#call(['lbldis', fname, '0', oname])
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run'
		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_this version assumes that rads is from one instance
		tmp = (rads - apriori_rads) / sensitivity['ref']
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T
		jacobian['ref'] = tmp 

		if clean:
			[os.unlink(ff) for ff in [fname, '{0:s}.cdf'.format(oname)]]

		if pipe is not None:
			t= [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)

	############################################################################
	#_state_mods_lblrtm_########################################################
	#_water_vapor_profile_#
	#######################
	elif 'rel_hum' in jacobian_vars:
		dbg('rel_hum {0}'.format(lidx))
		if lidx:
			dbg('exiting')
			sys.exit(0)

		from libtools import mkdir_p
		from lblrtm_utils import link_tape3, read_tape5, write_tape5
		from shutil import rmtree

		#_pull out shift
		delta = sensitivity['rel_hum']

		#_construct labels
		label = 'rel_hum_{0:5.3f}_{1:07d}'.format(delta, pid)

		#_create new LBLRTM run directory 
		dir_lblrtm_n = '{0}.{1}'.format(dir_lblrtm, label)

		#_read in old tape5
		prof, pres, t5kw = read_tape5(os.path.join(dir_lblrtm, 'TAPE5'))

		#_directory crap
		cwd = os.getcwd()
		mkdir_p(dir_lblrtm_n)
		os.chdir(dir_lblrtm_n)
		dbg("Changing directoris {0} {1}".format(cwd, dir_lblrtm_n))

		#_shift new tape5 and write
		varr = '{0:s}_{1:s}'.format(t5kw['profile_source'], 'relative_humidity')
		setattr(prof, varr, getattr(prof, varr) + delta)
		write_tape5(prof, pres, **t5kw)

		try:
			link_tape3(**kwargs)
		except OSError:
			dbg(('WARNING: tape3 link failure', dir_lblrtm, dir_lblrtm_n))

		#_run LBLRTM and move back into run directory
		dbg('running lblrtm')
		call(['lblrtm'])
		dbg('completed lblrtm')	
	
		#_lbldis inputs
		cloud = deepcopy(clddef)
		fname = path(dir_lblrtm_n, 'lbldis.input.' + label)
		oname = 'lbldis.out.{0:s}'.format(label)
		lbldis_write(dir_lblrtm_n, filename=fname, clddef=cloud, **kwargs)
 
	#	call(['lbldis', fname, '0', oname])
		dbg('calling lbldis')
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run'
		dbg('lbldis_complete')

		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_find third column of jacobians
		#_this version assumes that rads is from on instance
		tmp = (rads - apriori_rads) / delta 
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T

		jacobian['rel_hum'] = tmp 
	
		if clean:
			dbg('Think about how to clean up LBLRTM sections')
			[os.unlink(ff) for ff in [fname, '{0}.cdf'.format(oname)]]

		#_WRS 2015.09.14 moves this from below lblrtm call... unsure why
		# it was up there.
		os.chdir(cwd)

		#_delete tmp rel_hum directory
		if clean:
			rmtree(dir_lblrtm_n)
			
		if pipe is not None:
			t = [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)
	
	#_move back into calling directory
	os.chdir(dir_pwd)

	'''
	#_put jacobian in  
	if output_array and type(jacobian_vars) == list:
		#_loop over jacobian_vars and put into a 2D array
		for i, sv in enumerate(jacobian_vars):
			if not i:
				tmp = jacobian[sv][:,None]
			else:
				tmp = append(tmp, jacobian[sv][:,None], axis=1)

		#_update return value
		jacobian = tmp
	'''
	return jacobian, matrix(apriori_rads).T, wave


def clear_sky_radiance(dir_lblrtm, clddef=None, lbldis_clear='lbldis.clear.cdf',
	**kwargs):
	'''
	dummy run to drop aerosol / cloud layers from retrieval

	Produces a file with the expected name of lbldis_output.clear.cdf
	that is commonly used in other places as a clear sky bounding for plots.
	
	lbldis_clear	str,	name of clear file output
	'''
	from scipy.io import netcdf
	from shis import convolve2shis
	from lblrtm_utils import write_lbldis_parameters as lbldis_write
	from numpy import zeros

	if 'LB_LIBRARY_PATH' not in os.environ:
		os.environ['LD_LIBRARY_PATH'] = '/opt/netcdf-4.1.2/lib/'

	cwd = os.getcwd()
	os.chdir(dir_lblrtm)

	#_get a label from lbldis clear
	label = '.'.join(lbldis_clear.split('.')[1:-1])

	#_create a name for this run based upon state variables tested
	fname = os.path.join(dir_lblrtm, 'lbldis_input.' + label)
	oname = 'lbldis_output.{0}'.format(label)
	lbldis_write(dir_lblrtm, filename=fname, clddef=[], **kwargs)

	#call(['lbldis', fname, '0', oname])
	if call(['lbldis', fname, '0', oname]):
		raise RuntimeError, 'lbldis failed to run'
	ncdf = netcdf.netcdf_file(oname+'.cdf','r')
	rads = ncdf.variables['radiance'][:] #_[nchan, ninstances]
	wave = ncdf.variables['wnum'][:]

	#_if not using microwindows, convolve to sensor
	if kwargs.get('dv') > 0:
		wave_tmp = wave.copy()
		wave, tmp0 = convolve2shis(wave_tmp, rads[:,0])

		n_rads = zeros((wave.size, 2))
		n_rads[:,0] = tmp0.copy().squeeze()
		rads = n_rads.copy()

	os.chdir(cwd)
	return rads, wave 


def lbldis_final(dir_lblrtm, clddef=None, clean=True,
	lbldis_output='lbldis_output.final.cdf', **kwargs):
	'''
	dummy run to create final retrieval cdf file

	Produces an file of the lbldis_output.final.cdf that is
	expected elsewhere.  This should probably be updated to 
	include out_label.
	'''
	from lblrtm_utils import write_lbldis_parameters as lbldis_write

	if 'LB_LIBRARY_PATH' not in os.environ:
		os.environ['LD_LIBRARY_PATH'] = '/opt/netcdf-4.1.2/lib/'

	pid = os.getpid()
	cwd = os.getcwd()
	os.chdir(dir_lblrtm)

	#_extract label
	label = '.'.join(lbldis_output.split('.')[1:-1])

	#_create a name for this run based upon state variables tested
	fname = os.path.join(dir_lblrtm, 'lbldis_input.{0}'.format(label))
	oname = 'lbldis_output.{0}'.format(label)
	lbldis_write(dir_lblrtm, filename=fname, clddef=clddef, **kwargs)

	if call(['lbldis', fname, '0', oname]):
		raise RuntimeError, 'lbldis failed to run'

	os.chdir(cwd)


def increment_thickness(dir_lbl, z=None, z_top=None, **kwargs):
	'''
	Expands passed layer the smallest geometric amount based
	on atmosphere in lblrtm directory

	dir_lbl	str,	path to lblrtm output
	z		flt,	layer base
	z_top	flt,	layer top
	'''
	from tape7 import tape7
	from numpy import append, diff

	atm = tape7(os.path.join(dir_lbl, 'TAPE7'))
	
	zmin = atm.zlevel1.min()
	zmax = atm.zlevel2.max()

	#_get max index
	imax = atm.zlevel1.size - 1

	#_get thickness and midlevel points
	depths = append(diff(atm.zlevel1), diff(atm.zlevel2)[-1])	
	midlev = atm.zlevel1 + depths/2

	#_get current level indices, idx0 will correspond to depth/midlev, idx1 not
	idx = range(atm.zlevel1.size)
	idx0 = abs(z - atm.zlevel1).argmin()
	idx1 = abs(z_top - atm.zlevel2).argmin()
	idx0 = idx0 if z > atm.zlevel1[idx0] else idx0-1
	idx1 = idx1 if z_top <= atm.zlevel2[idx1] else idx1+1

	#_figure out if we're going up or down
	#_First, are we at the edge already?
	if idx0 == 0 and idx1 == imax: 
		return z, z_top, atm.zlevel2[-1] - atm.zlevel1[0]
	if idx0 == 0:
		direction = 'up'
	elif idx1 == imax:
		direction = 'down'
	else: #_expand in shortest direction
		direction = 'up' if depths[idx0-1] > depths[idx1+1] else 'down'

	#_expand
	if direction == 'up':
		idx1 += 1
		z_top = midlev[idx1] #atm.zlevel2[idx1+1]
		z = midlev[idx0]
	elif direction == 'down':
		idx0 -= 1
		z_top = midlev[idx1]
		z = midlev[idx0] #atm.zlevel1[idx0-1]

	#_get new geometric depth
	dz = atm.zlevel2[idx1] - atm.zlevel1[idx0]

	return z, z_top, dz


def jacobian_forker(dir_lblrtm, jacobian_vars=['tau', 'ref'],  **kwargs):
	'''
	Treads each column of jacobian matrix
	dir_lblrtm		str,	location of parent lblrtm output for LBL-DIS
	jacobian_vars	list,	list of state variables for columns of jacobian

	returns
	K		dict,	jacobians		
	'''
	from libtools import setup_groups
	from multiprocessing import Process, Queue, Pipe
	from numpy import append
	dbg(jacobian_vars)

	#_make list of state vars without tau
	state = jacobian_vars[:]
	try:
		del state[state.index('tau')]
	except ValueError:
		pass

	#_Ki{sv}[lidx] || Ki[lidx]{sv}
	nlayer = len(kwargs.get('clddef'))

	#_start process for each column of jacobian
	args	= (dir_lblrtm,)
	nv		= len(state) 
	thread	= [None]*nv
	pipe	= [None]*nv
	pope	= [None]*nv

	#_initialize output for change, a priori values, and wavenumbers
	K    = [{} for i in range(nlayer)]
	rads = [{} for i in range(nlayer)]
	wave = [{} for i in range(nlayer)]

	#_need to get apriori_rads from first tau case, launch that first
	#_dict, matrix, array
	for lidx in range(nlayer):
		dbg('Getting tau a apriori radiances: nlayer {0}'.format(nlayer))
		K[lidx], rads[lidx]['tau'], wave[lidx]['tau'] =	\
			jacobians(dir_lblrtm, lidx=lidx, jacobian_vars=['tau'], **kwargs)

	apriori = rads[0]['tau'].copy()

	#_loop over each layer
	for lidx in range(nlayer):

		#_setup fork groups
		groups = setup_groups(state, **kwargs)
		for group in groups:

		  #_launch threads for non-Tau cases
 		  for i, sv in enumerate(group):

			#_skip first layer's tau
			#if (lidx and sv[:4] == 'surf'):
			if (not lidx and sv == 'tau')	\
			or (lidx and (sv[:4] == 'surf' or sv == 'rel_hum')):
				continue

			#_open communication pipe
			pipe[i], pope[i] = Pipe()

			#_launch subprocess
			kwargs.update({	'pipe'			: pope[i],
							'lidx'			: lidx, 
							'jacobian_vars'	: [sv], 
							'apriori_rads'	: apriori.A.squeeze() })
			thread[i] = Process(target=jacobians, args=args, kwargs=kwargs)
			thread[i].start()
	
		  #_collect output
 		  for i, sv in enumerate(group):

			#_skip first layer's tau
			if (not lidx and sv == 'tau') \
			or (lidx and (sv[:4] == 'surf' or sv == 'rel_hum')):
				continue

			#_get data
			Ktmp, rads[lidx][sv], wave[lidx][sv] = pipe[i].recv()
			K[lidx].update(Ktmp)
	
			#_wait for it to complete
			thread[i].join()

	#_set return vals... to just the last one in the list?  seems bad.
	wavenum = wave[0]['tau']
	apriori = rads[0]['tau']

##	wavenum = wave[0][sv]
##	apriori = rads[0][sv]

	#_put in desired format
	return K, apriori, wavenum
	
	
################################################################################
#_plotting_#####################################################################
################################################################################

	
def get_all_residuals(out_label, **kwargs):
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

	#	#_build this for checking what is missing
	#	final = 'lbldis_output.{0}.final.cdf'.format(out_label)
	#	clear = 'lbldis_output.{0}.clear.cdf'.format(out_label)
	#	kwargs.update({	'lbldis_output'	: final, 
	#						'lbldis_clear'	: clear, })
	#	kwargs.update(values)
	#	if miss == total:
	#		dbg('No data for {0}'.format(dtg))
	#		continue
	
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
	##	cpl_aod = append(col_pbl[:,None], col_ele[:,None], axis=1).max(1)

		#_get residual difference between CPL and S-HIS retrievals
		residuals = append(residuals, data.tau - cpl_aod)
	
	return residuals


def plot_real_summary(out_label='', residuals=None, **kwargs):
	'''
	fname	str,	path to output file of optimal estimation
	
	do full campaign comparison	
	'''
	import matplotlib.pyplot as plt
	from numpy import append, vstack, empty, linspace, tile 
	from hs3_utils import Flight_segment as F
	from libtools import epoch2iso, shrink_ticks
	from numpy import random
	from libtools import shrink_ticks
	import locale
	import matplotlib
	locale.setlocale(locale.LC_ALL, 'en_US')
	 
	#_initialize residual
	residuals = get_all_residuals(out_label, **kwargs) if residuals is None \
		else residuals 

	#_initialize artist objects
	fig, ax = plt.subplots(2,1) #_2 rows, 1 column, histogram and??

	#_setup histogram bins
	xmin, xmax = -1, 1
	bins = linspace(xmin, xmax, 101)

	#_gen
	rand_ok0 = random.normal(size=5000, loc=0.0, scale=0.025)
	rand_ok1 = random.normal(size=2000, loc=0.1, scale=0.050)
	rand_ok  = append(rand_ok0, rand_ok1)
	total_ok = append(residuals, rand_ok)
	total_nk = tile(residuals, 3)

	nstr = locale.format("%d", total_ok.size, grouping=True)
##	garbage = random.normal(size=10858, loc=0, scale=total_ok.std())	
	#_create histogram
	tmp = ax[0].hist(total_ok, bins=bins, alpha=0.5, label='w/ cloud')
	tmp = ax[0].hist(total_nk, bins=bins, alpha=0.5, label='w/o cloud')
##	tmp = ax[0].hist(garbage, bins=bins, alpha=0.5, label='garbage')
	ax[0].legend(loc='upper right', fontsize=8)
	ax[0].set_xlabel('S-HIS AOT minus CPL AOT', size='xx-small')
	ax[0].set_ylabel('Field of Views', size='xx-small')
	ax[0].set_xticks(linspace(xmin, xmax, 21))
	ax[0].text(-0.9, 200, 'n={0}'.format(nstr), fontsize=6)
##	ax[0].text(-0.9, 200, 'n={0}\n$\sigma_blue={1}\n\sigma_green={2}$'.format(nstr, total_ok.std(), total_nk.std()), fontsize=6)

	#_uncertainties
	def read_jeff(fname='/home/wsessions/REID_UNCERT_ALL.csv'):
		from numpy import array
		f = open(fname, 'r')
		vars = [v.strip() for v in f.readline().split(',')]
		data = {}
		[data.update({var : []}) for var in vars]
		for line in f.readlines():
			vals = line.split(',')
			for i, var in enumerate(vars):
				data[var].append(float(vals[i]))
		f.close()
		for var in vars:
			data[var] = array(data[var])
		rel_hum = data['temperature'] * 1.24		
		data['RH'] = rel_hum
		return data 

	#_read uncertainties
	uncert = read_jeff()
	x = uncert['TRUE_AOD']
	del uncert['TRUE_AOD']
	vars = uncert.keys()
	for i, var in enumerate(vars):
		ax[1].plot(x, uncert[var], label=var.upper())

	ax[1].set_xscale('log')		
	ax[1].set_xlim(0.05, 2.0) 
	ax[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	ax[1].set_xticks([0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 2])
#	ax[1].set_xticks(x)

	ax[1].set_xlabel('AOT', size='xx-small')
	ax[1].set_ylabel('Fractional Uncertainty', size='xx-small')
	ax[1].legend(loc='upper right', fontsize=8)

	#_shrink	
	[shrink_ticks(a) for a in ax]

	#_name file and save
	pname = 'histogram.shisoe-cpl.summary.{0}.png'.format(out_label)
	plt.savefig(pname)
	dbg(pname)


def plot_jeff_qsub(dtg=None, fidx=None, out_label='', dz_cpl=25, smooth=True,
	**kwargs):
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

	#_generate input file name
	flight = F(dtg=dtg, **kwargs) 
	fname = flight.FLIGHT_file.replace('.nc',
		'.{0}.retrieved'.format(out_label))

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
	ax0.legend()
	ax0.set_xticklabels(data.fov)
	ax0.grid(True)
	ax0.get_xaxis().set_visible(False)
	ax0.set_ylabel('Optical Thickness (532 nm)', size='xx-small')
	ax0.set_xlim([0,x.size])
	dbg(('test', x.size))
	shrink_ticks(ax0)

	#_CPL nonsense
	cpl_plt, cpl_x, cpl_y = flight.map_cpl(**kwargs) 
	cpl_max = cpl_plt.shape[1]
	cpl_nx, cpl_ny = cpl_plt.shape
	CTP = masked_where(flight.SHIS_CTP <= 0, flight.SHIS_CTP)
	CTZ = p2z(CTP)
	CTZ[CTZ <= 0] = -9999
	cb = ax1.pcolormesh(cpl_x, cpl_y[:450], cpl_plt.T[:450,:], vmin=0, vmax=1e-4,
			cmap=spectral, zorder=0)
	ax1.scatter(cpl_x, CTZ, marker='x', linewidth=0.5, s=4, color='yellow',
			zorder=1)
	xlim = ax1.xaxis.get_data_interval()
	ax1.set_xlim(xlim)
	ax1.set_ylim(0, cpl_y[:450][-1])
	ax1.set_yticks(linspace(0, cpl_y[:450][-1], 10))
	ax1.set_yticklabels(['{0:4.1f}'.format(vv) for vv in 
		linspace(0, 9, 10)])
###	ax1.set_yticks(linspace(0, cpl_y[-1], 11))
###	ax1.set_yticklabels(['{0:4.1f}'.format(vv) for vv in 
###		linspace(0, 18, 11)])
##	ax1.set_yticklabels(['{0:4.1f}'.format(vv) for vv in
##		linspace(0, cpl_y[-1], 11)*dz_cpl/1e3])
	ax1.set_xticks(arange(0, cpl_nx, cpl_nx/5))
	ax1.set_xticklabels([epoch2iso(ttt)[-8:] for
		ttt in flight.CPL_epoch[::cpl_nx/5]])
	shrink_ticks(ax1)
	ax1.set_ylabel('Height (km)', size='xx-small')
	ax1.set_xlabel('532 nm Backscatter', size='xx-small')
##	ax1.set_ylabel('HEIGHTS TOO LOW ~1KM', size='xx-small')
	ax1.plot([f0, f0], [0, cpl_y[-1]], 'r-')
	ax1.plot([f1, f1], [0, cpl_y[-1]], 'y-')
##	ax1.plot([f0, f0], [0, cpl_ny], 'r-')
##	ax1.plot([f1, f1], [0, cpl_ny], 'y-')
			
	fig.suptitle('Retrieval Comparison', size='xx-small')

	#_name file and save
	pname = fname.replace('.retrieved','.png')
	plt.savefig(pname)
	dbg(pname)



def plot_real_qsub(dtg=None, fidx=None, out_label='', dz_cpl=25, smooth=True,
	**kwargs):
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
	pname = fname.replace('.retrieved','.png')
	plt.savefig(pname)
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


def plot(r_obs, n_obs, r_lbl, n_lbl, fig=None, clddef={}, surf_temp=0, dv=0,
	pname='default.png', label='', last=False, r_clr=None, n_clr=None, sup=None,
	**kwargs):
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

	if clear:
		if r_clr.ndim != 1:
			r_clr = array(r_clr).squeeze()

		t_clr = planck_inv(r_clr/1e3, n_clr*100, domain='wavenumber')

	#_find which axes are empty
	i = 0
	for i, ax in enumerate(fig.axes):
		#_save final plot for solution
		if last:
			ax = fig.axes[-1]

			#_set supertitle if not set
			if fig is not None and sup is not None:
				fig.suptitle(sup, size='xx-small')

			break	
		elif not len(ax.lines) and i != len(fig.axes)-1: #_don't fill last
			#_plot empty, get out of loop
			break
		else:
			#_line not empty, continue
			continue
	else:
##		#_all plots full, get out
##		if not os.path.exists(pname):
##			fig.tight_layout()
##			plt.savefig(pname)
		
		return 

	#_plot new brightness temperatures and difference between last
	if dv > 0:
		arg = { 'linewidth' : 0.3 }
		ax.plot(n_obs, t_obs, color='r', **arg)
		ax.plot(n_lbl, t_lbl, color='g', **arg)
		if clear:
			ax.plot(n_clr, t_clr, color='#CCCCFF', **arg)

	elif dv < 0: #_only plot for microwindows
		arg = { 'linewidth' : 0.3, 'marker' : 'x' }
		ax.plot(n_obs, t_obs, color='r', **arg)
		ax.plot(n_lbl, t_lbl, color='g', **arg)
		ax.scatter(n_obs, t_obs, color='r', **arg)
		ax.scatter(n_lbl, t_lbl, color='g', **arg)
		if clear:
			ax.plot(n_clr, t_clr, color='#CCCCFF', **arg)

##	ax.plot(n_obs, t, color='k', **arg)
	#_if going to do differences, match them up
		
	#_set labels
	shrink_ticks(ax) 
	ax.grid(True)
	ax.set_xlim(n_lbl[0]-20, n_lbl[-1]+20)
	ax.set_title(tit, size='xx-small')


def plot_simulated_histogram(out_label='', plot_var='tau',
	dir_hs3=os.path.join(os.environ['PRODUCTS'],'hs3'), **kwargs):
	'''
	produce some plots for this god damn stuff

	out_label	str,	experiment is not used so that it can
						denote specific input files as opposed 
						to output.  in_label should have been used.
	'''
	import matplotlib.pyplot as plt
	from libtools import mkdir_p, shrink_ticks
	from numpy import array, append

	#_read in retrieval data
	fname = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.retrieved'.format(out_label)
	fname = os.path.join(dir_hs3, fname)
	try:
		truth, retr, uncrt = read_simulated_retrieval(fname, **kwargs)
	except IOError:
		dbg(('Did you complete this retrieval yet?', out_label))
		os._exit(23)

	#_what are we working with?
	state_vars = retr.keys() # truth.keys()
	dbg(state_vars)

	#_kick back if not available
	if plot_var not in state_vars:
		dbg(('Variable not retrieved:', plot_var))
		return

	#_color dictionary
	order = ['low', 'mid', 'high', 'very_high']
	ranges = {
		'low'	: {
			'color'	: 'blue',
			'tau'	: (-9999., 0.5),
			'ref'	: (0.5, 0.675),
		},
		'mid'	: {
			'color'	: 'green',
			'tau'	: (0.5, 1.0),
			'ref'	: (0.625, 0.75),
		},
		'high'	: {
			'color'	: 'orange',
			'tau'	: (1.0, 1.5),
			'ref'	: (0.75, 0.875),
		},
		'very_high'	: {
			'color'	: 'red',
			'tau'	: (1.5, 2.0),
			'ref'	: (0.875, 1.0),
		},
		}

	#_output directory for plots
	ppath = 'PLOTS_{0}'.format(out_label)
	if not os.path.exists(ppath):
		mkdir_p(ppath)

	#_produce for each state variable
	for sv in state_vars:
		#_get residual
		res = truth[sv] - retr[sv]

		#_initialize figure and axes
		fig = plt.figure()
		axes = []
		figtext = []

		ymax, xmax, xmin = 0, 0, 0
		nbin = 100
		axes.append(plt.subplot2grid((1,1), (0,0)))

		#_plot histogram of residuals
		(n, bins, patches) = axes[0].hist(res, bins=nbin,
				facecolor='b', edgecolor='b')
		xmax = max((xmax, bins.max()))
		xmin = min((xmin, bins.min()))
		ymax = max((ymax, n.max()))
			
		axes[0].plot([0,0],[0,n.max()], 'k--')

		axes[0].set_xlabel('diff in tau, {0}'.format(sv), size='x-small')
		axes[0].set_title('retrieval residuals'.format(out_label), 
							size='xx-small')

	##	#_DELETE WRS
	##	if sv == 'tau':
	##	  with open('REID_HIST.csv','w') as f:
	##		for nn, bb in zip(n,bins):
	##			fmt = '{0:7.5f}, {1:7.5f}\n'.format(bb,nn)
	##			f.write(fmt)

		#_get how high we can go and still be 95 % certain of the AOD
		m95 = truth[sv][certainty_garbage(res)]
		
		axes[0].set_ylabel('residuals, truth - oe_retrieval', size='x-small')

		#_max symmetric
		xbnd = max((abs(xmin), abs(xmax)))

		#_shrink...ticks...
		[ax.set_xlim(-xbnd, xbnd) for ax in axes]
		[ax.set_ylim(top=100) for ax in axes]
		[shrink_ticks(ax) for ax in axes]
	
		#_create image name
		pname = 'hist_{0}_{1}.png'.format(out_label, sv)
		pname = os.path.join(ppath, pname)
		plt.savefig(pname)
		plt.close()		
		
		dbg(pname)	


def plot_simulated_retrieval_by_var(var, out_label='', plot_var='tau',
	fsim=None, dir_hs3=os.path.join(os.environ['PRODUCTS'],'hs3'), **kwargs):
	'''
	produce some plots for this god damn stuff

	out_label	str,	experiment is not used so that it can
						denote specific input files as opposed 
						to output.  in_label should have been used.
	'''
	import matplotlib.pyplot as plt
	from libtools import mkdir_p, shrink_ticks
	from numpy import array, append
	import re

	res = re.search('SIMULATED.(.*?).nc', fsim)
	sim_experiment = res.group(1)

	#_read in retrieval data
	fname = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.{1}.retrieved'.format(
		sim_experiment, out_label)
	fname = os.path.join(dir_hs3, fname)
	try:
		truth, retr, uncrt = read_simulated_retrieval(fname, **kwargs)
	except IOError:
		dbg(('Did you complete this retrieval yet?', out_label))
		os._exit(23)

	#_what are we working with?
	state_vars = retr.keys() # truth.keys()
	dbg(state_vars)

	#_kick back if not available
	if plot_var not in state_vars:
		dbg(('Variable not retrieved:', plot_var))
		return

	#_color dictionary
	order = ['low', 'mid', 'high', 'very_high']
	ranges = {
		'low'	: {
			'color'	: 'blue',
			'tau'	: (-9999., 0.5),
			'ref'	: (0.5, 0.675),
		},
		'mid'	: {
			'color'	: 'green',
			'tau'	: (0.5, 1.0),
			'ref'	: (0.625, 0.75),
		},
		'high'	: {
			'color'	: 'orange',
			'tau'	: (1.0, 1.5),
			'ref'	: (0.75, 0.875),
		},
		'very_high'	: {
			'color'	: 'red',
			'tau'	: (1.5, 9999),
			'ref'	: (0.875, 1.0),
		},
		}

	#_output directory for plots
	ppath = 'PLOTS_{0}'.format(out_label)
	if not os.path.exists(ppath):
		mkdir_p(ppath)

	#_get potential settings for value
	values = list(set(truth[var]))
	nval = len(values)
	values.sort()

	''' make loop over pages, 6 per page '''
	nrow = 6
	npag = nval / nrow #_do we need this? 
	page = 0

	#_produce for each state variable
	for sv in state_vars:

		for j, value in enumerate(values):
			dbg((sv, var, value))

			#_initialize new figure
			if not j % nrow:
				page += 1
				fig = plt.figure()
				k = 0

			#_get indices of limit
			idx = truth[var] == value

			#_get residual
			res = truth[sv][idx] - retr[sv][idx]
	
			#_initialize figure and axes
			axis = fig.add_subplot(nrow, 1, k+1)
	
			ymax, xmax, xmin = 0, 0, 0
			nbin = 20
	
			#_plot uncertainty by truth
			fract = uncrt[sv][idx] / truth[sv][idx]
			x_truth, y_fract = [], []
			for value_truth in set(truth[sv]):
				x_truth.append(value_truth)
				y_fract.append(fract[truth[sv][idx] == value_truth])	
	
			##	axis.scatter(truth[sv][idx], fract, marker='x', color=col)
			axis.boxplot(y_fract, positions=x_truth, widths=0.05)
	
			axis.set_xticks(list(set(truth[sv])))
	
			#_get how high we can go and still be 95 % certain of the AOD
			m95 = truth[sv][certainty_garbage(res)]
			
			#_various things for uncertainty
			axis.set_ylim(bottom=0, top=2)
			axis.set_xlim(truth[sv].min()-0.1, truth[sv].max()+0.1)
			axis.set_ylabel('{0} = {1}'.format(var, value), size='x-small')
##			axis.set_xlabel('{0}, max_95 {1:5.2f}'.format(sv,m95), size='x-small')
			axis.grid(False)
			
			#_shrink...ticks...
			shrink_ticks(axis)

			#_save when page full	
			
			if k == (nrow-1) or j == nval-1: 
				#_create image name
				pname = 'hist_uncert.{0}.{1}.{2}.by_{3}.p{4:02d}.png'.format(
					sim_experiment, out_label, sv, var, page)
				pname = os.path.join(ppath, pname)
				plt.savefig(pname)
				plt.close()		
		
				dbg(pname)

			k += 1


def plot_simulated_retrieval(out_label='', plot_var='tau', fsim=None,
	dir_hs3=os.path.join(os.environ['PRODUCTS'],'hs3'), **kwargs):
	'''
	produce some plots for this god damn stuff

	out_label	str,	experiment is not used so that it can
						denote specific input files as opposed 
						to output.  in_label should have been used.
	'''
	import matplotlib.pyplot as plt
	from libtools import mkdir_p, shrink_ticks
	from numpy import array, append
	import re

	res = re.search('SIMULATED.(.*?).nc', fsim)
	sim_experiment = res.group(1)

	#_read in retrieval data
	fname = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.{1}.retrieved'.format(
		sim_experiment, out_label)
	fname = os.path.join(dir_hs3, fname)
	try:
		truth, retr, uncrt = read_simulated_retrieval(fname, **kwargs)
	except IOError:
		dbg(('Did you complete this retrieval yet?', out_label))
		os._exit(23)

	#_what are we working with?
	state_vars = retr.keys() # truth.keys()
	dbg(state_vars)

	#_kick back if not available
	if plot_var not in state_vars:
		dbg(('Variable not retrieved:', plot_var))
		return

	#_color dictionary
	order = ['low', 'mid', 'high', 'very_high']
	ranges = {
		'low'	: {
			'color'	: 'blue',
			'tau'	: (-9999., 0.5),
			'ref'	: (0.5, 0.675),
		},
		'mid'	: {
			'color'	: 'green',
			'tau'	: (0.5, 1.0),
			'ref'	: (0.625, 0.75),
		},
		'high'	: {
			'color'	: 'orange',
			'tau'	: (1.0, 1.5),
			'ref'	: (0.75, 0.875),
		},
		'very_high'	: {
			'color'	: 'red',
			'tau'	: (1.5, 9999),
			'ref'	: (0.875, 1.0),
		},
		}

	#_output directory for plots
	ppath = 'PLOTS_{0}'.format(out_label)
	if not os.path.exists(ppath):
		mkdir_p(ppath)

	#_produce for each state variable
	for sv in state_vars:
		#_get residual
		res = truth[sv] - retr[sv]

		#_get indices for ranges
	 	for label, opts in ranges.iteritems():
			min_val, max_val = opts[sv]
		#	idx0 = truth[sv] >= min_val 
		#	idx1 = truth[sv] < max_val
		#	idx = append(idx0[:,None], idx1[:,None], axis=1).all(axis=1)
			idx = (truth[sv] >= min_val) * (truth[sv] < max_val)
			ranges[label]['idx'] = idx

		#_initialize figure and axes
		fig = plt.figure()
		axes = []
		axis = plt.subplot2grid((2,4), (1,0), colspan=len(ranges.keys()))
		figtext = []

		ymax, xmax, xmin = 0, 0, 0
		nbin = 20
		for i, label in enumerate(order): #ranges):
			min_val, max_val = ranges[label][sv]

			idx = ranges[label]['idx']
			col = ranges[label]['color']
			if idx.sum() == 0: 
				continue	

			axes.append(plt.subplot2grid((2,4), (0,i)))

			#_plot histogram of residuals
			(n, bins, patches) = axes[i].hist(res[idx], bins=nbin,
											facecolor=col, edgecolor=col)
			xmax = max((xmax, bins.max()))
			xmin = min((xmin, bins.min()))
			ymax = max((ymax, n.max()))
			
			#_add label
			figtext.append('{0}\nmin {1:5.1f}\nmax {2:5.1f}'.format(label,
								min_val, max_val))

			axes[i].set_xlabel('difference in {0}'.format(sv), size='x-small')
			
			#_plot uncertainty by truth
			fract = uncrt[sv][idx] / truth[sv][idx]
			x_truth, y_fract = [], []
			for value in set(truth[sv][idx]):
				x_truth.append(value)
				y_fract.append(fract[truth[sv][idx] == value])	

		##	axis.scatter(truth[sv][idx], fract, marker='x', color=col)
			axis.boxplot(y_fract, positions=x_truth, widths=0.05)
		##	axis.plot(truth[sv][idx], fract, marker='x', color=col)
		##	axis.scatter(truth[sv][idx], uncrt[sv][idx], marker='x', color=col)

		axis.set_xticks(list(set(truth[sv])))
		axes[0].set_title('{0} (experiment label)'.format(out_label), 
							size='xx-small')

		#_get how high we can go and still be 95 % certain of the AOD
		m95 = truth[sv][certainty_garbage(res)]
		
		#_various things for uncertainty
		axis.set_ylim(bottom=0)
		axis.set_xlim(truth[sv].min()-0.1, truth[sv].max()+0.1)
		axis.set_ylabel('posterior uncercertainty', size='x-small')
		axis.set_xlabel('{0}, max_95 {1:5.2f}'.format(sv,m95), size='x-small')
		axis.grid(False)
		
		axes[0].set_ylabel('residuals, truth - oe_retrieval', size='x-small')

		#_max symmetric
		xbnd = max((abs(xmin), abs(xmax)))
		for i, label in enumerate(order):	
			axes[i].text(0, ymax-ymax/5., figtext[i], size='xx-small')

		#_shrink...ticks...
		[ax.set_xlim(-xbnd, xbnd) for ax in axes]
		[ax.set_ylim(top=ymax) for ax in axes]
		[shrink_ticks(ax) for ax in axes]
		shrink_ticks(axis)
	
		#_create image name
		pname = 'hist_uncert.{0}.{1}.{2}.png'.format(sim_experiment,
			out_label, sv)
		pname = os.path.join(ppath, pname)
		plt.savefig(pname)
		plt.close()		
		
		dbg(pname)


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


def read_simulated_retrieval(fname, **kwargs):
	'''
	Read in output file produced by optimal estimation test
	on simulated data. File should be space deliminated and 
	contain true values for state variables, retrieved values,
	and the uncertainty (unlabled columns in same order as sv)

	fname	str,	full path to file
	'''
	from numpy import array
	
	dbg(('opening', fname))
	with open(fname, 'r') as f:
		lines = [l.strip() for l in f.readlines()]

		#_pull out state variables
		header_vars = [l.strip() for l in lines[0].split(',')]
		state_vars = []
		for v in header_vars:
			if v not in state_vars:
				state_vars.append(v)
		nsv = len(state_vars)
		retr_vars = header_vars[nsv:]
		nrv = len(retr_vars) 
		del lines[0]

		#_pull out data
		truth, retrieved, uncertainty = {}, {}, {}
		[truth.update({sv:[]}) for sv in state_vars]
		[retrieved.update({sv:[]}) for sv in retr_vars]
		[uncertainty.update({sv:[]}) for sv in retr_vars]

		#_loop over lines and fill arrays
		for line in lines:
			cols = line.split(',')

			#_loop over state var columns
			for i, variable in enumerate(state_vars):
				truth[variable].append(float(cols[i]))

			#_loop over retrieve var columns (values and uncertainties)
			for i, variable in enumerate(retr_vars):
				retrieved[variable].append(float(cols[i+nsv]))
				uncertainty[variable].append(float(cols[i+nsv+nrv]))

	for sv in state_vars:
		truth[sv] = array(truth[sv])

	for rv in retr_vars:
		retrieved[rv] = array(retrieved[rv])
		uncertainty[rv] = array(uncertainty[rv])

	return truth, retrieved, uncertainty


################################################################################
#_tools_########################################################################
################################################################################


def check_state_boundaries(x, state_vars, sensitivity, style='max', 
	clddef=[], **kwargs):
	'''
	Looks through current state values against permitted
	values in the SSP file and other constraints such
	as observation height.
	'''
	from lblrtm_utils import check_reff_v_ssp as check_ssp
	from tape7 import tape7

	nlayer = len(clddef)
	nstate = len(state_vars)

	#_do check for each layer
	for lidx in range(nlayer):

		#_check optical depth
		if 'tau' in state_vars:
			tidx = state_vars.index('tau') + lidx * nstate
			x[tidx,-1] = 0.0001 if x[tidx,-1] < 0 else x[tidx,-1]
	
		#_check effective radius
		if 'ref' in state_vars:
			ridx = state_vars.index('ref') + lidx * nstate	
 
			x[ridx,-1], style = check_ssp(x[ridx,-1],
				dbnum=clddef[lidx]['dbnum'], **kwargs) 
			if style == 'max':
				x[ridx,-1] -= sensitivity['ref']
	
		#_check height of layer top and bottom
		if 'z' in state_vars:
			max_z = tape7(os.path.join(dir_lblrtm, 
					'TAPE7')).H1-sensitivity['z']*1.5
	
			zidx = state_vars.index('z') + lidx * nstate
			x[zidx,-1] = 0.01 if x[zidx,-1] < 0 else x[zidx,-1]
			x[zidx,-1] = max_z if x[zidx,-1] > max_z else x[zidx,-1]
	
		#_make sure bimodal ssp dist fraction is not outside [0,1]
		if 'ssp_dist' in state_vars:
			tidx = state_vars.index('ssp_dist')
			if x[tidx,-1] < 0:
				x[tidx,-1] = 0
			elif x[tidx,-1] > 1:
				x[tidx,-1] = 1

	return x


def read_retrieval(fidx=None, **kwargs): 
	'''
	read in values passed to  LBL-DIS
	'''
	pass


def certainty_garbage(data, level=95, tolerance=0.1, **kwargs):
	'''
	Find the values between we have a certainty of being within
	this tolerance value
	really garbagey way of testing this...

	And what if we're chopping off the ones that are within tolerance?
	There is no sorting being done.

	'''
	from numpy import arange
	d = data.copy()
	certainty = float((abs(d) < tolerance).sum()) / d.size 
	level = level / 100. if level > 1 else level

	#_generate a list of indices to return
	idx = arange(d.size)

	#_keep chopping off end until we hit it
	while certainty < level:
		d = d[:-1]
		certainty = float((abs(d) < tolerance).sum()) / d.size 

	return d.size - 1  


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

		if not err:
			print scream
		else:
			raise RuntimeError, scream


def gen_outlabel(dv=None, apriori=None, dynamic_CTP=None, out_label=None,
	cld=False, **kwargs):
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
		ctp = 'CTP-{0}'.format('-'.join(dynamic_CTP))
	except:
		ctp = 'NOCTP' 

	#_update out_label
	arg = (out_label, n_habit, abs(dv), cld, ctp, zkm)
	return '{0}_M{1}_DV{2}_{3}_{4}_{5}'.format(*arg)


def check_missing(dtg, fidx=None, dir_lblrtm_fmt=None, lbldis_output='',
	**kwargs):
	from hs3_utils import find_shiscplgdas as find_scg
	import re

	#_get 
	fname = find_scg(dtg, **kwargs)
	redtg = re.search('SHIS.CPL.GDAS.COLLOC.(\d+).(\d+).nc', fname)

	#_pull out start and end dtg
	dtg0 = '20{0}'.format(redtg.group(1))
	dtg1 = '20{0}'.format(redtg.group(2))

	total = 0

	#_loop over fidx to get
	missing_count = 0
	for i in fidx:
		#_build file path
		dir_lblrtm = dir_lblrtm_fmt.format(dtg0, dtg1, i, experiment)
		lbldis = os.path.join(dir_lblrtm, lbldis_output)
	
		#_check if missing
		if not os.path.exists(lbldis):
			missing_count += 1	

		#_wait what	
		total += 1

	#_report
	dbg('{0} missing {1:>6d} / {2:>6d}'.format(dtg, missing_count, len(fidx)))
	return missing_count, len(fidx)


################################################################################
#_SHELL_########################################################################
################################################################################


if __name__ == '__main__':
	from hs3_utils import Flight_segment
	from qsubmissions import lblrtm_hs3_lim, lbldis_hs3_lim
	from hs3_utils import find_shiscplgdas
	import re

	############################################################################
	#_simulated_################################################################
	############################################################################
	if 0:	#_single
		simulated_test(**namelist)
		exit()

	if 0:	#_qsub loop over namelist-defined file
		simulated_test_qsub(**namelist)
		exit()
	if 0:	#_write output from single run
		write_simulated_qsub(**namelist)
		exit()

	if 0:	#_qsub loop over vars/files listed here 
		#_loop over everything
		fnames = ['SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.nc',
				'SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.NOISE.0004485.nc', 
				'SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.NOISE.0021783.nc',
				'SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.NOISE.0030252.nc']
		vnames = ['ref', 'z', 'thickness' , 'surf_emis', 'rel_hum', 
					'instrument', 'surf_temp']
		base_label = namelist['out_label'][:]

		for fname in fnames:

			for i, current_var in enumerate(vnames):

				#_for 'total'
				if i != 0: #_only do total once
					continue
				label = 'ERR_{0}'.format('total')
				label = label.replace('.','_')
				for vname in vnames: 
					namelist['uncertainty'][vname] = 1 

			##	#_for individual	
			##	label = 'ERR_{0}'.format(current_var)
			##	for vname in vnames: 
			##		namelist['uncertainty'][vname] = int(vname == current_var)

				#_finalize label
				if label[-1] == '_':
					label = label[:-1]
				namelist['out_label'] = label
				namelist['fsim'] = fname
	
				#_never have both launch together. stq must complete before wsq
				dbg(label)
				simulated_test_qsub(**namelist)
			##	write_simulated_qsub(**namelist)

		exit()

	#_plot statistics from runs
	if 0:
		plot_simulated_retrieval_by_var('z', **namelist)
	##	plot_simulated_retrieval(**namelist)
		exit()


	#_yank out label
##	out_label = namelist['out_label']
	out_label = gen_outlabel(**namelist)
	##___SEPT 11th MADE THIS TO PLOT OLD THINGS FOR JEFF
##	dbg('NOT GENERATED OUT LABEL')
##	out_label = 'GDAS-2M-DV26-CLD-NOCTP-3KM'
	namelist.update({'out_label' : out_label})


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
		clear = 'lbldis_output.{0}.clear.cdf'.format(out_label)
		namelist.update({	'lbldis_output'	: final, 
							'lbldis_clear'	: clear, })
		namelist.update(values)
	
		#_run retrieval over real flight segment	
		file_out = real_case_qsub(**namelist)


	#_PLOTTING ONLY
	if 0:	#_WITH RETRIEVALS
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
	if 0:	#_GENERATE SUMMARY
	  for dtg, values in segments.iteritems():
		namelist['dtg'] = dtg
		final = 'lbldis_output.{0}.final.cdf'.format(out_label)
		clear = 'lbldis_output.{0}.clear.cdf'.format(out_label)
		namelist.update({	'lbldis_output'	: final, 
							'lbldis_clear'	: clear, })
		namelist.update(values)
	
		dbg(''); dbg(dtg)
		try:
			write_retrieved_qsub(**namelist)
		#	plot_jeff_qsub(**namelist)
	#		plot_real_qsub(**namelist)
		except:
			dbg((dtg, 'FAILED TO PLOT'))


	#_FULL CAMPAIGN SUMMARIES
	if 0:
		#_send to method that reads in all retrieval data
		# with indiv methods to handle plotting	
		plot_real_summary(**namelist)


	if 1: #_SEE WHAT'S COMPLETE
	  dbg(out_label)
	  for dtg, values in segments.iteritems():
		final = 'lbldis_output.{0}.final.cdf'.format(out_label)
		clear = 'lbldis_output.{0}.clear.cdf'.format(out_label)
		namelist.update({	'lbldis_output'	: final, 
							'lbldis_clear'	: clear, })
		namelist.update(values)
		check_missing(dtg, **namelist)		
			

	if 0:
	  for dtg, values in segments.iteritems():
		#_plot comparisons between OD from SHIS and CPL 
		# SEPARATE! CANNOT RUN WITH real_case_qsub()
		plot_real(**namelist)
