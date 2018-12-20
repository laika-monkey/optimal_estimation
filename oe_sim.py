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
from numpy import matrix, linspace
from numpy.linalg import inv, svd
from scipy.linalg import sqrtm
from scipy.io import netcdf
import os
import math
import time
import sys
from optimal_estimation import check_settings

'''
GYPonKAL_27		Done
GYPonKAL_26		Done
KALonGYP_27		Done
KALonGYP_26		Done

KALonKAL_27     Done
GYPonGYP_27     Done 
KALonKAL_26     Done 
GYPonGYP_26     Done 

GYPonKAL_27_tau-ref		THESE WERE ALL SUBMITTED ON 12.2
GYPonKAL_26_tau-ref
KALonGYP_27_tau-ref
KALonGYP_26_tau-ref

BOTHonKAL_26_tau-ref	RESUB 12.7
BOTHonKAL_27_tau-ref		THESE WERE ALL SUBMITTED ON 12.2
BOTHonGYP_26_tau-ref
BOTHonGYP_27_tau-ref


RERUN KALonKAL26 to see if it can retrieve tau/ref

'''

#_pull out desired experiment
experiment = 'HS3'
#segments   = experiments[experiment] 

# DIR_LOG  = os.path.expanduser('~/qsub_logs/')
DIR_LOG  = os.environ['LOG']
DIR_PROD = os.environ['PRODUCTS']
DIR_LBL  = 'LBL-RTM_simulated' #_generally want to use just LBL-RTM
DIR_TMP  = '/data/wsessions/TMP'
DEBUG    = 1

#_pick standard atmosphere to use, set to False for user defined
std_atmos = False # 'tropical'
model = {	False			: 0,	#_user defined prof 
			'tropical'		: 1,
			'midlat_summer'	: 2,
			'midlat_winter'	: 3,
			'subarc_summer'	: 4,
			'subarc_winter'	: 5,
			'us1976'		: 6 }[std_atmos]
dir_lblrtm	= '/data/wsessions/{0}/test_fov'.format(DIR_LBL) if not std_atmos \
			else '/data/wsessions/{0}/std_{1}'.format(DIR_LBL, std_atmos)

#_used instead of config_file
namelist = {
	#_if directory exists, skip
	'rerun'			: True,  

	#_set lblatm profile option
	'model'			: model,

	#_ok, see if this breaks everything
	'experiment'	: experiment,		#_used for LBL I/O naming 
#	'out_label'		: 'GYPonKAL_27',	#_for sim, name what is being retreived?
#	'out_label'		: 'GYPonKAL_26',
#	'out_label'		: 'KALonGYP_27',
#	'out_label'		: 'KALonGYP_26',

#	'out_label'		: 'GYPonGYP_27',	
#	'out_label'		: 'KALonKAL_27',	
#	'out_label'		: 'GYPonGYP_26',	
#	'out_label'		: 'KALonKAL_26_tau-ref',	
#	'out_label'		: 'KALonGYP_26_tau-ref',
#	'out_label'		: 'KALonGYP_27_tau-ref',
#	'out_label'		: 'GYPonKAL_26_tau-ref',
#	'out_label'		: 'GYPonKAL_27_tau-ref',

#	'out_label'		: 'BOTHonGYP_26_tau-ref',
#	'out_label'		: 'BOTHonGYP_27_tau-ref',
	'out_label'		: 'BOTHonKAL_26_tau-ref',
#	'out_label'		: 'BOTHonKAL_27_tau-ref',
#	'out_label'		: 'TESTING-TWO',

	#_these are separate so that uncertainty
	# can be looped via jv without attempts to retrieve
	# and covar matrix generated
	'state_vars'	: ['tau', 'ref'],	#_vars to retrieve
	'uncert_vars'	: ['z','thickness','surf_emis','surf_temp','rel_hum'],

	#_location of input radiances and output directory for fov plots
	'dir_shis'		: os.path.join(DIR_PROD, experiment.lower()), 
	'dir_out'		: os.path.join(DIR_PROD, experiment.lower()), 

	#_input sources (incomplete)
	'surf_temp'			: 'GDAS',	#_GDAS || LBL-ATM || LBL-ATM_old
	'profile_source'	: 'GDAS',	#_GDAS || SHIS || ECMWF || (only gdas works)
									# used for rel_hum jacobian
	'surf_temp_var'		: 'skt',	#_sst || skt, only used for ECMWF

	#_For testing of sensitivity to biases in apriori.
	# Dictionary should be popularted with with variable to arbitrarily
	# shift around the base apriori.
	'arbitrary_fudge'	: { 
		'surf_temp'	: linspace(-5, 5, 11),
		'z'			: linspace(-3, 3, 11),
		},

	#_directory of simulated run static case
	'dir_lblrtm'	: dir_lblrtm, 

	# 2. VARIABLES
	'dv'			: -26,			#_microwindow (window, no O3)
	'obs_angle'		: 0,			# 180 (zenith); 0 (nadir) 

	#_variables that are in JV but not SV are sources of uncertainty
	# that we are not attempting to back out

	#_SIMULATED
	############################################################################
	############################################################################
	#_label used for experiment (image output, flight file input for sim case)
	# fsim overwrites sim_experiment if both passed. If only sim_exp, fsim is
	# generated in methods as needed.
##	'sim_experiment'	: 'tau_ref',
#	'fsim'	: 'GYPSUM_SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.nc',
	'fsim'	: 'KAOLINITE_SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.nc',
##	'fsim'	: 'SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.nc',

	#_actual state in test (0.5, 0.5, 3, 1.0)
	#_a priori values for state vars
	'apriori' : [
					{	'type'		: 'aerosol', #_only use with simulated test
						'tau'		: 0.01, 
						'ref'		: 0.5,	
						'z'			: 3.0,
						'z_top'		: 3.01,
						'thickness' : 0.5,	
						'rel_hum'	: 1.0,		#_1 == apriori is GDAS
						'ssp_db' : 'ssp_db.mie_kaolinite.lognormal_sigma_0p699',
					#	'ssp_db' : 'ssp_db.shettle_dust.gamma_sigma_0p100',
					#	'ssp_db' : 'ssp_db.mie_gypsum.lognormal_sigma_0p699',
					#	'ssp_db' : 'ssp_db.mie_quartz.lognormal_sigma_0p699',
						'ref_wn'	: 900,	#_to match dummy_profiles.py	
						'surf_emis'	: [[100, 3000], [0.985, 0.985]],
						},
					{	'type'		: 'aerosol', #_only use with simulated test
						'tau'		: 0.01, 
						'ref'		: 0.5,	
						'z'			: 3.0,
						'z_top'		: 3.01,
						'thickness' : 0.5,	
						'ssp_db' : 'ssp_db.mie_gypsum.lognormal_sigma_0p699',
						'ref_wn'	: 900,	#_to match dummy_profiles.py	
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
						'instrument': 1,	#_0 for off, 1 for on
						'ssp_dist'	: 0.5,
						'model'		: 0,	#_LEAVE AT ZERO	
						},

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

	'dir_plot' : os.path.join(DIR_PROD, 'plots'),

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

#_sanity check
check_settings(**namelist)

#############################################################################_80
#_real/simulated_cases_#########################################################
################################################################################


def simulated_test_qsub(fsim=None, rerun=True, dir_plot='.', **kwargs):
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
	ppath = '{0}.{1}'.format(sim_experiment, out_label)
	ppath = os.path.join(dir_plot, ppath, 'FOV')
	kwargs.update({'dir_plot' : ppath})
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
			#_skip if being retrieved
			if variable in state_vars:
				continue

			for iii, layer in enumerate(kwgs['apriori']):
				#_else put in true value
				layer[variable] = value[i]
	
				if variable == 'z':
					layer['z_top'] = value[i]+.01
	
		#	#_else put in true value
		#	kwgs['apriori'][0][variable] = value[i]
	#
	#		if variable == 'z':
	#			kwgs['apriori'][0]['z_top'] = value[i]+.01
	
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
	header = [re.search('(.*?)-[\d.]+', t).group(1) \
		for t in true_vals[0].split(',')]
	fmt = []
	for j, v in enumerate(header):
		fmt.append('{{{0:d}:10.5f}}'.format(j))
	buf = len(header)	#_length of true vars (probably tau, ref, z == 3)
	header.extend(state_vars)	#_get known state, then add retrieved
	header.extend(['unc_{0}'.format(v) for v in state_vars])

	#_this is supposed to get the real values up front of the lines...
	lines = []
	for truth in true_vals:
		values = [re.search('.*?-([\d.]+)', t).group(1) \
			 for t in truth.split(',')]
		lines.append([float('{0}'.format(v)) for v in values])

	#_output format for results
	nv = len(state_vars)


	#_collect output
	line_out = [None] * flight.size
	q = 0
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
			input.merge_layers()	#_put same habit together
			x_std = load(open(pname, 'rb'))

			#_initialize an order to layers to make sure it's consistent
			if i == 0:
				layer_order_static = [s.split('/')[-1] for s in input.ssp_db]	
				nlayer = len(input.layers)
				#_later, use whatever the CURRENT layer order is to be keyed
				q = buf 
				for jj in range(nlayer):
				  for j, v in enumerate(state_vars):	#_uncertainty format (?)
					for k, v in enumerate(state_vars):  #_meeh
					#  fmt.append('{{{0:d}:10.5f}}'.format(buf + k + j*nv +jj*nv))
					  fmt.append('{{{0:d}:10.5f}}'.format(q))
					  q += 1

				fmt = ','.join(fmt) 
				fmt += '\n'
			
			#_initialize received tuple tau_00 ref_00 unc_t00 unc_r00
			tup_retr = [None] * (nv*2) * nlayer 

			#_pull out retrieved layer
			for j, x in enumerate(input.layers):
				
				#_format retrieved values for output
				for ssp_db in layer_order_static:
	
					#_pull out index within static ssp order	
					dbnum = x['dbnum']	
					current_ssp = input.ssp_db[dbnum].split('/')[-1]
					idx = layer_order_static.index(current_ssp)

					#_loop over state vars, put into tup_retr
					for ii, v in enumerate(state_vars):
						layer_idx = idx*nv*2
		
						#_add retrieval values to tuple
						if v == 'tau':
							tup_retr[layer_idx + ii] = x[v][0]
						else:	
							tup_retr[layer_idx + ii] = x[v]
	
						#_add uncertainty values to tuple
						tup_retr[layer_idx + nv + ii] = x_std[ssp_db][v] 

			#	tup_retr = [x[v] for v in state_vars]			#_ret values
			#	tup_retr.extend([x_std[i][v] for v in state_vars])	#_unc values

		#	x = input.layers[0]		#_pull out first layer (dict = {'tau' 'ref'}
		#	x['tau'] = x['tau'][0]
		#	x_std = load(open(pname, 'rb'))

		#	#_format retrieved values for output
		#	tup_retr = [x[v] for v in state_vars]			#_retrieved values
		#	tup_retr.extend([x_std[v] for v in state_vars]) #_uncert values
		#	lines[i].extend(tup_retr)

		except IOError:
			dbg('WARNING: missing {0}'.format(fname))
			tup_retr = [-9999 for v in state_vars]
			tup_retr.extend([-9999 for v in state_vars])
			lines[i].extend(tup_retr)
		
	#	if len(input.layers) != 1:
	#		raise RuntimeError, 'not yet implemented for more than one layer'
		lines[i].extend(tup_retr)
		line_out[i] = fmt.format(*lines[i])

	print fmt

	#_if there is more than one layer, extend header
	for i in range(nlayer - 1):
		header.extend(state_vars)	#_get known state, then add retrieved
		header.extend(['unc_{0}'.format(v) for v in state_vars])
	
	#_add ssp file order to header
	try:
		[header.append(' ' + ssp) for ssp in input.ssp_db]
	except:
		pass

	#_write output as we go
	dbg(file_out)

	#_initialize output (fmt produces header with state variable names)
	with open(file_out, 'w') as f:
		#_write header
		fmt =','.join(['{{{0:d}:>10s}}'.format(i) for i,v in enumerate(header)])
		f.write(fmt.format(*header) + '\n')

		#_write retrieved and true values
		[f.write(l) for l in line_out] 


def simulated_sensitivity_qsub(fsim=None, fidx=655, rerun=True, dir_plot='.',
	aods=[.24,.83,1.42,2.0], arbitrary_fudge={}, ref=1.1, z=3., **kwargs):
	'''
	fov		int,	location for desired base case.
					655 = tau==1.03, ref==1.1, z==3

	Loop over arbitrary biases in various apriori variables and get
	to find the change in optical depth in the retrieval as a result
	'''
	import matplotlib.pyplot as plt
	from hs3_utils import Flight_segment as F
	from numpy import diag, matrix, array, arange, append
	from lblrtm_utils import microwindow_average, microwindow_stdev
	from time import sleep
	from libtools import shrink_ticks, setup_groups, mkdir_p
	from os import fork, waitpid
	from multiprocessing import Process, Pipe
	from pickle import dump
	from libtools import combination

	#_build out label base, add values in loop
	order = arbitrary_fudge.keys()
	out_base = '_'.join(['SHIFT', '-'.join(['tau'] + order)])

	#_announce
	pid = os.getpid()
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
	ppath = '{0}.{1}'.format(sim_experiment, out_base)
	ppath = os.path.join(dir_plot, ppath, 'FOV')
	kwargs.update({'dir_plot' : ppath})
	mkdir_p(ppath)

	#_build potential combinations for variables	
	vars = [arbitrary_fudge[v] for v in order]
	combos = combination(vars)

	#_for easier transition from copied code
	for tau in aods:
	  #_find matching index
	  ridx = array(true_values['ref']) == ref
	  tidx = array(true_values['tau']) == tau 
	  zidx = array(true_values['z']) == z 
	  fidx = arange(flight.size)[ridx*tidx*zidx][0]	
	  fov = flight[fidx]
	  for vals in combos:
		#_add values to end of label
		out_vals = ['{0:05.2f}'.format(tau)]
		[out_vals.append('{0:05.2f}'.format(v)) for v in vals]
		out_label = out_base + '_' + '-'.join(out_vals)
	
		#_initialize plotting, for first ten instances only
		kwargs.update({	'sup'	: meh[fidx] + ' ' + out_label,
						'comment' : meh[fidx] + ' ' + out_label	})

		#_name of final output file
		lbldis_output = 'lbldis_output.SIMULATED.' + \
			'{0:04d}.{1}.{2}.final.cdf'.format(fidx, sim_experiment, out_label)
		lbldis_clear = 'lbldis_output.SIMULATED.' + \
			'{0:04d}.{1}.{2}.clear.cdf'.format(fidx, sim_experiment, out_label)
		posterior_pickle = 'posterior.SIMULATED.' + \
			'{0:04d}.{1}.{2}.pk'.format(fidx, sim_experiment, out_label)

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
			#_skip ones wen're actually retrieving.
			if variable in state_vars:
				continue

			#_else put in true value
			kwgs['apriori'][0][variable] = value[fidx]
	
			if variable == 'z':
				kwgs['apriori'][0]['z_top'] = value[fidx] + .01
	
		#_change apriori values for test
		kwgs['apriori'][0]['tau'] = tau
		for i, value in enumerate(vals):
			if variable in state_vars:
				raise RuntimeError, "What are you doing?"
	
			varname = order[i]

			#_shift var
			if varname == 'surf_temp':
				kwgs.update({'shift_sfct' : value })
			else:
				kwgs['apriori'][0][varname] += value

			#_adjust layer top if variable is geometric height
			if order[i] == 'z':
				kwgs['apriori'][0]['z_top'] = kwgs['apriori'][0]['z'] + 0.01 

		#_write kwarg pickle file
		fkwgs = 'kwargs.f{0:04d}.{1}.{2}.{3}.pk'.format(i, sim_experiment,
			out_label, pid)
		fkwgs = os.path.join(DIR_TMP, fkwgs)
		dump(kwgs, open(fkwgs, 'wb'))
		#_submit job
		dbg(('TESTING FOV NUMBER', i, out_label))
		cmd = ' '.join((qsub, scpt, str(i), fkwgs))
		dbg(cmd); dbg('')
		os.system(cmd)


def read_sensitivity_qsub(sv='tau', fsim=None, ref=1.1, z=3., 
	dir_plot='.', aods=[.24,.83,1.42,2.], arbitrary_fudge={}, **kwargs):
	'''
	fov		int,	location for desired base case.
					655 = tau==1.03, ref==1.1, z==3

	Loop over arbitrary biases in various apriori variables and get
	to find the change in optical depth in the retrieval as a result

	Read in files generated by simulated_s_q
	'''
	import matplotlib.pyplot as plt
	from lblrtm_utils import Lbldis_input
	from hs3_utils import Flight_segment as F
	from numpy import diag, matrix, recarray, array, append, arange
	from lblrtm_utils import microwindow_average, microwindow_stdev
	from time import sleep
	from libtools import shrink_ticks, setup_groups, mkdir_p
	from os import fork, waitpid
	from multiprocessing import Process, Pipe
	from pickle import dump
	from libtools import combination

	#_build out label base, add values in loop
#	order = arbitrary_fudge.keys()
#	out_base = '_'.join(['SHIFT', '-'.join(order)])
	order = arbitrary_fudge.keys()
	out_base = '_'.join(['SHIFT', '-'.join(['tau'] + order)])

	#_announce
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

	#_create directory for plots
	ppath = '{0}.{1}'.format(sim_experiment, out_base)
	ppath = os.path.join(dir_plot, ppath, 'FOV')
	kwargs.update({'dir_plot' : ppath})
	mkdir_p(ppath)

	#_build potential combinations for variables	
	vars = [arbitrary_fudge[v] for v in order]
	combos = combination(vars)
	ntau = len(aods)

	#_intialize a recarray
	dtype = [(varname, 'f8') for varname in order]
	dtype.append((sv, 'f8'))
	dtype.append(('retr_{0}'.format(sv), 'f8'))
	data = recarray((len(combos)*ntau), dtype=dtype)
	
	fmt = 'lbldis_input.SIMULATED.{0:04d}.tau_ref_z.ERR_total.final'
	dir_input = os.path.join(os.environ['PRODUCTS'], 'LBL-RTM_simulated', 
			'UNCERTAINTY_TESTS_2015.09.29', 'inputs')

	#_for easier transition from copied code
	k =0
	for tau in aods:
	  #_find matching index
	  ridx = array(true_values['ref']) == ref
	  tidx = array(true_values['tau']) == tau 
	  zidx = array(true_values['z']) == z 
	  fidx = arange(flight.size)[ridx*tidx*zidx][0]	
	  for i, vals in enumerate(combos):

		#_add values to end of label
		out_vals = ['{0:05.2f}'.format(tau)]
		[out_vals.append('{0:05.2f}'.format(v)) for v in vals]
		out_label = out_base + '_' + '-'.join(out_vals)

		#_read in original thing
		if vals == [0., 0.]:
			  #_read in non-shifted cases
			  fname = fmt.format(fidx)
			  fname = os.path.join(dir_input, fname)
	
		else:
			#_name of final output file
			lbldis_input = 'lbldis_input.SIMULATED.' + \
				'{0:04d}.{1}.{2}.final'.format(fidx, sim_experiment, out_label)

			#_read in shifted final input
			fname = os.path.join(kwargs.get('dir_lblrtm'), lbldis_input)

		if not os.path.exists(fname):
			raise RuntimeError, 'TESTING'
		input = Lbldis_input(fname)

		#_add current shift values to recarray
		[setattr(data[k], varname, vals[j]) for j, varname in enumerate(order)]

		#_add current aod value to recarray
		setattr(data[k], 'tau', tau)
		setattr(data[k], 'retr_tau', input.layers[0]['tau'][0]) 
		k += 1
##
##	  #_add current shift values to recarray
##	  [setattr(data[k], varname, 0) for j, varname in enumerate(order)]
##
##	  #_add current aod value to recarray
##	  setattr(data[k], 'tau', tau)
##	  setattr(data[k], 'retr_tau', input.layers[0]['tau'][0]) 
	  
	return data


################################################################################
#_plotting_#####################################################################
################################################################################

	
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


def plot_simulated_histogram(out_label='', plot_var='tau', dir_plot='.',
	dir_hs3=os.path.join(os.environ['PRODUCTS'],'hs3'), files=[], **kwargs):
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
	dbg(("TESTING"))
	dbg(files)
	if len(files) == 0:
		fname = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.retrieved'.format(out_label)
		fname = os.path.join(dir_hs3, fname)
		try:
			truth, retr, uncrt = read_simulated_retrieval(fname, **kwargs)
		except IOError:
			dbg(('Did you complete this retrieval yet?', out_label))
			os._exit(23)
	else:
		for i, f in enumerate(files):
			try:
				if i == 0:	#_initialize
					truth, retr, uncrt = read_simulated_retrieval(f, **kwargs)
				else:		#_append
					t, r, u = read_simulated_retrieval(f, **kwargs)
					[truth[sv].append(t[sv]) for sv in truth]
					[retr[sv].append(r[sv]) for sv in retr]
					[uncrt[sv].append(u[sv]) for sv in uncrt]
			except IOError:
				dbg(('Did you complete this retrieval yet?', f))
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
	ppath = '{0}'.format(out_label)
	ppath = os.path.join(dir_plot, ppath)
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
	fsim=None, dir_hs3=os.path.join(os.environ['PRODUCTS'],'hs3'),
	dir_plot='.', files=[], **kwargs):
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

	if len(files) == 0:
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
	else:
		for i, f in enumerate(files):
			res = re.search('SIMULATED.(.*?).nc', f)
			sim_experiment = res.group(1)
			f = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.{1}.retrieved'.format(
				sim_experiment, out_label)
			f = os.path.join(dir_hs3, f)
			try:
				if i == 0:	#_initialize
					truth, retr, uncrt = read_simulated_retrieval(f, **kwargs)
				else:		#_append
					t, r, u = read_simulated_retrieval(f, **kwargs)
					for sv in truth:
						truth[sv] = append(truth[sv], t[sv]) 
					for sv in retr:
						retr[sv] = append(retr[sv], r[sv]) 
					for sv in uncrt:
						uncrt[sv] = append(uncrt[sv], u[sv]) 
			except IOError:
				dbg(('Did you complete this retrieval yet?', f))
				os._exit(23)

		res = re.search('SIMULATED.(.*?).nc', fsim)
		sim_experiment = res.group(1).split('.')[0]
		out_label = 'summary'

	#_what are we working with?
	state_vars = retr.keys() # truth.keys()
	dbg('');dbg(state_vars)

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
	ppath = '{1}.{0}'.format(out_label, sim_experiment)
	ppath = os.path.join(dir_plot, ppath)
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
				fig.suptitle('Fractional Uncertainty')

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

			#_pull out values across all over vars, at truth value
			for value_truth in set(truth[sv]):
				x_truth.append(value_truth)
				y_fract.append(fract[truth[sv][idx] == value_truth])	

			#_pull out number of retrievals going into these
			y_bar = [len(y) for y in y_fract]

			##	axis.scatter(truth[sv][idx], fract, marker='x', color=col)
			axis.boxplot(y_fract, positions=x_truth, widths=0.05)
			axis.set_xticks(list(set(truth[sv])))

			abar = axis.twinx()
			bar_width = (max(x_truth) - min(x_truth)) / float(len(x_truth))\
				 - 0.05
			
			abar.bar(x_truth-bar_width/2., y_bar, bar_width, 
				alpha=0.24, color='brown')	
			#_get how high we can go and still be 95 % certain of the AOD
		##	m95 = truth[sv][certainty_garbage(res)]
			
			#_various things for uncertainty
		#	axis.text(truth[sv].max()-0.2, 1.6, 'n={0}'.format(len(y_fract))
			axis.set_ylim(bottom=0, top=2)
			axis.set_xlim(truth[sv].min()-0.1, truth[sv].max()+0.1)
			axis.set_ylabel('{0} = {1}'.format(var, value), size='xx-small')
			axis.set_ylabel('{0} = {1}'.format(var, value), size='x-small')
##			axis.set_xlabel('{0}, max_95 {1:5.2f}'.format(sv,m95), size='x-small')
			axis.grid(False)
			
			#_shrink...ticks...
			shrink_ticks(axis)
			shrink_ticks(abar)

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
	dir_hs3=os.path.join(os.environ['PRODUCTS'],'hs3'), dir_plot='.', **kwargs):
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

	raise RuntimeError, 'Not sure what this was doing, but needs updating.'

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
	ppath = '{1}.{0}'.format(out_label, sim_experiment)
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
	##	m95 = truth[sv][certainty_garbage(res)]
		
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


def read_simulated_retrieval(fname, **kwargs):
	'''
	Read in output file produced by optimal estimation test
	on simulated data. File should be space deliminated and 
	contain true values for state variables, retrieved values,
	and the uncertainty (unlabled columns in same order as sv)

	fname	str,	full path to file
	'''
	from numpy import array
	import re
	
	dbg(('opening', fname))
	with open(fname, 'r') as f:
		#_get rid of extraneous whitespace
		lines = [l.strip() for l in f.readlines()]

		#_first find out many ssp files are listed at the end of the header
		header_vars = [l.strip() for l in lines[0].split(',')]
		ssps = []
		for v in header_vars:
			if re.search('/', v):
				ssps.append(v)

		nlayer = len(ssps)

		#_pull out state variables
		state_vars = []
		retr_vars = []
		for v in header_vars[:-nlayer]:
			if v not in state_vars and not re.search('unc_', v):
				state_vars.append(v)
		nsv = len(state_vars)

		for v in header_vars[nsv:]:
			if v not in retr_vars and not re.search('unc_', v) \
				and not re.search('/', v):
				retr_vars.append(v)

		state_vars = list(set(state_vars))
		nrv = len(retr_vars) 
		del lines[0]

		#_pull out data
		truth, retrieved, uncertainty = {}, {}, {}
		[truth.update({sv:[]}) for sv in state_vars]
		for ssp in ssps:
			retrieved.update({ ssp : {} })
			uncertainty.update({ ssp : {} })

			for sv in retr_vars:
				#_retrieved[<layer_#>][<state_var>] = array
				retrieved[ssp].update({sv:[]})
				uncertainty[ssp].update({sv:[]})

		#_loop over lines and fill arrays
		for line in lines:
			cols = line.split(',')

			#_loop over state var columns
			for i, variable in enumerate(state_vars):
				truth[variable].append(float(cols[i]))

			#_loop over retrieve var columns (values and uncertainties)
			for i, variable in enumerate(retr_vars):
				ibuf = i + nsv
				for ii, ssp in enumerate(ssps):
					retrieved[ssp][variable].append(float(cols[ibuf + ii*(nrv*2)]))
					uncertainty[ssp][variable].append(float(cols[ibuf + nrv + ii*(nrv*2)]))

	for sv in state_vars:
		truth[sv] = array(truth[sv])

	for rv in retr_vars:
		for ssp in ssps: 
			retrieved[ssp][rv] = array(retrieved[ssp][rv])
			uncertainty[ssp][rv] = array(uncertainty[ssp][rv])

	return truth, retrieved, uncertainty


################################################################################
#_tools_########################################################################
################################################################################


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


def generate_delta_sfc_temp(delta=10, dt=0.25, dir_lblrtm='.', **namelist):
	''' stupid way to copy test_fov diretory for sfc temp tests '''
	from numpy import abs, arange
	from shutil import copytree as copy

	bot = delta/2  - delta
	top = bot + delta + dt 
	temps = arange(bot, top, dt)

	for t in temps:
		if t == 0: continue
		dir_name = '{0}_{1:4.2f}'.format(dir_lblrtm, t)
		dbg(dir_name) 
		copy(dir_lblrtm, dir_name)


def gen_outlabel(dv=None, apriori=None, dynamic_ATP=None, dynamic_CTP=None, 
	out_label=None, cld=False, **kwargs):
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

	try:
		atp = 'ATP-{0}'.format('-'.join(dynamic_ATP))
	except:
		atp = 'NOATP' 

	#_update out_label
	arg = (out_label, n_habit, abs(dv), cld, ctp, zkm, atp)
	return '{0}_M{1}_DV{2}_{3}_{4}_{6}_{5}'.format(*arg)


def check_missing(dir_lblrtm=None, fsim=None, dir_shis=None, 
	out_label=None, **kwargs):
	from hs3_utils import Flight_segment 
	import re

	file_seg = os.path.join(dir_shis, fsim)
	flight = Flight_segment(file_seg=file_seg)
	fidx = range(flight.size)

	#_pull out thing defining what the simulated file is varied across
	svars = re.search('SIMULATED.(.*?).nc$', fsim).group(1)

	#_get 
	total = 0
	tmp_outlabel = '{0}.{1}'.format(svars, out_label)

	#_loop over fidx to get
	missing_count = 0
	for i in fidx:
		#_get output file name
		args = (i, tmp_outlabel)
		lbldis_output = 'lbldis_output.SIMULATED.{0:04d}.{1}.final.cdf'.format(*args)
			
		#_build file path
		lbldis = os.path.join(dir_lblrtm, lbldis_output)
	
		#_check if missing
		if not os.path.exists(lbldis):
			missing_count += 1	

		#_wait what	
		total += 1

	#_report
	dbg('');dbg('{2} missing {0:>6d} / {1:>6d}'.format(missing_count, len(fidx), tmp_outlabel))
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
	#_FOR SINGLE, (specify fsim) 
	if 0:	#_qsub loop over namelist-defined file (MUST ALSO RUN WRITE_SIM)
		simulated_test_qsub(**namelist)
		exit()
	if 1:	#_write output from single run 
	  for a in ['BOTHonKAL_26_tau-ref',
		'BOTHonKAL_27_tau-ref',
		'BOTHonGYP_26_tau-ref',
		'BOTHonGYP_27_tau-ref']:
		namelist['out_label'] = a
		write_simulated_qsub(**namelist)
	  exit()
	if 0: #_SEE WHAT'S COMPLETE
	  dbg(namelist['out_label'])
	  check_missing(**namelist)		
	  exit()
	if 0: #_plot statistics from runs
		for sv in namelist['uncert_vars']:
			try:
				plot_simulated_retrieval_by_var(sv, **namelist)
			except KeyError:
				dbg('crash out of non-truth set variables: {0}'.format(sv))

		#_loop over all uncertainty vars and plot AOT variance by them
		plot_simulated_retrieval(**namelist)
		exit()

	############################################################################
	#_FOR EVERYTHING.  Ignores fsim above, USE FOR UNCERTAINTY TESTING, NOTHING
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

				################################################################
				if True:	#_for 'total', ie. full retrieval 
					if i != 0: #_only do total once
						continue
					label = 'ERR_{0}'.format('total')
					label = label.replace('.','_')
					for vname in vnames: 
						namelist['uncertainty'][vname] = 1 

				else:	#_for individual
					label = 'ERR_{0}'.format(current_var)
					for vname in vnames: 
						namelist['uncertainty'][vname] = int(vname==current_var)
				################################################################

				#_finalize label
				if label[-1] == '_':
					label = label[:-1]
				namelist['out_label'] = label
				namelist['fsim'] = fname
	
				#_never have both launch together. stq must complete before wsq
				dbg(label)
			##	check_missing(**namelist)			#_see which are done 
			##	simulated_test_qsub(**namelist)		#_perform pseudo-obs retr
			##	write_simulated_qsub(**namelist)	#_after stq, write .retr
		
			##	for sv in namelist['uncert_vars']:	#_plot everything up
			##		try:
			##			plot_simulated_retrieval_by_var(sv, **namelist)
			##		except KeyError:
			##			dbg('crash out of non-truth set vars: {0}'.format(sv))

		#_do all for particular type
		for sv in namelist['uncert_vars']:
			namelist['out_label'] = 'ERR_total'
			try:
				plot_simulated_retrieval_by_var(sv, files=fnames, **namelist)
			except KeyError:
				dbg('crash out of non-truth set variables: {0}'.format(sv))
		exit()


	#_PLOTTING ONLY
	if 0:	#_WITH RETRIEVALS
	  print "NOTHING IN SEGMENTS FOR SIMULATED CASES. HALT"
	  raise RuntimeError, '23'
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


	#_GENERATE SURF TEMP TEST CASES
	from numpy import arange
	dir_lblrtm_base = namelist['dir_lblrtm'] 
	if 0:
##	for t in arange(-5, 5.25, 1):
##	for t in arange(-4.5, 5.25, .5):
#	for t in arange(-4.75, 5.25, .25):
#	for t in [-5]:
	  if t == 0:
		continue #_skip normal situation

	  dir_lblrtm = '{0}_{1:4.2f}'.format(dir_lblrtm_base, t)
	  namelist.update({
			'dir_lblrtm' : dir_lblrtm,
			'out_label'  : 'SFC_TEMP_{0:4.2f}'.format(t),
			'shift_sfct' : t,
			})
	  if 0:	#_make copies of test directory
		generate_delta_sfc_temp(**namelist)
	  if 0:	#_qsub loop over namelist-defined file
		fnames = ['SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.nc',
				'SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.NOISE.0004485.nc', 
				'SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.NOISE.0021783.nc',
				'SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.NOISE.0030252.nc']
		for fsim in fnames:
			namelist.update({'fsim' : fsim})
		#	simulated_test_qsub(**namelist)
			write_simulated_qsub(**namelist)
		#	check_missing(**namelist)		
	  if 0: #_plot statistics from runs
		for sv in namelist['uncert_vars']:
			try:
				plot_simulated_retrieval_by_var(sv, **namelist)
			except KeyError:
				dbg('crash out of non-truth set variables: {0}'.format(sv))
		#_loop over all uncertainty vars and plot AOT variance by them
	##	plot_simulated_retrieval(**namelist)
		exit()


	#_Sensitivity matrices
	if 0:
		#_setup output label for testing

		#_submit jobs
#		simulated_sensitivity_qsub(**namelist) 
		read_sensitivity_qsub(**namelist)
		#write_sensitivity_qsub(**namelist)
		#

