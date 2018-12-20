#!/usr/bin/env python
# 2014.12.07
################################################################################
# Produce simulated S-HIS data to retrieve with OE                             #
################################################################################
from numpy import arange, array, linspace
import os
import sys

kwargs = {
	#_number of processors to use
	'nproc' : 20,

##	'vars'	: ['tau', 'ref', 'surf_temp'],
##	'vars'	: ['tau', 'ref', 'z'],
##	'vars'	: ['tau', 'ref'],
	'vars'	: ['tau'],

	'ranges': {
		'tau'		: linspace(0.05, 2, 40), #_4)
		'ref'		: linspace(0.5, 1, 21),	#_4
		'z'			: linspace(1, 4, 4), 
		'thickness'	: linspace(0.5, 3., 6)
##		'surf_temp'	: arange(-5, 6),
##		'surf_temp'	: linspace(-286.8, 296.8, 11),
		},

	#_add noise?
	'noisey'	: False,

	#_some constant values
	'clddef'	: {
		'z'		: 3,	#_keep stable for now
		'ref'	: 0.5,	#_same
		'dbnum'	: 0,
		'ref_wn': 900,
		},

	#_formatting
	'fmt'	: {
		'tau'		: '10.5f',
		'ref'		: '10.5f',
		'z'			: '10.5f',
		'surf_temp'	: '10.5f',
		},

	#_lblrtm input directory
	# going to need to rerun with higher vertical resolution (500 m?)
	'dir_lblrtm'	: '/data/wsessions/LBL-RTM_simulated/std_tropical',

	#_used... for... HBB?
	'real_shis'		: 'SHIS_rdr20130916T210344end20130917T020353' + 
						'sdr20140110T030108_rad_pcfilt.nc',
	
	#_name of raw forward model output
	'fname_lbl'		: 'lblrtm.SIMULATED.cdf',

	#_location of HS3 output
	'dir_hs3'		: os.path.join(os.environ['PRODUCTS'], 'hs3'),
	'dir_shis'		: os.path.join(os.environ['PRODUCTS'], 'hs3'),

	}

#_generate fname_shis from settings
a = ('_'.join(kwargs['vars']), '.NOISE'*kwargs['noisey'])
f = 'SHIS.CPL.COLLOC.SIMULATED.{0}{1}.hdf'.format(*a)
kwargs.update({'fname_shis' : f}) 


#############################################################################_80
#############################################################################_80
#############################################################################_80


def run_main(fname_shis='simulated.cdf', fmt={'tau':'10.5f'}, **kwargs):
	'''
	atmospheric profile is arbitrary from gdas/lblrtm, rest is garbage
	'''
	from hs3_utils import write_collocated
	from libtools import combination as combo
	from libtools import dtg2epoch as d2e

	#_temporarily pull out for size since profile has to come first
	vars, ranges = kwargs.get('vars'), kwargs.get('ranges')
	combos = combo([ranges[v] for v in vars])
	size = len(combos)

	#_create name of file to store real values
	fname_txt = fname_shis.replace('.hdf', '.txt').replace('COLLOC',
				 'GDAS.COLLOC')
	fname_txt = os.path.join(kwargs.get('dir_hs3'), fname_txt)

	#_write combos in order to output
	fmts = ['{{{0:d}:{1:s}}}'.format(i,key) for i, key in enumerate(fmt)]
	fmt = ''.join(fmts)
	with open(fname_txt, 'w') as f:
		fmt0 =''.join(['{{{0:d}:>8s}}'.format(i) for i,v in enumerate(vars)])
		fmt1 =''.join(['{{{0:d}:8.5f}}'.format(i) for i,v in enumerate(vars)])
		f.write(fmt0.format(*vars) + '\n')
		[f.write(fmt1.format(*a)+'\n') for a in combos]
		print fname_txt, 'written...'

	#_use single profile
	prof = fake_prof(size, **kwargs)

	#_make fake epochs to to maintain order 
	epoch_start = d2e('2014010100')

	#_generate fake S-HIS data
	shis, notes = fake_shis_sge(latitude=prof.latitude[0], 
		longitude=prof.longitude[0], **kwargs)

	#_add individual seconds to fake time data
	epoch = arange(shis.size) + epoch_start
	shis.epoch[:] = epoch

	#_read in fov profile from GDAS (PUT IN METADATA)
	cpl	= fake_cpl(size, lat=prof.latitude[0], lon=prof.longitude, **kwargs)

	#_write to file
	write_collocated(shis, cpl, prof, fname_shis, notes=notes, **kwargs)
	

################################################################################
################################################################################
################################################################################


def fake_cpl(size=1, lat=0, lon=-45., **kwargs):
	''' generate fake cpl object '''

	#_pass fake SHIS profile, and then expand a blank recarray from there
	from numpy import ndarray, recarray, zeros
	from libtools import julian2epoch as j2e

	#_read in cpl data
	dtype = [   ('ext_532', ndarray),
	            ('tau_532', ndarray),
	            ('tau_type', ndarray),
	            ('layertop', ndarray),
	            ('layerbot', ndarray),
	            ('epoch', 'f8'),
	            ('latitude','f8'),
	            ('longitude','f8')]
	cpl = recarray((size,), dtype=dtype)
	cpl.__setattr__('fname', 'SIMULATED')

	#_calculate epoch from Year and Julian float
	year = 1982
	jday = 54 
	time = j2e(year, jday)

	#_add data to recarray
	cpl.epoch[:]        = time
	cpl.latitude[:]     = lat
	cpl.longitude[:]    = lon 
	for t in xrange(size):
	    cpl.ext_532[t]  = zeros((10)) 
	    cpl.tau_532[t]  = zeros((10))
	    cpl.tau_type[t] = zeros((10))

	return cpl


def fake_prof(size=1, nz=47, **kwargs):
	''' generate fake profile '''
	from hs3_utils import read_profile_gdas
	from numpy import ndarray, recarray, linspace, zeros
	
	dtype_dict = {	'latitude' : 'f8', 'longitude' : 'f8', 'epoch' :'f8',
					'sfc_temperature' : 'f4', 'ozone_mixing_ratio' : ndarray,
					'relative_humidity' : ndarray, 'geopotential_height' : ndarray,
					'temperature' : ndarray }
	dtype = [(n, t) for n, t in dtype_dict.iteritems()]
	
	names = ('latitude', 'longitude', 'epoch', 'sfc_temperature',
            'ozone_mixing_ratio', 'relative_humidity', 'geopotential_height',
            'temperature')

	#_initialize fake profile recarray
	prof = recarray((size,), dtype=dtype)

	#_read in a fake clear sky profile, copy it over #size times
	for attr in names:
		if dtype_dict[attr] == ndarray:
			for i in range(size):
				prof.__getattribute__(attr)[i] = zeros((nz))
		else:
			prof.__getattribute__(attr)[:] = 0.

	#_add pressure
	setattr(prof, 'pressure', linspace(1000., 1, nz))
	setattr(prof, 'fname', 'SIMULATED')

	return prof


def fake_prof_old(size=1, idx=[2],
	real_shis='SHIS_rdr20130916T210344end20130917T020353' + 
				'sdr20130917T131448_rad.nc', **kwargs):
	''' generate fake profile '''
	from hs3_utils import read_profile_gdas
	from numpy import array, recarray

	#_read in real profile to copy
	real_prof = read_profile_gdas(real_shis, array(idx), dummy=True)

	#_initialize fake profile recarray
	prof = recarray((size,), dtype=real_prof.dtype)

	#_read in a fake clear sky profile, copy it over #size times
	for attr in prof.dtype.names:
		prof.__getattribute__(attr)[:] = real_prof.__getattribute__(attr)

	#_add pressure
	setattr(prof, 'pressure', real_prof.pressure)
	setattr(prof, 'fname', 'SIMULATED')

	return prof


def fake_shis_sge(vars=['tau'], ranges={'tau':[0,]}, dir_lblrtm='.', dz=0.2,
	latitude=0, longitude=0, clddef={}, noisey=False, real_shis=None, **kwargs):
	''' generate fake shis retrieval '''
	from numpy import recarray, array, ndarray, tile
	from libtools import combination as combo
	from lblrtm_utils import write_lbldis_parameters as lbldis_write
	from subprocess import call
	from shis import convolve2shis as c2shis
	from scipy.io import netcdf
	from scipy.interpolate import splev, splrep, interp1d
	from hs3_utils import read_shis
	from numpy.random import random
	from multiprocessing import Process, Queue, Pipe
	from libtools import setup_groups
	from time import sleep

	#_read in real thing for HBB_NESR
	dir_hs3 = kwargs.get('dir_hs3')
##	ugh = '{0}/SHIS_rdr20130816T033458end20130816T080506sdr20130817T165537_rad_pcfilt.nc'.format(dir_hs3)
##	real_shis = read_shis(ugh, idx=array([0])) 
	real_shis = read_shis(real_shis, idx=array([0])) 
	hbb = None

    #_lbldis gets persnickety about this
	if 'LB_LIBRARY_PATH' not in os.environ:
		os.environ['LD_LIBRARY_PATH'] = '/opt/netcdf-4.1.2/lib/'

	#_run lbl-rtm with parameters
	combos = combo([ranges[v] for v in vars])

	#_initialize recarray to return
	dtype = [   ('radiances', ndarray),
                ('relative_humidity', ndarray),
                ('temperature', ndarray),
                ('sfc_temperature', 'f4'),
                ('ozone_mixing_ratio', ndarray),
                ('epoch', 'f8'),
                ('latitude','f8'),
                ('longitude','f8'),
                ('altitude', 'f4')  ]
	shis = recarray((len(combos)), dtype)
	shis.__setattr__('fname', 'SIMULATED')

	#_some fake profiles to shove into SHIS
	fk_tmp = tile(-9999, 101)
	fk_rel = tile(-9999, 101)
	fk_ozn = tile(-9999, 101)

	#_add fake float values to some
	shis.sfc_temperature[:] = -9999. 
	shis.epoch [:] = -9999
	shis.latitude[:] = latitude 
	shis.longitude[:] = longitude 
	shis.altitude[:] = 18.

	#_initialize list for metadat
	notes = []

	#_split combinations in processing groups
	groups = setup_groups(combos, **kwargs)

	#_name of script to launch
	dir_log = '/data/wsessions/logs'
	env = ','.join(['='.join((var, os.environ[var])) \
        for var in ['PYTHONPATH', 'PRODUCTS', 'PATH', 'LOG']])
	qsub = ' '.join((   'qsub -v', env, '-o', dir_log, '-e', dir_log,
                        '-cwd -S /opt/ShellB3/bin/python'))
	srpt = os.path.expanduser('~/lib/sge_std_dummy_profiles.py')

	#_ THIS COULD BE FORKED
	#_build clddef dict
	for j, opts in enumerate(combos):
		#_get pid
		pid = os.getpid()

		#_launch lbldis
		cloud = [clddef.copy()]
		[cloud[0].update({v : opts[q]}) for q, v in enumerate(vars)]
		cloud[0]['z_top'] = cloud[0]['z'] + dz

		#_build prefix for labeling
		fmt	= ['{0:s}-{1:05.2f}'.format(v, opts[q]) for q,v in enumerate(vars)]
		notes.append(','.join(fmt))
		fname = '{0:s}/{1:s}.dummy.{2:07d}.namelist'.format(dir_lblrtm,
														'_'.join(fmt), pid)
		oname = '{0:s}/{1:s}.dummy.{2:07d}.output'.format(dir_lblrtm,
														'_'.join(fmt), pid)
		lname = '{0:s}/{1:s}.dummy.{2:07d}.lock'.format(dir_lblrtm,
														'_'.join(fmt), pid)
		#_lname == lockfile name
		with open(lname, 'wb') as f: pass

		#_write input paramater file
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kwargs)

		#_launch
		cmd = ' '.join((qsub, srpt, fname, oname, lname)) 
		os.system(cmd)

	for j, opts in enumerate(combos):
		#_generate lockfile name
		fmt	= ['{0:s}-{1:05.2f}'.format(v, opts[q]) for q,v in enumerate(vars)]
		notes.append(','.join(fmt))
		oname = '{0:s}/{1:s}.dummy.{2:07d}.output'.format(dir_lblrtm,
				'_'.join(fmt), pid)
		lname = '{0:s}/{1:s}.dummy.{2:07d}.lock'.format(dir_lblrtm,
				'_'.join(fmt), pid)

		#_loop until ready
		while os.path.exists(lname):
			print 'Waiting on {0}'.format(lname)
			sleep(10)
			continue

		#_get thread
		rads, wave = read_lbldis_fork(oname)

		#_convolve to shis profile
		rads = rads / 1000.
		wn, rd = c2shis(wave, rads, **kwargs)
		rd *= 1000.

		#_put into recarray
		shis[j].radiances			= rd 
		shis[j].relative_humidity	= fk_rel 
		shis[j].temperature			= fk_tmp 
		shis[j].ozone_mixing_ratio	= fk_ozn 

	#_interpolate blackbody reference
	f1d = splrep(real_shis.wavenumber, real_shis.hbb_nesr, k=3)
	hbb = splev(wn, f1d)

	#_if testing with noise, generate
	if noisey:
		from random import gauss
		for j in xrange(shis.size):

			noise = array([gauss(0, sigma) for sigma in hbb]) 
			shis[j].radiances = shis[j].radiances.A.squeeze() + noise

	#_add wavenumber to object
	setattr(shis, 'wavenumber', wn)
	setattr(shis, 'pressure', tile(-9999, 101))
	setattr(shis, 'hbb_nesr', hbb)

	#_build string for metadata
	# don't bother making array so write_collocated can remain stable
	notes = '\n'.join(notes) 

	#_return
	return shis, notes


def run_lbldis_fork(fname, oname):
	'''	launch instance of lbldis from SGE call '''
	from scipy.io import netcdf
	from subprocess import call

	if 'LB_LIBRARY_PATH' not in os.environ:
		os.environ['LD_LIBRARY_PATH'] = '/opt/netcdf-4.1.2/lib/'

	#_launch
	call(['lbldis', fname, '0', oname])


def read_lbldis_fork(oname):
	''' read lbldis output from main call '''
	from scipy.io import netcdf

	#_read in radiances
	ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
	rads = ncdf.variables['radiance'][:].squeeze() #_nchan, ninstances
	wave = ncdf.variables['wnum'][:] 
	ncdf.close()

	return rads, wave


################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
	try:
		#_path of submission to SGE when passed kwarg filename
		from pickle import load

		#_get kwargs from file
		fname, oname, lname = sys.argv[1:]
		run_lbldis_fork(fname, oname)
		os.unlink(lname)

	except ValueError:
		#_generate kwarg from settings
		run_main(**kwargs)
