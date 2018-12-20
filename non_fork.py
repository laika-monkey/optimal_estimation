#!/usr/bin/env python
# 2014.12.07
################################################################################
# Produce simulated S-HIS data to retrieve with OE                             #
################################################################################
from numpy import arange, linspace


kwargs = {
##	'vars'	: ['tau', 'ref', 'surf_temp'],
##	'vars'	: ['tau', 'ref', 'z'],
	'vars'	: ['tau', 'ref'],
##	'vars'	: ['tau'],
	'ranges': {
		'tau'	: linspace(0.05, 2, 40), #_4)
		'ref'	: linspace(0.5, 1, 21),	#_4
		'z'		: linspace(1, 4, 4), 
##		'tau'	: linspace(0.1, 1, 10),
##		'ref'	: linspace(0.1, 3, 11),
##		'surf_temp'	: arange(-5, 6),
##		'surf_temp'	: linspace(-286.8, 296.8, 11),
		},

	#_add noise?
	'noisey'	: True,

	#_some constant values
	'clddef'	: {
		'z'		: 3,	#_keep stable for now
		'ref'	: 0.5,	#_same
		'dbnum'	: 0,
		'ref_wn': 900,
		},

	#_formatting
	'fmt'	: {
		'tau'	: '10.5f',
		'ref'	: '10.5f',
		'z'		: '10.5f',
		'surf_temp'	: '10.5f',
		},

	#_name of simulated S-HIS data file
	#'fname_shis'	: 'flight198202230000.198202230215.cdf',

##	#_tau+ref+surf_temp
##	'fname_shis'	: 'SHIS.CPL.COLLOC.SIMULATED.TAU_REF_NOISE.hdf',

	#_tau+ref
	'fname_shis'	: 'SHIS.CPL.COLLOC.SIMULATED.TAU_REF_NOISE.hdf',

	#_tau+ref+z
##	'fname_shis'	: 'SHIS.CPL.COLLOC.SIMULATED.TAU_REF_Z.hdf',
##	'fname_shis'	: 'SHIS.CPL.COLLOC.SIMULATED.198202260000.198202260215.hdf',

##	#_tau+ref
##	'fname_shis'	: 'SHIS.CPL.COLLOC.SIMULATED.198202240000.198202240215.hdf',

##	#_tau only
##	'fname_shis'	: 'SHIS.CPL.COLLOC.SIMULATED.198202250000.198202250215.hdf',

	'real_shis'		: 'SHIS_rdr20130916T210344end20130917T020353' + 
						'sdr20130917T131448_rad.nc',
	
	#_name of raw forward model output
	'fname_lbl'		: 'lblrtm_SIMULATED.cdf',

	#_lblrtm input directory
	# going to need to rerun with higher vertical resolution (500 m?)
	'dir_lblrtm'	: '/data/wsessions/LBL-RTM/test_fov',

	}


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

	#_write combos in order to output
	fmts = ['{{{0:d}:{1:s}}}'.format(i,key) for i, key in enumerate(fmt)]
	fmt = ','.join(fmts)
##	fmt = '{0:10.5f}, {1:10.5f}, {2:10.5f}\n'
	with open('simulation.txt', 'w') as f:
		fmt0 ='_'.join(['{{{0:d}:10s}},'.format(i) for i,v in enumerate(vars)])
		fmt1 =' '.join(['{{{0:d}:8.5f}},'.format(i) for i,v in enumerate(vars)])
		f.write(fmt0.format(*vars) + '\n')
		[f.write(fmt1.format(*a)+'\n') for a in combos]

	#_use single profile
	prof = fake_prof(size, **kwargs)

	#_make fake epochs to to maintain order 
	epoch_start = d2e('2014010100')

	#_generate fake S-HIS data
	shis, notes = fake_shis(latitude=prof.latitude[0], 
							longitude=prof.longitude[0], **kwargs)

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


def fake_prof(size=1, idx=[2],
	real_shis='SHIS_rdr20130916T210344end20130917T020353' + 
				'sdr20130917T131448_rad.nc', **kwargs):
	''' generate fake profile '''
	from hs3_utils import read_profile_gdas
	from numpy import array, recarray

	#_read in real profile to copy
	real_prof = read_profile_gdas(real_shis, array(idx))

	#_initialize fake profile recarray
	prof = recarray((size,), dtype=real_prof.dtype)

	#_read in a fake clear sky profile, copy it over #size times
	for attr in prof.dtype.names:
		prof.__getattribute__(attr)[:] = real_prof.__getattribute__(attr)

	#_add pressure
	setattr(prof, 'pressure', real_prof.pressure)
	setattr(prof, 'fname', 'SIMULATED')

	return prof


def fake_shis(vars=['tau'], ranges={'tau':[0,]}, dir_lblrtm='.', dz=0.2,
	latitude=0, longitude=0, clddef={}, noisey=False, **kwargs):
	''' generate fake shis retrieval '''
	from numpy import recarray, array, ndarray, tile
	from libtools import combination as combo
	from lblrtm_utils import write_lbldis_parameters as lbldis_write
	from subprocess import call
	from shis import convolve2shis as c2shis
	from scipy.io import netcdf
	from scipy.interpolate import splev, splrep, interp1d
	from hs3_utils import read_shis
	import os
	from numpy.random import random

	#_read in real thing for HBB_NESR
	ugh = '/data/wsessions/hs3/SHIS_rdr20130816T033458end20130816T080506' + \
			'sdr20130817T165537_rad_pcfilt.nc'
	real_shis = read_shis(ugh, idx=array([0])) 
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

	#_ THIS COULD BE FORKED
	#_build clddef dict
	for j, opts in enumerate(combos):
		cloud = [clddef.copy()] 
		[cloud[0].update({v : opts[i]}) for i, v in enumerate(vars)]
		cloud[0]['z_top'] = cloud[0]['z'] + dz

		#_build prefix for labeling
		fmt	= ['{0:s}-{1:05.2f}'.format(v, opts[i]) for i,v in enumerate(vars)]
		notes.append(','.join(fmt))
		fname = '{0:s}/{1:s}.namelist'.format(dir_lblrtm, '_'.join(fmt))
		oname = '{0:s}/{1:s}.output'.format(dir_lblrtm, '_'.join(fmt))

		#_write input paramater file
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kwargs)
		call(['lbldis', fname, '0', oname])

		#_read in radiances
		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze() #_nchan, ninstances
		wave = ncdf.variables['wnum'][:] 
		ncdf.close()

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
		for j in xrange(shis.size):
			posneg = random(hbb.size)
			posneg[posneg < 0] = -1
			posneg[posneg >= 0] = 1
			noise = hbb * posneg
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


if __name__ == '__main__':
	run_main(**kwargs)
