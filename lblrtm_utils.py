#!/usr/bin/env python
# encoding: utf-8
'''
AUTHOR:	Walter R. Sessions
	April 2014
'''
import os, sys, logging, numpy, inspect

from numpy.ma import masked_array, masked_where
from numpy import nan

import libtools as lt

import matplotlib
if 'DISPLAY' not in os.environ: matplotlib.use('Agg')	


LOG = logging.getLogger(__name__)
#_look into logging

#_add julian day to keyword arguments
debug = 1
DIR_PROD = os.environ['PRODUCTS']
JOIN = os.path.join

#############################################################################_80
#_MAIN_######################################################################_80
#############################################################################_80


def run_main():
	return 0


def run_lblrtm_profile(prof,
	dir_lblrtm_out=os.path.join(DIR_PROD, 'LBL-RTM'),
	tape3=os.path.join(DIR_PROD, 'TAPE3_default'),
	**kwargs):
	'''
	run lbl-rtm spectra for given profile and options

	INPUT:
	fname		str		path to file with GDAS data
	tape3		str		path to TAPE3 input file

	OUTPUT:
	A directory under dir_lblrtm_out/fname/NNNNN with the
	desired output type as described by kwargs

	'''
	from libtools import mkdir_p
	from os import symlink, chdir
	
	#_create outdir
	dbg(dir_lblrtm_out, 5)
	mkdir_p(dir_lblrtm_out)
	chdir(dir_lblrtm_out)

	#_build header for tape5 comment line
	header = kwargs.get('out_label') 
				
	#_write_TAPE5
	if write_tape5(prof, header=header, **kwargs):
		pass		#_generate input
	else:
		raise RuntimeError, 'TAPE5 not written'

	#_copy lines file
	if not os.path.exists('TAPE3'):
		symlink(tape3, 'TAPE3')

	dbg(('running lblrtm'))
	if os.system('lblrtm'):
		dbg('LBL-RTM failed')


def run_lblrtm_file(fname,
	dir_gdas=os.path.join(DIR_PROD, 'gdas'),
	dir_lblrtm_out=os.path.join(DIR_PROD, 'LBL-RTM'),
	tape3=os.path.join(DIR_PROD, 'TAPE3_default'),
	cloud_free=True, ocean_only=True, **kwargs):
	'''
	run lbl-rtm spectra for given profile and options

	INPUT:
	fname		str		path to file with GDAS data
	tape3		str		path to TAPE3 input file
	dir_gdas	str		path to collocated input files with gdas/caliop
	cloud_free	bool	require that caliop cloud top be... negative?
	ocean_only	bool	require GDAS land fraction to be zero

	OUTPUT:
	A directory under dir_lblrtm_out/fname/NNNNN with the
	desired output type as described by kwargs

	'''
	from libtools import mkdir_p, setup_groups
	from shutil import copy
	from os import symlink, chdir, fork, waitpid
	from netCDF4 import Dataset
	from numpy import arange
	
	#_create outdir
	mkdir_p(dir_lblrtm_out)

	#_open hdf and read in GDAS and CALIOP data
	with Dataset(fname, 'r', format='HDF4') as hdf:
		gdas = read_gdas(hdf, **kwargs)
		caliop = read_caliop(hdf, **kwargs)	

	#_remove non-clear sky cases
	if cloud_free:
		cidx 	= caliop.cloud <= 0 # -1e-9		#_look at cloud top height
		gdas 	= gdas[cidx]
##		gdas 	= [gdas[i] for i in cidx]
		caliop 	= caliop[cidx]
		
	#_limit to ocean
	if ocean_only:
		idx 	= gdas.land_fraction == 0
		gdas	= gdas[idx]
##		gdas	= [gdas[i] for i in idx]
		caliop	= caliop[idx]

	if not gdas.size:
		dbg('no suitable profiles')
		return	
	
	#_get file name to prefix output
	prefix = '.'.join(fname.split('/')[-1].split('.')[:-1])

	#_break into groups to fork
	groups = setup_groups(range(gdas.size), **kwargs)
	
	#_for use in a really overly obtuse line used to build the header		
	h_fields = ['lat=', 'lon=', 'time=', 'cloud=']

	#_loop over every location, produce LBL-RTM output
	for group in groups:
	  children = []
	  for n in group:
		pid = fork()
		if pid != 0:
			children.append(pid)
		elif pid == 0:
			prof = gdas[n]
			cloud = caliop.cloud[n]

			#_create run directory
			dir_lblrtm_run = os.path.join(	dir_lblrtm_out, prefix,
											str(n).zfill(5))
			dbg(dir_lblrtm_run, 5)
			mkdir_p(dir_lblrtm_run)
			chdir(dir_lblrtm_run)

			#_initialize Profile object
##			prof = Profile(	scene.temperature.pressure,
##							scene.temperature,
##							scene.sfc_temperature,
##							scene.water_vapor,
##							scene.altitude,
##							scene.ozone	)

			#_build header for tape5 comment line
			h_vals = [prof.latitude, prof.longitude, prof.time, cloud]
			h_vals = [str(h) for h in h_vals] 
			header = ','.join([''.join(v) for v in zip(h_fields, h_vals)]) 
				
			#_write_TAPE5
			if write_tape5(prof, header=header, **kwargs):
				pass		#_generate input
			else:
				continue	#_write_tape5 returns false if something is wrong

			#_copy lines file
			if not os.path.exists('TAPE3'):
				symlink(tape3, 'TAPE3')

			dbg(('running lblrtm', str(n)))
			if os.system('lblrtm'):
##				raise RuntimeError, 'LBL-RTM failed ' + str(res)
				dbg(('LBL-RTM failed', n))

			os._exit(0)		

	  for kid in children: waitpid(kid, 0)


def run_lbldis(	dir_lblrtm, dir_lbldis_out=None, out_lbldis=None, 
				verbose='1', label=None, **kwargs):
	'''
	When passed an output directory of lblrtm with
	appropriate TAPE5 and OD files, attempts to run LBL-DIS

	Options for parameter file should be in kwargs. See
	write_lbldis_parameters() for more information.

	dir_lbldis_out	str		if set, output netcdf will be written there
							as 'lbldis_out.nc.'  This will take overrule
							out_lbldis values
	out_lbldis		str		if dir_lbldis_out is None, output will be written
							to out_lbldis + '.cdf'

	returns lbldis output status
	'''
	#_open TAPE5 from lblrtm run and get wavenumber range
	file_tape5 = os.path.join(dir_lblrtm, 'TAPE5')
	with open(file_tape5, 'r') as f:
		lines = f.readlines()
		v1, v2 = lines[2].split()[:2]	#_assumes this line is always wn
		kwargs.update({'v1' : float(v1), 'v2' : float(v2)})

	#_write lbldis import file to output dir
	file_lbldis			= os.path.join(dir_lblrtm, 'lbldis_input')
	file_microwindow	= os.path.join(dir_lblrtm, 'microwindows.txt')
	if label is not None:
		file_lbldis			= '.'.join((file_lbldis, label))
		file_microwindow	= '.'.join((file_microwindow, label))
	kwargs.update({
		'file_microwindow'	: file_microwindow,  
		'filename'			: file_lbldis	})
	write_lbldis_parameters(dir_lblrtm, **kwargs)

	#_if not defined, put output in same directory
	if dir_lbldis_out is None and out_lbldis is None:
		dir_lbldis_out = dir_lblrtm
		out_lbldis = os.path.join(dir_lbldis_out, 'lbldis_out')
	elif out_lbldis is None:
		out_lbldis = os.path.join(dir_lbldis_out, 'lbldis_out')

	if label is not None:
		out_lbldis = '.'.join((out_lbldis, label))

	cmd = ' '.join(('LD_LIBRARY_PATH=/opt/netcdf-4.1.2/lib lbldis',
						file_lbldis, verbose, out_lbldis))
	return os.system(cmd)


def init_layers(clddef, ssp_dist=1, dbnum=0, ref_wn=900, **kwargs):
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

	DEPRECATED.  Terrible idea indeed.
	'''
	#_loop over each layer and perform certain tests
	for i, cld in enumerate(clddef):
		#_make sure tau is a list
		if type(cld['tau']) != list:
			cld['tau'] = [cld['tau']]

		#_first add database number and reference wavenumber
		#_loop over db files and define
		db = dbnum if type(dbnum) != list else dbnum[i]
		rw = ref_wn if type(ref_wn) != list else ref_wn[i]
		cld.update({	'dbnum'		: db, 
						'ref_wn'	: rw	})

		#_then break up ODs 
		if ssp_dist != 1:

		##	#_add basis values
		##	cld.update(clddef)

			#_linearly fractionalize bimodal distribution
			cld['tau'] *= ssp_dist
	
			#_locked into bimodal
			ssp_dist = 1 - ssp_dist
	
			if i > 1:
				raise RuntimeError, 'Only currently supporting bimodal dists.'

	
##	return cld
##	for i, cld in enumerate(clddef):
##
##		#_if using bimodal dist, split linearly
##		if ssp_dist != 1:
##			clddef[i]['tau'][0] = ssp_dist * clddef[i]['tau'][0]
##			ssp_dist = 1 - ssp_dist
##		
##			if i > 1:
##				raise RuntimeError, 'Not currently setup to handle trimodal.'


def clddef_bounds(cld, ftape7, min_tau=0.01):
	''' check cld for cross boundary clouds, adjust '''
	import tape7
	from numpy import array, append

	#_read in tape7
	t7 = tape7.tape7(ftape7)

	#_loop over clouds
	new_cld = []
	for cloud in cld[:]:

		#_calc if within one level
		z0, z1 = cloud['z'], cloud['z_top']

		#_wait, no, they should be equal to a level
##		idx0 = append(t7.zlevel1, z0).argsort().argmax()
##		idx1 = append(t7.zlevel1, z1).argsort().argmax()
		idx0 = abs(t7.zlevel1 - z0).argmin()
		idx1 = abs(t7.zlevel2 - z1).argmin()

		#_same level, move on
		if idx0 == idx1:
			new_cld.append(cloud)
			continue
		elif idx0 == idx1 + 1:
			#_close enough
			new_cld.append(cloud)
	
		#_create new cloud dictionary
		new_cloud = cloud.copy() 

		#_for my sanity	
		idx1 = abs(t7.zlevel1 - z1).argmin()

		#_figure out how many whole levels it is spread across
		nlev = idx1 - idx0 - 1	

		#_how large is layer
		dz = z1 - z0
		tau = new_cloud['tau'][:]

		#_linearly break up tau
		for k in xrange(nlev+1)+idx0:
			ztop = t7.zlevel2[k]
			zbot = t7.zlevel1[k]
			frc = (ztop - zbot) / dz
			new_tau = [ t*frc for t in tau ]	

##			#_only allow optical depths above tau_min
##			tmp	= array(new_tau)
##			tmp[tmp < tau_min] = tau_min
##			new_tau = tmp.tolist()

			#_build cloud entry
			new_cloud = cloud.copy() 
			new_cloud.update({'tau' : new_tau, 'z' : zbot, 'z_top' : ztop})

			#_add to new clddef
			new_cld.append(new_cloud)

	return new_cld


def write_lbldis_parameters(
	#_location of source lblrtm directory
	dir_lblrtm, 

	#_start, ending wavenumber and spectral resolution
	v1=600, v2=3000, dv=0.5,

	#_solar position, zen < 0 indicates below hor
	solar_zen=0., solar_azi=0., solar_dst=1.,

	#_view angle of scene.  0=nadir, 180=zenith 
	obs_angle=0.,

	nstreams=16,
	use_HG=False,
	scattering_off=False,
	surf_temp=-1,
	enable_flux_profile=False, 
	filename='lbldis_input', 

	#_see below
	microwindows=None,
	clddef=[],

	comment='LBLDIS file auto-generated by ' + 
			'lbldis_utils.write_lbldis_parameters()',

	#_surface emissivity. first contains wavnum, second is emissivity. between
	# points, values are linearly interpolated.  Put as many as you'd like
	surf_emis=[[100, 3000], [0.985, 0.985]],

	#_this is something i've done nothing with...
	solar_datafile=JOIN(DIR_PROD, 'lbldis_inputs',
						'solar.kurucz.rad.1cm-1binned.full_disk.asc'),

	#_single scattering property file list
	ssp_db_files=[	JOIN(DIR_PROD, 'lbldis_inputs',
						'ssp_db.shettle_dust.gamma_sigma_0p100'),
					JOIN(DIR_PROD, 'lbldis_inputs',
						'ssp_db.shettle_sand.gamma_sigma_0p100')],

	#_bimodal distribution 
	ssp_dist=1,

	**kwargs):
	'''
	Produce LBLDIS input file. If layer definitions are passed in clddef,
	then an LBL-RTM run directory must be passed as well so that the profile
	can be checked and adjusted for cross boundary levels.

	create input text file for LBLDIS. inputs (via keywords), as follows:

	required input keywords:
	dir_lblrtm	str,	path for input profile
	v1	flt,	starting wavenumber
	v2	flt,	ending wavenumber
	dv	flt,	wavenumber increment for spectral radiance calculation

	Note on v1, v2, dv	: these can contain other values which 
	    indicate microwindows, microwindow files, etc. See LBLDIS 
	    instructions for more details. Note, if the run() method is 
	    used, the microwindow input will he handled correctly.

	microwindows	set to either name of file to use for microwindows
					or to a list of tuples(2) containing the start and end
					of channels. If 'file_microwindow' is in kwargs, the
					file will overwrite that file. 

	ssp_db_files	list(str),	list of string path+filenames for ssp_db files
	solar_datafile	str,		path+filename to solar spectrum file
	clddef	list,	or tuple of dictionaries, containing the following 
					fields for each simulated cloud layer:

	    dbnum	int,	index into ssp_db_files list, specifying the ssp 
						database that should be used for scattering properties
	    ref		flt,	effective radius of integrated particle mixture [um]
	    tau		list,	floats of optical depth of cloud layer.  Can
						be multiple values to have multiple instancds of 
						same layer
	    z		flt,	height of cloud [km]; In this case, the cloud will
						"fill" whatever layer encloses this particular altitude.
						MAY NEED TO ADD VERY SMALL NUMBER TO MAKE SURE
						CORRECT LAYER IS BEING FILLED.  THIS IS A BOUNDARY!!
	    ref_wn	flt,	reference wavenumber for tau calculation; pick a
						non-positive number to assume the short wavelength,	
						geometric optics limit (Q_ext = 2)

	optional input keywords (defaults noted where applicable):

	use_HG		bool,	(default false), on whether to use HG phase functions 
						for all scattering calculations.
	filename	str,	defaults to 'lbldis_input' in the current directory.
	solar_zen	flt,	solar zenith angle; value < 0 implies sun below
						horizon [deg]
	solar_azi	flt,	solar relative azimuth angle [deg]
	solar_dst	flt,	solar distance (in AU, 1=Earth)
	obs_angle	flt,	observer zenith angle [deg]; 0==nadir, 180 zenith.
	nstreams	int,	number of DISORT streams; default is 16
	lblrtm_dir	str,	path to LBLRTM optical depth run. 
						default is './lblrtm_output'
	surf_temp	flt,	surface skin temperature; defaults to -1, 
						which implies the temperature at the bottom
						level of the LBLRTM run is reused.
	surf_emis	list[n,2] shaped array, containing wavenumber (first column), 
						and emissivity (second column); this emissivity 
						spectrum is linearly interpolated to the simulated
						wavenumber range. default is a gray emis=0.985 surface;
						the array is defined as:[ [100, 0.985], [3000, 0.985] ]

	enable_flux_profile bool,	to true to compute the per-level
								spectral flux profile.

	If any clouds were input with thicknesses (e.g., ztop/zbottom), then the 
	levels keyword is required; this keyword is normally just sent to LBLRTM

	ssp_dist	float,	fractional (<=1), setting this causes two, coincident
						layers to be modeled with the first having ssp_dist
						of the OD

	'''	
	import numpy as np
	ilines = []
	# to shorten syntax.
	l = ilines.append

	#_kludge
	for i, cld in enumerate(clddef):
		#_make sure tau is a list
		if type(cld['tau']) != list:
			clddef[i]['tau'] = [cld['tau']]

	#_if clddef passes layers that crosses prof boundaries,
	# split up and adjust clddef (requires TAPE7)
	clddef = clddef_bounds(clddef, os.path.join(dir_lblrtm, 'TAPE7'))

	try:
		numtau = len(clddef[0]['tau'])
	except IndexError:
		numtau = 0

	numcld = len(clddef)
	numssp = len(ssp_db_files)

	#_line 1: commment
##	l('LBLDIS file auto-generated by lbldis_utils.write_lbldis_parameters()')
	l(comment)

	#_line 2: number of streams
	l('{0:<16d}Number of streams'.format(nstreams))

	#_line 3: solar information
	l('{0:4.1f} {1:4.1f} {2:4.1f} '.format(solar_zen, solar_azi, solar_dst) + \
	       'Solar zenith angle (deg), relative azimuth (deg), ' + \
	       'solar distance (a.u.)' )

	#_line 4: zenith angle
	l('{0:<16.1f}'.format(obs_angle) + \
	       'Zenith angle (degrees): 0 -> upwelling, 180 -> downwelling' )

	#_line 5: wavelength range/resolution
	if type(microwindows) is str:
		#_file already generated, print name
	    l('-1 0 0 ' +  microwindows)
	elif hasattr(microwindows, '__iter__'):
		#_write microwindow file
	    l('-1 0 0 ' + write_microwindow(microwindows, **kwargs))
	else:
		#_else use prescribed range, resolution
	    l(	'{0:7.2f} {1:7.2f} {2:6.3f} '.format(v1, v2, dv) + 
			'v1, v2, and dv [cm-1]' )

	if scattering_off:
	    dbg('WARNING: Turning scattering off!')
	    cloud_param_spec = -numtau
	else:
	    cloud_param_spec = numtau

	#_line 6: cloud parameters
	l( '{0:<16d}'.format(cloud_param_spec) +
	   'Cloud parameter option flag: 0: ' +
	   'reff and numdens, >=1:  reff and tau')

	#_line 7: number of cloud layers
	l('{0:<16d}Number of cloud layers'.format(numcld))

	#_line 7.1-7.ncloud: cloud layer details
	for n in range(numcld):

		dbnum = clddef[n]['dbnum']	#_which ssp database to use (0 == first)
		if use_HG:					# and -1 == Henyey-Greensteinfielder
			dbnum = dbnum + 1
			dbnum = -1 * dbnum
		
		#_some things may send z... low..
		tmp_z = max([clddef[n]['z'], 1e-4]) # clddef[n]['z'], 
		p1 = '{0:1d} {1:8.4f} {2:8.5f} {3:5.0f} '.format(dbnum,
														tmp_z, 
														clddef[n]['ref'],
														clddef[n]['ref_wn'])
		p2 = " ".join(['{0:10.5f}'.format(float(t))
	                   for t in clddef[n]['tau']])
		l(p1 + p2)

	#_line 8: lblrtm output path 
	l(dir_lblrtm)

	#_line 9: ascii solar source function file path
	l(solar_datafile)

	#_line 10: number of scattering property databases
	l('{0:<16d}Number of scattering property databases'.format(numssp))

	#_line 10.1 - 10.n_databases
	for ssp_db_file in ssp_db_files:
	    l(ssp_db_file)

	#_line 11: surface temperature
	l( '{0:<12.2f}'.format(surf_temp) +
	   'Surface temperature (specifying a negative value takes ' +
	   'the value from profile)' )

	emis	= np.c_[surf_emis].T
	numemis = emis.shape[0]

	#_line 12: number of lines in the surface emiss spectrum
	l( '{0:<12d}Number of surface spectral '.format(numemis) +
	   'emissivity lines (wnum, emis)')

	#_line 12.1-12.n_lines
	for n in range(numemis):
	    l( '{0:8.3f} {1:8.6f}'.format(emis[n,0],emis[n,1]) )

	l('{0:d}'.format(enable_flux_profile))

	with open(filename, 'w') as f:
		for line in ilines:
			f.write(line + '\n')

	return filename


def write_microwindow(windows, file_microwindow='microwindow.txt', **kwargs):
	'''
	writes microwindow file described by third method in lbldis.instructions

	windows				iter(Nx2)	Array containing microwindow pairs
	file_microwindow	str			path to microwindow file
	
	returns:
		file_microwindow

	N windows
	vstart[0]	vend[0]
	vstart[1]	vend[1]
	...
	vstart[N]	vend[N]
	'''
	with open(file_microwindow, 'w') as f:
		nwin = len(windows)
		f.write(str(nwin) + '\n')

		for window in windows:
			if len(window) != 2:
				raise ValueError, 'windows should contain two elements'

			v0, v1 = window
			if v0 > v1:
				raise ValueError, 'window values should increase'

			chan = '{0:10.5f} {1:10.5f}\n'.format(v0, v1)
			f.write(chan)

	return file_microwindow


def microwindow_average(rads, wave, mw, error=False,
	mw_file='/home/wsessions/src/lbldis/microwindows.h', **kwargs):
	'''
	Averages radiance values to pre-specified channel values found in LBL-DIS

	rads	np.ndarray,	radiance values
	wave	np.ndarray,	corresponding wavenumber values
	mw		-int,		used to read in channels form mw_file
			array,		if array is passed, uses those instead of mw header
						WAIT NO NEED CHANNEL WIDTH NOT JUST WNUM
						LEAVING IN SO TO PREVENT THIS FROM COMING BACK UP
		
	error	bool,		if averaging the RMSE, adjust average by number
						of channels, sqrt

	returns radiance values for channels and new wavenumbers for sanity,
	but don't use the new wavenumbers for anything. They should be discarded.

	This does not perfectly match the way microwindows will come out of
	lbldis, but the errors were on the order of 1/100 Kelvin and are being
	ignored.
	'''
	from numpy import zeros, sqrt, isnan, array
	from numpy.ma import masked_where, mean

	#_read in channel bounds (nchan, 2) 
	channels = read_microwindows_header(mw_file)[abs(mw+1)] 
	nchan, t = channels.shape

	#_initialize microwindow wavenumbers and channels
##	new_rads = zeros((nchan)) - 9999.
##	new_wave = zeros((nchan)) - 9999.

	new_rads = []
	new_wave = []

	#_loop over chan
	for n in xrange(nchan):

		#_find nearest indices for the wavenumbers
		idx0 = abs(wave - channels[n,0]).argmin() 
		idx1 = abs(wave - channels[n,1]).argmin() #+ 1 
	
		#_pull out channel values
		rad_chan = rads[idx0:idx1]

		#_if missing region, skip
		if len(rad_chan) == 0: continue

		rad_chan = masked_where(rad_chan < 0, rad_chan)
		rad_chan = masked_where(isnan(rad_chan), rad_chan)

		#_if adjusting error
		factor = 1 if not error else sqrt(idx1-idx0)

		new_rads.append(mean(rad_chan) / factor)
		new_wave.append(mean([wave[idx0], wave[idx1]]))
	###	new_rads[n] = mean(rad_chan) / factor
	###	new_wave[n] = mean([wave[idx0], wave[idx1]])

	##	new_wave[n] = wave[idx0:idx1].mean()

	#_return new averages
	new_rads = array(new_rads)
	new_wave = array(new_wave)	

	return new_rads, new_wave


def microwindow_stdev(rads, wave, mw,
	mw_file='/home/wsessions/src/lbldis/microwindows.h', **kwargs):
	'''
	Returns standard deviation values to pre-specified
	channel values found in LBL-DIS

	rads	np.ndarray,	radiance values
	wave	np.ndarray,	corresponding wavenumber values
	mw		-int,		used to read in channels form mw_file

	returns radiance values for channels and new wavenumbers for sanity
	'''
	from numpy import array, zeros
	from numpy.ma import masked_where

	#_read in channel bounds (nchan, 2) 
	channels = read_microwindows_header(mw_file)[abs(mw+1)]
	nchan, t = channels.shape

	#_initialize microwindow wavenumbers and channels
##	new_rads = zeros((nchan)) - 9999.
##	new_wave = zeros((nchan)) - 9999.
	new_rads = []
	new_wave = []

	#_loop over chan
	for n in xrange(nchan):

		#_find nearest indices for the wavenumbers
		idx0 = abs(wave - channels[n,0]).argmin() 
		idx1 = abs(wave - channels[n,1]).argmin() 

		#_pull out channel values
		rad_chan = rads[idx0:idx1]

		#_skip in empty regions
		if len(rad_chan) == 0: continue

		rad_chan = masked_where(rad_chan < 0, rad_chan)

		new_rads.append(rad_chan.std())
		new_wave.append(wave[idx0:idx1].mean())
	#	new_rads[n] = rad_chan.std()
	#	new_wave[n] = wave[idx0:idx1].mean()

	new_rads = array(new_rads)
	new_wave = array(new_wave)

	#_return new averages
	return new_rads, new_wave


def plot_mw():
	''' read in a random spectrum and plot microwindows over it '''
	import matplotlib.pyplot as plt
	from hs3_utils import Flight_segment
	from numpy import array, append
	from libtools import shrink_ticks
	from libgeo import planck_inv
	from oe_utils import read_airs

	mw = read_microwindows_header()	
	
	fig = plt.figure()
	axL = fig.add_subplot(111)
	axL.set_title('microwindows available')

##	f = Flight_segment(dtg='20130820190000')
#3	rads = f.SHIS_radiances[0]
##	wave = f.SHIS_wavenumber
	f = '/odyssey/isis/users/wsessions/sips/airs/AIRS.2013.08.24.004.L1B.AIRS_Rad.v5.0.21.0.G13238114023.hdf' 
	rads, wave = read_airs(f, 0, 0)
	rads = rads.AIRS_radiances[0]
	print rads.size, wave.size
	Tb = planck_inv(rads, wave*100, domain='wavenumber')	
##	Tb = planck_inv(f.SHIS_radiances[0], f.SHIS_wavenumber*100,
##		domain='wavenumber')	
	axL.plot(wave, Tb, lw='0.1', color='gray')	

	axR = axL.twinx()
	minx = 1000
	for mw_num, mw_chans in mw.iteritems():	
		nchan, nbound = mw_chans.shape

		#_add an extra layer of Nones
		mw_plot = append(mw_chans, array([None]*nchan)[:,None], 1).flatten()

		#_plot line
##		axR.plot(mw_plot, [mw_num]*nchan*(nbound+1))
		axR.hlines([mw_num]*nchan, xmin=mw_chans[:,0], xmax=mw_chans[:,1],lw=10)

		#_get left bound
		minx = min(minx, mw_chans[0][0])

	#_set bounds
	axR.grid(True)
	axR.set_ylim(-1, mw_num + 1)
	axL.set_ylim(200,310)
	axL.set_xlim(left=minx-100)
	
	#_turn on markers for all options
	axR.set_yticks(mw.keys())
	axR.set_yticklabels([nn*-1-1 for nn in mw.keys()])
	axR.set_xlabel('wavenumber (cm-1)', size='xx-small')
	shrink_ticks(axL)
	
	plt.savefig('microwindows.png')


def read_microwindows_header(fname='/home/wsessions/src/lbldis/microwindows.h'):
	''' read in microwindows.h file and return dictionary of channel bounds '''
	import re
	from numpy import array
	dbg(fname, 5)

	#_if you don't understand regular expressions, this is a bad place to start
	re_chan = re.compile('if\(select_microwindows\s*==\s*(\d*)\).*?{(.*?)}',
			re.DOTALL)
	re_size = re.compile('num_microwindows\s*=\s*(\d+)')
	re_vals = re.compile('microwindows\[\d+]\[\d+]\s*=\s*([\d.]+);')
	re_comm = re.compile('/\*.*\*/')

	#_initialize returned dict
	mw_options = {}

	with open(fname, 'r') as f:
		#_read in everything, ditch comments
		lines = f.read()

		#_create list of matching blocks of channel numbers
		blocks = re_chan.findall(lines)

		for block_tup in blocks:
			#_option id and block strings are initially in tuple
			mw_key, mw_block = block_tup
		
			#_bulk of block text
			nchan = int(re_size.findall(mw_block)[0])	

			#_pull out microwindow strings
			bounds = re_vals.findall(mw_block)

			if len(bounds) != nchan*2:
				dbg((nchan, bounds))
				raise ValueError, 'incorrect number of channels'

			#_convert to float and rearrange
			bounds = array(bounds, dtype=float).reshape(nchan, 2)

			#_add to dictionary
			mw_options.update({int(mw_key) : bounds})	

	return mw_options


def read_lbldis_out(fname, fname_input=None, **kwargs):
	''' read in netcdf file produced by lbldis '''
	from netCDF4 import Dataset
	from numpy import float64
	import re

	dbg(fname)
	cdf = Dataset(fname, 'r')
	label = '.'.join(fname.split('/')[-1].split('.')[1:-1])
	dir_lbl = cdf.Path_to_LBLRTM_gas_ODs.replace('-RAW','')
	#_...

	wnum = cdf.variables['wnum'][:].squeeze()
	rads = cdf.variables['radiance'][:].squeeze()
	rads /= 1e3	#_convert mW->W (GO THROUGH AND MAKE SURE CONSISTENT...toomany)

	#_get resolution from input file
	if fname_input is None:
		fname_input = '{0:s}/lbldis_input.{1}'.format(dir_lbl, label)

	#_moving to iris messed up a bunch of shit.  why am I even using tmp dirs?
	if re.search('/scratch/wsessions/tmp//', fname_input):
		tmp = fname_input.split('/')[-3:]
		tmp = '/'.join(tmp)
		dir_kludge = os.path.join(DIR_PROD, 'LBL-RTM_hs3')
		fname_input = os.path.join(dir_kludge, tmp)

		kludge = os.path.dirname(fname_input) 
	else:
		kludge = False 
	lbldis_input = Lbldis_input(fname_input)

	#_get metadata from LBL-RTM files
##	with open(os.path.join(cdf.Path_to_LBLRTM_gas_ODs.replace('-RAW',''),'TAPE5')) as f:
##	with open(os.path.join(cdf.Path_to_LBLRTM_gas_ODs,'TAPE5')) as f:
	if kludge:
		f = open(os.path.join(kludge, 'TAPE5'))
	else:
		
		try:
			f = open(os.path.join(cdf.Path_to_LBLRTM_gas_ODs,'TAPE5'))

		except IOError:
			f = open('/odyssey/isis/users/wsessions/LBL-RTM_simulated/std_tropical/TAPE5')
 
	#_get first line
	header = f.readline()

	lat		= float64(re.search('lat=([\d.-]+)', header).group(1))
	lon		= float64(re.search('lon=([\d.-]+)', header).group(1))
	time	= float64(re.search('time=([\d.-]+)', header).group(1))
	try:
		cloud = float64(re.search('cloud=([\d.-]+)', header).group(1))
	except AttributeError:
		cloud = None
	f.close()

	
	#_if multiple instances, pass multiple spectrums.
	# MAKE THIS LESS STUPID
	n_instances = len(cdf.dimensions['n_instances'])
	if n_instances > 1:
		spectrums = []
		for i in xrange(n_instances):
			spectrums.append(Spectrum(	radiances=rads[:,i], wavenumber=wnum,
										label=label, latitude=lat,
										longitude=lon, time=time, 
										cloud=cloud, **kwargs	))
			#_set resolution value on spectrum attribute
			setattr(spectrums[i], 'resolution', lbldis_input.resolution)
	else:
		spectrums = Spectrum(	radiances=rads, wavenumber=wnum, label=label,
								latitude=lat, longitude=lon, time=time,
								cloud=cloud, **kwargs)
		#_set resolution value on spectrum attribute
		setattr(spectrums, 'resolution', lbldis_input.resolution)

	cdf.close()	

	return spectrums

##	currently ignored
##	flux_up, flux_down, flux_beam, cld_height, cld_phase, cld_tau
##	cld_tau_wnum, cld_reff, sfc_temperature, sfc_emissivity, pwv, height
##	pressure, temperature, mixing_ratio


class Lbldis_input(object):
	def __init__(self, fname):
		'''
		read in lbldis parameters file and create object.
		currently only managing what I need to pull out... DEADLINES 2.28.15

		ONLY WORKS FOR SINGLE LAYER AT THE MOMENT
	
		I forget if some of these lines above the layer items can move, check
		later.
		'''
		from numpy import arange
		with open(fname,'r') as f:
			dbg(fname)
			lines = f.readlines()

			#_add resolution parameters
			wn0, wn1, res = lines[4].split()[:3]
			setattr(self, 'resolution', float(res)) 

			#_read in cloud
			nlay = int(lines[6].split()[0])
			cld = []
			idx = 7 
			for i in range(nlay):
				#_for each layer, read in optical depth
				cols = lines[idx].split()
				dbnum, z = int(cols[0]), float(cols[1])
				ref, ref_wvnum = float(cols[2]), float(cols[3]) 
				tau = [float(tt) for tt in cols[4:]]

				#_update cloud dictionary
				cld.append({'tau' : tau, 'dbnum' : dbnum, 'z' : z, 
					'ref' : ref, 'ref_wvnum' : ref_wvnum})

				idx += 1 

			#_idx for line after layers
			# i+8 == directory for run
			# i+9 == solar
			#_read in list of ssp databases
			idx += 2; nssp = int(lines[idx].split()[0]); idx += 1
			ssp_db = []
			[ssp_db.append(lines[i].strip()) for i in arange(nssp) + idx] 

			#_add layer parameters
			setattr(self, 'layers', cld)
			setattr(self, 'ssp_db', ssp_db)

	def merge_layers(self):
		''' combine optical depths of layers using same ssp linearly '''
		new_layers = {}	#_keyed by db_num
		for layer in self.layers:
			dbnum	= layer['dbnum']

			#_update
			try:
				L			= new_layers[dbnum]
				L['z']		= min(L['z'], layer['z']) 
				L['z_top']	= max(L['z'], layer['z']) 
				L['tau'][0] += layer['tau'][0]
				
			#_initialize	
			except KeyError:
				new_layers.update({ dbnum : {} })
				L = new_layers[dbnum]
				L['tau']		= layer['tau']
				L['ref']		= layer['ref']
				L['ref_wvnum']	= layer['ref_wvnum']
				L['dbnum']		= dbnum	
				L['z']			= layer['z']
				L['z_top']		= layer['z'] + 0.01

		#layers_list = [None] * len(self.ssp_db)
		#for i, ssp in enumerate(self.ssp_db):
		#	try:
		#		layers_list[i] = new_layers[i]

		#	except KeyError:
		#		print 'WHAT" NOW'
		#		print layers_list
		#		print new_layers

		#		raise KeyError, 'DERP {0}'.format(i)	
	
		layers_list = []	
		for key, values in new_layers.iteritems():
			layers_list.append(values)

		self.layers = layers_list[:] 

	def sum_tau(self, layer_type='aerosol'):
		'''
		sum all lbldis optical depth outputs. MOVE THIS METHOD TO CLASS

		layer_type should be either aerosol or cloud
		'''

		#_keywords in aerosol ssp discrimination
		if layer_type == 'aerosol':
			ssps = ['dust','sand','gypsum','kaolinite','quartz']

		elif layer_type == 'cloud':
			ssps = ['wat', 'ice']

		else:
			raise RuntimeError, 'Invalid layer type {0}'.format(layer_type)

		#_initialize tau
		tau = 0

		#_pull out ssp labels
		ssp_labels = [ssp.split('.')[1] for ssp in self.ssp_db]
		ssp_labels = [ssp.split('_')[1] for ssp in ssp_labels]

		#_loop over layers
		for i, layer in enumerate(self.layers):

			#_check if aerosol ssp
			if ssp_labels[layer['dbnum']] in ssps:
			    tau += layer['tau'][0]

		return tau


def write_kwargs(name='full', sensor='nir_to_co2',
	path_out=os.path.join(DIR_PROD, 'LBL-RTM'), 
	opts={	'ref' : [0.],
			'tau' : [0.],
			'zkm' : [[0,0]]},
	order = ['tau','ref','zkm'], 
	dummy_levs = 0.5,	#_no greater spacing than 500m
	**kwargs):
	'''
	create pickles for kwargs options to pass to run_lbl/dis runs 
	name	str,	'full' or 'mw', sets runs to either full spectrum
					lbldis runs or a specific mw channel, currently
					defaulted to option -14.  Name is then used to name
					output file
	sensor	str,	references name in sensors.py module defining 
					starting and ending wavenumbers as well as spectral
					resolution ('full', 'airs', 'nir_to_co2', etc)
	path_out	str,	output directory 

	opts	dict,	contains the sensitivity axes for the lbldis runs.
					Don't go too nuts with trying to feed this things
					that aren't in a similar form to the above just yet.	
	'''
	from libtools import combination
	import sensors
	import pickle

	#_options to loop over (for lbldis) 
	options = [ opts[values] for values in order ]
	runs	= combination(options, dtype=str)

	#_expand
	dv = 0.5 if name == 'full' else -24

	#_intialize with clear-sky case
	label_fmt = '_'.join((name, 'r{0:4.2f}_od{1:3.1f}_z{2:3.1f}-{3:3.1f}'))
	default_opts = sensors.__getattribute__(sensor)
	default_opts.update({
			'dv'    : dv,

			#_lbldis options
			'merge' : 1,
			'iemit' : 0,
			'iatmm' : 1,
			'ipunch': 1,
			'ioddd' : 0,
			'cntnm' : 5,
			'nproc' : 1,

			#_<experiment>_<dtg0>.<dtg1>_fov<000000>
			'dir_lblrtm_fmt' : os.path.join(path_out, '{3}_{1}.{2}_fov{0:06d}'),
			'label' : label_fmt.format(0., 0., 0., 0.), })

	#_LBLDIS
	kw_lbldis = [default_opts]
	for ref, aod, zkm in runs:
		#_build output kwargs dictionary
		kw = default_opts.copy()
		kw.update({
			#_this is hideous.  Make it a dicitonary of lists
			'clddef': [{
				'z'     : float(zkm[0]),    #_height of layer 
				'z_top' : float(zkm[1]),    #_height of layer top 
				'ref'	: float(ref),       #_effective radius
				'tau'   : [float(aod)],     #_optical depth
				'ref_wn': 900.0, #_11.1um	# at this wnum
				'dbnum' : 0,    }],         #_which ssp file to use

			'label' : label_fmt.format(ref, aod, zkm[0], zkm[1]),   })

		kw_lbldis.append(kw)

	#_LBLRTM
	kw_lblrtm = {   'clddef' : [{'z' : 2.0, 'z_top' : 4.0},
                                {'z' : 3.0, 'z_top' : 5.0},
                                {'z' : 4.0, 'z_top' : 6.0}],
					'dummy_levs' : dummy_levs,	}
	kw_lblrtm.update(default_opts)

	#_write lblrtm run options to file
	file_kwarg = os.path.join(DIR_PROD, 'LBL-RTM', 'kwarg_lblrtm.pk')
	pickle.dump(kw_lblrtm, open(file_kwarg, 'wb'))
	dbg(file_kwarg)

	#_write lbldis run options to file
	file_kwarg = os.path.join(DIR_PROD, 'LBL-RTM', 
				'{0}_kwarg_lbldis.pk'.format(name))
	pickle.dump(kw_lbldis, open(file_kwarg, 'wb'))
	dbg(file_kwarg)


class Ssp(object):
	import re
	def __new__(data, fname=os.path.join(os.environ['PRODUCTS'], 
				'lbldis_inputs/ssp_db.shettle_dust.gamma_sigma_0p100'), 
				**kwargs):
		from numpy import array, ndarray, recarray

		with open(fname, 'r') as f:
			lines = f.readlines()

			header = {
				'wavenumber' : {
					'idx'	: 1,
					'units'	: 'cm-1',
					'type'	: 'f8',	},	
				'wavelength' : { 
					'idx'	: 0,
					'units'	: 'um',
					'type'	: 'f8',	},	
				'reff' : {
					'idx'	: 2,
					'units'	: 'um',
					'type'	: 'f8',	},	
				'extinction': {
					'idx'	: 3,
					'units'	: 'um^2',
					'type'	: 'f8',	},	
				'scattering': {
					'idx'	: 4,
					'units'	: 'um^2',
					'type'	: 'f8',	},	
				'absorption': {
					'idx'	: 5,
					'units'	: 'um^2',
					'type'	: 'f8',	},	
				'w0': {
					'idx'		: 6,
					'long_name' : 'single_scattering_albedo',
					'type'		: 'f8',	},	
				'g': {
					'idx'		: 7,
					'long_name'	: 'asymmetry parameter',
					'type'		: 'f8',	},	
				'Qext': {
					'idx'		: 8,
					'long_name'	: '',
					'type'		: 'f8',	},	
				'Qabs': {
					'idx'		: 9,
					'long_name'	: '',
					'type'		: 'f8',	},
				'Qsca': {
					'idx'		: 10,
					'long_name'	: '',
					'type'		: 'f8',	},	
				'volume': {
					'idx'	: 11,
					'units'	: 'um3',
					'type'	: 'f8',		},	
				'proj_area': {
					'idx'	: 12,
					'units'	: 'um2',
					'type'	: 'f8',		},
				'phase_function': {
					'units'	: 'radians',
					'type'	: ndarray,	},	
			}

			#_build dtype for recarray
			dtype = []
			for col, opts in header.iteritems():
				dtype.append((col, opts['type']))

			#_initialize recarra
			data = recarray((0,), dtype)
		
			#_assumes that these data always start header lines	
			data.__setattr__('gamma', float(lines[1].split()[-1]))
			data.__setattr__('nline', int(lines[2].split()[0]))
			data.__setattr__('nangle', int(lines[3].split()[0]))
			data.__setattr__('fname', fname.split('/')[-1])
	
			#_read in phases	
			angles = array([float(i) for i in lines[5].split()])
			data.__setattr__('angles', angles)
			if data.angles.size != data.nangle:
				raise ValueError, 'metadata does not match phase angles'	

			#_loop over ssp data and store phase function values
			for line in lines[9:]:
				#_split columns up
				cols = line.split()

				data.resize(data.size + 1)
				for col, opts in header.iteritems():
					try:
						data[-1].__setattr__(col, cols[opts['idx']]) 
					except KeyError:
						data[-1].__setattr__('phase_function', array(cols[13:]))

##					#_set units and long_names
##					if 'units' in opts:
##						data[-1].__getattribute__(col).__setattr__('units',
##																opts['units'])

			def plot():
				''' create plot of ssp information '''
				import matplotlib.pyplot as plt
				from numpy import pi
	
				ax0 = plt.subplot2grid((3,2), (0,0))
				ax0.plot(data.angles)
				ax0.set_xlabel('radians')
				ax0.set_title('phase angles')
				ax0.grid(True)

				ax1 = plt.subplot2grid((3,2), (0,1))
				ax1.set_title('wavenumbers, gamma: ' + str(data.gamma))
				wnum = list(set(data.wavenumber))
				wnum.sort()
				ax1.plot(wnum)
				ax1.grid(True)

				ax2 = ax1.twinx()
				ref = list(set(data.reff))
				ref.sort()
				ax2.plot(ref)

				ax3 = plt.subplot2grid(	(3,2), (1,0), polar=True, 
										colspan=2, rowspan=2)
				ax3.plot(pi*data.angles/180., data.phase_function[0])	
				plt.savefig('ssp.png')

			data.__setattr__('plot', plot)
	
			return data


def get_ssp_file(ssp_db_files=None, dbnum=0, **kwargs):
	''' returns path to database being used... don't ask '''
	return ssp_db_files[dbnum]


def check_reff_v_ssp(r_eff, delta=0.1 ,ssp_file=None, dbnum=0, **kwargs):
	''' 
	Checks to make sure the effective radius is within
	the range of the ssp database provided

	r_eff		float,	effective radius value
	ssp_db_file	 list,	ssp file being used for lbldis run	
	'''
	from ssp_db import read_ssp_db

	#_pull out db path
	if ssp_file is None:
		ssp_db_files	= kwargs.get('ssp_db_files')
		ssp_file		= ssp_db_files[dbnum] 

	#_read in effective radius range
	ssp = read_ssp_db(ssp_file)['r_eff']
	dbg((ssp.max(), ssp.min(), r_eff))	
	#_check if in bounds
	if r_eff >= ssp.max() - delta:
		return ssp.max()- delta, 'max'
	elif r_eff < ssp.min():
		return ssp.min(), 'min'
	else:
		return r_eff, None


def _gen_profile():
	from numpy import recarray, linspace
	#_when doing this for real, make class that uses dictionary
	# to generate the dtype, units, et cetera
	size	= 25
	dtype 	= [	('pres','f4'),('tdry','f4'),('w','f4'),('alt','f4'),
				('oz','f4'),('co','f4')] 
	d 		= recarray((size,),dtype=dtype)
	d.pres	= linspace(1000.,100,size)
	d.tdry	= linspace(310,250,size)
	d.w		= linspace(1,0,size)
	d.alt	= linspace(0,100,size)
	d.oz	= linspace(0,100,size)
	d.co	= linspace(300,0,size)
	
	return d


def write_in_lblrtm_cld(path_lblrtm_run='.', **kwargs):
	'''
	aersl = 5 in TAPE5 requires in_lblrtm_cld file to define
	aerosol layers.  This writes one.
	'''
	fname = os.join.path(path_lblrtm_run, 'in_lblrtm_cld')

	with open(fname, 'w') as f:
		pass
		#_record 1: number of frequency points

		#_record 2: n_freq frequencies at which layer spectral OD are prov

		#_record 3: number of layers for which the OD are provided

		#_repeat 4.? for each layer
		#_record 4.1: layer index, layer pressur

		#_record 4.2: cloudodlayer(n_freq) n_layer records of OD 
		#				with n_freq values in each record


def read_tape5(fname, profile_source='GDAS', **kwargs):
	'''
	read in TAPE5 file back into a profile object that can be 
	handed back to write_tape5

	fname   str,    relative or absolute path to TAPE5 file
	'''
	import re
	from numpy import array, mean, diff, round
	from collections import Counter

	#_open old tape5 file
	f = open(fname, 'r')
	
	#_initialize return dictionary
	t5kw = {}

	#_rec 1.1, header comment
	t5kw.update({'header' : f.readline()[1:].strip()})

	#_rec 1.2, options
	#_COULD JUST SPLIT
	rel2 = re.compile('HI=(\d)\s*F4=(\d)\s*CN=(\d)\s*AE=(\d)\s*EM=(\d)\s*'+
	                    'SC=(\d)\s*FI=(\d)\s*PL=(\d)\s*TS=(\d)\s*AM=(\d)\s*'+
	                    'MG=(\d)\s*LA=(\d)\s*(?:OD|MS)=(\d)\s*XS=(\d)\s*(\d)'+
	                    '+\s+(\d+)')
	opts = ['hirac','lblf4','cntnm','aersl','iemit','iscan','filtr','iplot',
	        'itest','iatmm','merge','laser','ioddd','xsect','mpts','npts']
	res = rel2.search(f.readline())
	[t5kw.update({opt : int(res.group(i+1))}) for i, opt in enumerate(opts)]

	#_read rec 1.2a
	if t5kw['cntnm'] == 6:
	    raise RuntimeError, 'not yet implemented'

	#_read rec 1.2.1
	if t5kw['iemit'] == 2:
	    opts = ['inflag','otflag','julday']
	    [ t5kw.update({opt[i] : float(val)}) for i, val \
	        in enumerate(f.readline().split()) ]

	#_rec 1.3: wavenumber bounds, etc
	if t5kw['hirac'] > 0 or t5kw['aersl'] > 0 \
	or t5kw['iatmm'] == 1 or t5kw['laser'] > 0:
	    t5kw['v1'], t5kw['v2'] = [float(d) for d in f.readline().split()]

	#_rec 1.3.a+b
	if 0: #_nmol_scal > 0??
	    raise RuntimeError, 'not yet implemented'

	#_rec 1.4: temperature and emissivity boundary
	if t5kw['iemit'] == 1 or (t5kw['iemit'] ==2 and t5kw['otflag'] == 2):
	    opts = ['surf_temp','bndem0','bndem1','bndem2',
	           'bndrf0','bndrf1','bndrf2','sfcrfc']
	    [ t5kw.update({opts[i] : float(val)}) for i, val \
	            in enumerate(f.readline().split()) ]

	#_rec 1.6a:
	if t5kw['merge'] in [35,36,40,41,45,46]:
	    raise RuntimeError, 'not yet implemented'

	#_rec...3.1?
	opts = ['model','itype','ibmax','nozero','noprnt','nmol',
	        'ipunch','ifxtyp','munits','re','hspace','vbar']
	[t5kw.update({opts[i] : float(val)}) for i, val \
	        in enumerate(f.readline().split())]

	#_rec 3.2: path parameters
	opts = ['h1','h2','angle']
	[t5kw.update({opts[i] : float(val)}) for i, val \
	        in enumerate(f.readline().split())]
	
	#_rec 3.3b: altitudes of layer boundaries, normally 8 per line
	# We calculate these from the pressure levels, so don't read them
	# in for now
	dummy_levs = 1.0
	levs = []
	while 1:
		vals = f.readline().split()

		try:
			levs.extend([float(v) for v in vals])
		except ValueError:
			#_hitting 'profile'
			break

	    #_exit at line with "## levels in user defined profile"
		if vals[-1] == 'profile':
			break

	#_this might break things
	count = Counter(diff(levs))
	dummy_levs = count.most_common(1)[0][0]
	t5kw['dummy_levs'] = dummy_levs

	#_rec 3.?: LBL-RTM profile levels
	nlev = int(vals[0])
	opts = [['altitude','pressure','temperature','char','jchar'],
	        ['water_vapor','co2amt','ozone','dum00','dum01','dum02']]
	[ t5kw.update({opt : []}) for opt in opts[0] ]
	[ t5kw.update({opt : []}) for opt in opts[1] ]
	for i in range(nlev*2):
	    optz = opts[i%2]
	    line = f.readline().split()
	    for j, val in enumerate(line):
	        try:
	            t5kw[optz[j]].append(float(val))
	        except ValueError: #_allow strings to pass
	            t5kw[optz[j]].append(val)

	#_hrmph
	letter = t5kw['jchar'][0][2]
	try:
	    source = { 'C' : 'GDAS', 'A' : 'SHIS' }[letter]
	except KeyError:
	    dbg('WARNING: UNKNOWN UNITS')
	    source = 'UNKNOWN'

	t5kw['profile_source'] = source	

	#_put all previous in numpy arrays
	for opt in opts[0]:
	    t5kw[opt] = array(t5kw[opt])
	for opt in opts[1]:
	    t5kw[opt] = array(t5kw[opt])

	#_kludge
	t5kw['co2amt'] = mean(t5kw['co2amt'])

###	#_rec 3.7: selected cross sections
###	opts = ['ixmols','iprfl','ixsbin']
###	line = f.readline().split()
###	[ t5kw.update({opt : int(line[i])}) for i, opt \
###	        in enumerate(opts) ]

	#_rest is static in current write_tap5

	#_put into fake profile attributes
	class dummy(object):
	    def __init__(self, profile_source='GDAS', surf_temp=None,
						water_vapor=None, h1=None, header=None,
	                    temperature=None, ozone=None, **kwargs):
			import re

			source = profile_source
			fmt = '{0:s}_{1:s}'
			surf_temp = temperature[0] if surf_temp is None else surf_temp
			setattr(self, fmt.format(source, 'temperature'), temperature)
			setattr(self,	fmt.format(source,	'ozone_mixing_ratio'),	ozone)
			setattr(self, fmt.format(source, 'sfc_temperature'), surf_temp)
			setattr(self, fmt.format(source, 'relative_humidity'), water_vapor)
			self.SHIS_altitude = h1	
	
			reg = re.compile('lat=([-\d.]+),lon=([-\d.]+),time=([\d.]+)')	
			res = reg.search(header)
			if res:
				self.SHIS_epoch = float(res.group(3))
				self.SHIS_longitude = float(res.group(2))
				self.SHIS_latitude = float(res.group(1))

	prof = dummy(**t5kw)

	#_close old tape 5
	f.close()

	#_pull out pressure
	pressure = t5kw['pressure']
	del t5kw['pressure']

	return prof, pressure, t5kw


def write_tape5(prof, 

	#_vertical pressure levels
	pressure,
	
	#_dummy levels.  If not None, arbitrary levels will
	# be inserted for LBL-ATM at value provided in km
	dummy_levs=None,	#_0.5 == 500 m 

	#_'GDAS' or 'SHIS'
	profile_source='GDAS',

	#_Location info
	latitude=0.0,
	longitude=0.0,
	epoch = 0.,
	
	#_ppm? no clue.
	co2amt=395.,

	label='', 
	#_DEFAULTS BELOW!!!
	#_SEE lblrtm_instructions.html#RECORD_1.2a FOR MORE DETAIL_#
	#_80 Character String, identify run
	header='Template TAPE5 File',

	#_HIRAC Options
	# 0 Not activated
	# 1 Voight Profile
	# 2/3 Not in LBLRTM
	# 4 Non Local Thermodynamic Equilibrium (TAPE4 requirements)
	# 9 Central line contribution omitted (functions 1-3)
	hirac=1, 

	#_LBLF4 Flag
	# 0 not activated, bound is 64 halfwidths
	# 1 lbl bound is 25cm-1 for all layer pressures (altitudes)
	# 2 lbl bound is "      for pres > .5mb
	#            and 5cm-1 for pres < .5mb 
	#	(saves comp time,small accuracy loss)
	lblf4=1,

	#_Continuum Flag
	# 0 No continuum calculated 
	# 1 All continua calc, include rayleigh
	# 2 H2O self not calc, all others calc
	# 3 H20 foreign not calc, all other calc
	# 4 H2O self and foreign not calc, all others calc
	# 5 Rayleigh ext not calc, all others calc
	# 6 Individual contiuum scale factors input (requires 1.2a rec) 
	cntnm=1,

	#_Aerosol flag (LOWTRN)
	# 0 no aerosols used
	# 1 internal LOWTRAN aerosol models (DEPRECATED)
	# 5 spectral optical depths by layer from file 'in_lblrtm_cld' 
	# 7 user defined aerosol models (DEPRECATED)
	# 9 precalculated aerosols (TAPE20 from previous aerosol run)
	aersl=0,

	#_SOMETHING
	# 0 optical depth only
	# 1 raidance and transmittance  
	# 2 solar radiance (requires prev calc OD or trasn and bin sr)
	# 3 radiance analytic Jacobian/derivative (requirements in doc)
	iemit=1,

	#_Scanning function
	# 0 none
	# 1 scanning function 
	# 2 interpolating procedure 
	# 3 FFT scan 
	iscan=0,

	#_0 off, 1 on
	filtr=0,	#_filter
	iplot=0,	#_plot lblb (PLTLBL)
	itest=0,	#_Test?
	iatmm=1,	#_LBLATM

	#_Merge options
	# 0-9, 12-18, 22-28, 32, 35-36, 40-43, 45-46
	# 0 KFILE rewound after each layer
	# > 1 kfile not rewound, see above link
	merge=0,

	#_Laser options
	# 0 no, 1 single laser frq, 2 multiple laser freq
	laser=0,

	#_Optical Depth layer controls
	# 0 normal, 1 uses exct dv for each layer and interpolates
	# 2 uses ecact calc dv, (IMRG must be 1), 3, 4
	ioddd=0,

	#_Cross section flag
	# No cross sections included in calc, cross sections in calc
	xsect=0,

	#_Number of OD values printed for beg and end of each panel
	# as the result of convolution for the current layer
	mpts=0,

	#_Number of values printed for the beg and end of each panel
	# as a result of merge of current layer withi previous layers
	npts=0,	

	#_cntnm=6 options
	xself=1.,xfrgn=1.,xco2c=1.,xo3cn=1.,xo2cn=1.,xn2cn=1.,xrayl=1.,

	#_input flag for solar rad calcs
	# 0 prev calc R&T from TAPE12 (default), 1 prev calc OD from T12
	# 2 upwell R&T from T12, downwell R&T from SOL.PATH.T2 and
	#	solar reflectance from SOL.REFLECTANCE
	# 3 input previously calc trans from T12
	inflag=0,

	#_ouput flag for sol rad calc (to TAPE13)
	# 0 output attenuated solar rad (default)
	# 1 output total rad
	# 2 output total rad, including effects of downwelling and thrm
	# reflection
	otflag=0,

	v1=649.,		#_beginning wavenumber
	v2=1614.,		#_end wavenumber (defaults to AIRS windows) 
	sample=4, 		#_number of sample points per mean halfwidth (1-4)
	dvset=0.,
	alfal0=0.04, 	#_average collision broadened halfwidth(cm-1/atm
	avmass=36, 		#_average molecular masss
	dptmin=2e-4,	#_minimum molecular OD below which lines reject
	dptfac=1e-3, 	#_mult factor something something
	ilnflg=0, 		#_binary record of line rejection info
					# 0 not rec, 1 write line reject to REJ1, REJ4
					# 2 read line rejection info from REJ1, REJ4
	dvout=0,		#_selected dv grid for OD something something 

	nmol_scal=0,	#_enables scaling of atm prof for selected species.
	nmol=7, 
	hmol=1, 
	xmol=0,
	# 1:  H2O  2:  CO2  3:     O3  4:  N2O  5:   CO  6:  CH4  7:    O2 8: NO  
	# 9:  SO2 10:  NO2 11:    NH3 12: HNO3 13:   OH 14:   HF 15:  HCL 16: HBR 
	#17:   HI 18:  CLO 19:    OCS 20: H2CO 21: HOCL 22:   N2 23:  HCN 24: CH3CL
	#25: H2O2 26: C2H2 27:   C2H6 28:  PH3 29: COF2 30:  SF6 31:  H2S 32: HCOOH
    #33:  HO2 34:    O 35: CLONO2 36:  NO+ 37: HOBR 38: C2H4 39: CH3OH

	#_emissivity parameters
	bndem0=.985,   	#_coefficients
	bndem1=0.,   	# emissivity = bndem0+bndem1*V+bndem2*V**2
	bndem2=0.,
	bndrf0=0.,   	#_reflectivity
	bndrf1=0.,
	bndrf2=0.,
	sfcrfc='s',		#_'s' for specular, 'l' for lambertian

	#_pathname for optical depth files (var+##)
	pthodl='ODtest_',
	nodlay=1000,

	#_LBLATM params
	# 0 user supplied atm prof
	# 1 tropical
	# 2 midlat summer
	# 3 midlat winter
	# 4 subarctic summer
	# 5 subarctic winter
	# 6 US standard 1976
	model=0,		#_user supplied atmospheric profile
	aflag=None,		#_no idea, see Holz code for why??

	#_type of path
	# 1 horizontal path, contant p/T, uses rec 3.2H
	# 2 slant path from H1->H2, uses rec 3.2
	# 3 slant path from H1 to HSPACE, use rec 3.2
	itype=2,

	#_layering for LBLRTM
	# 0 	layers are generated internally
	# >0	is the number of layer boundaries to use
	# <0	set to pressure levs
	ibmax=0,

	#_zeroes out absorber amounts which are less than 0.1% total
	nozero=1,	#_suppress zeroing
	noprnt=1,	#_0 full out, 1 selects short printout
	ipunch=1,	#_0 layer data not written, 1 write to ITAPE7, 2 to AJ
	ifxtyp=0,	#_leave these blank on tape 7
	munits=0,	#_0 write molecular column amounts or 1 mixing rats to TAPE7
	re=0,		#_select radius of earth 
	hspace=100,	#_altitude def for space (km)
	ref_lat=0,	#_latitude of lcoation of calc (degrees)

	h1=100,		#_height of observer (ibmax neg, mbar)
	h2=0,		#_for itype2, end point altitude (km)
	angle=180.,	#_zenirh angle at h1. 0 for uplooking, 180 for downlooking

	ixmols=3,	#_number of x-section molecules to be inputed (max 35)
	iprfl=0,	#_0 user input profile, 1 standard profile from LBLATM
	ixsbin=0,	#_to deselect p convol of x-sect. 0 xsec con, 1 xsect not

	clddef=None,	#_if this lbldis cloud layer option is passed, adjust levs

	**kwargs ):
	'''
	lblrtm_instructions.html#RECORD_1.2a FOR MORE DETAIL_#
	generate tape 5 file in accordance with whatever
	'''
	from numpy import append, arange, array
	from libtools import strictly_increasing
	from libtools import unique, epoch2dtg
	from libgeo import p2z

	def add_levels(arr, dv):
		''' make sure no space larger than dv exists in arr '''
		from numpy import append, arange, array, diff
		arr = array(arr)
    
		#_loop until all gaps are filled
		while (diff(arr) > dv+1e-8).any():

			#_get index of first gap that needs filling
			tmp = arr[append(diff(arr) > dv+1e-8, [False])][0]
        
			#_loop over sections that have gap
			#_add value to end, sort
			arr = append(arr, tmp+dv)
			arr.sort()
        
		return arr

	try:
		#_set which profile to use
		prof_t = getattr(prof, '{0}_temperature'.format(profile_source))
		prof_rh = getattr(prof, '{0}_relative_humidity'.format(profile_source))
		ozone = getattr(prof, '{0}_ozone_mixing_ratio'.format(profile_source))
		sfc_temp = getattr(prof, '{0:s}_sfc_temperature'.format(profile_source))
		jchar = {'GDAS' : 'C', 'SHIS' : 'A', 'ECMWF' : 'C'}	#_ECMWF SFC ONLY
		o3_jchar = jchar[profile_source]
	except AttributeError:
		raise RuntimeError, 'implement'
		#_check if using standard atmospheric p

	#_do a really ugly conversation of pressure to altitude
	altitude = p2z(pressure) / 1e3
	if altitude.size != prof_t.size:
		dbg('vertical coordinate and profile mismatch')
		return False

	#_options dependency check
	##	check_opts(**kwargs)
	## need to add all keywords to kwargs.var

	#_open TAPE5 file
	fname = 'TAPE5'
	dbg(('writing', fname), 5)
	fid = open(fname,'w')

	#_rec 1.1:
	header = 'time={0}, lat={1:5.2f}, lon={2:5.2f}'.format(
			epoch2dtg(epoch, full=True), latitude, longitude) 
	fid.write('$%-79s\n' % (header))			#_line 1: comment

	#_rec 1.2
	l = ' HI=%1i F4=%1i CN=%1i AE=%1i EM=%1i SC=%1i FI=%1i PL=%1i ' % \
			(hirac,lblf4,cntnm,aersl,iemit,iscan,filtr,iplot)
	l += 'TS=%1i AM=%1i MG%2s LA=%1i OD=%1i XS=%1i %4i %4i\n' % \
			(itest,iatmm,str(merge).rjust(2,'='),laser,ioddd,xsect,mpts,npts)
	fid.write(l)								#_line 2: settings

	#_rec 1.2a: only matters if cntnm=6, continuum multipliers
	if cntnm == 6:
		fmt = '%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n'
		cont_mult = {	2:	[1,1,0,0,0,0,0],		#_wco
						3:	[0,0,0,0,0,0,0],		#_wnc
						4:	[1,1,0,1,0,0,0],		#_wvo
						5:	[0,0,0,1,0,0,0],		#_ozo
						6:	[0,0,1,0,1,1,1]	}		#_dry
		line = fmt % cont_mult[cntnm]
		fid.write(line)								#_line 3: cont multipliers
	##		f.write('%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n' 
	##		% (xself,xfrgn,xco2c,xo3cn,xo2cn,xn2cn,xrayl)) 

	#_rec 1.2.1: bunch of notes on this for TOTAL upwelling/downwelling
	if iemit == 2:
		fid.write('%i5%i5  %i3\n' % (inflag,otflag,julday))

	#_rec 1.3: wavenumber bounds, et cetera
	if hirac>0 or aersl>0 or iemit==1 or iatmm==1 or laser>0:
		l = '%10.3f%10.3f\n' % (v1,v2)

##		l = '%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f' % \
##			(v1,v2,sample,dvset,alfal0,avmass,dptmin,dptfac)
##		l += '    %1i     %10.3f   %2i\n' % (ilnflg,dvout,nmol_scal)
		fid.write(l)

	#_rec 1.3.a+b (not in holz)
	if nmol_scal > 0:
		fid.write(str(hmol)+'\n')
		fid.write('%15.7f'*8+'\n' % xmol)
        
    #_rec 1.4: temperature and emissivity boundary
	if iemit == 1 or (iemit == 2 and otflag == 2):
		l = '%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f    %1s\n' % \
		(	sfc_temp, bndem0, bndem1, bndem2, 
			bndrf0, bndrf1, bndrf2, sfcrfc)
		fid.write(l)    

	#_direct input of emissivity from EMISSIVITY file if <0
	if bndem0 < 0: write_EMISSIVITY(v1, v2, dv, nlim)
	#_direct input of boundary REFLECTIVITY file if <0
	if bndrf0 < 0: write_REFLECTIVITY(**kwargs)

	#_rec 1.5: Analytic Jacobian calculation
	if iemit == 3 and merge in [40,41,42,43]: err('Jacobian not implemented!')

	#_rec 1.6a: not in holz
	#_These options requires separate optical depth files for each layer
	if merge in [35,36,40,41,45,46]:
		f.write('%55a %4i' % (pthodl, nodlay))

	#_rec 2.1: not in holz, not implemented
	if iatmm == 0: err('Not implemented (RECORD_2.1)')

	#_rec 3.1: LBLATM
	vbar = (v1+v2)/2 				#_get mean wavenumber

	#_make temporary arrays of LBLATM quantities
	temperature = prof_t 
	water_vapor = prof_rh 

	#_limit profile to where temperature is not missing
	idx			= temperature != -9999
	temperature = temperature[idx]
	pressure	= pressure[idx]
	water_vapor = water_vapor[idx]
	ozone		= ozone[idx]

	#_figure out what's going on with altitude values
	if type(altitude) is int:
		from libgeo import p2z #_SUPER SLOPPY RIGHT NOW
		altitude = p2z(pressure)
		units = -1
	else:
		altitude = altitude.copy()[idx]
		units = 1 
	
	#_limit to below flight level
	idx			= altitude <= h1 
	temperature = temperature[idx]
	pressure	= pressure[idx]
	water_vapor = water_vapor[idx]
	ozone		= ozone[idx]
	altitude	= altitude[idx]
	paltitude	= altitude.copy()	#_keep a copy of levels that we have 
									# column data for separate from all
									# levels being computed
	nlev		= len(altitude)

	dbg('INCORPORATE THE ABILITY TO USE NEG IBMAX TO DENOTE PRESSURE LEVS')
	dbg('BY FIXING LBLDIS THING BELOW')

	#_check for lbldis levels, LBL-RTM will interpolate between these
	# based upon the LBL-ATM values below.  IBMAX does not have to equal
	# nlev
	if clddef is not None and dummy_levs is None:
		add_levs = []
		for cld in clddef:
			add_levs.append(cld['z'])
			add_levs.append(cld['z_top'])

		#_add these to level boundaries
		altitude = append(altitude, add_levs)
		altitude = array(unique(altitude))

	#_when not using clddef layer to define where to split levs
	# possible to use arbitrary breakdown
	#_if desired, never let space between leveels be greater than dummy_levs
	elif dummy_levs is not None:
		altitude = add_levels(altitude, dummy_levs)

	if len(altitude) > 120:
		from time import sleep
		dbg('WARNING! Max LBLDIS vertical levels set to 120')
		sleep(10)

	#_number of boundaries	
	ibmax = len(altitude)

	l = '%5i%5i%5i%5i%5i%5i' % (model, itype, ibmax, nozero, noprnt, nmol)
	l+= '%5i%2i %2i%10.3f%10.3f%10.3f\n' % (ipunch,ifxtyp,munits,re,hspace,vbar)
	fid.write(l)

	#_rec 3.2: path parameters
	l = '%10.3f%10.3f%10.3f\n' % (h1, h2, angle)
	fid.write(l)

	#_rec 3.3b: write altitudes of layer boundaries, eight per line
	for group in setup_groups(altitude, 8):
		fmt		= '%10.3f'*group.size
		fmt		+= '\n'
		line 	= fmt % tuple(group)
		fid.write(line) 
	else:
		pass

	#_rec 3.4: atm prof boundaries (start with model==0, user defined)
	if model == 0:
		l = '%5i levels in user defined profile\n' % nlev
		fid.write(l)

		#_calculate standard profile to use for molecules not provided
		# 1 tropical
		# 2 midlat summer
		# 3 midlat winter
		# 4 subarctic summer
		# 5 subarctic winter
		# 6 US standard 1976
		month = int(epoch2dtg(epoch)[4:6])
		if abs(latitude) < 20:
			aflag = 1
		elif latitude <= -20 and latitude >= -60:#_SH MIDLAT
			aflag = 2 if month > 9 or month <= 3 else 3
		elif latitude < -60:							#_SH SUBARC
			aflag = 4 if month > 9 or month <= 3 else 5
		elif latitude >= 20 and latitude <= 60:	#_NH MIDLAT
			aflag = 3 if month > 9 or month <= 3 else 2
		elif latitude > 60:							#_NH SUBARC
			aflag = 5 if month > 9 or month <= 3 else 4
	
		#_rec 3.5: level definitions
		afgstr = [str(aflag)] * 8
		
		#_check that profile is in ascending order
		if not strictly_increasing(paltitude):
			dbg('WARNING: Profile not strictly increasing')
			return False
		
		#_rec 3.6.1-N
		for i in arange(nlev):		
			afgstr = [str(aflag)] * 8

			#_[H20, CO2, O3, N20, CO, CH4, O2]	
			# afgstr[0] == H20, et cetera
			if cntnm == 1:
				if hasattr(prof, 'co'):	
					cntnmvar = [water_vapor[i], co2amt, 0.0, 
								prof.carbon_monoxide[i], 0., ' ']
					afgstr[4:6]	= ['A']*2
					coflag = 1
				else:
					cntnmvar = [water_vapor[i], co2amt, ozone[i],
								0., 0., ' ']
				afgstr[0] = 'C'	#_ozone units, C=>g/kg

			elif cntnm == 2:
				cntnmvar = [water_vapor[i], 0., 0., 0., 0., ' ']
				afgstr[4:6]	= ['A']*2
			
			elif cntnm == 3:
				cntnmvar = [water_vapor[i], 0., 0., 0., 0., ' ']
				afgstr[4:6]	= ['A']*2
			
			elif cntnm == 4:
				cntnmvar = [water_vapor[i], 0.0, ozone[i], 0., 0., ' ']
				afgstr[4:6]	= ['A']*2
				afgstr[0]	= 'A'
			
			elif cntnm == 5:
				cntnmvar = [water_vapor[i], co2amt, ozone[i], 0., 0., ' ']
				afgstr[4:6]	= ['A']*2
				afgstr[0]	= 'C'	#_ozone units, C==g/kg
			
			elif cntnm == 6:
				cntnmvar = [0., co2amt, 0., 0., 0., ' ']
				afgstr[4:6] = ['A']*2
			'''
			HAR = 1-6         - default to value for specified model atmosphere
              = " ",A         - volume mixing ratio (ppmv):
              = B             - number density (cm-3)
              = C             - mass mixing ratio (gm/kg)
              = D             - mass density (gm m-3)
              = E             - partial pressure (mb)
              = F             - dew point temp (K) *H2O only*
              = G             - dew point temp (C) *H2O only*
              = H             - relative humidity (percent) *H2O only*
              = I             - available for user definition
			GDAS == C*1000, SHIS == A
			'''
			o3_jchar = o3_jchar if ozone[i] >= 0 else '1'
			afgstr = 'HA{0:1s}111AA11'.format(o3_jchar)
##			afgstr = 'CA'+''.join(afgstr)	if prof.jchar is None	\
##											else prof.jchar
			val	= (paltitude[i], pressure[i], temperature[i],
					'     AA   ', afgstr)
			l	= '%10.3f%10.4f%10.3f%10s%8s\n' % tuple(val)
			fid.write(l)

			l 	= '%10.3f%10.3f%10.3f%10.3f%10.3f%30s\n' % tuple(cntnmvar)
			fid.write(l)
###			fid.write('%40s\n'%' ')		#_???
	
	else:	#_use standard model atmosphere.  Do I need anything else here?
		pass
	
	#_rec 3.7: selected cross sections
	l = '%5i%5i%5i' % (ixmols,iprfl,ixsbin)
	fid.write(l)
	
	#_THIS DOES NOT MATCH THE DOCUMENTATION
	#_rec 3.7.1: cross section molecules to be used (this is where the 3 is frm)
	fid.write('CCL4      F11       F12 \n')	#_xsname
	fid.write('%5i%5i  \n' % (2,0))			#_???
	
	#_rec 3.???: 'CFC11 and 12 are half of the official values
	fid.write('%10.3f     AAA\n' % altitude.min())
	fid.write('%10.3e%10.3e%10.3e\n' % (1.105e-04,1.343e-04,2.527e-04))
	fid.write('%10.3f     AAA\n' % altitude.max())
	fid.write('%10.3e%10.3e%10.3e\n' % (1.105e-04,1.343e-04,2.527e-04))
	fid.write('-1.\n')

##	check_opts(prof, angle=angle)

	#_eof
	fid.write('%'*7)

	#_close file handle
	fid.close()
	
	return fname
    
'''
 ---------------------------------------------------------------------------------------------------
|                                FORMAT FOR 'EMISSION' and 'REFLECTION' FILES                       |
|                                --------------------------------------------                       |
|                                                                                                   |
|      Record 1:                                                                                    |
|                V1,      V2,      DV,      NLIM                                                    |
|                                                                                                   |
|              1-10,   11-20,   21-30,     35-40                                                    |
|                                                                                                   |
|             E10.3,   E10.3,   E10.3,     5X,I5                                                    |
|                                                                                                   |
|                                                                                                   |
|             V1   beginning wavenumber of input (should be less than V1 on RECORD 1.3 from TAPE5)  |
|                                                                                                   |
|             V2   ending wavenumber of input (should be greater than V2 on RECORD 1.3 from TAPE5)  |
|                                                                                                   |
|             DV   spectral spacing of input (cannot be less than 0.5 cm-1)                         |
|                                                                                                   |
|           NLIM   number of points in input (must be less than or equal to 4040)                   |
|                                                                                                   |
|                                                                                                   |
|      Record 2(i), i=1,NLIM:                                                                       |
|                                                                                                   |
|               ZDATA                                                                               |
|                                                                                                   |
|                1-10                                                                               |
|                                                                                                   |
|               E10.3                                                                               |
|                                                                                                   |
|                                                                                                   |
|          ZDATA   Spectrally dependent boundary emission or reflection value                       |
|                         (one spectral value per card)                                             |
|                                                                                                   |
 ---------------------------------------------------------------------------------------------------
'''


def write_emissivity(v1=None,v2=None,dv=None,nlim=None,**kwargs):
	''' write emissivity file... '''
	import numpy as np

	if None not in (v1,v2,dv,nlim):
	    pass
	elif v2 == None and None not in (v2,dv,nlim):
	    v2 = (np.arange(nlim)*dv + v1)[-1]
	elif dv == None and None not in (v1,v2,nlim):
	    dv = (v2-v1)/(nlim-1.) #np.linspace(v1,v2,nlim)
	elif nlim == None and None not in (v1,v2,dv):
	    nlim = np.arange(v1,v2,dv).size
	else:
	    err('Requires at least three of v1,v2,dv,nlim to be defined')
    
	wavelengths = np.arange(nlim)*dv + v1
    
	f = open('EMISSIVITY','w')
	f.write('%10.3e%10.3e%10.3e     %i5' % (v1,v2,dv,nlim))

	for n in xrange(nlim):
	    l = wavelengths[n]
	    emis = '111.'
	    f.write('%10.3f%10.1f'%(emis,l))
	    print 'WRITING DUMMY EMISSIVITY FILE'
	f.close()

def write_reflectivity(): pass


def link_tape3(tape3_loc=os.path.join(DIR_PROD, 'TAPE3_default'), dest='.',
	**kwargs):
	''' if we ever get bigger '''
	sym_file = os.path.join(dest,'TAPE3')
	dbg(('linking', tape3_loc, sym_file)) 
	os.symlink(tape3_loc, sym_file)


def read_tape6(**kwargs):
	'''
	LBLRTM Output File
	'''
	pass
	
def read_tape7(**kwargs):
	'''
	Molecular column amounts from LBLATM (IATM=1)
	'''
	pass
	
def read_tape9(**kwargs):
	'''
	File of effective lines for LBLF4 created by LINF4
	'''
	pass
	
def read_tape10(**kwargs):
	'''
	Optical depth results from LBL calculations
	IMRG=0,	last layer
	IMRG=1,	layer by layer 
	'''
	pass
	
def read_tape11(fname='TAPE11', **kwargs):
	'''
	Spectral results from SCANFN and INTRPL
	JEMIT=-1:	absorption
	JEMIT= 0:	transmittance
	JEMIT= 1:	radiance
	'''
	from libtools import rec_read, merge
	from libgeo import wavenumber2frequency, planck_inv
	from numpy import linspace, arange, hstack, pi, array, append
	
	fid = open(fname, 'rb')		#_open binary

	#_get header? Should be a bunch of stuff about endian, but LATER
	try:
		n = 352					#_uncertain, buff for binary header	
		fmtt = 'f4'
		test = rec_read(fid, [('test', str(n)+fmtt)])
	except MemoryError:
		fid.seek(0)
		n = 132					#_uncertain, buff for binary header	
		fmtt = 'f8'
		test = rec_read(fid, [('test', str(n)+fmtt)])

	dtype = [('v1', fmtt),('v2',fmtt),('dv','f4'),('npt','i4')]

	rads = array([]); wavs = array([]); tran = array([])
	eof = False								#_loop over binary reads
	while not eof:
		v1,v2,dv,npt = rec_read(fid, dtype)	#_get band information

		#_read brightness temperature and flux
		radiances	= rec_read(fid,[('data', str(npt)+'f4')])
		trans		= rec_read(fid,[('data', str(npt)+'f4')])
		wavenum 	= arange(v1,v2,dv)

		rads = append(rads,radiances['data'])
		tran = append(tran,trans['data'])
		wavs = append(wavs,wavenum)

		if npt != 2400: eof = True			#_is there not a more elegant way?

	fid.close()								#_close file handle
	
	#_convert to W/m2/sr/cm-1
	rads *= 1e4
	
	wavs = (wavs[1:] + wavs[:-1])/2 		#_get midpoint
	data = Spectrum(radiances=rads[:wavs.size], wavenumber=wavs,
		transmissivity=tran[:wavs.size], label='LBL-RTM', **kwargs)
##	data.radiances *= 10**7					#_convert W->mW, m-2->cm-2
	return data

def read_tape12(**kwargs):
	'''
	monochromatic results
	IEMIT=0:	optical depth
	IEMIT=1:	radiance/transmittance (includes aerosols for IAERSL=1)
	'''
	pass


def airs_channels(fname='/'.join((	os.environ['PRODUCTS'], 
									'airs/AIRS_channel.dat'))):
	''' read in AIRS channel wavenumber data, ascending '''
	from numpy import genfromtxt
	return genfromtxt(fname)


def read_dir(dir_in='/'.join((os.environ['PRODUCTS'],'COL')),**kwargs):
	''' reads collocated data for initial testing '''
	from glob import glob
	from libtools import merge
	from netCDF4 import Dataset
	from numpy import recarray, ndarray
	
	dtype 	= [	('cloud','f8'),('latitude','f8'),('longitude','f8')]
	CALIOP 	= recarray((0,),dtype=dtype)
	dtype 	= [	('latitude','f8'),		('longitude','f8'),
				('temperature',ndarray),('water_vapor',ndarray),
				('ozone',ndarray),		('time','f8'),
				('altitude',ndarray),	('sfc_temperature','f8')	]
	GDAS	= recarray((0,),dtype=dtype)
	dtype	= [	('spectrum', Spectrum),('latitude','f8'),
				('longitude','f8'),('time','f8')]
	AIRS	= recarray((0,),dtype=dtype)
	
	#_loop over available files
	for file in glob(dir_in+'/COL*AIRS*hdf')[30:40]:
#		dbg(file)
		hdf 	= Dataset(file,'r',format='HDF4')
		CALIOP 	= merge((CALIOP, read_CALIOP(hdf)))
		AIRS 	= merge((AIRS, read_AIRS(hdf)))
		GDAS 	= merge((GDAS, read_GDAS(hdf)))
##		MODIS 	= merge((CALIOP, read_MODIS(hdf)))
		hdf.close()		#_close filehandle

	return AIRS, CALIOP, GDAS #,MODIS
			
def read_airs(hdf, **kwargs):
	''' 
	read all AIRS data from collocated data files in a directory
	dir_in,	str,	path to input
	
	RETURN
	data,	np.recarray,	contains AIRS radiances, location, time
	'''
	from glob import glob
	from netCDF4 import Dataset
	from numpy import recarray, ndarray, arange, tile
	from numpy.ma import masked_where
	
	dtyp = [('spectrum',	Spectrum),
			('latitude',	'f8'),
			('longitude',	'f8'),
			('time',		'f8')	]
	data = recarray((0,),dtype=dtyp)
	dbg(hdf)
	if type(hdf) != Dataset: 					#_allow hand pass
		hdf = Dataset(hdf, 'r', format='HDF4')	#_open filehandle
		close = True
	else:
		close = False
	
	#_get AIRS wavenumbers
	wave = airs_channels()
	mask = AIRS_ancillary().l2_ignore
	
	rad = hdf.variables['AIRS_Radiances'][:].squeeze()	#_mW/m2/cm-1/sr (conver)
	lat = hdf.variables['AIRS_Latitude'][:].squeeze()
	lon = hdf.variables['AIRS_Longitude'][:].squeeze()
	tim = hdf.variables['AIRS_Time'][:].squeeze()
	
##	if rad.ndim == 3:
##		rad = rad.squeeze()
##		lat = lat.squeeze()
##		lon = lon.squeeze()
##		tim = tim.squeeze()

	#_if only one view in the file, then squeeze kind of kills it
##	if rad.ndim == 1: rad = sing(rad)

	#_never mind, if only one view, assume borked
	if rad.ndim == 1: return False

	#_convert to W/m2/cm-1/sr
	idx = rad != -9999
	rad[idx] /= 1e3
	
	#_mask where rad err too great
	mask = tile(mask, (rad.shape[0], 1))	#_masked_where won't broadcast
	rad	= masked_where(mask, rad)
	
	size = lat.size
	data.resize(size)					#_expand record array
	arg = { 'wavenumber' : wave, 'label' : 'AIRS' }
	try:
		for n in arange(size):
			arg.update({'radiances' : rad[n,:],
						'lat'		: lat[n],
						'lon'		: lon[n],
						'time'		: tim[n]	})
			data.spectrum[n-size] = Spectrum(**arg)
	except IndexError:
		for n in arange(size):
			arg.update({'radiances' : rad[:],
						'lat'		: lat[n],
						'lon'		: lon[n],
						'time'		: tim[n]	})
			data.spectrum[n-size] = Spectrum(**arg)
			
	data.latitude[:]	= lat
	data.longitude[:] 	= lon
	data.time[:] 		= tim
	
	if close: hdf.close()					#_close filehandle
	return data								#_return recarray with data

def read_caliop(hdf, cloud_var='CALIPSO_Layer_Top_Mean_5km', **kwargs):
	''' 
	read all CALIOP data from collocated data files in a directory
	dir_in,	str,	path to input
	
	RETURN
	data,	np.recarray,	contains 
	'''
	from glob import glob
	from netCDF4 import Dataset
	from numpy import recarray, ndarray, arange
	from numpy.ma import masked_outside, masked_where, masked_array
	
	dtyp = [('cloud','f8'),('latitude','f8'),('longitude','f8')]
	data = recarray((0,),dtype=dtyp)

	if not type(hdf) == Dataset: 				#_allow handle to be passed
		hdf = Dataset(hdf,'r',format='HDF4')	#_open filehandle
		close = True
	else:
		close = False
	lat = hdf.variables['CALIPSO_Latitude_1km'][:].squeeze()
	lon = hdf.variables['CALIPSO_Longitude_1km'][:].squeeze()
	clf = hdf.variables[cloud_var][:].squeeze()

	#_CALIPSO_Fraction_Cloudy_1km 5km
	###_CALISPO_Layer_Top_Mean_5km
	###_CALIPSO_Freature_Optical_Depth_532_Mean_5km

##	(1) cloud optical thickness of top layer =0 or 
##		Column_Optical_Depth_Cloud_532 = 0 or/and
## 	(2) cloud  layer top altitude equal to -999 (i.e., 
##		cloud top layer altitude is not retrieved) or/and
##	(3) CAD score = negative (are cloud free or no cloud feature detected) .
##		CAD score = -70 to -100 is aerosol and CAD score = 70 to 100 ensure 
##		cloud.
	
	size = clf.size
	data.resize(data.size+size)				#_expand record array
	data.cloud[:] 		= clf
	data.latitude[:] 	= lat.max(axis=0)	#_just don't
	data.longitude[:] 	= lon.max(axis=0)	#_just don't use these

	#_mask FillValues
##	data.cloud = masked_where(data.cloud == -8888, data.cloud)
##	data.cloud = masked_where(data.cloud == -101, data.cloud)

	if close: hdf.close()					#_close filehandle
	return data								#_return recarray with data


def read_gdas_raw(fname, variable='BLAG', **kwargs):
	pass
	import pygrib as grb


def read_gdas(hdf, **kwargs):
	''' 
	read all GDAS data from collocated data files in a directory
	dir_in,	str,	path to input
	
	RETURN
	data,	np.recarray,	contains GDAS radiances, location, time
	'''
	from glob import glob
	from netCDF4 import Dataset
	from numpy import recarray, ndarray, array, arange
	from libgeo import rh2wv, kg2ppmv, Z2z
	from libtools import dtg2epoch	
##	dtyp = [	('latitude','f8'),		('longitude','f8'),
##				('temperature',ndarray),('water_vapor',ndarray),
##				('ozone',ndarray),		('time','f8'),
##				('altitude',ndarray),	('sfc_temperature','f8'),
##				('land_fraction', 'f4')	]
##	data = recarray((0,),dtype=dtyp)
	data = Profile2() 
	
	if not type(hdf) == Dataset: 				#_allow handle to be passed
		hdf = Dataset(hdf,'r',format='HDF4')	#_open filehandle
		close = True
	else:
		close = False
	
	lat = hdf.variables['GDAS_Latitude'][:].squeeze()
	lon = hdf.variables['GDAS_Longitude'][:].squeeze()
	pre = hdf.variables['GDAS_Pressure'][:].squeeze()
	tem = hdf.variables['GDAS_Temperature'][:].squeeze()
	wvc = hdf.variables['GDAS_Relative_Humidity'][:].squeeze()
	ozn = hdf.variables['GDAS_Ozone_Mixing_Ratio'][:].squeeze()
	gpt = hdf.variables['GDAS_Geopotential_Height'][:].squeeze()
	sfc = hdf.variables['GDAS_Surface_Temperature'][:].squeeze()
	frc = hdf.variables['GDAS_Land_Fraction'][:].squeeze()
	tim = hdf.variables['AIRS_Time'][:].squeeze()	#_is this the same...
	epo = tim + dtg2epoch('1993010100')

	size = lat.size
	data.resize(size)			#_expand record array
	arg = {'attrv':[pre], 'attrn':['pressure']}
	for n in arange(size):
		try:
			data[n].temperature	= var(tem[:,n], **arg)
			data[n].water_vapor	= var(rh2wv(wvc[:,n], tem[:,n], pre*100), **arg)
			data[n].ozone		= var(ozn[:,n]*1e3, **arg)
			data[n].altitude 	= Z2z(gpt[:,n]) / 1000.
		except IndexError:
			data[n].temperature = tem[:]
			data[n].water_vapor = rh2wv(wvc[:], tem[:], pre*100)
			data[n].ozone 		= ozn[:]*1e3
			data[n].altitude 	= Z2z(gpt[:]) / 1000.

		data[n].pressure	= pre
		data[n].latitude	= lat[n]
		data[n].longitude	= lon[n]
		data[n].time		= tim[n],
		data[n].epoch		= epo[n]
		data[n].land_fraction	= frc[n]
		data[n].sfc_temperature	= sfc[n]

	if close: hdf.close()			#_close filehandle
	return data						#_return recarray with data


def match_lbldis_out(data, dir_gdas=os.path.join(os.environ['PRODUCTS'],'gdas'),
	**kwargs):
	''' return airs Spectrum object for this lbldis output '''

	file_gdas = os.path.join(dir_gdas, data.label + '.hdf')
	airs = read_airs(file_gdas)
	
	idx = abs(airs.time - data.time).argmin()
	return airs[idx].spectrum

	
def convolve2airs(data, f0=None, f1=None, interp_type='linear', **kwargs):
	'''
	do convolutions from tabulated spectral response functions
	
	rads_in		ndarray		Mx1 radiance data, (LBL?)
	freq_in		ndarray		Mx1 frequencies for radiances (also LBL?)
	channels	ndarray		Kx1	channel numbers for output (AIRS?) (freq?)
	
	kind		string		type of interpolation, see scipy.interpolate
	interp_type	string		type of interpolation for srfvalues
	'''
	from numpy import dot, arange, ones, zeros, diff, array, ceil, floor, append
	from scipy.interpolate import interp1d
	dbg('convolving...', 3)
	rads_in, freq_in = data.radiances, data.wavenumber
	
	#_check shape
	if rads_in.ndim == 1: rads_in = sing(rads_in).T

	#_get spectral response information
	srf = Srf(**kwargs)	

	#_subset to range specified in channels
##	f 		= interp1d(arange(srf.chanid.size), srf.chanid, kind='nearest')
##	chanid 	= srf.chanid[channels]
##	freq	= srf.freq[channels]
##	srfval	= srf.srfval[channels,:]
##	width	= srf.width[channels]
	nchan_airs, nfwhm_airs 	= srf.srfval.shape
	nchan_lbl, nscene 		= rads_in.shape
	
	#_create a nchan_airs x nfwhm_airs array that contains
	# all the frequencies that are influenced by each channel
	#_so srfreq[0,:] will be a small range of freqs that incoming
	# radiation on chanid[0] will be spread along
	srfreq = dot(sing(srf.width).T, sing(srf.fwgrid)) 
	srfreq += dot(sing(srf.freq).T, ones((1, nfwhm_airs)))

	#_check that they are unformly spaced and increasing
	dfin 	= diff(freq_in)
	dv 		= dfin[0]
	
	f0 		= freq_in[0] if f0 is None else f0
	f1 		= freq_in[-1] if f1 is None else f1
	nfin 	= len(freq_in)
	
	rads_out = zeros((nchan_airs, nscene))
	freq_out = srf.freq.copy()
	
	#_truncate what we can do srfreq(nchan_airs, 471)
	#_The ones dropped off are where the SRF is either entirely outside
	# of the channels within freq_in, so they do not influence results
	loop_idx = arange(nchan_airs)
	tmp = append(	sing(srfreq[:,0] - f0 >= 0).T,
					sing(f1 - srfreq[:,-1] >= 0).T, axis=1).all(axis=1)
	loop_idx = loop_idx[tmp]
	
	#_any spots where the SRF goes outside limits, fill (probably works??)
	srfreq[srfreq[:, 0] < f0-dv,  0] = f0
	srfreq[srfreq[:,-1] > f1+dv, -1] = f1
	
	#_apply SRF to each channel 
	for j in loop_idx:
		#_get current wavenumber span for this channel
		v0 = srfreq[j,0]
		v1 = srfreq[j,-1]

		# find the indices of the current SRF in freq_in 
		v0ind = abs(freq_in - v0).argmin()
		v1ind = abs(freq_in - v1).argmin()
		if freq_in[v0ind] < v0: v0ind += 1	#_keep them in range of srfreq
		if freq_in[v1ind] > v1: v1ind -= 1
		
		#_locations 
		vidx = arange(v0ind, v1ind)
		
		# interpolate the SRF to a subinterval of the fin grid
		# s1 = interp1(srfreq(j,:), srfval(j,:), fin(v1ind:v2ind), 'spline');
		f 	= interp1d(srfreq[j,:], srf.srfval[j,:], interp_type) 
		s1 	= f(freq_in[vidx]) 	#_get how much each of the desired freqs are
								# influenced 

		# normalize the SRF
		s1 /= s1.sum()			#_cannot influence more than 100%, dawg

		# apply the SRF rin[vrnge] is the array of radiances 
		# that will influence this channel (j) in rads_out with
		# s1 describing how much
		rads_out[j,:] = dot(s1, rads_in[vidx]).copy()
		
	label = 'airs-conv_' + data.label

	rads_out = rads_out.squeeze()

	return Spectrum(radiances=rads_out, wavenumber=freq_out, label=label,
					inherit=data, **kwargs)


def interp_index(x, y, x_new, **kwargs):
	'''
	NOT DONE!  Make axis splined interpolation to speed up above
	
	create function to interpolate over a single index of 
	a multi-dimensional array
	
	'''
	from scipy.interpolate._fitpack import _bspleval
	from scipy.interpolate import interp1d
	
	f = interp1d(x, y, axis=-1)

def for_bob_loop(dir_in='/'.join((os.environ['PRODUCTS'],'COL')), **kwargs):
	from glob import glob
	from netCDF4 import Dataset
	from os import fork, waitpid, chdir
	from libtools import setup_groups
	from sys import platform
	from random import randint
	from shutil import copy
	
	groups = setup_groups(glob(dir_in+'/COL*AIRS*hdf'), **kwargs)
	osx = True if platform == 'darwin' else False
	
	for group in groups:
		children = []
		pid = randint(1,10000) if osx else fork()
		
		if not osx and pid != 0: children.append(pid)
		elif osx or pid == 0:
			for fname in group:
				#_skip if file already done
				fname_out = '/'.join((os.environ['PRODUCTS'], 'spectra',
					fname.split('/')[-1].replace('.hdf', '.nc')))

				#_skip the finished
				if os.path.exists(fname_out):
					dbg((fname_out, 'complete!'))
					continue

				dbg(fname)
				hdf 	= Dataset(fname, 'r', format='HDF4')
				CALIOP 	= read_CALIOP(hdf)
				AIRS 	= read_AIRS(hdf)
				GDAS 	= read_GDAS(hdf)
				hdf.close()

				for_bob(AIRS, CALIOP, GDAS, fname=fname, pid=pid, **kwargs)
	
				if not osx: os._exit(0)
	
		if not osx:
			for kid in children: waitpid(kid, 0)
		
		
#############################################################################_80
#_CLASSES_######################################################################
################################################################################


class Spectrum(object):
	inherit_v = ['radiances', 'wavenumber', 'frequency', 'Tb', 'transmissivity',
				'flux', 'label', 'time', 'latitude', 'longitude', 'cloud']
	def __init__(self, radiances=None, wavenumber=None, frequency=None, Tb=None,
		transmissivity=None, flux=None, label='s', time=-9999, latitude=-9999, 
		longitude=-9999, spectrum_file=False, cloud=-9999, inherit=False,
		**kwargs):
		from numpy import recarray, diff, where
		from libgeo import wavenumber2frequency, planck_inv
		from numpy.ma import masked_outside

		rad_range = [-1e-6,9.e7] 		#_outside these limits, mask data

		if radiances != None:
			self.radiances	= masked_outside(	radiances, rad_range[0], 
												rad_range[1])
			self.wavenumber = wavenumber				#_cm-1
			self.Tb			= planck_inv(self.radiances, wavenumber*100,
								domain='wavenumber')	#_K
			self.frequency	= wavenumber2frequency(wavenumber)
			self.label		= label 
			self.flux		= flux
			self.latitude	= latitude
			self.longitude	= longitude
			self.time		= time
			self.cloud		= cloud
			self.transmissivity = transmissivity
		else:
			self.radiances 	= None
			self.Tb			= None
			self.wavenumber = None
			self.frequency 	= None
			self.label 		= None
			self.flux 		= None
			self.latitude	= None
			self.longitude	= None
			self.time		= None
			self.cloud		= None
			self.transmissivity = None
		
		if inherit:
			for v in Spectrum.inherit_v:
				if hasattr(self.__getattribute__(v), '__iter__'):
					pass
				elif self.__getattribute__(v) in [None, -9999]:
					self.__setattr__(v, inherit.__getattribute__(v))

		#_for matching... or something
		self.spectrum_file	= spectrum_file
		
	def plot(self, figure=None, title='', show=False, Tb_diff=None, **kwargs):
		import matplotlib.pyplot as plt
		from numpy import where, diff, zeros
		'''
		plots current spectra to figure.axes.  Does nothing to
		persistently save it currently.
		'''

		dbg('PUT IN OPTION TO PLOT ABSORPTION LINES')
		#_if new, craete figure
		if figure == None:
			self.figure = plt.figure()
			self.figure.add_subplot(311)	#_radiance
			self.figure.add_subplot(312)	#_brightness temp
			self.figure.add_subplot(313)	#_difference
			self.figure.axes[0].set_title(title)
			ls = 'b-'
		else:
			self.figure = figure
			ls = 'r-'
			
		if Tb_diff == None:
			Tb_diff = zeros((self.wavenumber.size))
		
		#_loop over indices near each other
		arg = {'label':self.label, 'linewidth':.3}
		start = 0
		for end in where(diff(self.wavenumber) > 50)[0]:
			#_radiances
			self.figure.axes[0].plot(	self.wavenumber[start:end],
										self.radiances[start:end], ls, **arg)
			#_brightness temp
			self.figure.axes[1].plot(	self.wavenumber[start:end],
										self.Tb[start:end], ls, **arg)	
			self.figure.axes[2].plot(	self.wavenumber[start:end],
										Tb_diff[start:end], ls, **arg)
				
			start = end+1
		else:
			self.figure.axes[0].plot(	self.wavenumber[start:],
										self.radiances[start:], ls, **arg)
			self.figure.axes[1].plot(	self.wavenumber[start:],
										self.Tb[start:], ls, **arg)
			self.figure.axes[2].plot(	self.wavenumber[start:],
										Tb_diff[start:], ls, **arg)
		
		#_turnon grids
		[xxx.grid(True) for xxx in self.figure.axes]

		#_plot differences
	
		cymin, cymax = self.figure.axes[0].yaxis.get_data_interval()
		tymin, tymax = self.figure.axes[1].yaxis.get_data_interval()

		xlim0 = self.figure.axes[0].xaxis.get_data_interval()
		xlim1 = self.figure.axes[1].xaxis.get_data_interval()
		xlim2 = self.figure.axes[2].xaxis.get_data_interval()
		self.figure.axes[0].set_xlim(xlim0)
		self.figure.axes[1].set_xlim(xlim0)
		self.figure.axes[2].set_xlim(xlim0)

		self.figure.axes[0]
		self.figure.canvas.draw()
		
		if show: plt.show()
		
	def write(self, fname=None, 
		dir_out='/'.join((os.environ['PRODUCTS'], 'spectra')), **kwargs):
		'''write data to label.nc'''
		from netCDF4 import Dataset
		from libtools import mkdir_p
		
		fname = fname if fname != None \
			else '.'.join((self.label, 'nc'))
		pname = '/'.join((dir_out, fname.split('/')[-1]))
		
		#_append if file already present
		if not os.path.exists(pname):
			mkdir_p(dir_out)
			dbg(('writing', pname))
			cdf = Dataset(pname, 'w', format='NETCDF3_CLASSIC')
			
			#_initialize dimensions
			cdf.createDimension('time', None)
			cdf.createDimension('wavenumber', self.radiances.size)
			
			#_initialize variables
			cdf.createVariable('latitude', 	'f8', ('time',), fill_value=-9999)
			cdf.createVariable('longitude', 'f8', ('time',), fill_value=-9999)
			cdf.createVariable('time', 		'f8', ('time',), fill_value=-9999)
			cdf.createVariable(	'wavenumber', 'f8', ('wavenumber'),
								fill_value=-9999)
			cdf.createVariable( 'cloud', 'i4', ('time',), fill_value=-9999)
			cdf.createVariable(	'radiances', 'f8', ('time', 'wavenumber',),
								fill_value=-9999)
			cdf.createVariable(	'Tb', 		'f8', ('time', 'wavenumber'), 
								fill_value=-9999)
			cdf.createVariable(	'wavenumber','f8', ('wavenumber'),
								fill_value=-9999)
			cdf.createVariable(	'frequency', 'f8', ('wavenumber'),
								fill_value=-9999)
			dbg((self.frequency.shape, self.radiances.size))	
			cdf.variables['wavenumber'][:]		= self.wavenumber
			cdf.variables['frequency'][:]		= self.frequency
			
			cdf.variables['wavenumber'].units 	= 'cm-1'
			cdf.variables['Tb'].units 			= 'K'
			cdf.variables['radiances'].units 	= 'W/m2/um/sr'
			
			cdf.setncattr('label', self.label)
		else:
			dbg(('appending', pname))
			cdf	= Dataset(pname, 'a', format='NETCDF3_CLASSIC')
		
		#_write spectrum data
		idx = cdf.variables['time'].size	#_last index of UNLIMITED
		for var in ['radiances', 'Tb']:
			cdf.variables[var][idx] 	= self.__getattribute__(var).squeeze()
		cdf.variables['latitude'][idx] 	= self.latitude
		cdf.variables['longitude'][idx] = self.longitude
		cdf.variables['time'][idx] 		= self.time
		cdf.variables['cloud'][idx]		= self.cloud
		
		cdf.close()
	
	
class Srf(object):
	def __init__(self, srffile=os.environ['PRODUCTS']+
		'/airs/AIRS_srftables_51118v4.hdf', **kwargs):
		'''
		u'To convert fwgrid to a frequency grid, 
		freqgrid = fwgrid*width + freq
		
		Do a linear interpolation of srfval and fwgrid 
		(or freqgrid) to calculate the SRF on any disired 
		frequency grid; Do not extrapolate the SRF beyond 
		the endpoints of fwgrid 
		'''
		from netCDF4 import Dataset
		dbg(('reading', srffile), 5)
		hdf = Dataset(srffile, mode='r', format='HDF4')
	
		#_store variables
		for variable, values in hdf.variables.iteritems():
			
			#_store their attributes
			for attr in hdf.variables[variable].ncattrs():
				val = hdf.variables[variable].__getattribute__(attr)
				dbg((variable, attr, val.strip()), 5)
##				tmp  = self.__getattribute__(variable)
##				tmp.__setattr__(attr, attr, attr)
				
			self.__setattr__(variable, values[:].squeeze())
				
			
		for attr in hdf.ncattrs():
			self.__setattr__(attr, hdf.__getattribute__(attr))
		
		hdf.close()

	
class AIRS_ancillary(object):
	'''
	return AIRS ancillary data
	
	Paolo provided L2.chan_prop.2003.11.19.v8.1.0.tobin.anc which may
	be reduncanty with AIRS_srftables_51118v4.hdf, but I have not checked
	as of May 2014.  Primarily loading this for RAD QUAL/L2_Ingore columns
	'''
	def __init__(self, nchan=2378, AIRS_ancfile=os.environ['PRODUCTS']+
		'/airs/L2.chan_prop.2003.11.19.v8.1.0.tobin.anc', **kwargs):
		#_from numpy import genfromtxt (not fixed ncol, so can't use)
		from numpy import tile
		
		self.chanid			= tile(nan, nchan)
		self.l2_ignore 		= tile(nan, nchan)
		self.srf_centroid 	= tile(nan, nchan)
		
		#_read in data
		with open(AIRS_ancfile) as f:
			n = 0
			for line in f.readlines()[123:]:
				cols = line.split()
				self.chanid			= int(cols[0])
				self.srf_centroid[n]= float(cols[1])
				self.l2_ignore[n] 	= int(cols[12])

				n += 1
		
#############################################################################_80
################################################################################
################################################################################

def match_collocated_spectrum(s, **kwargs):
	'''
	for whatever stupid reason, I didn't save the AIRS spectra with
	the collocated data, so now I have to read in and find the appropriate
	scenes via times
	'''
	
	raise RuntimeError, 'DO NOT USE'
	from libtools import epoch2date, dtg2epoch
	from glob import glob
	
	date = epoch2date(s + dtg2epoch('1993010100'))
	year, mon, day = date[:4], date[5:7], date[8:10]
	
	#_loop over files until this time is found
	for fname in glob(dir_col +'/COL_AIRS'):
		pass
	
	return fname


def get_surf_temp2(fov, surf_temp_src='LBL-ATM', surf_temp_var='skt',
	shift_sfct=False, ecmwf_file=JOIN(DIR_PROD,'era-interim','hs3_2013.nc'),
	dir_lblrtm='.', **kwargs):
	'''
	Surface temperature source seems to be a big hangout
	fov				hs3_class,	Field of view from Flight_segment() class
	surf_temp_src	str,		Select surface temperature source
								GDAS || SHIS || ECMWF || LBL-ATM
								Not all are currently implemented
								Selecting LBL-ATM results in profile
								used from LBL-RTM run (-1)
	surf_temp_var	str,		Name of variable to use in input file 

	LBL-ATM, GDAS, and SHIS are natively associated with the fov object,
	however the ECMWF and any other are external and take this kludging. 
	'''
	from netCDF4 import Dataset
	from libtools import ecmwf_day2dtg as d2d
	from libtools import dtg2epoch as d2e
	from libtools import epoch2dtg as e2d
	from numpy import array
	from tape7 import tape7

	#_alert if shifting arbitrarily
	if shift_sfct:
	#	raise RuntimerError, "STOP USING THIS"
		dbg('Shifting surface temperature by {0} K!'.format(str(shift_sfct)))

	#_only do this if ECMWF 
	if surf_temp_src in ['GDAS', 'SHIS']:
		return getattr(fov, '{0}_sfc_temperature'.format(surf_temp_src))+shift_sfct 
	elif surf_temp_src == 'LBL-ATM':
	    return -1
	elif surf_temp_src == 'LBL-ATM_old':
		t7 = tape7(os.path.join(dir_lblrtm, 'TAPE7'))
		return t7.tlevel1[t7.zlevel1.argmin()]
	elif type(surf_temp_src) != str:
		#_arbitrary value to be passed
		return surf_temp_src
	else:
		#_move onto kludgey cases
	    pass

	#_open netcdf file
	ncdf = Dataset(ecmwf_file, 'r')

	#_get indices
	yidx = abs(fov.SHIS_latitude - ncdf.variables['latitude']).argmin()
	xidx = abs(fov.SHIS_longitude - ncdf.variables['longitude']).argmin()
	epoch = array([d2e(d2d(t)) for t in ncdf.variables['time']])
	tidx = abs(fov.SHIS_epoch - epoch).argmin()

	dbg((d2d(ncdf.variables['time'][tidx]), ncdf.variables['latitude'][yidx], \
	    ncdf.variables['longitude'][xidx]))
	dbg((e2d(fov.SHIS_epoch), fov.SHIS_latitude, fov.SHIS_longitude))

	surf_temp = ncdf.variables[surf_temp_var][tidx,yidx,xidx]
	ncdf.close()

	return surf_temp

def get_surf_temp(fov, surf_temp_src='LBL-ATM', surf_temp_var='skt',
	shift_sfct=False, ecmwf_file=JOIN(DIR_PROD,'era-interim','hs3_2013.nc'),
	dir_lblrtm='.', **kwargs):
	'''
	Surface temperature source seems to be a big hangout
	fov				hs3_class,	Field of view from Flight_segment() class
	surf_temp_src	str,		Select surface temperature source
								GDAS || SHIS || ECMWF || LBL-ATM
								Not all are currently implemented
								Selecting LBL-ATM results in profile
								used from LBL-RTM run (-1)
	surf_temp_var	str,		Name of variable to use in input file 

	LBL-ATM, GDAS, and SHIS are natively associated with the fov object,
	however the ECMWF and any other are external and take this kludging. 
	'''
	from netCDF4 import Dataset
	from libtools import ecmwf_day2dtg as d2d
	from libtools import dtg2epoch as d2e
	from libtools import epoch2dtg as e2d
	from numpy import array
	from tape7 import tape7

	#_alert if shifting arbitrarily
	if shift_sfct:
	#	raise RuntimerError, "STOP USING THIS"
		dbg('Shifting surface temperature by {0} K!'.format(str(shift_sfct)))

	#_only do this if ECMWF 
	if surf_temp_src in ['GDAS', 'SHIS']:
		return getattr(fov, '{0}_sfc_temperature'.format(surf_temp_src))+shift_sfct 
	elif surf_temp_src == 'LBL-ATM':
	    return -1
	elif surf_temp_src == 'LBL-ATM_old':
		t7 = tape7(os.path.join(dir_lblrtm, 'TAPE7'))
		return t7.tlevel1[t7.zlevel1.argmin()]
	elif type(surf_temp_src) != str:
		#_arbitrary value to be passed
		return surf_temp_src
	else:
		#_move onto kludgey cases
	    pass

	#_open netcdf file
	ncdf = Dataset(ecmwf_file, 'r')

	#_get indices
	yidx = abs(fov.SHIS_latitude - ncdf.variables['latitude']).argmin()
	xidx = abs(fov.SHIS_longitude - ncdf.variables['longitude']).argmin()
	epoch = array([d2e(d2d(t)) for t in ncdf.variables['time']])
	tidx = abs(fov.SHIS_epoch - epoch).argmin()

	dbg((d2d(ncdf.variables['time'][tidx]), ncdf.variables['latitude'][yidx], \
	    ncdf.variables['longitude'][xidx]))
	dbg((e2d(fov.SHIS_epoch), fov.SHIS_latitude, fov.SHIS_longitude))

	surf_temp = ncdf.variables[surf_temp_var][tidx,yidx,xidx]
	ncdf.close()

	return surf_temp


def convert_airs_time(s):
	'''covert from seconds from 1 jan 1993 to time ('1970.01.01_00:00:00')'''
	from libtools import dtg2epoch, epoch2date
	return epoch2date(s+dtg2epoch('1993010100'))


def sing(x):
	'''generate singleton (Nx1)'''
	from numpy import array
	return array([x])


def check_opts(prof,angle=180.,**kwargs):
	''' look through options for tape5 and make sure no dependencies fail'''
	from numpy import cos, pi
	
	if (prof.pressure[1:]-prof.pressure[:-1]).max()*cos(angle*pi/180.) < -1e-8:
		raise RuntimeError, 'pressure and angle mismatch'
		#_ cos(angle) should match sign of difference
		
###	if kwargs.get('laser') > 1 and kwargs.get('mpts') > 0: 
###		err('Cannot run MPTS > 0 and LASER > 1')
###	if kwargs.get('laser') > 1 and kwargs.get('npts') > 0: 
###		err('Cannot run NPTS > 0 and LASER > 1')
###
###	mrg = kwargs.get('merge')
###	iod = kwargs.get('ioddd')
###	if iod == 2 and mrg != 1: err('IOD=2 requires MERGE=1')
###	if iod == 3 and not (mrg == 1 or (mrg >= 40 and mrg <= 43)):
###		err('IOD=3 requires MERGE={1,40-43}')
###
###	if kwargs.get('iemit') == 3 and not (mrg >= 40 and mrg <= 43):
###		err('IEMIT=3 requires MERGE={40-43}')
###
###	infg = kwargs.get('inflag'); otfg = kwargs.get('otflag')
###	if otfg == 1 and infg != 0: err('OTFLAG=1 requires INFLAG=0')
###	if otfg == 2 and infg != 2: err('OTFLAG=2 requires INFLAG=2')
###
###	if kwargs.get('v2') - kwargs.get('v1') > 2020:
###		err('Wavenumber bandwidth cannot exceed 2020cm-1')


def err(msg): raise RuntimeError, msg


def setup_groups(objects,nproc=1,**kwargs):
	'''
	USAGE:	GROUPS 			= setup_groups( <ITER_OBJ>, <INT> )
		[[1,2,3],[4,5,6]]	= setup_groups( range(6), 3 )

	Break iterable object into groups of NP size to be forked.
	objects	: iterable object 
	nproc	: number of processors, int
	'''
	n_obs	= len( objects )
	ng		= int( n_obs / nproc ) 			#_number of full groups
	extra	= 1 if ( n_obs % nproc ) else 0	#_check for residual 
	groups 	= ['none'] * (ng + extra)		#_setup 2D list 
	for i in range( ng ): groups[i] = objects[ (nproc*i) : nproc*(i+1) ]	
	if extra != 0: groups[ng] = objects[nproc*ng:]
##	dbg(( objects, groups ), l=9 )
	return groups

	
class var(masked_array):
	''' Setup to have basic expected information about 3d variables '''	
	def __new__(cls,input,attrv=None,attrn=['wavenumber']):
		from numpy import where, array, arange, empty, isnan
		from numpy.ma import masked_outside
		
		rad_range = [-1e-6,9.e7] 		#_outside these limits, mask data

		#_initialize object
		obj = masked_outside(input,rad_range[0],rad_range[1])

		#_check number of dimensions, then limit attributes to be made
		ndim 	= input.ndim
		attrn 	= attrn[-ndim:]

		#_set dimensional attributes, e.g., var.attrn[0] => attrv[0]
		for i in arange(ndim):
			attr 	= str(attrn[i])
			value 	= array(attrv[i])
			obj.__setattr__(attr,value)
			if input.shape[i] != value.size:
				err = 	'\nERROR: Dim size mismatch!\n' +\
						'ERROR: ' + attr + ' ' + str(i) +'\n'+ \
						'ERROR_VSIZE: ' + str(input.shape[i]) + '\n' \
						'ERROR_DSIZE: ' + str(value.size) + '\n' 
				raise AttributeError(err) 

		return obj
		
	def __array_finalize__(self,obj): 
		if obj is None: return	

	
def dbg(msg, l=1, err=False):
	''' 
	if global debug is set to true, be more verbose 
	msg	: str, Message to be :rinted
	l	: int, Debug level of message.  Set higher for lower level 
			messages.  As debug increases, noisiness should also.
	'''
	msg = to_string(msg)
	if hasattr(msg,'__iter__'): msg = ' '.join(msg)

	if debug >= l:
		curf 	= inspect.currentframe()
		calf 	= inspect.getouterframes(curf,2)
		file, line, method = calf[1][1:4]
		file	 = '.'.join(file.split('/')[-1].split('.')[:-1])
		scream 	= '[%s.%s.%i] %s' % (file,method,line,msg)
	
		if not err:
			print scream
		else: 
			raise RuntimeError, scream

		
def to_string(a):
	if hasattr(a,'__iter__'): return [ str(b) for b in a ]
	else: return str(a)

	
if __name__ == '__main__':
	sys.exit(run_main())
