#!/usr/bin/env python



'''
NOTES
Speed up Spectrum/plank_inv.  Figure out if yield can be used.

Run LBL-RTM/DIS with aerosol layer at 3km with varying AOD values
adjust reference wavelength

Setup bimodal aerosol distribution (layer at same level, 2 r_eff)
Would I need or want to change the OD?
'''


import sys, os, imp
if 'DISPLAY' not in os.environ:
	import matplotlib
	matplotlib.use('Agg')
	

JOIN = os.path.join
DIR_PROD = os.environ['PRODUCTS']
DIR_HS3 = JOIN(os.environ['PRODUCTS'], 'hs3')
DIR_LBL = JOIN(os.environ['PRODUCTS'], 'LBL-RTM')
DIR_COL = JOIN(os.environ['PRODUCTS'], 'colloc')
DIR_LOG = os.environ['LOG']
# DIR_LOG = os.path.expanduser('~/qsub_logs/')


DEBUG = 1
#_based on SHIS segmented files #_TODO, lookup for CPL


################################################################################
#_MAIN_#########################################################################
################################################################################
#_a few common reminder functions


def process_period(fidx=None, kidx=[0], nproc=10, require_fname_dtg=True):
	'''
	given a start and end dtg, process all flight segments

	includes creation of collocate files from shiscpl,
	collocated netcdf file using both

	plots of flights
	'''
	from libtools import dtg2epoch
	from lblrtm_utils import write_kwargs
	import matplotlib.pyplot as plt
	from glob import glob
	import re

	fidx = [0]

	#_employ a first check on dtg using the filename, but
	# double check with actual metadata
	re_date = re.compile('(\d{12}).(\d{12}).nc')
	files_seg = glob(DIR_HS3 + '/SHIS.CPL.GDAS*nc')

	#_create collocated data netcdf files between dates
	dtg0 = '20130820193743'
	dtg1 = '20130820193744'

##	#_write kwarg files that define the types of runs to be done on eaach fov
##	# edit those in lblrtm_utils()
##	kwarg_names = ['full', 'mw'] #_microwindow makes it -14
##	for name in kwarg_names:
##		write_kwargs(name=name)

##	TO GET LBLRTM OUTPUT, JUST MAKE TEST CASES THROUGH RUN_LBLRTM FOR NOW

	#_plot only the flight data?
	for file_seg in files_seg:
		#_first swipe check
		res = re_date.search(file_seg)

		#_not all files have dtg in name, but it's faster to search
		# if within a date range without actually reading in the data
		try:
			if dtg0 > '20' + res.group(2) or dtg1 < '20' + res.group(1):
				continue

		except AttributeError:
			if require_fname_dtg:
				continue
			else:
				pass

		#_read in flight segment
		flight = Flight_segment(file_seg=file_seg)

		#_skip outside window, not bothering to limit because check is fast
		if dtg0 > flight.dtg1 or flight.dtg0 > dtg1:
			continue

#		try:
			#_plot flight segment without LBL-RTM output
			#flight.plot_flight(fidx=[330], nproc=nproc)
#			flight.plot(fidx=fidx, kidx=kidx, nproc=nproc, plot_full=False)	

#		except:
#			dbg(('failed', file_seg))


def simulate_flights(require_fname_dtg=True):
	'''
	run lblrtm and lbldis on flight
	
	either put in something that pauses until LBLRTM completes, or run twice
	otherwise the qsubmissions will potentially get picked up at the same time
	without the implementation of a lockfile or something.
	'''
	from qsubmissions import lblrtm_hs3_lim, lbldis_hs3_lim
	import re
	from glob import glob

	kwarg_names = ['full', 'mw'] #_microwindow makes it -14

	#_employ a first check on dtg using the filename, but
	# double check with actual metadata
	re_date = re.compile('(\d{12}).(\d{12}).nc')
	files_seg = glob(DIR_HS3 + '/SHIS.CPL.GDAS*nc')

	#_clear sky segments
	dtg0 = '20130820193743'
	dtg1 = '20130820193744'

	for file_seg in files_seg:
		#_first swipe check
		res = re_date.search(file_seg)

		#_not all files have dtg in name, but it's faster to search
		# if within a date range without actually reading in the data
		try:
			if dtg0 > '20' + res.group(2) or dtg1 < '20' + res.group(1):
				continue
		except AttributeError:
			if require_fname_dtg:
				continue
			else:
				pass

		flight = Flight_segment(file_seg=file_seg)

		#_skip outside window, not bothering to limit because check is fast
		if dtg0 > flight.dtg1 or flight.dtg0 > dtg1:
			continue

		#_submit jobs to run for all kidx/fidx for these flights
		opts = {'file_flight' : file_seg}

		#_submit lblrtm job
		lblrtm_hs3_lim(**opts)

###		for name in ['mw']: #kwarg_names:
###			#_all above jobs need to be done before anything else
###			opts.update({'file_kw' : '{0:s}_kwarg_lbldis.pk'.format(name)})
###			lbldis_hs3_lim(kidx=[0], **opts)


def compare_simulations():
	pass


################################################################################
#_SUB_##########################################################################
################################################################################

'''
to produce flight segment files
create_collocation_files(dir_shis=<path_to_SHIS_rdr>, dir_cpl=<path_to_cpl>)
create_segment_files()
'''

def create_collocation_files(dtg0=0, dtg1=1e9, dir_shis=DIR_HS3, pcfilt=True,
	dir_cpl=JOIN(DIR_PROD, 'cpl'), dir_col=DIR_COL, **kwargs):
	''' 
	From dir_shis, attempts to create collocationed segment files in
	dir_col (default $PRODUCTS/colloc) from the shiscpl script.
	The COLLOC files are then used to create segment files.

	THIS PART IS OUT OF DATE, FLIGHTS FILE DEPRECATED

	dtg0	str,	set if earliest date desired
	dtg1	str,	set if end date desired, default to last in dir_shis
	save_old bool,	true to skip already generated files	

	dir_col	str,	NOT CURRENTLY USED!  Check shiscpl script.
	'''
	import re, glob
	from libtools import dtg2iso, newdtg2
	from libtools import dtg2epoch as d2e
	from libtools import epoch2dtg as e2d 
	from numpy import array

	#_get first dtg for flights
	#_read in all available hours (what is last time?)
	re_dtg = re.compile('SHIS.CPL.COLLOC.(\d{12}).(\d{12}).hdf')
	re_cpl = re.compile('OP_\d{5}\w?_(\d{6})_(\d{6})_(\d{6}).nc')
	re_shs = re.compile('rdr(\d{8})T(\d{6})end(\d{8})T(\d{6})sdr\d{8}T\d{6}')
	re_xdr = re.compile('OP_(\d{5}\w?)_')
	
	#_collect start and end times for available SHIS files
	estart = []
	files_shis	= glob.glob(dir_shis + '/*rad{0}.nc'.format('_pcfilt'*pcfilt))
	files_cpl	= array(glob.glob(dir_cpl + '/nc/*'))
	files_xdr	= array(glob.glob(dir_cpl + '/xdr/*'))

	#_build dictionary of xdr files using julian day
	xdr = {}
	for file_xdr in files_xdr:
		res = re_xdr.search(file_xdr)
		key = res.group(1)
	
		#_add key and file to dict 
		xdr.update({key : file_xdr})

	#_create list of CPL start and end times
	cpl_bounds = {} 
	for file_cpl in files_cpl:
		res = re_cpl.search(file_cpl)

		if not res:
			continue

		dtg0 = '20{0:6s}{1:6s}'.format(res.group(1), res.group(2))
		dtg1 = '20{0:6s}{1:6s}'.format(res.group(1), res.group(3))

		#_search for julian date key
		key = re_xdr.search(file_cpl).group(1)

		start	= d2e(dtg0, full=True)
		end		= d2e(dtg1, full=True)

		#_meeeeeh, day crossover. This whole section is bad,
		# but so is the naming convention.
		if start > end:
			end += 86400

		try:
		 cpl_bounds.update({key : {'bounds' : [start, end], 'file' : xdr[key]}})
		except KeyError:
		 continue

	#_loop over SHIS files to find overlapping periods
	collocation = {}
	for file_shis in files_shis:
		cpl_list = []
		
		#_pull out SHIS start and end times
		res	= re_shs.search(file_shis)

		#_file didn't match, move on
		if not res: 
			continue

		#_convert to epoch time
		str = d2e(res.group(1) + res.group(2), full=True)
		end = d2e(res.group(3) + res.group(4), full=True)

		#_loop over cpl start and end times
		for key, values in cpl_bounds.iteritems():
 
			#_is ANY part of the file overlapping?
			cpl_str, cpl_end = values['bounds']
			if	(cpl_str >= str and cpl_str <= end) or	\
				(cpl_end >= str and cpl_end <= end) or	\
				(cpl_str <= str and cpl_end >= end):		#_cpl encompasses
				cpl_list.append(values['file']) 

		#_add shis file to collocation dictionary	
		if len(cpl_list):
			collocation.update({file_shis : cpl_list})

	#_loop over shis files, collocate and collect names of collocation files
	files_col = []
	err_files = []
	for file_shis, files_cpl in collocation.iteritems():
		for file_cpl in files_cpl:
			dbg(file_shis) 
			dbg(file_cpl)
		#	if file_shis != '/data/wsessions/hs3/SHIS_rdr20130830T042645end20130830T091836sdr20140109T122247_rad_pcfilt.nc':
		#		continue
			try:
				#_run collocation code, get filename
				cmd			= ' '.join(('shiscpl', file_shis, file_cpl))
				p			= os.popen(cmd, 'r')
				stdio		= p.readlines()
				file_col	= stdio[-1].rstrip() 
			except IndexError:
				dbg(('no go', cmd))
				continue

			#_check for collocation file
			dbg((file_col, '\n\n'))
			res = re_dtg.search(file_col)
			if not res:

				dbg(('problem with collocation of', file_col))
				err_files.append([file_shis, file_cpl])
				continue
			
	#_write files that failed to error file
	err_log = os.path.join(dir_col, 'failed_files')
	with open(err_log, 'w') as f:
		[f.write('{0} {1}\n'.format(s, c)) for s, c in err_files]


def create_segment_files_qsub(dir_col=DIR_COL, queue='qsub', **kwargs):
	'''
	submit job to queue for generation

	Trying to split based on qsub or sbatch, 2016 06 27

	'''
	from glob import glob
	from pickle import dump
	from libtools import dummy_sbatch as write_sbatch
    #_qsubmission prefix        
	scpt = os.path.expanduser('~/lib/run_create_collocation.py')

	#_get list of collocation files
	files_col = glob(dir_col + '/SHIS.CPL.COLLOC.*hdf')
	for fname in files_col:

		#_create pickle file
		fn = fname.split('/')[-1]
		fkw = JOIN(dir_col, 'kwargs.{0:d}.{1}.pk'.format(os.getpid(), fn))
		dump(kwargs, open(fkw, 'wb'))

		if queue == 'qsub':
			env = ','.join(['='.join((var, os.environ[var])) \
				for var in ['PYTHONPATH', 'PRODUCTS', 'PATH', 'LOG']])
			qsub = ' '.join((   'qsub -v', env, '-o', DIR_LOG, '-e', DIR_LOG,
				'-cwd -S /opt/ShellB3/bin/python'))
	
			#_put it all together and submit
			cmd = ' '.join((qsub, scpt, fname, fkw))
			dbg(cmd) 
			os.system(cmd)

		elif queue == 'sbatch':
			unique_id = fname.split('/')[-1][:-4]
			sname = write_sbatch(' '.join((scpt, fname, fkw)),
				unique_id=unique_id)

			#_submit job
			cmd = ' '.join(('sbatch', sname))
			os.system(cmd)

		else:
			raise RuntimeError, 'Invalid queue option: {0}'.format(queue)

def create_segment_files(dir_col=DIR_COL, **kwargs):
	'''
	reads in colloc files produced by create_collocation_files()
	from shiscpl and creates single matched file

	This method is simply a simple way of regenerating based 
	upon what is in the dir_col directory because I kept forgetting how.
	
	dir_col	str,	path to SHIS.CPL.COLLOC files	
	'''
	from glob import glob
	from libtools import setup_groups

	files_col = glob(dir_col + '/SHIS.CPL.COLLOC.*hdf')

	#_loop over all files in collocation directory
	for file_col in files_col:
		#_read in sensor data
		shis, cpl, prof = read_collocated(file_col, **kwargs)

		#_if there is no shis data, skip file
		if shis == False: 
		#if not shis:
			continue

		#_write COLLOC file
		write_collocated(shis, cpl, prof, file_col, **kwargs)


def read_collocated(fname, **kwargs):
	''' read in colloc index file data associated with shiscpl output '''
	from netCDF4 import Dataset
	import re
	from numpy import diff, append

	try:
	  with Dataset(fname, 'r') as hdf:
		sidx = hdf.variables['SHIS_Index'][:] - 1#_the count seems to start at 1
		cidx = hdf.variables['CPL_Index'][:] - 1
		w = hdf.variables['Weights'][:]

		#_reformat as strictly strings
		f_shis	= '{0}'.format(hdf.SHIS_File)
		f_cpl	= '{0}'.format(hdf.CPL_File)	

	except RuntimeError:
		#_HDF on IRIS
	  from pyhdf.SD import SD, SDC

	  hdf = SD(fname, SDC.READ)
	  sidx = hdf.select('SHIS_Index')[:] - 1
	  cidx = hdf.select('CPL_Index')[:] - 1
	  w = hdf.select('Weights')[:]

	  #_reformat as strings
	  f_shis = hdf.attributes()['SHIS_File']
	  f_cpl = hdf.attributes()['CPL_File']

	#_swap xdr file with nc file
	f_cpl = cpl_xdr2nc(f_cpl)

	#_read SHIS and CPL indices
	shis	= read_shis(f_shis, sidx, **kwargs)
	cpl		= read_cpl(f_cpl, cidx[:,0], **kwargs)
	prof	= read_profile_gdas(f_shis, sidx, **kwargs)

	#_check for jumps. Collocation code only looks at location,
	# does not seem to account for backtracking
	print len(shis), cidx, cidx.shape 
	idx_jump = abs(shis.epoch - cpl.epoch) < 150 
##	idx_jump = diff(cpl.epoch) > 0
##	idx_jump = append(True, idx_jump)
	s_dict = shis.__dict__.copy()
	c_dict = cpl.__dict__.copy()
	g_dict = prof.__dict__.copy()
	shis = shis[idx_jump]
	cpl = cpl[idx_jump]
##	cpl = cpl[idx_jump,:] 2016.06.30 for JUMP
	prof = prof[idx_jump]
	shis.__dict__ = s_dict
	cpl.__dict__ = c_dict
	prof.__dict__ = g_dict


	if prof == False:
		dbg(('WARNING! Missing data for', fname))
		return False, False, False

	return shis, cpl, prof



def collocate_modis(fname, thresh_time=3600*6, modis_meta=None, **kwargs):
	from glob import glob
##	modis = read_collocated_modis(fname, **kwargs)
##	write_collocated_modis(fname, modis, **kwargs)

	#_get list of files to collocate with,
	files = glob('{0}/SHIS.CPL.GDAS.COLLOC.13*nc'.format(os.path.join(DIR_PROD,
		'hs3')))

	#_loop over files and match locations, then write it out
	for file in files:
		#_read in flight file to get times
		flight = Flight_segment(file_seg=file)

		#_read in all modis for that period
		e0 = flight.SHIS_epoch[0] - thresh_time 
		e1 = flight.SHIS_epoch[-1] + thresh_time 
		modis_meta = read_modis_location(e0, e1, **kwargs)
		modis = read_collocated_modis2(file, modis_meta=modis_meta, **kwargs)
		write_collocated_modis(file, modis, **kwargs)


def write_collocated_modis(fname, modis, **kwargs):
	''' write out collocated modis information '''
	from netCDF4 import Dataset
	import os

	fout = fname.replace('CPL.GDAS','MODIS')
	if os.path.exists(fout): return

	cdf = Dataset(fout, 'w', format='NETCDF3_CLASSIC')	

	#_initialize dimension for field of views
	cdf.createDimension('fov', modis.size, **kwargs) 

	#_loop over recarray dtype and add to file
	for i, varname in enumerate(modis.dtype.names):
		vartype = modis.dtype[i]
		vcdf = cdf.createVariable(varname, vartype, 'fov')
		vcdf[:] = getattr(modis, varname)
	##	setattr(vcdf, varname, getattr(modis, varname))

	cdf.close()


def read_collocated_modis(fname, dir_mod=os.path.join(DIR_PROD, 'MODIS'),
	**kwargs):
    '''based on location and time in fname, find closest modis overpass data '''
    import re
    from libtools import julian2epoch, julian2dtg, epoch2dtg
    from glob import glob
    from numpy import array, recarray
    from modis import read_modis_aod_latlon as read_modis

	#_get list of available MODIS Aerosol data
    files = glob('{0}/M*D04_3K*hdf'.format(dir_mod))
    files.sort()

	#_setup regex to pull out dates from filenames
    re_tm = re.compile('M.D04_3K.A(\d{4})(\d{3}).(\d{2})(\d{2})')

    epoch = []
    for f in files:
        reg = re_tm.search(f)
        year = int(reg.group(1))
        jday = int(reg.group(2))
        hour = int(reg.group(3))
        min = int(reg.group(4))

        epoch.append(julian2epoch(year, jday) + hour*3600 + min*60) 

	#_put files into dictionary keyed by time, then sort times list
    scantimes = dict(zip(epoch, files))
    epoch.sort()
    epoch = array(epoch)

	#_get desired times
    flight = Flight_segment(file_seg=fname)

    #_init array to store optical depths
    dtype = [('aod', 'f4'), ('epoch', 'f8'), ('distance', 'f4'), 
            ('lat', 'f4'), ('lon', 'f4')]
    modis_aod = recarray((flight.size), dtype=dtype)

    for i, shis_epoch in enumerate(flight.SHIS_epoch):
        idx0 = abs(shis_epoch-86400 - epoch).argmin()
        idx1 = abs(shis_epoch+86400 - epoch).argmin()

		#_look at file before and after
        lat = flight.SHIS_latitude[i]
        lon = flight.SHIS_longitude[i]
        dist = 999999999999.
        time = 999999999999.
        aod = -9999. 
        for idx in range(idx0, idx1+1):
            fmod = scantimes[epoch[idx]]	
            maod, mdist, mtime = read_modis(fmod, lat, lon, shis_epoch) 

            if mdist < dist:
                time = mtime
                dist = mdist
                aod = maod

        #_put data into recarray
        modis_aod.aod = aod
        modis_aod.distance = dist
        modis_aod.epoch = time
        modis_aod.lat = lat 
        modis_aod.lon = lon 
 
    return modis_aod	


def read_collocated_modis2(fname, dir_mod=os.path.join(DIR_PROD, 'MODIS'),
	modis_meta=None, thresh_space=111e3, thresh_time=3600*6, **kwargs):
	'''based on location and time in fname, find closest modis overpass data '''
	import re
	from libtools import julian2epoch, julian2dtg, epoch2dtg
	from glob import glob
	from numpy import array, recarray, sqrt, cos, mean
	from modis import read_modis_aod_xy as read_modis
	from libgeo import great_circle
	from netCDF4 import Dataset

	#_get desired times
	flight = Flight_segment(file_seg=fname)

    #_init array to store optical depths
	dtype = [('aod', 'f4'), ('epoch', 'f8'), ('distance', 'f4'), 
            ('lat', 'f4'), ('lon', 'f4')]
	modis_aod = recarray((flight.size), dtype=dtype)

	#_initialize as missing
	modis_aod.aod[:] = -9999.

	#_loop over flight
	for i, fov in enumerate(flight):

		#_look at file before and after
		shis_lat = fov.SHIS_latitude
		shis_lon = fov.SHIS_longitude
		shis_tim = fov.SHIS_epoch

		#_find where scans are within time window
		tidx = abs(shis_tim - modis_meta.epoch) < thresh_time 
		if tidx.sum() == 0:
			continue

		#_search those for ones close enough to the location
		modis_time = modis_meta[tidx]
		req = (shis_lat, shis_lon); d = []
		
		#_this first with bad measurements
		delta_y = 111e3 * (modis_time.lat - shis_lat)
		delta_x = 111e3 * cos(shis_lat) * (modis_time.lon - shis_lon)
		r = sqrt(delta_x**2 + delta_y**2)
		didx = r < thresh_space
		if didx.sum() == 0:
			continue

		modis_dist = modis_time[didx] 
		[d.append(great_circle(req, mod)) for mod \
			in zip(modis_dist.lat, modis_dist.lon)]

		#_get index of shortest distance
		ind = array(d).argmin()
		modis_final = modis_dist[ind]

		#_open file and get aod
		x, y = modis_final.x, modis_final.y
		with Dataset(modis_final.file, 'r') as cdf:
		##	aod = cdf.variables['Effective_Optical_Depth_Average_Ocean'][1,x,y]
			aod = cdf.variables['Optical_Depth_Land_And_Ocean'][x,y]

		#_pull out closest point and check that it is close enough
		modis_aod[i].lat = modis_final.lat 
		modis_aod[i].lon = modis_final.lon 
		modis_aod[i].epoch = modis_final.epoch
		modis_aod[i].distance = d[ind]
		modis_aod[i].aod = aod
 
	return modis_aod	


def read_modis_location(epoch0=0, epoch1=9e9, 
	dir_mod=os.path.join(DIR_PROD, 'MODIS'),
	**kwargs):
	'''
	read in all lat/lon/time data for modis and put into a dict
	epoch#	start and end times in unix time
	'''
	from glob import glob
	from netCDF4 import Dataset
	from numpy import recarray, tile
	from libtools import modis2epoch, epoch2dtg

	files = glob('{0}/M*D04*hdf'.format(dir_mod))

	file = []
	lats = []
	lons = []
	epoc = []
	xxxx = []
	yyyy = []
	
	for p, fname in enumerate(files):
		cdf = Dataset(fname, 'r')
		epoch = modis2epoch(cdf.variables['Scan_Start_Time'][:].flatten())

		#_the extra subsetting for min is to avoid some that are
		# set to 1993 01 01
		if epoch[epoch > 1e9].min() > epoch1 or epoch.max() < epoch0:
			cdf.close()
			continue

		dbg(fname)
		epoc.extend(epoch)
		n = len(cdf.variables['Latitude'][:].flatten())
		nx, ny = cdf.variables['Latitude'][:].shape		
		lats.extend(cdf.variables['Latitude'][:].flatten())
		lons.extend(cdf.variables['Longitude'][:].flatten())
		file.extend([fname] * n)

		xxxx.extend(tile(arange(nx), (ny,1)).T.flatten())
		yyyy.extend(tile(arange(ny), (nx,1)).flatten())

		cdf.close()

	ns = len(files[0])
	dtype = [('file', 'a{0:d}'.format(ns)), ('lat', 'f4'), ('lon', 'f4'),
		('epoch', 'f8'), ('x', 'i4'), ('y', 'i4')]
	out = recarray((len(file)), dtype=dtype)
	out.file[:] = file
	out.lat[:] = lats
	out.lon[:] = lons
	out.epoch[:] = epoc
	out.x[:] = xxxx
	out.y[:] = yyyy
	return out


def write_collocated(shis, cpl, prof, fname, dir_shis=DIR_HS3, notes=None,
	**kwargs):
	''' write collocated flight segment file from read_collocated() '''
	from netCDF4 import Dataset
	from numpy import array, vstack#, #str
	import re

	#_generate new flight segment file name
	fname = fname.split('/')[-1]
	fname = re.sub('.hdf$', '.nc', fname)
	fname = re.sub('.CPL.', '.CPL.GDAS.', fname)
	fname = JOIN(dir_shis, fname)	
	dbg(fname)

	#_pull out start and end times or default to arbitrary date after 1970
	try:
		re_dtg 	= re.compile('(\d{12}).(\d{12}).nc$')
		res		= re_dtg.search(fname)	
		dtg0, dtg1 = '20' + res.group(1), '20' + res.group(2)
	except AttributeError:
		dtg0, dtg1 = '198202230000', '198202230216'

	with Dataset(fname, 'w') as cdf:
		#_create dimensions
		nfov	= shis.size
		nwvn	= shis.wavenumber.size 
		nz_cpl	= cpl.ext_532[0].size
		nz_shis = shis.pressure.size
		nz_gdas = prof.pressure.size 
		nlay	= 10 

		cdf.createDimension('wavenumber', nwvn)
		cdf.createDimension('fov', nfov)
		cdf.createDimension('layer', nlay)
		cdf.createDimension('nz_cpl', nz_cpl) 
		cdf.createDimension('nz_gdas', nz_gdas) 
		cdf.createDimension('nz_shis', nz_shis) 

		#_write global attributes
		cdf.file_cpl	= cpl.fname	
		cdf.file_shis	= shis.fname
		cdf.file_gdas	= prof.fname
		cdf.dtg0		= dtg0
		cdf.dtg1		= dtg1
 
		#_make sure everything is in temporal order
		idx	= shis.epoch.argsort()
		rad_shis = vstack(shis.radiances[idx])
		tem_gdas = vstack(prof.temperature[idx])
		tem_shis = vstack(shis.temperature[idx])
		ozo_gdas = vstack(prof.ozone_mixing_ratio[idx])
		ozo_shis = vstack(shis.ozone_mixing_ratio[idx])
		rel_gdas = vstack(prof.relative_humidity[idx])
		rel_shis = vstack(shis.relative_humidity[idx])
		gpt_gdas = vstack(prof.geopotential_height[idx])
		
		idx = cpl.epoch.argsort()
		tau = vstack(cpl.tau_532[idx])			
		ext = vstack(cpl.ext_532[idx])			
		typ = vstack(cpl.tau_type[idx])			

		#_create variables
		#_write variables
		v = {'fill_value' : -9999.}

		c		= cdf.createVariable('CPL_tau_532', 'f4', ('fov', 'layer'), **v)
		c[:]	= tau
		c.units = 'unitless'	

		c		= cdf.createVariable('CPL_tau_type', 'i4', ('fov', 'layer'),**v)
		c[:]	= typ
		c.units = '1=PBL, 2=AEROSOL, 3=CLOUD'

		c		= cdf.createVariable('CPL_ext_532', 'f8', ('fov', 'nz_cpl'),**v)
		c[:]	= ext
		c.units = '???'

		c		= cdf.createVariable('CPL_latitude', 'f8', ('fov',), **v)
		c[:]	= cpl.latitude
		c.units = 'degrees_north'

		c		= cdf.createVariable('CPL_longitude', 'f8', ('fov',), **v)
		c[:]	= cpl.longitude
		c.units = 'degrees_east'

		c		= cdf.createVariable('CPL_epoch', 'f8', ('fov',), **v)
		c[:]	= cpl.epoch
		c.units = 'sec since Jan 1 1970'	

		c		= cdf.createVariable('SHIS_radiances','f8',
									('fov','wavenumber'),**v)
		c[:]	= rad_shis
		c.units = 'W/m2/sr/cm-1'
	
		c		= cdf.createVariable('SHIS_wavenumber','f4',('wavenumber'), **v)
		c[:]	= shis.wavenumber
		c.units = 'cm-1'	
	
		c		= cdf.createVariable('SHIS_HBB_NESR','f4',('wavenumber'), **v)
		c[:]	= shis.hbb_nesr 

		c		= cdf.createVariable('SHIS_latitude', 'f8', ('fov',), **v)
		c[:]	= shis.latitude
		c.units = 'degrees_north'

		c		= cdf.createVariable('SHIS_longitude', 'f8', ('fov',), **v)
		c[:]	= shis.longitude
		c.units = 'degrees_east'
	
		c		= cdf.createVariable('SHIS_epoch', 'f8', ('fov',), **v)
		c[:]	= shis.epoch
		c.units = 'sec since Jan 1 1970'	

		c		= cdf.createVariable('SHIS_altitude', 'f4', ('fov',), **v)
		c[:]	= shis.altitude
		c.units = 'm'
		c.long_name = 'sensor altitude'

		c		= cdf.createVariable('SHIS_pressure', 'f4', ('nz_shis',), **v)
		c[:]	= shis.pressure 
		c.units = 'hPa'

		c		= cdf.createVariable('SHIS_temperature', 'f4',
									('fov','nz_shis',), **v)
		c[:]	= tem_shis 
		c.units = 'K'

		try:
			c		= cdf.createVariable('SHIS_CTT', 'f4', ('fov',), **v)
			c[:]	= shis.cloud_top_temperature 
			c.units = 'K'

			c		= cdf.createVariable('SHIS_CTP', 'f4', ('fov',), **v)
			c[:]	= shis.cloud_top_pressure
			c.units = 'hPa'

			c		= cdf.createVariable('SHIS_COT', 'f4', ('fov',), **v)
			c[:]	= shis.cloud_optical_depth
			c.units = ''
		except AttributeError:
			pass #_meeeh

		c		= cdf.createVariable('SHIS_relative_humidity', 'f4', 
				('fov','nz_shis',), **v)
		c[:]	= rel_shis 
		c.units = 'percent'

		c		= cdf.createVariable('SHIS_ozone_mixing_ratio', 'f4',
				('fov','nz_shis',), **v)
		c[:]	= ozo_shis 
		c.units = 'ppmv'

		c		= cdf.createVariable('SHIS_surface_temperature', 'f4', 
				('fov',), **v)
		c[:]	= shis.sfc_temperature
		c.units = 'K'

		c		= cdf.createVariable('GDAS_latitude', 'f8', ('fov',), **v)
		c[:]	= prof.latitude
		c.units = 'degrees_north'

		c		= cdf.createVariable('GDAS_longitude', 'f8', ('fov',), **v)
		c[:]	= prof.longitude
		c.units = 'degrees_east'

		c		= cdf.createVariable('GDAS_pressure', 'f4', ('nz_gdas',), **v)
		c[:]	= prof.pressure 
		c.units = 'hPa'

		c		= cdf.createVariable('GDAS_temperature', 'f4',
				('fov','nz_gdas',), **v)
		c[:]	= tem_gdas
		c.units = 'K'

		c		= cdf.createVariable('GDAS_relative_humidity', 'f4', 
				('fov','nz_gdas',), **v)
		c[:]	= rel_gdas 
		c.units = 'percent'

		c		= cdf.createVariable('GDAS_ozone_mixing_ratio', 'f4',
				('fov','nz_gdas',), **v)
		c[:]	= ozo_gdas 
		c.units = 'kg/kg (converted to g/kg in reading)'

		c		= cdf.createVariable('GDAS_surface_temperature', 'f4',
				('fov',), **v)
		c[:]	= prof.sfc_temperature 
		c.units = 'K'

		#_quick check
	#	if	abs(shis.latitude - prof.latitude).max() > 0.5	\
	#		or abs(shis.longitude - prof.longitude).max() > 0.5:
	#
	#		print shis.latitude - prof.latitude).max()
	#		raise ValueError, 'Collocation failed!'

		for plat, slat in zip(prof.latitude, shis.latitude):
			if abs(plat - slat) > 0.5 and plat != -8888.:
				print 'PLAT, SLAT', plat, slat
				raise ValueError, 'Collocation failed!: SLAT,PLAT'
		for plon, slon in zip(prof.longitude, shis.longitude):
			if abs(plon - slon) > 0.5 and plon != -8888.:
				raise ValueError, 'Collocation failed!: SLON,PLON'

		#_add arbitrary comments
		if notes is not None:
			cdf.notes = notes


def read_cpl_layers(fname, **kwargs):
	'''
	2016 June 30
	Use to manually read in and process CPL layers
	
	fname	str,	Collocation file.  No more flight file garbage 
	'''
	#_read in collocated data
	s, c, p = read_collocated(fname, **kwargs)

		



def read_all_retrievals(out_labels=[], **kwargs):
	from libtools import merge

	for i, out_label in enumerate(out_labels):
		kwargs.update({'out_label' : out_label})
		if i == 0:
			data = read_retrieved_campaign(**kwargs)
		else:
			data = merge((data, read_retrieved_campaign(**kwargs)))

	return data


def read_retrieved_campaign(**kwargs):
	'''
	Used with read_retrieved_segment to read 
	entire flight campaign during hs3
	'''
	from libtools import merge
	from flight_namelists import experiments
	
	for i, (dtg, values) in enumerate(experiments['HS3'].iteritems()):

		kwargs.update(values)

		if i == 0:
			data = read_retrieved_segment(dtg=dtg, **kwargs)	

		else:
			tmp = read_retrieved_segment(dtg=dtg, **kwargs)
			print 'TESTING', tmp.size
			data = merge((data, tmp)) 

	return data
	

def read_retrieved_segment(dtg=None, out_label='TMP',
	experiment='HS3', modis=False, cpl=False, **kwargs):
	'''
	Read in the final solution for a segment from optical estimation crap

	out_label	str, needed for input file names
	experiment	str, needed for input file directory names

	'''
	from lblrtm_utils import Lbldis_input
	from oe_real_hs3 import check_missing, find_BT_diff
	from numpy import recarray

	#_open flight segment handle
	flight = Flight_segment(dtg=dtg, **kwargs)
	fidx = range(flight.size)

	#_see if anything is available
	missing_no, total = check_missing(dtg, **kwargs)
	
	if missing_no == total:
		dbg('No retrievals for {0}'.format(dtg))
		return False

	#_initialize output array
	dtype = [('tau', 'f4'), ('ref', 'f4'), ('ssp', 'a20'),
			('epoch', 'f8'), ('latitude', 'f4'), ('longitude', 'f4'),
			('z_bot', 'f4'), ('z_top', 'f4'), ('out_label', 'a100'), 
			('sensor', 'a10'), ('cmask', 'i2')]
##	data = recarray((flight.size,), dtype)
	data = recarray((0,), dtype)

	#_pull out formating specification
	dir_lblrtm_fmt = kwargs.get('dir_lblrtm_fmt')

	#_if modis wanted
	if modis:
		from netCDF4 import Dataset
		fmodis = flight.FLIGHT_file.replace('CPL.GDAS', 'MODIS')
		cdf_modis = Dataset(fmodis, 'r') 

	#_get cloud mask for threshold
	mask_dbt = find_BT_diff(flight, bt_thresh=7., **kwargs)

	#_loop over flight and add
	for idx, i in enumerate(fidx):

		#_generate input filename
		dtgfov = (flight.dtg0, flight.dtg1, i, experiment, flight.dtg0[2:])
		dir_lblrtm = dir_lblrtm_fmt.format(*dtgfov)
		lbldis_out = 'lbldis_input.{0}.final'.format(out_label)
		lbldis_out = os.path.join(dir_lblrtm, lbldis_out)

		#_open output
		try:
			lbldis = Lbldis_input(lbldis_out)
		except IOError:
			print 'Missing {0}'.format(lbldis_out)
			continue

		#_merge layers with same habit
		lbldis.merge_layers()

		#_get location information
		lat = flight[i].SHIS_latitude
		lon = flight[i].SHIS_longitude
		epoch = flight[i].SHIS_epoch

		#_build lists based on ssp db
		for layer in lbldis.layers:
			ssp = lbldis.ssp_db[layer['dbnum']]
			habit = ssp.split('.')[1]

			tau = layer['tau'][0]
			ref = layer['ref']

			z_b = layer['z']
			z_t = layer['z_top']

			#_add to recarray
			data.resize(data.size + 1)
			data[-1] = (tau, ref, habit, epoch, lat, lon, z_b, z_t, out_label,
						'S-HIS', mask_dbt[i])

		#_if modis desired, read in as well
		if modis:
			tau = cdf_modis.variables['aod'][i]

			#_if too far, mask
			if cdf_modis.variables['distance'][i] > 3e4 or tau < 0:
				tau = -9999.

			#_matched with SHIS epoch
			data.resize(data.size + 1)
			data[-1] = (tau, -9999, 'aerosol', epoch, lat, lon, -9999, -9999,
					'MODIS', 'MODIS', False)

		if cpl:
			from numpy import append
			taus = flight[i].CPL_tau_532[:5]
			type = flight[i].CPL_tau_type[::2]

			t2 = append((type == 1)[:,None], (type == 2)[:,None], 1)
			aerosol = sum(taus[(taus > 0) * t2.any(1)])
			cloud = sum(taus[(taus > 0) * (type == 3)])

			#_using SHIS epoch	
			data.resize(data.size + 2)
			data[-1] = (aerosol, -9999, 'aerosol', epoch, lat, lon, -9999,-9999,
					'CPL', 'CPL', False)
			data[-2] = (cloud, -9999, 'cloud', epoch, lat, lon, -9999,-9999,
					'CPL', 'CPL', False)

	return data
	

def read_cpl_flighttrack(fname, idx=None, dir_cpl=JOIN(DIR_PROD, 'cpl', 'nc'),
	**kwargs):
	''' read flight track only '''
	from netCDF4 import Dataset
	from numpy import ndarray, recarray, arange
	from libtools import julian2epoch as j2e 
	from numpy.ma import masked_where

	#_find filename for this dtg
	fname = JOIN(dir_cpl, fname)
	dbg(fname)
	h_cpl = Dataset(fname, 'r', format='NETCDF3_CLASSIC')

	#_read in cpl data
	dtype = [	('epoch', 'f8'),
				('latitude','f8'),
				('longitude','f8')]
	cpl = recarray((0,), dtype=dtype)
	cpl.__setattr__('fname', fname)

	if idx is None:
		idx = arange(h_cpl.variables['Time'][:].size)

	#_calculate epoch from Year and Julian float
	year = int(h_cpl.variables['Year'][:])
	jday = h_cpl.variables['Time'][idx]
	time = j2e(year, jday)

	#_add data to recarray
	cpl.resize(time.size)
	cpl.epoch[:]		= time 
	cpl.latitude[:]		= h_cpl.variables['Lat'][idx]
	cpl.longitude[:]	= h_cpl.variables['Lon'][idx]
	return cpl 



def read_cpl(fname, idx=None, dir_cpl=JOIN(DIR_PROD, 'cpl', 'nc'),
	max_layers=10, manual_tau=True, cad_thresh=0.2, **kwargs):
	'''
	fname	str,		filename of cpl data in netcdf format in dir_cpl
	dir_cpl	str,		path to cpl data files in netcdf format
	idx		ndarray,	collocated indices with SHIS for flight segment
	'''
	from netCDF4 import Dataset
	from numpy import ndarray, recarray, arange
	from libtools import julian2epoch as j2e 
	from numpy.ma import masked_where
	from scipy.integrate import trapz
	import re

	#_find filename for this dtg
	fname = JOIN(dir_cpl, fname)
	dbg(fname)
	h_cpl = Dataset(fname, 'r', format='NETCDF3_CLASSIC')
		
	#_CPL channels
	wavelengths = ['355', '532', '1064']

	#_read in cpl data
	dtype = [	
				('ext_355', ndarray),
				('tau_355', ndarray),
				('ext_532', ndarray),
				('tau_532', ndarray),
				('ext_1064', ndarray),
				('tau_1064', ndarray),
				('tau_type', ndarray),
				('depolarization', ndarray),
				('layertop', ndarray),
				('layerbot', ndarray),
				('epoch', 'f8'),
				('latitude','f8'),
				('longitude','f8')]
	if manual_tau:
		dtype.extend([('tau_ssec', ndarray),
					('layertop_ssec', ndarray),
					('layerbot_ssec',ndarray),
					('tau_type_ssec', ndarray)])
	cpl = recarray((0,), dtype=dtype)
	cpl.__setattr__('fname', fname)

	if idx is None:
		idx = arange(h_cpl.variables['Time'][:].size)

	#_calculate epoch from Year and Julian float
	year = int(h_cpl.variables['Year'][:])
	print h_cpl.variables['Time'][:].shape
	jday = h_cpl.variables['Time'][idx]
	time = j2e(year, jday)

	#_add data to recarray
	cpl.resize(time.size)
	cpl.epoch[:]		= time 
	cpl.latitude[:]		= h_cpl.variables['Lat'][idx]
	cpl.longitude[:]	= h_cpl.variables['Lon'][idx]
	for t in xrange(idx.size):
		i = idx[t]
		cpl.ext_355[t]	= h_cpl.variables['Extprofile_sig_ch1'][i][::-1]
		cpl.tau_355[t]	= h_cpl.variables['Layertau_ch1'][i][:]
		cpl.ext_532[t]	= h_cpl.variables['Extprofile_sig_ch2'][i][::-1]
		cpl.tau_532[t]	= h_cpl.variables['Layertau_ch2'][i][:]
		cpl.ext_1064[t]	= h_cpl.variables['Extprofile_sig_ch3'][i][::-1]
		cpl.tau_1064[t]	= h_cpl.variables['Layertau_ch3'][i][:]

		#_delete this
		if 0:
			import matplotlib
			matplotlib.use('Agg')
			from libtools import shrink_ticks
			import matplotlib.pyplot as plt
			res = re.search('(\d+.\d+).nc$', fname)
			pname = 'PROF_FROM-NC_f{0:04d}_{1}.png'.format(t, res.group(1))
			pname = os.path.join(os.environ['PRODUCTS'], 'plots', 'CAD', pname)	
			fig, ax = plt.subplots(1, 3)
			y = range(900)
	
			ax[0].plot(h_cpl.variables['Extprofile_sig_ch1'][i][::-1], y, 'k--',
				linewidth=0.3)
			ax[0].plot(cpl.ext_355[t], 'k-', linewidth=0.3)
			ax[1].plot(h_cpl.variables['Extprofile_sig_ch2'][i][::-1], y, 'k--',
				linewidth=0.3)
			ax[1].plot(cpl.ext_532[t], 'k-', linewidth=0.3)
			ax[2].plot(h_cpl.variables['Extprofile_sig_ch3'][i][::-1], y, 'k--',
				linewidth=0.3)
			ax[2].plot(cpl.ext_1064[t], 'k-', linewidth=0.3)
	
			[a.set_xlim(0, 5e-5) for a in ax]	
			[shrink_ticks(a) for a in ax]	
			plt.savefig(pname)
			print 'plotting', pname
			plt.close()

		cpl.tau_type[t]	= h_cpl.variables['Layertype'][i][:]

		cpl.depolarization[t] = h_cpl.variables['Depolorization'][i][::-1] #_sic

	#_calculate optical depth from ext
	if manual_tau:
		from scipy.signal import argrelextrema as extrema
		from numpy import greater, less, vstack, zeros

		#_loop over each field of view to sort this out			
		for i, fov in enumerate(cpl):
			imin, imax = {}, {} 

			#_get local maximum indices
			for wl in wavelengths:
				imax[wl] = extrema(getattr(fov,'ext_{0}'.format(wl)),greater)[0]
				imin[wl] = extrema(getattr(fov,'ext_{0}'.format(wl)),less)[0]

			#_loop over all wavelengths and their layers, find top and bottom
			# of layers
			idx_top, idx_bot = {}, {} 
			for key_wl in wavelengths:
			  idx_top.update({key_wl : []})
			  idx_bot.update({key_wl : []})
			  init = 0

			  #_go through peaks to associate them with a layer
			  for j, peak in enumerate(imax[key_wl]):

				#_if peak is within last layer, skip
				if peak < init:
					continue

				thresh = 1e-8
				extinction = getattr(fov, 'ext_{0}'.format(key_wl))

				#_find below the peak where it gets closest to clear
				noise = extinction[init:peak] < thresh 
				ibot = arange(extinction.size)[noise].max() + init
					
				#_find above peack where it gets closest to clear
				noise = extinction[peak:] < thresh 
				itop = arange(extinction[peak:].size)[noise].min() + peak	

				if j and (ibot - idx_top[key_wl][-1]) <= 20:
					#_if close together, merge layers
					idx_top[key_wl][-1] = itop

				else:

					idx_top[key_wl].append(itop)
					idx_bot[key_wl].append(ibot)

				#_reset initial index
				init = idx_top[key_wl][-1] 

			#_loop over just defined layers and calculate 
			# optical thickness and try to discriminate layer type
			for j, wl in enumerate(['532']):
		##	for j, wl in enumerate(wavelengths):
			
				#_initialize recarray with size n array
				extinction = getattr(fov, 'ext_{0}'.format(wl))
				cpl.tau_ssec[i] = zeros((max_layers)) - 9999. 
				cpl.tau_type_ssec[i] = [''] * max_layers
	
				#_loop over each level				
				for z, (ibot, itop) in enumerate(zip(idx_bot[wl], idx_top[wl])):

					#_integrate over extinction and put into recarray	
					tmp = extinction[ibot:itop]
					tmp = tmp[tmp > -1e-8]

					#_discriminate between cloud and aerosol
					depo = cpl[i].depolarization[ibot:itop]
					idx = (depo < 1.) * (depo > -1e-8)
					depo = depo[idx] #_filter weird crap
					cad = depo.std() < cad_thresh 
				##	cad = depo.mean() > cad_thresh 
				
					type = 'cloud' if cad else 'aerosol'

					cpl.tau_ssec[i][z] = trapz(tmp) * 30
					cpl.tau_type_ssec[i][z] = type

			#_diagnostic plot
			if 1:
				from libtools import epoch2iso
				import matplotlib
				matplotlib.use('Agg')
				import matplotlib.pyplot as plt
				from libtools import shrink_ticks
				
				#_for optical depth on figure
				props = dict(facecolor='white', alpha=0.99, edgecolor='none')

			##	fig_shape = (5, 4)
				fig_shape = (3, 4)
				fig = plt.figure()
				axes = []	
				fig.suptitle(epoch2iso(fov.epoch))
				for j, wl in enumerate(wavelengths):

					ax = plt.subplot2grid(fig_shape, loc=(0,j), rowspan=2)
					axes.append(ax)

					#_init draw area and values
					x = getattr(fov, 'ext_{0}'.format(wl))
					y = arange(x.size)
		
					#_plot line and extremes
					ax.plot(x, y, linewidth=0.3)
					ax.scatter(x[imin[wl]], y[imin[wl]], marker='x')
					ax.scatter(x[imax[wl]], y[imax[wl]], marker='d')

					#_draw diagnosed layer lines
					for z, (ibot, itop) in enumerate(zip(idx_bot[wl],
														idx_top[wl])):
		
						#_plot a green line at bottom, red line at top of layer
						ax.plot([0, 1e-3], [itop, itop], 'r-', linewidth=0.5,
							zorder=0)
						ax.plot([0, 1e-3], [ibot, ibot], 'g--', linewidth=0.5,
							zorder=1)

						#_we calc'd above, use that for 532	
						if wl != '532':	
							tmp = x[ibot:itop]
							tmp = tmp[tmp > -1e-8]
							tau = trapz(tmp) * 30
							txt = '{0:5.3f}'.format(tau)
							xxt = .65e-4
						else: 
							tau = cpl[i].tau_ssec[z]
							typ = cpl[i].tau_type_ssec[z]
							txt = '{0:5.3f},{1}'.format(tau, typ)
							xxt = .25e-4

						#_label layer
						ax.text(xxt, (itop+ibot)/2., txt, size='xx-small',
							bbox=props)

					from numpy import linspace
					ax.set_title('EXT {0}'.format(wl))	
					ax.set_ylim(0, y.size)
					ax.set_xlim(0, 1e-4)
					xticks = linspace(0, 1e-4, 5)
					xlabs = ['{0:3.2e}'.format(a) for a in xticks] 
					ax.set_xticks(xticks)
					ax.set_xticklabels(xlabs, rotation=15)
					
					if j:
						ax.yaxis.set_visible(False)

				#_what does depol look like?
				ax = plt.subplot2grid(fig_shape, loc=(0,3), rowspan=2)
				for ibot, itop in zip(idx_bot['532'], idx_top['532']):
					ax.plot([0, 1e-3], [itop, itop], 'r-', linewidth=0.5,
						zorder=0)
					ax.plot([0, 1e-3], [ibot, ibot], 'g--', linewidth=0.5,
						zorder=1)

					tmp = fov.depolarization[ibot:itop]
					idx = (tmp < 1.) * (tmp > -1e-8)
					print tmp.shape
					print idx.shape
					print 'TESTING, WALTER'
					tmp = tmp[idx]
					dep = tmp.std()
					mnn = tmp.mean()
				##	dep = tmp.mean()
					ax.text(0.5, (itop+ibot)/2., '{0:5.3f}, {1:5.3f}'.format(dep, mnn),
						size='xx-small', bbox=props)

				ax.plot(fov.depolarization, y, linewidth=0.3)
				ax.set_title('DEPOL')
				ax.set_xlim(0, 1)
				axes.append(ax)

				#_plot cpl
				from matplotlib.cm import spectral, jet
			##	for c, wl in enumerate(wavelengths):
				for c, wl in enumerate(['532']):
					ax = plt.subplot2grid(fig_shape, loc=(2+c,0), colspan=4)
		
					extinct = vstack(getattr(cpl, 'ext_{0}'.format(wl)))		
					ax.pcolormesh(extinct.T, vmin=4e-7, vmax=1e-4, cmap=spectral)
					ax.plot([i,i], [0,900], 'r-', linewidth=1)
					ax.set_ylabel('{0} EXT'.format(wl), size='xx-small')
					ax.set_xlim(0, cpl.size)
					axes.append(ax)	

				#_finish up	
				[shrink_ticks(ax) for ax in axes]
	
				res = re.search('(\d+.\d+).nc$', fname)
				pname = 'cadstd_algo_{0}_f{1:04d}.png'.format(res.group(1),i)
				dir_plot = os.path.join(os.environ['PRODUCTS'], 'plots', 'CAD')
				pname = os.path.join(dir_plot, pname)	
				if not os.path.exists(dir_plot):
					from libtools import mkdir_p
					mkdir_p(dir_plot)
				
			##	if os.path.exists(pname): os.unlink(pname)
				plt.savefig(pname)
				dbg(pname)
				plt.close()

	return cpl 


def read_profile_gdas(fname_shis, idx=None, dir_shis=DIR_HS3, pcfilt=True,
	path_gdas=os.path.join(os.environ['PRODUCTS'], 'gdas_prof'), dummy=False,
	**kwargs):
	'''
	This one uses output from Greg Quinn's gdas2hdf script

	Read in atmospheric profile data from HS3 SHIS retrievals
	and GDAS fields, collocated separately.

	path_gdas	str,		path to gdas1 files output from gdas2hdf
	idx			ndarray,	index array for matched fovs with SHIS

	The time matching in this is pretty stupid.  Looks for roughly the
	nearest file within a 6 hour window.
	'''
	from numpy.ma import masked_where
	from libgeo import ppmv2kg, p2z, rh2wv
	from libclass import var
	from numpy import arange, array, recarray, zeros, append, tile
	import re
#	import h5py
	from pyhdf.SD import SDC, SD
	from netCDF4 import Dataset
	from libtools import newdtg2, dtg2epoch, epoch2dtg

	#_open hd5 file
	dbg(fname_shis)
	tmp			= re.search('_rdr(\d{8})T(\d{6})', fname_shis)
	shis_dtg	= tmp.group(1) + tmp.group(2)

	#_dtg for today and tomorrow for gdas files
	dtg0, dtg1 = shis_dtg, newdtg2(shis_dtg, 24, full=True)
	dtgs_gdas = [ '{0}{1:02d}'.format(dtg_loop[:8], hr_loop)	\
				for dtg_loop in [dtg0, dtg1]					\
				for hr_loop in [0,6,12,18]]

	#_kludge for dummy profs
	if dummy:
		dtgs_gdas = [dtgs_gdas[0]]

	epochs_gdas = [ dtg2epoch(dtg_loop) for dtg_loop in dtgs_gdas ]

	#_open atm.h5 file and get flight segment indices
	fname_pro = re.sub('rad{0}.nc'.format('_pcfilt'*pcfilt),
		'atm_prof_rtv_bc.h5', fname_shis)
	fname_pro = JOIN(dir_shis, fname_pro)
	fname_shis = JOIN(dir_shis, fname_shis)
	dbg(fname_shis)
	dbg(fname_pro)
	if os.path.exists(fname_pro):
	#	with h5py.File(fname_pro, 'r') as hdf:
		with SD(fname_pro, SDC.READ) as hdf:
		#	sidx = hdf.get('Flight_Segment')[:] != 0
			sidx = hdf.select('Flight_Segment')[:] != 0
	else:
		with Dataset(fname_shis, 'r') as h_rad:
			#_if not specified, use all
			sidx = tile(True, (len(h_rad.dimensions['time']))) 

	#_if no idx specified for matching, load all
	idx = arange(sidx.sum()) if idx is None else idx
		
	#_do a stupid conversion because there's a bug in pyhdf
	from numpy import arange
	sidx = arange(sidx.size)[sidx]

	#_initialize returned recarray
	dtype = [	('latitude', 'f8'),
				('longitude', 'f8'),
				('epoch', 'f8'),
				('sfc_temperature', 'f4'),
				('ozone_mixing_ratio', ndarray),
				('relative_humidity', ndarray),
				('geopotential_height', ndarray),
				('temperature', ndarray)	]
	data = recarray((idx.size,), dtype)

	#_get time of profiles from radiance file
	with Dataset(fname_shis, 'r') as h_rad:
		#_...
		date = [dtg2epoch('20{0:06.0f}00'.format(a)) 
				for a in h_rad.variables['date'][idx]]
		time = h_rad.variables['refTimeSec'][idx] 
		epochs_shis	= array(date) + time

	#_open all potential gdas files
	hdfs_gdas = []
	sfc, lat, lon, ozn, tmp, hum, pot = [], [], [], [], [], [], []
	for gdas_dtg in dtgs_gdas:
		arg = (shis_dtg, gdas_dtg) 	
		fname_gdas = JOIN(path_gdas,'gdas_prof{0:14s}.{1:10s}.hdf'.format(*arg))
		dbg(fname_gdas)

		#_open handle and add to list
		handle = SD(fname_gdas, SDC.READ)
		temp = handle.select('Surface_Temperature')[:]
		sfc.append(temp[sidx])
		temp = handle.select('Latitude')[:]
		lat.append(temp[sidx])
		temp = handle.select('Longitude')[:]
		lon.append(temp[sidx])
		temp = handle.select('Ozone_Mixing_Ratio')[:]
		ozn.append(temp[sidx,:])
		temp = handle.select('Temperature')[:]
		tmp.append(temp[sidx,:])
		temp = handle.select('Relative_Humidity')[:]
		hum.append(temp[sidx,:])
		temp = handle.select('Geopotential_Height')[:]
		pot.append(temp[sidx,:])
		hdfs_gdas.append(handle)
	#	handle = Dataset(fname_gdas, 'r')
	#	sfc.append(handle.variables['Surface_Temperature'][sidx])
	#	lat.append(handle.variables['Latitude'][sidx])
	#	lon.append(handle.variables['Longitude'][sidx])
	#	ozn.append(handle.variables['Ozone_Mixing_Ratio'][sidx,:])
	#	tmp.append(handle.variables['Temperature'][sidx,:])
	#	hum.append(handle.variables['Relative_Humidity'][sidx,:])
	#	pot.append(handle.variables['Geopotential_Height'][sidx,:])
	#	hdfs_gdas.append(handle)

	#_loop over every shis time, use associated GDAS profile
	for j, epoch_shis in enumerate(epochs_shis):

		#_pull out appropriate collocated index
		i = idx[j]

		#_find closest timestamp between observation and GFS
		sec_diff = abs(epoch_shis - epochs_gdas)
		eidx = sec_diff.argmin()
		if sec_diff[eidx] > 3600 * 4. and not dummy:
			raise RuntimeError, 'Time between GFS and measurement > 4 hr'
		hdf = hdfs_gdas[eidx]

		#_add data into recarray
		data.sfc_temperature[j]		= sfc[eidx][i]
		data.latitude[j]			= lat[eidx][i]
		data.longitude[j]			= lon[eidx][i] 
		data.epoch[j]				= epochs_gdas[eidx]
	
		data.ozone_mixing_ratio[j]	= ozn[eidx][i,:]
		data.temperature[j]			= tmp[eidx][i,:]
		data.relative_humidity[j]	= hum[eidx][i,:] 
		data.geopotential_height[j] = pot[eidx][i,:]
	
	#_set pressure once since it is static
	setattr(data, 'fname', fname_gdas[0])
	setattr(data, 'pressure', hdf.select('Pressure')[:])
#	setattr(data, 'pressure', hdf.variables['Pressure'][:])

	#_close out filehandles
##	[hdf.close() for hdf in hdfs_gdas]
	return data


def read_cod(flight, dtg, fidx=None, out_label=None, thresh=0.3, 
	dir_lbl=os.path.join(DIR_PROD, 'LBL-RTM_hs3'), experiment='HS3', **kwargs):
	'''
	flight		Flight_segment object
	dtg 
	Read in cloud optical thicknesses for CPL and SHIS
	Plot 532 backscatter of CPL
	Have two scatter dots respresenting if CPL and SHIS identified cloud layers

	'''
	from lblrtm_utils import Lbldis_input
	import re
	from numpy import recarray

	#_get flight information
	d0, d1 = flight.dtg0, flight.dtg1

	#_build base date dir structure
	dir_dtg = '{3}/{0}_{1}.{2}_fov{{0:06d}}'.format(experiment, d0, d1, d0[2:])

	#_create regular expression for water cloud ssp file
	reg = re.compile('ssp_db.mie_wat.gamma_sigma_0p100')

	#_initialize recarray
	dtype = [('epoch', 'f8'), ('fidx', 'i4'), 
			('COD_SHIS', 'f4'), ('COD_CPL', 'f4')]

#	#_get full length of flight if no indices passed
	if fidx is None:
		nf = len(flight)
		fidx = range(nf)
	else:
		nf = len(fidx)
	data = recarray((nf,), dtype)

	def _get_cpl_cod(flight, fidx):
		from numpy.ma import masked_where
		from numpy import array

		''' pull out field of views cloud optical depth from CPL '''
		fov = flight[f]

		#_create string of text containing CPL aod values for layers
		tau_532     = fov.CPL_tau_532[:]
		tau_type    = fov.CPL_tau_type[::2]
		tau_type_str = [None, 'PBL', 'AEROSOL', 'CLOUD']
		tau_532 = masked_where(tau_532 <= 0, tau_532)

		#_go down where not masked.... I guesss?
		cod = 0.
		lll = min(10-tau_532.mask.sum(), 5) 
		for ii in xrange(lll):
			try:
			#	print fov.CPL_tau_532
			#	print fov.CPL_tau_type 
			#	print tau_type[ii], tau_532[ii]
				if tau_type_str[tau_type[ii]] == 'CLOUD':
					cod += tau_532[ii] 
			except IndexError:
				print 'FAILED ON THIS'
				print ii
				print tau_532
				print tau_type
				print fov.CPL_tau_532
				print fov.CPL_tau_type
				raise IndexError, 'FAILED ON THIS'

		return cod

	#_loop over fields of view
	for i, f in enumerate(fidx):
		#_pull out flight time
		epoch = flight[f].CPL_epoch

		#_read in lbldis_input for this file
		fname = 'lbldis_input.{0}.final'.format(out_label)	
		lbl_input = os.path.join(dir_lbl, dir_dtg.format(f), fname) 

		try:
			#_read in layers
			lbl = Lbldis_input(lbl_input)

			cod = 0.
			for j, layer in enumerate(lbl.layers):
				#_pull out current layers ssp
				ssp = lbl.ssp_db[layer['dbnum']]
			
				#_if it's the water cloud, include its OD in total
				if reg.search(ssp):
					cod += sum(layer['tau'])

		except IOError:
			#_missing input
			cod = -9999.	

		#_add values for cod to record array for SHIS and CPL
		data[i] = (epoch, f, cod, _get_cpl_cod(flight, f))	

	return data

	
def read_shis(fname, idx=None, dir_shis=DIR_HS3, average=False, pcfilt=True,
	**kwargs):
	'''
	fname	str,		filename in dir_shis for SHIS_rdr file
	idx		ndarray,	index array for collocated SHIS/CPL data
	average	bool,		average SHIS retrievals

	Read in raw SHIS_rdr files.  If idx is not defined, all read in.		
	'''
	from netCDF4 import Dataset
	from numpy.ma import masked_where
	from numpy import ndarray, recarray, arange, array, isnan, tile
	from libtools import dtg2epoch as d2e
	import re
	from time import sleep
	from pyhdf.SD import SD, SDC

	#_open file handle
	fname = JOIN(dir_shis, fname) if not os.path.exists(fname) else fname
	fname_pro = re.sub('rad{0}.nc'.format('_pcfilt'*pcfilt),
		'atm_prof_rtv_bc.h5', fname)
	#_principle component filtered and bias corrected

	dbg(fname)
	dbg(fname_pro)
	h_rad = Dataset(fname, 'r', format='NETCDF3_CLASSIC')

	if os.path.exists(fname_pro):
		n_pro = SD(fname_pro, SDC.READ)

		#_pull out profile segments from flight
		sidx = n_pro.select('Flight_Segment')[:] != 0
	else:
		dbg("No profile from SHIS for this segment")
		n_pro = False

		#_should be entire data fields
		sidx = arange(h_rad.variables['date'][:].size)

	#_get location data
	dtype = [	('radiances', ndarray),
				('relative_humidity', ndarray),
				('temperature', ndarray),
				('sfc_temperature', 'f4'),
				('ozone_mixing_ratio', ndarray),
				('epoch', 'f8'), 
				('latitude','f8'),
				('longitude','f8'),
				('altitude', 'f4'),
				('cloud_top_temperature', 'f4'),
				('cloud_top_pressure', 'f4'),
				('cloud_optical_depth', 'f4')	]
	shis = recarray((0,), dtype)
	shis.__setattr__('fname', fname)

	idx = arange(h_rad.variables['date'][:].size) if idx is None else idx

	#_read in time
	date	= [d2e('20{0:06.0f}00'.format(a)) 
				for a in h_rad.variables['date'][idx]]
	time	= h_rad.variables['refTimeSec'][idx] 
	epoch	= array(date) + time
	shis.resize(epoch.size)

	#_field of view angle from nadir, left is neg
	setattr(shis, 'wavenumber', h_rad.variables['wavenumber'][:])
	setattr(shis, 'hbb_nesr', h_rad.variables['HBB_NESR'][:])
	setattr(shis, 'angle', h_rad.variables['FOVangle'][:])
	fov = h_rad.variables['FOVangle'][:]

	shis.epoch[:]		= epoch
	shis.latitude[:]	= h_rad.variables['Latitude'][idx]
	shis.longitude[:]	= h_rad.variables['Longitude'][idx]
	shis.altitude[:]	= h_rad.variables['Altitude'][idx] / 1e3 

	#_dont' always have n_pro
	if n_pro:
		setattr(shis, 'pressure', n_pro.select('Plevs')[:])
		shis.sfc_temperature[:]			= n_pro.select('TSurf')[sidx][idx]
		shis.cloud_top_temperature[:]	= n_pro.select('CTT')[sidx][idx]
		shis.cloud_top_pressure[:]		= n_pro.select('CTP')[sidx][idx]
		shis.cloud_optical_depth[:]		= n_pro.select('COT')[sidx][idx]
		tem = n_pro.select('TAir')[:,sidx]
		rel = n_pro.select('RelHum')[:,sidx]
		ozn = n_pro.select('O3VMR')[:,sidx]
	else:
		dbg('ASSUMING PRESSURE LEVELS!!! NO ATM FILE!!!')
		plevs = array([0.005, 0.0161, 0.0384, 0.0769, 0.137, 0.2244, 0.3454,
		0.5064, 0.714, 0.9753, 1.2972, 1.6872, 2.1526, 2.7009, 3.3398, 4.077,
		4.9204, 5.8776, 6.9567, 8.1655, 9.5119, 11.0038, 12.6492, 14.4559,
		16.4318, 18.5847, 20.9224, 23.4526, 26.1829, 29.121, 32.2744, 35.6505,
		39.2566, 43.1001, 47.1882, 51.5278, 56.126, 60.9895, 66.1253, 71.5398,
		77.2396, 83.231, 89.5204, 96.1138, 103.0172, 110.2366, 117.7775,
		125.6456, 133.8462, 142.3848, 151.2664, 160.4959, 170.0784, 180.0183,
		190.3203, 200.9887, 212.0277, 223.4415, 235.2338, 247.4085, 259.9691,
		272.9191, 286.2617, 300, 314.1369, 328.6753, 343.6176, 358.9665,
		374.7241, 390.8926, 407.4738, 424.4698, 441.8819, 459.7118, 477.9607,
		496.6298, 515.72, 535.2322, 555.1669, 575.5248, 596.3062, 617.5112,
		639.1398, 661.192, 683.6673, 706.5654, 729.8857, 753.6275, 777.7897,
		802.3714, 827.3713, 852.788, 878.6201, 904.8659, 931.5236, 958.5911,
		986.0666, 1013.948, 1042.232, 1070.917, 1100.])
		sleep(4)
		setattr(shis, 'pressure', plevs)
		shis.sfc_temperature[:]			= -9999 
		shis.cloud_top_temperature[:]	= -9999 
		shis.cloud_top_pressure[:]		= -9999
		shis.cloud_optical_depth[:]		= -9999
		tem = tile(-9999, (plevs.size, sidx.size))
		rel = tile(-9999, (plevs.size, sidx.size))
		ozn = tile(-9999, (plevs.size, sidx.size))
	
	rad = h_rad.variables['radiance'][:]

	#_welcome to the pain zone
	if average:
		rad = shis_average(rad, fov, idx=idx)

	for n in xrange(idx.size):
		i = idx[n]
		shis.radiances[n]			= rad[i]
		shis.temperature[n]			= tem[:,i]		
		shis.relative_humidity[n]	= rel[:,i]		
		shis.ozone_mixing_ratio[n]	= ozn[:,i]		

	return shis


def shis_average(rad, fov, idx=None, nlen=5, nwid=5, **kwargs): 
	'''
	nlen		int		along track size (inc. nadir)
	nwid		int		along sweet size (inc. nadir)
	collocate	bool	average at collocation points

	creates an average spectral radiance from surrounding
	measurements, from an area NLENxNWID

	WARNING!  Original version would drop ones that were
				too close to the segment boundary.  This
				check is no longer performed.
	'''
	from libtools import strictly_decreasing as mono
	from numpy import array, append, recarray, arange, zeros, vstack

	rad_new	= rad.copy()
	nf, nw	= rad.shape

	#_calc how far from nadir in each direction to go
	nalong	= int(nlen) / 2	#_intentionally rounding
	nwing	= int(nwid) / 2

	#_this should be constant at 14, but meh
	s, e = 0, 2 
	while(mono(fov[:e])):
		e += 1
	else:
		nangle = e - 1

	#_total number of sweeps
	nsweep = fov.size / nangle 

	#_find nadir and sanity check
	fov		= vstack(fov).reshape(nsweep, nangle)
	nadir	= abs(fov).argmin(axis=1)
	if nadir.std() > 1e-5:
		raise RuntimeError, 'nadir should be at same index throughout'
	n = nadir[0]

	#_get list of indices
	idx = idx if idx is not None else arange(n, nf+1, nangle)

	#_loop over each idx
	for i in idx:
		#_index range for each sweep, 5 sweeps
		idxn = array([i + n*nangle for n in arange(-nalong, nalong+1)]) 

		#_make total array in indices
		idxn = [ arange(idxn[n]-nwing, idxn[n]+nwing+1) for n in arange(nwid) ]
		idxn = array(idxn)
		idxn = idxn[idxn - nwing > 0]
		idxn = idxn[idxn + nwing + 1 < nf - 1]

		#_pull out sweeps + side lobs and put together
		r = rad[idxn]

		#_put into recarray to return 
		rad_new[i] = r.mean(axis=0)

	return rad_new


#############################################################################_80
#_METHODS_###################################################################_80
#############################################################################_80

##_commented WRS 2015.05.13 (maybe breaks?)
def find_file_shis(dtg, full=False, pcfilt=True, file_type='_rad_pcfilt.nc',
	path_hs3=DIR_HS3):
	'''
	given an dtg+minutes, find containing HS3 file
	dtg	str	YEARMONTHDAYHOURMINSEC
	'''
	from libtools import dtg2epoch as d2e
	from libtools import epoch2dtg as e2d
	from libtools import unique
	import re, glob
	from numpy import array, append

	#_read in all available hours (what is last time?)
	reg = 'rdr(\d{8})T(\d{6})end(\d{8})T(\d{6})sdr(\d{8})T(\d{6})'
	re_hs3 = re.compile(reg)

	#_collect start and end times for available files
	estart = []
	eend = []
	eunk = []
	for fname in glob.glob(path_hs3 + '/*'):
		res = re_hs3.search(fname)

		#_if matches, add date/times
		if res:
			start = d2e(res.group(1) + res.group(2), full=True)
			end = d2e(res.group(3) + res.group(4), full=True)
			unk = d2e(res.group(5) + res.group(6), full=True)

			estart.append(start)
			eend.append(end)
			eunk.append(unk)

	#_convert to numpy arrays
	estart = array(estart)
	eend = array(eend)
	eunk = array(eunk)

	#_find if the desired time lies between any of these
	epoch = d2e(dtg, full=True)
	idx = (epoch >= estart).reshape(1, estart.size)
	idx = append(idx, (epoch <= eend).reshape(1, eend.size), axis=0)
	idx = idx.all(axis=0)

	if idx.sum() == 0:	#_time not available
		dbg(('time not available shis', dtg))
		return False

	#_pull out times
	start = unique(estart[idx], unique=True)		
	end = unique(eend[idx], unique=True)		
	unk = unique(eunk[idx], unique=True)		

	#_rebuild prefix name and return it, or use index on a globbed list
	dtgs = e2d(start, full=True)
	dtge = e2d(end, full=True)
	dtgu = e2d(unk, full=True)
	f = 'SHIS_rdr{0:0>8s}T{1:0>6s}end{2:0>8s}T{3:0>6s}sdr{4:0>8s}T{5:0>6s}'
	f = f.format(dtgs[:8], dtgs[8:], dtge[:8], dtge[8:], dtgu[:8], dtgu[8:])

	#_not really sure... asking for full filename to be returned?
	if full:
		f += '_rad{0}.nc'.format('_pcfilt'*pcfilt) # file_type
		f = JOIN(path_hs3, f)

	return f


def find_file_cpl(dtg, file_fmt='nc', path_cpl=JOIN(DIR_PROD, 'cpl')):
	'''
	given an dtg+minutes, find containing CPL file
	dtg	str	YEARMONTHDAYHOURMINSEC
	'''
	from libtools import dtg2epoch as d2e
	from libtools import epoch2dtg as e2d
	from libtools import unique
	import re, glob
	from numpy import array, append

	#_read in all available hours (what is last time?)
	reg = 'OP_\d{5}\w?_(\d{6})_(\d{6})_(\d{6}).nc'
	re_cpl = re.compile(reg)
	re_xdr = re.compile('OP_(\d{5}\w?)_')

	#_get lsits of avail files
	files_cdf = array(glob.glob(path_cpl + '/nc/*.nc'))

	#_convert desired time to epoch seconds
	epoch = d2e(dtg, full=True)
	
	#_collect start and end times for available files
	time_dict = {}
	for fname in files_cdf:
		res = re_cpl.search(fname)

		#_if matches, add date/times
		if res:
			dtg0 = '20{0:6s}{1:6s}'.format(res.group(1), res.group(2))
			dtg1 = '20{0:6s}{1:6s}'.format(res.group(1), res.group(3))
			e0 = d2e(dtg0, full=True)
			e1 = d2e(dtg1, full=True)

			#_dumb way to do a day crossover check
			if e1 < e0:
				e1 += 86400.

			if epoch >= e0 and epoch <= e1:
				break

	#_chnage to desired file format
	if file_fmt == 'xdr':
		fname = cpl_nc2xdr(fname)

	#_pull out file name
	return fname


def cpl_nc2xdr(fname, dir_cpl=JOIN(DIR_PROD, 'cpl'), **kwargs):
	'''
	the xdr files had a weird way of keeping track of legs
	so i'm using this stupid method to match them to netcdf
	'''
	from glob import glob
	import re

	fname = JOIN(dir_cpl, 'nc', fname)

	#_get lists of files
	files_cdf = glob(dir_cpl + '/nc/O*nc')
	files_xdr = glob(dir_cpl + '/xdr/O*xdr')

	re_xdr = re.compile('OP_(\d{5}\w?)_')

	#_get julian date key - if not there, should crash
	key = re_xdr.search(fname).group(1)

	#_build list associating cdf files to key
	for file_xdr in files_xdr:
		res = re_xdr.search(file_xdr)

		#_if matching key, return
		if res and res.group(1) == key:
			return file_xdr


def cpl_xdr2nc(fname, dir_cpl=JOIN(DIR_PROD, 'cpl'), **kwargs):
	'''
	the xdr files had a weird way of keeping track of legs
	so i'm using this stupid method to match them to netcdf
	'''
	from glob import glob
	import re

	fname = JOIN(dir_cpl, 'xdr', fname)

	#_get lists of files
	files_cdf = glob(dir_cpl + '/nc/O*nc')
	files_xdr = glob(dir_cpl + '/xdr/O*xdr')

	re_xdr = re.compile('OP_(\d{5}\w?)_')

	#_get julian date key - if not there, should crash
	key = re_xdr.search(fname).group(1)

	#_build list associating cdf files to key
	for file_cdf in files_cdf:
		res = re_xdr.search(file_cdf)

		#_if matching key, return
		if res and res.group(1) == key:
			return file_cdf
	

def find_shiscplgdas(dtg, dir_shis=DIR_HS3, **kwargs):
	'''
	intput
		dtg		str(14),	YYYYMMDDhhmmss
		dir_shis	str,		path to collocation files

	output
		fname	str,		path to file if one contains that dtg, else False

	Assumes no overlap and that file names are accurate
	'''
	from glob import glob
	from numpy import array
	import re
	from libtools import dtg2epoch as d2e

	#_get list of available files
	files_scg = glob(dir_shis + '/*')
	files_scg.sort()

	#_convert dtg to epoch time
	e = d2e(dtg, full=True)
	re_dtg = re.compile('SHIS.CPL.GDAS.COLLOC.(\d{12}).(\d{12}).nc')

	for fname in files_scg:
		res = re_dtg.search(fname)
		if res:
			dtg0, dtg1 = '20' + res.group(1), '20' + res.group(2)
			#_check if date in file range
			if d2e(dtg0, full=True) <= e and d2e(dtg1, full=True) >= e:
				return fname
			
	else:
		return False


#############################################################################_80
#_CLASSES_###################################################################_80
#############################################################################_80


from numpy import arange, recarray, ndarray
class Flight_segment(recarray):
	''' MAKE THIS ALL MORE GENERAL FOR PROFILING '''
	def __new__(cls, dtg=None, file_seg=None, **kwargs):
		from netCDF4 import Dataset
		from numpy import array, recarray#, masked_where

		#_find the file we want
		if dtg is None and file_seg is None:
			dbg('specify desired segment time or file')
			return False
		elif file_seg is None and dtg is not None:
			file_seg = find_shiscplgdas(dtg, **kwargs)
			dbg((dtg, file_seg))
			if not file_seg:
				return False 

		with Dataset(file_seg, 'r') as cdf:
			#_read in dimensional data
			nfov	= len(cdf.dimensions['fov'])
			dtype	= [	('CPL_tau_532', ndarray),
						('CPL_tau_type', ndarray),
						('CPL_ext_532', ndarray),
						('CPL_latitude', 'f8'),
						('CPL_longitude', 'f8'),
						('CPL_epoch', 'f8'),
						('SHIS_radiances', ndarray),
						('SHIS_latitude', 'f8'),
						('SHIS_longitude', 'f8'),
						('SHIS_epoch', 'f8'),
						('SHIS_altitude', 'f4'),
						('SHIS_CTP', 'f4'),
						('SHIS_CTT', 'f4'),
						('SHIS_COT', 'f4'),
						('SHIS_temperature', ndarray),
						('SHIS_relative_humidity', ndarray),
						('SHIS_ozone_mixing_ratio', ndarray),
						('SHIS_sfc_temperature', 'f4'),
						('GDAS_temperature', ndarray),
						('GDAS_relative_humidity', ndarray),
						('GDAS_ozone_mixing_ratio', ndarray),
						('GDAS_sfc_temperature', 'f4'),
						('GDAS_latitude', 'f8'),
						('GDAS_longitude', 'f8'),
						('SHIS_tau', 'f4'),
						('SHIS_ref', 'f4')	]
			obj	= recarray.__new__(cls, nfov, dtype=dtype)

			def strictly_dec(A):
				return all(x < y for x, y in zip(A, A[1:]))
			def strictly_inc(A):
				return all(x > y for x, y in zip(A, A[1:]))

			#_if pressure values are not decreasing, reverse values
			rg = 1 if strictly_inc(cdf.variables['GDAS_pressure'][:]) else -1
			rs = 1 if strictly_inc(cdf.variables['SHIS_pressure'][:]) else -1
			obj.CPL_latitude[:]		= cdf.variables['CPL_latitude'][:]
			obj.CPL_longitude[:]	= cdf.variables['CPL_longitude'][:]
			obj.CPL_epoch[:]		= cdf.variables['CPL_epoch'][:]
			obj.SHIS_latitude[:]	= cdf.variables['SHIS_latitude'][:]
			obj.SHIS_longitude[:]	= cdf.variables['SHIS_longitude'][:]
			obj.SHIS_epoch[:]		= cdf.variables['SHIS_epoch'][:]
			obj.SHIS_altitude[:]	= cdf.variables['SHIS_altitude'][:]
			obj.SHIS_sfc_temperature[:]	\
								= cdf.variables['SHIS_surface_temperature'][:]
			try:
				obj.SHIS_CTT[:] = cdf.variables['SHIS_CTT'][:]
				obj.SHIS_COT[:] = cdf.variables['SHIS_COT'][:]
				obj.SHIS_CTP[:] = cdf.variables['SHIS_CTP'][:]
			except KeyError:
				pass #_meh

			obj.GDAS_sfc_temperature[:]	\
								= cdf.variables['GDAS_surface_temperature'][:]
			obj.GDAS_latitude[:]	= cdf.variables['GDAS_latitude'][:]
			obj.GDAS_longitude[:]	= cdf.variables['GDAS_longitude'][:]
			
			#_default these to missing, fill with separate method
			obj.SHIS_tau[:]	= -9999.
			obj.SHIS_ref[:]	= -9999.
	
			'''
			wvn = cdf.variables['SHIS_wavenumber'][:]
			hbb = cdf.variables['SHIS_HBB_NESR'][:]
			spr = cdf.variables['SHIS_pressure'][::rs]
			gpr = cdf.variables['GDAS_pressure'][::rg]
			setattr(obj[n], 'SHIS_wavenumber', wvn)
			setattr(obj[n], 'SHIS_HBB_NESR', hbb)
			setattr(obj[n], 'SHIS_pressure', spr)
			setattr(obj[n], 'GDAS_pressure', gpr)
			'''
			for n in xrange(nfov):
				obj[n].CPL_tau_532	= cdf.variables['CPL_tau_532'][n,:]
				obj[n].CPL_tau_type = cdf.variables['CPL_tau_type'][n,:]
				obj[n].CPL_ext_532	= cdf.variables['CPL_ext_532'][n,:]

				obj[n].SHIS_radiances = cdf.variables['SHIS_radiances'][n,:]/1e3
				obj[n].SHIS_ozone_mixing_ratio \
						= cdf.variables['SHIS_ozone_mixing_ratio'][n,::rs]
				obj[n].SHIS_temperature \
						= cdf.variables['SHIS_temperature'][n,::rs]
				obj[n].SHIS_relative_humidity \
						= cdf.variables['SHIS_relative_humidity'][n,::rs]
				obj[n].GDAS_ozone_mixing_ratio \
						= cdf.variables['GDAS_ozone_mixing_ratio'][n,::rg]*1e3
				obj[n].GDAS_temperature \
						= cdf.variables['GDAS_temperature'][n,::rg]
				obj[n].GDAS_relative_humidity \
						= cdf.variables['GDAS_relative_humidity'][n,::rg]

			#_add last variables
			setattr(obj, 'SHIS_wavenumber', cdf.variables['SHIS_wavenumber'][:])
			setattr(obj, 'SHIS_HBB_NESR', cdf.variables['SHIS_HBB_NESR'][:])
			setattr(obj, 'SHIS_pressure', cdf.variables['SHIS_pressure'][::rs])
			setattr(obj, 'GDAS_pressure', cdf.variables['GDAS_pressure'][::rg])

			setattr(obj, 'GDAS_file', cdf.file_gdas)
			setattr(obj, 'SHIS_file', cdf.file_shis)
			setattr(obj, 'CPL_file', cdf.file_cpl)
			setattr(obj, 'FLIGHT_file', file_seg)

			setattr(obj, 'dtg0', cdf.dtg0)
			setattr(obj, 'dtg1', cdf.dtg1)

			setattr(obj, 'collocated', True)

			try:
				setattr(obj, 'notes', cdf.notes)
			except:
				pass

		return obj 


	def __array_finalize__(self, obj):
		#_new object construct, skip
		if obj is None:
			return
	
		#_otherwise, copy over default
		self.info = getattr(obj, 'info', None)


	def map_cpl(self, dz_cpl=25, maxz=20000., **kwargs):
		'''
		Convert CPL data from relative to sensor to AGL
		dz_cpl	int,	vertical resolution in meters
		maxz	float,	maximum height to bin...?	
		'''
		from numpy import arange, ceil, vstack, zeros

		#_create a plotting field running from 0 - 20k
		nz = ceil(2e4 / dz_cpl) + 1
		nc = self.size 
		cpl_fill = zeros((int(nc),int(nz)))		

		cpl_x = arange(nc)

		#_transform from sensor to height
		alt = dz_cpl*arange(nz)

		for k, shot in enumerate(self.CPL_ext_532):

			#_find surface in retrieval
			idx = ((self.SHIS_altitude[k]*1e3 - alt) > 0).sum() - 1 

			#_fill in zeros from bottom of aircraft to sfc
			cpl_fill[k,:idx] = shot[-idx:]

		return cpl_fill, cpl_x, alt 


	def plot(self, fidx=None, kidx=[0], **kwargs):
		'''
		plot entire flight with model data
		
		usage:
			Flight_data.plot(**FLIGHT[dtg])

		'''
		from numpy import arange, vstack
		from libtools import combination, setup_groups

		#_don't plot uncollocated data
		if not self.collocated:
			self.collocate()

		#_plot each field of view
		try:
			cpl, cpl_x, cpl_y = self.map_cpl(**kwargs)
		except ValueError:
			from numpy import arange, linspace
			dbg('WARNING: CPL mapping failed')
			cpl_x, cpl_y = arange(self.size), linspace(0, 22500, 900)
			cpl = arange(self.size*900).reshape(self.size, 900) 
		
		#_do all if none pass
		if fidx is None:
			fidx = range(self.size)
		elif type(fidx) != list:
			fidx = [fidx]

		#_get groups to fork
		runs = combination((fidx, kidx))
		groups = setup_groups(runs, **kwargs)

		for group in groups:
			children = []
			for f, k in group:
				pid = os.fork()
				if pid != 0:
					children.append(pid)
				elif pid == 0:
					self.plot_fov(f, kidx=k, cpl_plt=cpl, **kwargs)
					os._exit(0)

			for kid in children:
				os.waitpid(kid, 0)

		
	def plot_profiles(self, fidx=None, **kwargs):
		from numpy import arange, vstack
		from libtools import setup_groups
		import matplotlib.pyplot as plt

		#_don't plot uncollocated data
		if not self.collocated:
			self.collocate()

		for i in fidx:
			p = self.profile[i]
			fig = plt.figure()
			ax0 = fig.add_subplot(211); ax1 = fig.add_subplot(212)
			ax0.plot(p.temperature_shis, p.altitude_shis)
			ax0.plot(p.temperature_gdas, p.altitude_gdas)
			ax1.plot(p.water_vapor_shis, p.altitude_shis)
			ax1.plot(p.water_vapor_gdas, p.altitude_gdas)
			pname = 'profile_{0:10s}_{1:04d}.png'.format(self.dtg, i)
			dbg((pname, len(fidx)))
			plt.savefig(pname)
			plt.close()

	def plot_flight(self, fidx=None, **kwargs):
		''' plot only flight data, really just a shortcut for plot_model=0 '''
		from numpy import arange, vstack
		from libtools import setup_groups

		#_plot each field of view
		cpl, d0, d1 = self.map_cpl(**kwargs)

		fidx = [fidx] if type(fidx) == int else fidx
	
		#_get groups to fork
		fidx	= arange(self.size) if fidx is None else fidx
		groups	= setup_groups(fidx, **kwargs)

		for group in groups:
			children = []
			for i in group:
				pid = os.fork()
				if pid != 0:
					children.append(pid)
				elif pid == 0:
					self.plot_fov(i, cpl_plt=cpl, plot_model=False, **kwargs)
					os._exit(0)

			for kid in children: os.waitpid(kid, 0)


	def plot_clear(self, fidx=None, **kwargs):
		'''
		stupid thing to just get out some comparison plots for bob

		fix the plot() function if non-clear sky cases are desired
		'''
		from numpy import arange, vstack
		from libtools import setup_groups

		#_don't plot uncollocated data
		if not self.collocated:
			self.collocate()

		#_plot each field of view
		cpl = vstack(self.cpl.ext_532)
	
		#_get groups to fork
		fidx = arange(self.shis.size) if fidx is None else fidx	
		groups = setup_groups(fidx, **kwargs)
	
		for group in groups:
			children = []
			for i in group:
				pid = os.fork()
				if pid != 0:
					children.append(pid)
				elif pid == 0:
					self.plot_fov(i, cpl_plt=cpl, plot_clear=False, **kwargs)
					os._exit(0)

			for kid in children: os.waitpid(kid, 0)


	def plot_fov(self, fov, cpl_plt=None, experiment='hs3', 
		dir_lblrtm_fmt='{3}_{0}.{1}_fov{2:06d}', dir_lblrtm=None, 
		plot_tape7=False, lbldis_output=False,
		lbldis_clear='lbldis_output.clear.cdf', 
		retrieval_output=False, zoom=True, dir_plot='./plots', dz_cpl=25., 
		out_label='hs3', rerun=True, **kwargs):
		'''
		plot field of view dashboard
	
		INPUT
		fov			int,		index of field of view
								this controls the location, SHIS
		cpl			ndarray,	Raw cpl data, used for context.
								Pass to speed up.
		dir_lblrtm_fmt		string,		format string for LBL-RTM output
		plot_model	bool,		set to false if no LBL-DIS plots desired.
								useful for plotting flight data only.
		lbldis_full	False|str,	plot the full spectrum
		lbldis_mw	False|str,	plot the microwindows
		'''
		import matplotlib.pyplot as plt
		from libnva import draw_map
		import re, pickle
		from numpy import diff, meshgrid, round, vstack
		from numpy import where, linspace, meshgrid
		from libgeo import p2z, planck_inv
		from libtools import mkdir_p, shrink_ticks, dtg2iso, epoch2iso
		from tape7 import tape7
		from lblrtm_utils import microwindow_average
		from lblrtm_utils import read_lbldis_out, Lbldis_input
		from netCDF4 import Dataset
		from scipy.interpolate import interp1d 
		from matplotlib.cm import spectral, jet, gist_ncar
		from libcmap import rgbcmap	#_calipso
		from matplotlib.colors import LogNorm
		from libgeo import ppmv2kg, p2z, rh2wv
		from numpy.ma import masked_where

		pname = 'dashboard_{2}_{0}_f{1:04d}.png'.format(self.dtg0, fov, out_label)
		pname = JOIN(dir_plot, pname)
		if os.path.exists(pname) and not rerun:
				dbg(('already generated, skipping...'))
				return

		#_field of view label string for prefixing directories 
	#	fov_str = dir_lblrtm_fmt.format(self.dtg0, self.dtg1, fov, experiment)
		fov_str = dir_lblrtm_fmt.format(self.dtg0, self.dtg1, fov, experiment,
				self.dtg0[2:])

		#_current time
		date_fov = epoch2iso(self.SHIS_epoch[fov])

		#_gets overwritten if model run
		title = 'FOV={0:04d}, {1}'.format(fov, date_fov)
		path_fov = JOIN(DIR_PROD, 'LBL-RTM', fov_str) if dir_lblrtm is None \
				else dir_lblrtm	

		dbg(('testing', path_fov))

		#_find index for cpl_full in flight
		cpl_idx = abs(self.SHIS_epoch[fov] - self.CPL_epoch).argmin()
##		cpl_plt = vstack(self.CPL_ext_532) if cpl_plt is None else cpl_plt
		if cpl_plt is None:
			cpl_plt, cpl_x, cpl_y = self.map_cpl(dz_cpl=dz_cpl, **kwargs)

		else:
			#_UNTESTED
			cpl_x = arange(self.size)
			cpl_y = arange(cpl_plt.shape[1])*dz_cpl 	

		#_index SHIS data?
		fig = plt.figure()
		shape = (3, 11)	#_halve columns (3, 8) default)
		ax_rad = plt.subplot2grid(shape, (0,3), rowspan=2, colspan=8)
		ax_co2 = plt.subplot2grid(shape, (0,0), rowspan=3, colspan=3)
		ax_cpl = plt.subplot2grid(shape, (2,3), colspan=4)
		ax_pro = plt.subplot2grid(shape, (2,7), colspan=2)
		ax_map = plt.subplot2grid(shape, (2,9), colspan=2)
	
		#_get coords
		lats, lons	= self.CPL_latitude, self.CPL_longitude
		lat, lon	= lats[cpl_idx], lons[cpl_idx]
		time		= self.CPL_epoch[cpl_idx]
		ulat, llat, ulon, llon = lat+25, lat-25, lon+25, lon-25
		corn		= [llat, ulat, llon, ulon]
	
		#_draw location map
		lx	= [0,0,0,0]
		map = draw_map(grid=corn, corners=corn, delta=10, ax=ax_map, 
						labx=lx, laby=lx)
		map.scatter(lon, lat, color='r', marker='+', zorder=1)
		map.plot(lons, lats, linewidth=0.4, color='black')
		ax_map.set_title('lat={0:7.3f}\nlon={1:7.3f}'.format(lat, lon), 
						size='xx-small') 
		
		#_draw temp profile
		prof		= self[fov].copy()
		SHIS_alt	= p2z(self.SHIS_pressure)
		GDAS_alt	= p2z(self.GDAS_pressure)
		ax_pro.plot(prof.SHIS_temperature, SHIS_alt, linewidth=0.4, color='r')
		ax_pro.plot(prof.GDAS_temperature, GDAS_alt, linewidth=0.4, color='b')
		ax_pro.set_ylim(0, 20)
		ax_pro.set_yticks(GDAS_alt[::4])
		ax_pro.set_yticklabels(round(self.GDAS_pressure[::4], 2))
		ntmp = prof.GDAS_temperature.size
		ax_pro.set_xticks(linspace(150,300,4))
		ax_pro.set_xlim(150, 300)
		ax_pro.grid(True)
		ax_pro.set_title('gdas=blue, shis=red', size='xx-small')
		ax_pro.set_xlabel('temperature (K)', size='xx-small')

		#_if also want to plot the profile for lblrtm 
		if plot_tape7:
			prof7 = tape7(JOIN(path_fov, 'TAPE7'))
			ax_pro.scatter(prof7.tbar, prof7.zlevel1, color='k', marker='o')

		#_interpolate along temperature profile and place flight X
		CTP = masked_where(self.SHIS_CTP <= 0, self.SHIS_CTP)
		CTZ = p2z(CTP)
		ctz = p2z(CTP[fov]) if CTP[fov] >= 0 else 0 
		alt = prof.SHIS_altitude * 1e3
		f = interp1d(SHIS_alt, prof.SHIS_temperature) 
		try:
			ax_pro.scatter(f(alt), alt, marker='x', s=4, color='black')
			ax_pro.scatter(f(ctz), ctz, marker='x', s=4, color='green') 
		except ValueError:
			dbg('WARNING: failed to plot profile')		

		#_plot SHIS
		arg =	{'color' : 'k', 'linewidth' : 0.3, 'zorder' : 1} 
		shis_rads	= prof.SHIS_radiances.copy()
		shis_wnum	= self.SHIS_wavenumber.copy()
		shis_Tb		= planck_inv(shis_rads, shis_wnum*100, domain='wavenumber')
		ax_rad.plot(shis_wnum, shis_Tb, **arg)
		ax_co2.plot(shis_wnum, shis_Tb, **arg)

		#_plot LBLDIS
		if lbldis_output:
			label = '.'.join(lbldis_output.split('/')[-1].split('.')[1:-1])
			lbldis = read_lbldis_out(JOIN(path_fov, lbldis_output))
			lbldis_clr = read_lbldis_out(JOIN(path_fov, lbldis_clear))
			lbldis_inp = Lbldis_input(JOIN(path_fov,
				'lbldis_input.{0}'.format(label)))

			#_silly little kludge
			plot_type = 'full' if lbldis.wavenumber.size > 100 else 'mw'

			#_plot lbldis radiances
			start = 0
			arg0 = { 'linewidth':.3, 'zorder' : 3 }
			arg1 = { 'marker' : 'x', 'zorder' : 3 }

			if plot_type == 'full':
				ax_rad.plot(lbldis.wavenumber, lbldis.Tb, color='r', **arg0)

			else:
				#_average SHIS to microwindows
				mw_shis_rd, mw_shis_wv = microwindow_average(shis_rads,
							shis_wnum, int(lbldis.resolution), error=False)
				mw_shis_Tb = planck_inv(mw_shis_rd, mw_shis_wv*100,
									domain='wavenumber')

				#_plot lbldis microwindow output on top of SHIS
				ax_rad.scatter(lbldis.wavenumber, lbldis.Tb, color='r', **arg1)
				ax_rad.scatter(lbldis.wavenumber, lbldis_clr.Tb, 
					marker='o', zorder=2, color='#00CCFF', edgecolors='b', 
					linewidth=0.5)

				#_plot differences in MW channels
				ax_dif = ax_rad.twinx()
				ax_dif.set_ylim(-10, 10)
				ax_dif.scatter(lbldis.wavenumber, mw_shis_Tb - lbldis.Tb,
						color='magenta', marker='d')#, **arg1)
				ax_dif.plot([0, 1e9], [0,0], linewidth=0.5, color='magenta')
				ax_dif.set_ylabel('${\Delta}K$', size='xx-small')

				#_scatter plot SHIS at same windows 
				ax_rad.scatter(mw_shis_wv, mw_shis_Tb, marker='o', color='k')
				ax_co2.scatter(mw_shis_wv, mw_shis_Tb, marker='o', color='k')

			'''
			#_interp SHIS to LBL-DIS microwindows
			x = clear_mw.wavenumber[clear_mw.wavenumber.mask == False]
			f = interp1d(shis_wnum, shis_Tb)
			shis_mw = f(x)

			#_pull out only unmasked LBL-DIS microwindows
			dirt_mw = dirty_mw.Tb[clear_mw.wavenumber.mask == False]

			#_plot SHIS at microwindows
			ax_rad.scatter(x, f(x), marker='x', color='r')

			#_plot difference between microwindows and shis
			ax_mw = ax_rad.twinx()
			ax_mw.scatter(x, shis_mw-dirt_mw, marker='o',
						color='magenta', s=0.75) 
			ax_mw.set_ylim(-20, 20)
			ax_mw.set_ylabel('S-HIS minus LBL-RTM.aero (K)', size='xx-small')
			'''	

		#_limit to data area
		xlim0 = ax_rad.xaxis.get_data_interval()
		ylims = (200, 300) if not zoom else (280, 300)
		xlims = (xlim0[0], 2000) if not zoom else (750, 1300)
		ax_rad.set_xlim(*xlims) #xlim0[0], 2000)
		ax_rad.set_ylim(*ylims) #200, 300)
		ax_co2.set_xlim(680, 800)
		ax_co2.set_ylim(200, 300)

		#_turnon grids
		ax_rad.grid(True)
		ax_co2.grid(True)
		ax_rad.set_xlabel('wavenumber ($cm^{-1}$)', size='xx-small')
		ax_rad.set_ylabel('brightness temperature (K)', size='xx-small')
		
		#_set title
		ax_rad.set_title(title, size='small')

		#_generate cpl colormap
		calipso = rgbcmap('calipso')

		#_plot cpl data
		cpl_max = cpl_plt.shape[1]
		cpl_nx, cpl_ny = cpl_plt.shape
		cpl_xx, cpl_yy = meshgrid(cpl_x, cpl_y)
		CTZ[CTZ <= 0] = -9999
		cb = ax_cpl.pcolormesh(cpl_x, cpl_y, cpl_plt.T, vmin=-4e-7,
			vmax=1e-4, cmap=calipso)
		ax_cpl.scatter(cpl_x, CTZ, marker='x', linewidth=0.5, s=4,
			color='yellow') 
		ax_cpl.scatter(cpl_x[fov], CTZ[fov], marker='o', s=2, color='red') 
		ax_cpl.set_title('CPL BACKSCAT 532_EXT', size='xx-small')
		xlim = ax_cpl.xaxis.get_data_interval()
		ax_cpl.set_xlim(xlim)
		ax_cpl.set_ylim(0, cpl_y[-1]) #cpl_ny)
		ax_cpl.set_yticks(linspace(0, cpl_y[-1], 11))
	##	ax_cpl.set_yticklabels(['{0:4.1f}'.format(vv) for vv in 
	##		linspace(0, cpl_ny, 11)*dz_cpl/1e3]) 
		try:
			ax_cpl.set_xticks(arange(0, cpl_nx, cpl_nx/5))
			ax_cpl.set_xticklabels([epoch2iso(ttt)[-8:] for 
				ttt in self.CPL_epoch[::cpl_nx/5]]) 
		except ZeroDivisionError:
			dbg('not sure why this error is happening')
		ax_cpl.set_ylabel('HEIGHTS TOO LOW ~1KM', size='xx-small')
		ax_cpl.plot([cpl_idx, cpl_idx], [0, cpl_y[-1]], 'r-')
		
		#_create string of text containing CPL aod values for layers
		tau_532		= self.CPL_tau_532[cpl_idx][:]
		tau_type	= self.CPL_tau_type[cpl_idx][::2]
		tau_str		= []
		tau_type_str = [None, 'PBL', 'AEROSOL', 'CLOUD']
		tau_532 = masked_where(tau_532 <= 0, tau_532)

		if retrieval_output:
			pass

		for ii in xrange(10 - tau_532.mask.sum()):
			try:
				tau_str.append('{0:7.2f} {1}'.format(tau_532[ii],
								tau_type_str[tau_type[ii]]))
			except:
				pass

		#_put retrieved value at end if available
		if lbldis_output:
			tau = {}
			z = []
 
			#_sum up aod for each type of layer 
			for layer in lbldis_inp.layers:
				try:
					tau[layer['dbnum']] += layer['tau'][0] 
				except KeyError:
					tau[layer['dbnum']] = layer['tau'][0]
				
				#_make list of heights
				z.append(layer['z'])
		
			tau_str.append(' ')
			for db, value in tau.iteritems():
			#	tau_str.append('{0:7.2f}, SSP={1:d}'.format(value, db))
				tmp = lbldis_inp.ssp_db[db].split('.')[1].split('_')[-1][:3]
				tau_str.append('{0:7.2f}, SSP={1}'.format(value, tmp))	

			#_ADD MARKER FOR USED Z_LAYER
			z = [a * 1e3 for a in z]
			ax_cpl.scatter([cpl_idx]*len(z), z, marker='o', s=4,
				color='blue', zorder=10) 

		tau_str = '\n'.join(tau_str)

		#_label CPL aod values
		diff_nx = cpl_nx / 2.8
		xtextloc = cpl_idx-diff_nx if cpl_nx-cpl_idx < diff_nx else cpl_idx+20
		props = dict(boxstyle='round', facecolor='white', alpha=.5)
		ax_cpl.text(cpl_x[10], cpl_y[-1]-1000, tau_str, color='k',
					verticalalignment='top', size='xx-small', bbox=props)

		#_shrink labels
		[shrink_ticks(xx) for xx in fig.axes]
					
		plt.tight_layout()
		mkdir_p(dir_plot)
	##	pname = 'dashboard_{2}_{0}_f{1:04d}.png'.format(self.dtg0,fov,out_label)
	##	pname = JOIN(dir_plot, pname)
		dbg(pname)
		plt.savefig(pname)
		plt.close()
	

	def plot_simulated_fov(self, fov, kidx=0, cpl_plt=None,
		fov_fmt='hs3_{0:s}_fov{1:06d}',
		plot_mw=True, plot_full=True, plot_model=True, plot_tape7=False,
		plot_clear=True, plot_dirty=True, zoom=True, broke_mode=False,
		file_kw='mw_kwarg_lbldis.pk', dir_out=DIR_HS3, **kwargs):
		'''
		THIS IS THE OLDER VERSION THAT WAS NOT INTENDED TO BE USED FOR
		RETRIEVAL CODE.
	
		[um+full clear-sky vs. um]		[map]
		[						 ]		[Tprof]
		[cpl -/+ 50 idx]		
		
		plot field of view
	
		INPUT
		fov	int,		index of field of view
							this controls the location, SHIS
		kidx	int,		index of kwarg optsions in file_kw
							this controls the r_eff, aod, layer z
		cpl		ndarray,	Raw cpl data, used for context. Pass to speed up.
		fov_fmt	string,		format string for LBL-RTM output
		plot_model	bool,	set to false if no LBL-DIS plots desired.
							useful for plotting flight data only.
		plot_full	bool,	plot the full spectrum,	(default=False)
		plot_mw		bool,	plot the microwindows,	(default=True)	
		plot_model	bool,	
		'''
		if broke_mode:
			dbg('CURRENTLY PLOTTING SINGLE MODEL OUTPUT FOR SANITY')

		import matplotlib.pyplot as plt
		from libnva import draw_map
		import re, pickle
		from numpy import diff, meshgrid, round, vstack, where, linspace, ma
		from libgeo import p2z, planck_inv
		from libtools import mkdir_p, shrink_ticks, dtg2iso, epoch2iso
		from tape7 import tape7
		from lblrtm_utils import microwindow_average
		from lblrtm_utils import read_lbldis_out, Lbldis_input
		from netCDF4 import Dataset
		from scipy.interpolate import interp1d 
		from matplotlib.cm import spectral, jet
		from libgeo import ppmv2kg, p2z, rh2wv
		from numpy.ma import masked_where

		#_read in run options
		file_kwarg = os.path.join(DIR_PROD, 'LBL-RTM', file_kw)
		kw = pickle.load(open(file_kwarg, 'rb'))[kidx]
		try:
			fov_fmt = kw['dir_lblrtm_fmt']
		except:
			pass
		#_define location for output from model for field of view
		dbg((fov_fmt, fov, self.dtg0, self.dtg1))
		fov_str = fov_fmt.format(fov, self.dtg0, self.dtg1)
	
		#_current time
		date_fov = epoch2iso(self.SHIS_epoch[fov])

		#_gets overwritten if model run
##		title = 'FOV={0:04d}, {1:s}'.format(fov, dtg2iso(self.dtg, full=True))
		title = 'FOV={0:04d}, {1:s}'.format(fov, date_fov)
	
		if plot_model:
			path_fov = JOIN(DIR_PROD, 'LBL-RTM', fov_str)	

			#_format for outputs
			fmt_f ='lbldis_out.full_r{0:4.2f}_od{1:3.1f}_z{2:3.1f}-{3:3.1f}.cdf'
			fmt_m ='lbldis_out.mw_r{0:4.2f}_od{1:3.1f}_z{2:3.1f}-{3:3.1f}.cdf'
		
			#_pull out kwarg layer case
			k_pick = JOIN(DIR_PROD, 'LBL-RTM', file_kw) 
			kw = pickle.load(open(k_pick, 'rb'))
			try:
				cld = kw[kidx]['clddef'][0]
				ref, aod, z0, z1 =	cld['r_eff'],	cld['tau'][0],\
									cld['z'],		cld['z_top']
			except KeyError:
				ref, aod, z0, z1 = 0., 0., 0., 0.		

			title = 'R_EFF={0:4.1f}, AOD={1:3.1f}, Z={2:1.0f}-{3:1.0f}km,'\
					+' FOV={4:d} {5:s}'
			title = title.format(ref, aod, z0, z1, fov, date_fov)
	
			#_read in lbldis data
			if plot_mw:
				f_dirty_mw = fmt_m.format(ref, aod, z0, z1)
				path_fov = JOIN(DIR_PROD, 'LBL-RTM', fov_str)	
				if broke_mode:
					fov_tmp = fov_fmt.format(460, self.dtg0, self.dtg1)
					path_tmp = JOIN(DIR_PROD, 'LBL-RTM', fov_tmp)	
					dirty_mw = read_lbldis_out(JOIN(path_tmp, f_dirty_mw))

				else:
					dirty_mw = read_lbldis_out(JOIN(path_fov, f_dirty_mw))
	
			if plot_full:
				f_dirty_fl = fmt_f.format(ref, aod, z0, z1)		
				dirty_fl = read_lbldis_out(JOIN(path_fov, f_dirty_fl))
	
		#_find index for cpl_full in flight
		cpl_idx = abs(self.SHIS_epoch[fov] - self.CPL_epoch).argmin()
		cpl_plt = vstack(self.CPL_ext_532) if cpl_plt is None else cpl_plt
	
		#_index SHIS data?
		fig = plt.figure()
		shape = (3, 11)	#_halve columns (3, 8) default)
		ax_rad = plt.subplot2grid(shape, (0,3), rowspan=2, colspan=8)
		ax_co2 = plt.subplot2grid(shape, (0,0), rowspan=3, colspan=3)
		ax_cpl = plt.subplot2grid(shape, (2,3), colspan=4)
		ax_pro = plt.subplot2grid(shape, (2,7), colspan=2)
		ax_map = plt.subplot2grid(shape, (2,9), colspan=2)
	
		#_get coords
		lats, lons	= self.CPL_latitude, self.CPL_longitude
		lat, lon	= lats[cpl_idx], lons[cpl_idx]
		time		= self.CPL_epoch[cpl_idx]
		ulat, llat, ulon, llon = lat+25, lat-25, lon+25, lon-25
		corn		= [llat, ulat, llon, ulon]
	
		#_draw location map
		lx	= [0,0,0,0]
		map = draw_map(grid=corn, corners=corn, delta=10, ax=ax_map, 
						labx=lx, laby=lx)
		map.scatter(lon, lat, color='r', marker='+', zorder=1)
		map.plot(lons, lats, linewidth=0.4, color='black')
##		plat, plon = lats[max(0,cpl_idx-200):cpl_idx:10], \
##						lons[max(0,cpl_idx-200):cpl_idx:10] 
##		map.scatter(plon, plat, color='k', marker='.', s=.1, zorder=0)
		ax_map.set_title('lat={0:7.3f}\nlon={1:7.3f}'.format(lat, lon), 
						size='xx-small') 
		
		#_draw temp profile
		prof		= self[fov].copy()
		SHIS_alt	= p2z(self.SHIS_pressure)
		GDAS_alt	= p2z(self.GDAS_pressure)
		ax_pro.plot(prof.SHIS_temperature, SHIS_alt, linewidth=0.4, color='r')
##		ax_pro.plot(prof.GDAS_temperature, GDAS_alt, linewidth=0.4, color='b')
		ax_pro.set_ylim(0, 20)
		ax_pro.set_yticks(GDAS_alt[::4])
		ax_pro.set_yticklabels(round(self.GDAS_pressure[::4], 2))
		ntmp = prof.GDAS_temperature.size
		ax_pro.set_xticks(linspace(150,300,6))
		ax_pro.set_xlim(150, 300)
		ax_pro.grid(True)
##		ax_pro.set_title('gdas=blue, shis=red', size='xx-small')
		ax_pro.set_title('S-HIS PROF', size='xx-small')
		ax_pro.set_xlabel('temperature (K)', size='xx-small')

		#_if also want to plot the profile for lblrtm 
		if plot_model and plot_tape7:
			prof7 = tape7(JOIN(path_fov, 'TAPE7'))
			ax_pro.scatter(prof7.tbar, prof7.zlevel1, color='k', marker='o')

		#_interpolate along temperature profile and place flight X
		alt = prof.SHIS_altitude * 1e3
		f = interp1d(SHIS_alt, prof.SHIS_temperature) 
		ax_pro.scatter(f(alt), alt, marker='x', s=4, color='black')

		if plot_model:	
			#_plot lbldis radiances
			start = 0
			arg0 = { 'linewidth':.3 }
			arg1 = { 'marker' : 'x' }

			if plot_full:
				ax_rad.plot(dirty_fl.wavenumber, dirty_fl.Tb, 
								color='r', **arg0)

			if plot_mw:
				ax_rad.scatter(dirty_mw.wavenumber, dirty_mw.Tb,
								color='r', **arg1)
##				ax_co2.scatter(dirty_mw.wavenumber, dirty_mw.Tb,
##								color='r', **arg1)
		
		#_plot SHIS
		arg =	{'color' : 'k', 'linewidth' : 0.3} 
		shis_rads	= prof.SHIS_radiances.copy()
		shis_wnum	= self.SHIS_wavenumber.copy()
		shis_Tb		= planck_inv(shis_rads, shis_wnum*100, domain='wavenumber')
		ax_rad.plot(shis_wnum, shis_Tb, **arg)
		ax_co2.plot(shis_wnum, shis_Tb, **arg)

		#_if plotting model, plot the differences between shis and ???	
		#_NOT CHANNEL AVERAGE_#
		if plot_model:
			if plot_mw:
				#_get average over mw channels and invert. Is solution different
				# than averaging over shis_Tb?
				mw_shis_rd, mw_shis_wv = microwindow_average(shis_rads, shis_wnum,
											int(dirty_mw.resolution), error=False)
				mw_shis_Tb = planck_inv(mw_shis_rd, mw_shis_wv*100,
										domain='wavenumber')

				#_scatter plot it
				ax_rad.scatter(mw_shis_wv, mw_shis_Tb, marker='o', color='k')
				ax_co2.scatter(mw_shis_wv, mw_shis_Tb, marker='o', color='k')

			if plot_full:
				#_not implemented yet
				pass

				#_interp SHIS to LBL-DIS microwindows
				x = clear_mw.wavenumber[clear_mw.wavenumber.mask == False]
				f = interp1d(shis_wnum, shis_Tb)
				shis_mw = f(x)

				#_pull out only unmasked LBL-DIS microwindows
				dirt_mw = dirty_mw.Tb[clear_mw.wavenumber.mask == False]

				#_plot SHIS at microwindows
				ax_rad.scatter(x, f(x), marker='x', color='r')
	
				#_plot difference between microwindows and shis
				ax_mw = ax_rad.twinx()
				ax_mw.scatter(x, shis_mw-dirt_mw, marker='o',
							color='magenta', s=0.75) 
				ax_mw.set_ylim(-20, 20)
				ax_mw.set_ylabel('S-HIS minus LBL-RTM.aero (K)', size='xx-small')
	
		#_limit to data area
		xlim0 = ax_rad.xaxis.get_data_interval()
		ylims = (200, 300) if not zoom else (280, 300)
		xlims = (xlim0[0], 2000) if not zoom else (750, 1300)
		ax_rad.set_xlim(*xlims) #xlim0[0], 2000)
		ax_rad.set_ylim(*ylims) #200, 300)
		ax_co2.set_xlim(680, 800)
		ax_co2.set_ylim(200, 300)

		#_turnon grids
		ax_rad.grid(True)
		ax_co2.grid(True)
		ax_rad.set_xlabel('wavenumber ($cm^{-1}$)', size='xx-small')
		ax_rad.set_ylabel('brightness temperature (K)', size='xx-small')
		
		#_set title
		ax_rad.set_title(title, size='small')

		#_plot cpl data
		flight_idx = abs(self.SHIS_altitude[fov] - arange(900)*30).argmin()
		flight_idx = (arange(900)*30)[flight_idx]
		cpl_max = cpl_plt.shape[1]
		cpl_nx, cpl_ny = cpl_plt.shape
		cb = ax_cpl.pcolormesh(cpl_plt.T, vmin=0, vmax=1e-4, cmap=spectral)
		ax_cpl.set_title('CPL BACKSCAT 532_EXT', size='xx-small')
##		ax_cpl.plot([cpl_idx, cpl_idx], [-1, flight_idx], 'r-')
##		ax_cpl.plot([cpl_idx, cpl_idx], [-1, 1e3], 'r-')
		xlim = ax_cpl.xaxis.get_data_interval()
		ax_cpl.set_xlim(xlim)
	##	ax_cpl.set_ylim(150, cpl_ny)
	##	ax_cpl.set_ylim(0, cpl_ny)
		ax_cpl.get_yaxis().set_visible(False)
		ax_cpl.get_xaxis().set_visible(False)

		#_create string of text containing CPL aod values for layers
		tau_532		= self.CPL_tau_532[cpl_idx][:]
		tau_type	= self.CPL_tau_type[cpl_idx][::2]
		tau_str		= [' ']
		tau_type_str = [None, 'PBL', 'AEROSOL', 'CLOUD']
		tau_532 = masked_where(tau_532 <= 0, tau_532)
####	tau_type = masked_where(tau_type == 0, tau_type)
##		tau_532 = masked_where(tau_532 == -8.8, tau_532)
##		for ii in xrange(10 - tau_532.mask.sum()):
##		for ii in xrange(5 - tau_type.mask.sum()):
		for ii in xrange(10 - tau_532.mask.sum()):
			try:
				tau_str.append('{0:7.2f} {1:s}'.format(tau_532[ii],
								tau_type_str[tau_type[ii]]))
			except:
				pass
			##	dbg(('TAU ERROR', tau_532, tau_type, ii, fov))

		tau_str = '\n'.join(tau_str)

		#_label CPL aod values
		diff_nx = cpl_nx / 2.8
		xtextloc = cpl_idx-diff_nx if cpl_nx-cpl_idx < diff_nx else cpl_idx+20
		ax_cpl.text(xtextloc, cpl_max-10, tau_str, color='r',
					verticalalignment='top', size='xx-small')

		#_shrink labels
##		[shrink_ticks(xx) for xx in [ax_map, ax_cpl, ax_rad, ax_pro]]
		[shrink_ticks(xx) for xx in fig.axes]
					
		plt.tight_layout()
	
		dir_out = JOIN(dir_out, 'plots_{0:s}'.format(self.dtg0))
		mkdir_p(dir_out)
		if plot_model:	
			pname = 'hs3_{5:s}_f{4:04d}_r{0:03.1f}_'
			pname += 'aod{1:03.1f}_z{2:1.0f}-{3:1.0f}'
			pname = pname.format(ref, aod, z0, z1, fov, self.dtg0)
		else:
			pname = 'hs3_{0:s}_f{1:04d}'.format(self.dtg0, fov)
 	
		pname = '.'.join((pname, 'png'))
		pname = JOIN(dir_out, pname)
		dbg(pname)
		plt.savefig(pname)
		plt.close()
	
	
	def shis_average(self, nlen=5, nwid=5, collocate=True, **kwargs): 
		'''
		nlen		int		along track size (inc. nadir)
		nwid		int		along sweet size (inc. nadir)
		collocate	bool	average at collocation points

		WHICH CAME FIRST, THIS OR THE STAND ALONE METHOD? 2015.05??

		creates an average spectral radiance from surrounding
		measurements, from an area NLENxNWID
		'''
		from libtools import strictly_decreasing as mono
		from numpy import append, recarray, arange, zeros, vstack

		#_if already collocated, reject averaging
		if self.collocated:
			dbg('cannot average collocated data', 3)
			return
	
		fov = self.shis.angle.copy()

		#_calc how far from nadir in each direction to go
		nalong = int(nlen) / 2	#_intentionally rounding
		nwing = int(nwid) / 2
	
		#_this should be constant at 14, but meh
		s, e = 0, 2 
		while(mono(fov[:e])):
			e += 1
		else:
			nangle = e - 1

		#_total number of sweeps
		nsweep = fov.size / nangle 
	
		#_find nadir and sanity check
		fov = vstack(fov).reshape(nsweep, nangle)
		nadir = abs(fov).argmin(axis=1)
		if nadir.std() > 1e-5:
			raise RuntimeError, 'nadir should be at same index throughout'
		n = nadir[0]
	
		#_create new recarray to fill
		shis_avg = recarray((0,), dtype=self.shis.dtype)

		#_get list of indices
		idx = self.shis.idx if collocate \
							else arange(n, self.shis.size+1, nangle)
	
		#_currently doesn't check to see if 
		# any indices are within 2 of data boundary
		idx_final, idx_keep = [], []
		for i in idx:
			#_index range for each sweep, 5 sweeps
			idxn = [i + n*nangle for n in arange(-nalong,nalong+1)] 
	
			#_skip if outside data area
			if 	(arange(idxn[0]-2, idxn[0]+3) < 0).any() or \
				(arange(idxn[-1]-2, idxn[-1]+3) > self.shis.size-1).any():
				continue
	
			#_pull out sweeps + side lobs and put together
			rn = [	self.shis.radiances[idxn[n]-nwing:idxn[n]+nwing+1] 
					for n in arange(nwid)]
			r = append([], rn)

			#_put into recarray to return 
			shis_avg.resize(shis_avg.size+1)
			shis_avg[-1].radiances = self.shis.radiances[i].copy()
			shis_avg[-1].radiances[:] = r.mean(axis=0)
	
			#_store which indices were actually used
			idx_final.append(i)
			idx_keep.append(idx.tolist().index(i))	
	
		#_reduce rest of metadata to appropriate dimensionality
		idx = idx_final[:]
		shis_avg.epoch[:] = self.shis.epoch[idx]
		shis_avg.angle[:] = self.shis.angle[idx]
		shis_avg.latitude[:] = self.shis.latitude[idx]
		shis_avg.longitude[:] = self.shis.longitude[idx]
		shis_avg.altitude[:] = self.shis.altitude[idx]

		#_collocate CPL data using potentially degraded averaged indices
		if collocate:
			#_reduce cpl 
			fname = self.cpl.fname
			self.cpl = self.cpl[self.cpl.idx[idx_keep,0]]
			self.cpl.__setattr__('idx', None)
			self.cpl.__setattr__('fname', fname)
			self.collocated = True	
	
			#_reduce profile
			fname = self.profile.fname
			self.profile = self.profile[idx] 
			self.profile.__setattr__('fname', fname)
 
		#_update object	
		fname = self.shis.fname
		self.shis = shis_avg
		self.shis.__setattr__('idx', None)
		self.shis.__setattr__('fname', fname)

		self._strip_collocation()


	def require_profile(self):
		''' if called, reduces to only FOVs with retrieved temperature '''
		from numpy import append, array, vstack
		from numpy.ma import masked_where
 
		if not self.collocated:
			dbg('data must be collocated')
			return

		#_create single matrix to work with
		temp = vstack(self.profile.temperature_shis)
		temp = masked_where(temp == -9999, temp)

###		#_find where entire profile is missing
###		idx = temp.mask.all(axis=1) == False

		#_find where most of the profile is missing
		n0, n1 = temp.shape
		idx0 = temp.mask.sum(axis=1) < n1/2. 
		idx1 = self.profile.flag == 0
		idx = append(idx0[None,:], idx1[None,:], axis=0).all(axis=0)

		self._reduce(idx)

		#_temp
####		files = self._inherit()
###		f_shis = self.shis.fname
###		f_cpl = self.cpl.fname
###		f_prof = self.profile.fname

		#_reduce
	##	self.shis = self.shis[idx]
	##	self.cpl = self.cpl[idx]
	##	self.profile = self.profile[idx]

		#_set attrs
	#	self._inherit(files=files)
###		setattr(self.shis, 'fname', f_shis)
###		setattr(self.profile, 'fname', f_prof)
###		setattr(self.cpl, 'fname', f_cpl)
	#	setattr(self.shis, 'idx', None)
	#	setattr(self.cpl, 'idx', None)


	def read_shis_retrieval(self, fidx=None, fname_final='lbldis_input.final',
		dir_lblrtm_fmt=JOIN(DIR_LBL, '{3}_{0}.{1}_fov{2:06d}'),
		experiment='hs3', **kwargs):
		''' read in final retrievals for selected field of views '''
		from lblrtm_utils import Lbldis_input

		#_if non specified, try all
		fidx = range(self.size) if fidx is None else fidx

		#_loop over each
		for i in fidx:
			#_get output directory
			dir_lblrtm = dir_lblrtm_fmt.format(	self.dtg0, self.dtg1, i,
												experiment)

			#_get file name
			fname = JOIN(dir_lblrtm, fname_final)

			#_read input file
			input = Lbldis_input(fname)

			#_put optical depth into recarray
			try:
				tau = 0
				for layer in input.layers:
					tau += layer['tau'][1]	#_only true for files with 2
											# instances
				
##				self.SHIS_tau[i] = input.layers[0]['tau'][1] #_probably
				self.SHIS_tau[i] = tau 
			except:
				dbg(('no retrieval for fov', i))
				pass

	def cpl_layer(self, aerosol=True, pbl=True, cloud=True, clear=False,
		**kwargs):
		'''
		limit to cases where there are no layers in CPL
		aerosol, pbl, cloud	boolean,	leave in specified layer
		clear	boolean,				accept no layers

		'''
		if not self.collocated:
			dbg('data must be collocated')
			return

		from numpy import append, vstack 
		from numpy.ma import masked_where

		#_read in layer type (1=PBL, 2=AERO, 3=CLOUD)
		if clear:
####		idx = (vstack(self.cpl.tau_type) < 0).all(axis=1)
##			idx = (vstack(self.cpl.tau_type) == -8.8).all(axis=1)
			
			#_try just very low OD
			tmp = vstack(self.cpl.tau_532)
			tmp = masked_where(tmp < 0, tmp)
			tmp = tmp.sum(axis=1)
			idx = tmp == 0.0

		else:
			tmp = vstack(self.cpl.tau_type)
			idx = tile(True, (self.cpl.size, 1)) 

			if not pbl:
				tidx = (tmp != 1).any(axis=1)
				idx = append(idx, tidx[:,None], axis=1)

			if not aerosol:
				tidx = (tmp != 2).any(axis=1)
				idx = append(idx, tidx[:,None], axis=1)

			if not cloud:
				tidx = (tmp != 3).any(axis=1)
				idx = append(idx, tidx[:,None], axis=1)

			idx = idx.all(axis=1)
	
		self._reduce(idx)


	def _reduce(self, idx):
		files = self._inherit()
		self.shis = self.shis[idx]
		self.cpl = self.cpl[idx]
		self.profile = self.profile[idx]
		self._inherit(files=files)


	def _inherit(self, files=None):
		''' dumb dumb dumb '''
		if files is None:
			return self.shis.fname, self.cpl.fname, self.profile.fname
		else:
			setattr(self.shis, 'fname', files[0])
			setattr(self.profile, 'fname', files[1])
			setattr(self.cpl, 'fname', files[2])
			

	def to_ncdf(self, collocate=True, file_flight='test.nc', **kwargs):
		''' write collocated Flight_data to a netcdf file '''
		from netCDF4 import Dataset
		from numpy import array, vstack#, #str

		dbg(file_flight)
		with Dataset(file_flight, 'w') as cdf:
			#_create dimensions
			wave = self.shis.radiances[0].wavenumber
			nfov = self.shis.size
			nwvn = wave.size 
			nlay = self.cpl[0].tau_532.size
			ntyp = nlay / 2
			nrng = self.cpl[0].ext_532.size
			npre_shis = self.profile[0].temperature_shis.pressure.size
			npre_gdas = self.profile[0].temperature_gdas.pressure.size
			pressure_shis = self.profile[0].temperature_shis.pressure.copy()
			pressure_gdas = self.profile[0].temperature_gdas.pressure.copy()

			cdf.createDimension('wavenumber', nwvn)
			cdf.createDimension('fov', nfov)
			cdf.createDimension('layer', nlay)
			cdf.createDimension('layer_type', ntyp)
			cdf.createDimension('range', nrng) 
			cdf.createDimension('plevs_gdas', npre_gdas) 
			cdf.createDimension('plevs_shis', npre_shis) 

			#_write global attributes
			cdf.dtg = self.dtg
			cdf.file_cpl = self.cpl.fname	
			cdf.file_shis = self.shis.fname
			cdf.file_gdas = self.profile.fname

			#_make sure everything is in temporal order
			idx = self.shis.epoch.argsort()
			rad = vstack(self.shis.radiances[idx])
			idx = self.cpl.epoch.argsort()
			tau = vstack(self.cpl.tau_532[idx])			
			ext = vstack(self.cpl.ext_532[idx])			
			typ = vstack(self.cpl.tau_type[idx])			
			idx = self.profile.epoch.argsort()
			tem_shis = vstack(self.profile.temperature_shis[idx])
			ozo_shis = vstack(self.profile.ozone_shis[idx])
			wat_shis = vstack(self.profile.water_vapor_shis[idx])
			tem_gdas = vstack(self.profile.temperature_gdas[idx])
			ozo_gdas = vstack(self.profile.ozone_gdas[idx])
			wat_gdas = vstack(self.profile.water_vapor_gdas[idx])

			#_create variables
			#_write variables
			v = {'fill_value' : -9999.}
			c = cdf.createVariable('tau_532', 'f4', ('fov', 'layer'), **v)
			c[:] = tau
			c.units = 'unitless'	
			c = cdf.createVariable('tau_type', 'i4', ('fov', 'layer_type'), **v)
			c[:] = typ
			c.units = '1=PBL, 2=AEROSOL, 3=CLOUD'
			c = cdf.createVariable('ext_532', 'f8', ('fov', 'range'), **v)
			c[:] = ext
			c.units = '???'
			c = cdf.createVariable('radiances', 'f8', ('fov', 'wavenumber'),**v)
			c[:] = rad
			c.units = 'W/m2/sr/cm-1'	
			c = cdf.createVariable('wavenumber', 'f4', ('wavenumber'), **v)
			c[:] = wave
			c.units = 'cm-1'	
			c = cdf.createVariable('latitude', 'f8', ('fov',), **v)
			c[:] = self.shis.latitude
			c.units = 'degrees_north'
			c = cdf.createVariable('longitude', 'f8', ('fov',), **v)
			c[:] = self.shis.longitude
			c.units = 'degrees_east'
			c = cdf.createVariable('epoch_shis', 'f8', ('fov',), **v)
			c[:] = self.shis.epoch
			c.units = 'sec since Jan 1 1970'	
			c = cdf.createVariable('epoch_cpl', 'f8', ('fov',), **v)
			c[:] = self.cpl.epoch
			c.units = 'sec since Jan 1 1970'	
			c = cdf.createVariable('angle', 'f4', ('fov',), **v)
			c[:] = self.shis.angle
			c.units = 'degrees from nadir, left=neg'
			c = cdf.createVariable('altitude', 'f4', ('fov',), **v)
			c[:] = self.shis.altitude
			c.units = 'm'
			c.long_name = 'aircraft altitude AGL'
			c = cdf.createVariable('prof_alt_gdas', 'f4', ('plevs_gdas',), **v)
			c[:] = self.profile[0].altitude_gdas
			c.units = 'km'
			c.long_name = 'reallllly bad estimation of AGL of plevs'
			c = cdf.createVariable('prof_alt_shis', 'f4', ('plevs_shis',), **v)
			c[:] = self.profile[0].altitude_shis
			c.units = 'km'
			c.long_name = 'reallllly bad estimation of AGL of plevs'
			c = cdf.createVariable('pressure_gdas', 'f4', ('plevs_gdas',), **v)
			c[:] = pressure_gdas 
			c.units = 'hPa'
			c = cdf.createVariable('pressure_shis', 'f4', ('plevs_shis',), **v)
			c[:] = pressure_shis 
			c.units = 'hPa'
			c = cdf.createVariable('temperature_shis', 'f4', ('fov','plevs_shis',), **v)
			c[:] = tem_shis
			c.units = 'K'
			c = cdf.createVariable('temperature_gdas', 'f4', ('fov','plevs_gdas',),
				**v)
			c[:] = tem_gdas
			c.units = 'K'
			c = cdf.createVariable('relative_humidity_shis', 'f4', 
									('fov','plevs_shis',), **v)
			c[:] = wat_shis 
			c.units = 'percent'
			c = cdf.createVariable('relative_humidity_gdas', 'f4', 
									('fov','plevs_gdas',), **v)
			c[:] = wat_gdas 
			c.units = 'percent'
			c = cdf.createVariable('ozone_gdas', 'f4', ('fov','plevs_gdas',), **v)
			c[:] = ozo_gdas 
			c.units = 'ppmv or g/kg (check later)'
			c = cdf.createVariable('ozone_shis', 'f4', ('fov','plevs_shis',), **v)
			c[:] = ozo_shis 
			c.units = 'ppmv or g/kg (check later)'
			c = cdf.createVariable('sfc_temperature_shis', 'f4', ('fov',), **v)
			c[:] = self.profile.sfc_temperature_shis
			c.units = 'K'
			c = cdf.createVariable('sfc_temperature_gdas', 'f4', ('fov',), **v)
			c[:] = self.profile.sfc_temperature_gdas
			c.units = 'K'
			c = cdf.createVariable('jchar', str, ('fov',))
			for n in range(nfov):
				c[n] = self.profile.jchar[n]
			c.long_name = 'LBL-RTM species defs'


	def _strip_collocation(self):
		''' certain routines make collocation indices invalid '''
		[self.__getattribute__(s).__setattr__('idx', None) \
			for s in Flight_data.sensors]


from numpy import arange, ndarray
class Test(ndarray):
	def __new__(cls, shape, dtype=float, info=None, **kwargs):
		obj = ndarray.__new__(cls, shape, dtype, **kwargs)
		obj.info = info
		return obj 

	def __array_finalize__(self, obj):
		#_new object construct, skip
		if obj is None:
			return
	
		#_otherwise, copy over default
		self.info = getattr(obj, 'info', None)


class Spectrum(ndarray):
	def __new__(cls, radiances, wavenumber, frequency=None, Tb=None,
        transmissivity=None, flux=None, label='s', time=-9999, latitude=-9999,
        longitude=-9999, cloud=-9999, **kwargs):
		from numpy import recarray, diff, where
		from libgeo import wavenumber2frequency, planck_inv
		from numpy.ma import masked_outside

		rad_range = [-1e-6, 9.e7]	#_outside these limits, mask data
		rads = masked_outside(radiances, rad_range[0], rad_range[1])
		obj = ndarray.__new__(cls, rads.shape)
		obj[:] = rads[:]
		obj.wavenumber = wavenumber						#_cm-1
##		obj.Tb         = planck_inv(rads, wavenumber*100,
##									domain='wavenumber')    #_K
		obj.frequency  = wavenumber2frequency(wavenumber)
		obj.label      = label
		obj.flux       = flux
		obj.latitude   = latitude
		obj.longitude  = longitude
		obj.time       = time
		obj.cloud      = cloud
		obj.transmissivity = transmissivity

		return obj


	def __array_finalize__(self, obj):
		self.wavenumber = getattr(obj, 'wavenumber', None)
		self.Tb         = getattr(obj, 'Tb', None) 
		self.frequency  = getattr(obj, 'frequency', None) 
		self.label      = getattr(obj, 'label', None)
		self.flux       = getattr(obj, 'flux', None)
		self.latitude   = getattr(obj, 'latitude', None)
		self.longitude  = getattr(obj, 'longitude', None)
		self.time       = getattr(obj, 'time', None)
		self.cloud      = getattr(obj, 'cloud', None)
		self.transmissivity = getattr(obj, 'transmissivity', None)


	def plot(self, figure=None, ax=None, title='', show=False, Tb_diff=None,
		**kwargs):
		import matplotlib.pyplot as plt
		from numpy import where, diff, zeros
		'''
        plots current spectra to figure.axes.  Does nothing to
        persistently save it currently.
		'''

		#_if new, craete figure
		if figure is None and ax is None:
			self.figure = plt.figure()
			self.figure.add_subplot(311)    #_radiance
			self.figure.add_subplot(312)    #_brightness temp
			self.figure.add_subplot(313)    #_difference
			self.figure.axes[0].set_title(title)
			ls = 'b-'

			ax = self.figure.axes[0]

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
			self.figure.axes[0].plot(   self.wavenumber[start:end],
                                        self[start:end], ls, **arg)
            #_brightness temp
			self.figure.axes[1].plot(   self.wavenumber[start:end],
                                        self.Tb[start:end], ls, **arg)
			self.figure.axes[2].plot(   self.wavenumber[start:end],
                                        Tb_diff[start:end], ls, **arg)

			start = end+1
		else:
			self.figure.axes[0].plot(   self.wavenumber[start:],
                                        self[start:], ls, **arg)
			self.figure.axes[1].plot(   self.wavenumber[start:],
                                        self.Tb[start:], ls, **arg)
			self.figure.axes[2].plot(   self.wavenumber[start:],
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

		if show:
			plt.show()


	def write(self, fname=None, dir_out=JOIN(DIR_PROD, 'spectra'), **kwargs):
		'''write data to label.nc'''
		from netCDF4 import Dataset
		from libtools import mkdir_p

		fname = fname if fname != None \
            else '.'.join((self.label, 'nc'))
		pname = JOIN(dir_out, fname.split('/')[-1])

        #_append if file already present
		if not os.path.exists(pname):
			mkdir_p(dir_out)
			dbg(('writing', pname))
			cdf = Dataset(pname, 'w', format='NETCDF3_CLASSIC')

			#_initialize dimensions
			cdf.createDimension('time', None)
			cdf.createDimension('wavenumber', self.size)

			#_initialize variables
			cdf.createVariable('latitude',  'f8', ('time',), fill_value=-9999)
			cdf.createVariable('longitude', 'f8', ('time',), fill_value=-9999)
			cdf.createVariable('time',      'f8', ('time',), fill_value=-9999)
			cdf.createVariable( 'wavenumber', 'f8', ('wavenumber'),
                                fill_value=-9999)
			cdf.createVariable( 'cloud', 'i4', ('time',), fill_value=-9999)
			cdf.createVariable( 'radiances', 'f8', ('time', 'wavenumber',),
			                    fill_value=-9999)
			cdf.createVariable( 'Tb',       'f8', ('time', 'wavenumber'),
			                    fill_value=-9999)
			cdf.createVariable( 'wavenumber','f8', ('wavenumber'),
                                fill_value=-9999)
			cdf.createVariable( 'frequency', 'f8', ('wavenumber'),
                                fill_value=-9999)
			cdf.variables['wavenumber'][:] = self.wavenumber
			cdf.variables['frequency'][:] = self.frequency

			cdf.variables['wavenumber'].units = 'cm-1'
			cdf.variables['Tb'].units = 'K'
			cdf.variables['radiances'].units = 'W/m2/um/sr'

			cdf.setncattr('label', self.label)
		else:
			dbg(('appending', pname))
			cdf = Dataset(pname, 'a', format='NETCDF3_CLASSIC')

		#_write spectrum data
		idx = cdf.variables['time'].size    #_last index of UNLIMITED
		for var in ['radiances', 'Tb']:
		    cdf.variables[var][idx] = self.__getattribute__(var).squeeze()
		cdf.variables['latitude'][idx] = self.latitude
		cdf.variables['longitude'][idx] = self.longitude
		cdf.variables['time'][idx] = self.time
		cdf.variables['cloud'][idx] = self.cloud

		cdf.close()


def get_naaps_ctp(epoch, latitude, longitude, efold=1, diagnostic=False,  
	dir_naaps=os.path.join(os.environ['PRODUCTS'], 'NRL', 'NVA_CLIMO1misr'),
	specie='dust', conc_min=2.5e-8, **kwargs):
	'''
	return list of cloud top heights 
	
	dtg,	str{10},	date time group (yyyymmddhh)
	lat,	float,		latitude of location
	lon,	float,		longitude of location
	'''
	from libnva import find_nearest_dtg, read_conc
	from libtools import dtg2epoch, epoch2dtg
	from libgeo import sig2pres
	from numpy import array, arange, linspace, exp, greater, less
	from scipy.signal import argrelextrema as extrema
	from scipy.interpolate import interp1d 

	#_get nearest naaps file
	dtg_naaps = find_nearest_dtg(epoch)

	#_look for NAAPS file
	fname = '{0}_conc'.format(dtg_naaps)
	fname = os.path.join(dir_naaps, 'NAAPS', dtg_naaps[:6], fname) 
	
	#_check that file exists
	if not os.path.exists(fname):
		dbg('WARNING: NAAPS CTP not available {0}'.format(fname))
		return -9999., -9999.

	#_read in concentration file
	conc = read_conc(fname, **kwargs)

	#_pull out coordinates (this is bad, fix)
	varnames = conc.variable.tolist()
	idx_lat = varnames.index('lat')
	idx_lon = varnames.index('lon')
	idx_con = varnames.index('conc')
	idx_spe = varnames.index('species')
	idx_sga = varnames.index('sig_a')
	idx_sgb = varnames.index('sig_b')
	idx_sig = varnames.index('sigma')
	idx_sfc = varnames.index('sfc_pres')

	lon = conc.values[idx_lon]
	lat = conc.values[idx_lat]
	spe = conc.values[idx_spe].tolist()
	con = conc.values[idx_con][spe.index(specie)]	#_nz,ny,nx
	sga = conc.values[idx_sga]
	sgb = conc.values[idx_sgb] #_staggered levels (get midpoint? 'sigma')
	sig = conc.values[idx_sig]
	sfc = conc.values[idx_sfc]

	#_chec 
	idx_lat = abs(latitude - lat).argmin()
	idx_lon = abs(longitude - lon).argmin()
	con = con[:, idx_lat, idx_lon]

	#_get local maximum indices
	idx_max = extrema(con, greater)[0]
	idx_min = extrema(con, less)[0].tolist()

	#_if no max, return meh
	if len(idx_max) == 0:
		return -9999., -9999.

	#_set top pressure level
	nz = len(con)
	if nz == 25:
		ptop = 70. #_hPa
	elif nz == 30:
		ptop = 1.
	else:	
		raise ValueError, 'ERROR: Unknown sigma configuration'
 
	#_calculate cloud top pressure
	pressure = sig2pres(sig, psfc=sfc[idx_lat, idx_lon], ptop=ptop, **kwargs)

	#_find layer pressures
	layer_conc = con[idx_max] 
	layer_midp = pressure[idx_max]

	#_check direction of vertical coordinate	
	upward = 1 if all(x<y for x,y in zip(pressure, pressure[1:])) \
		else -1 

	#_loop over each layer, calc decay from max, find location ABOVE
	pres_top = []
	pres_bot = []
	conc_top = []
	conc_bot = []

	#_loop over local maximums
	for i, idx_gt in enumerate(idx_max):

		#_kludge around upward/downward issue
		try:
			idx_lt = idx_min[i] if upward else idx_min[i-upward]
		except IndexError:
			dbg(idx_min)
			dbg(idx_max)
			dbg('WARNING: skip this mess')
			continue

		idx_top = max([idx_gt, idx_lt]) + 1
		idx_bot = min([idx_gt, idx_lt])
		conc_sub = con[idx_bot:idx_top]
		pres_sub = pressure[idx_bot:idx_top]

		#_skip arbitrarily this layers
		if con[idx_bot] < conc_min:
			continue

		#_get decay concentration from max	
		conc_reduced = con[idx_gt] / exp(efold)

		#_look for nearest match
		f = interp1d(pres_sub[::-1], conc_sub[::-1])
		pres_hi = linspace(pres_sub[0], pres_sub[-1], 1000)
		conc_hi = f(pres_hi)
	
		idx_hi = abs(conc_hi - conc_reduced).argmin()
	
		#_midpoint
	#	pres_top.append(pressure[idx_gt])
	#	conc_decay.append(con[idx_gt])

		#_bounds
		pres_top.append(pres_hi[idx_hi])
		conc_top.append(conc_hi[idx_hi])
	
	#_find bottom of layer
	for i, idx_gt in enumerate(idx_max):

		idx_lt = ((idx_min-idx_gt) < 0).sum() - 1
		idx_lt = idx_min[idx_lt] if idx_lt >= 0 else 0

		#_kludge around upward/downward issue
	#	idx_lt = idx_min[i-1] # if upward else idx_min[i-upward]
		idx_top = max([idx_gt, idx_lt]) + 1
		idx_bot = min([idx_gt, idx_lt])
		conc_sub = con[idx_bot:idx_top]
		pres_sub = pressure[idx_bot:idx_top]

		#_skip arbitrarily this layers
		if con[idx_top] < conc_min:
			continue

		#_get decay concentration from max	
		conc_reduced = con[idx_gt] / exp(efold)

		#_look for nearest match
		f = interp1d(pres_sub[::-1], conc_sub[::-1])
		pres_hi = linspace(pres_sub[0], pres_sub[-1], 1000)
		conc_hi = f(pres_hi)
	
		idx_hi = abs(conc_hi - conc_reduced).argmin()
	
		#_bounds
		pres_bot.append(pres_hi[idx_hi])
		conc_bot.append(conc_hi[idx_hi])

	if len(pres_top) == 0 or len(pres_bot) == 0:
		return -9999., -9999.

	def plot_ctp(con, pressure, conc_top, pres_top, conc_bot, pres_bot, pname, 
		dir_plot=os.path.join(DIR_PROD, 'plots', 'CTP')):
		''' diagnostic, show naaps fields and what the results were '''
		from numpy import array
		import matplotlib.pyplot as plt
		from libtools import mkdir_p
		from libgeo import p2z

		pname = os.path.join(dir_plot, pname)
		if os.path.exists(pname):
			return

		fig, ax = plt.subplots(1)
		mkdir_p(dir_plot)
		ax.plot(con, pressure)
		ax.scatter(conc_top, pres_top)
		ax.scatter(conc_bot, pres_bot)
		ax.set_ylabel('pressure (hPa)')
		plt.gca().invert_yaxis()

		z_top, z_bot = p2z(pres_top[0]), p2z(pres_bot[0])
		plt.title('{0:7.1f} to {1:7.1f} m'.format(z_bot, z_top))

		plt.savefig(pname, dpi=40)
		plt.close()

	if 1:# diagnostic:
		pname = 'NAAPS_{0}T_{1:09.5f}LAT_{2:010.5f}LON.png'.format(dtg_naaps,
			latitude, longitude)
		plot_ctp(con, pressure, conc_top, pres_top, conc_bot, pres_bot, pname)

	#_get level
	return array(pres_top), array(pres_bot)


################################################################################
#_PLOT_#########################################################################
# most hinge on plot_fov. 
################################################################################


################################################################################
#_TOOLS_########################################################################
################################################################################


def datetime2epoch(date, time):
	''' convert goofiness to unix time '''
	from libtools import dtg2epoch as d2e
	if hasattr(time, '__iter__'):
		epoch = []
		for n in xrange(len(date)):
			d = str(int(date[n])).zfill(6)
			t = str(int(time[n])).zfill(6)
			epoch.append(d2e('20' + d + t, full=True))
	else:
		epoch = d2e('20' + str(int(date)).zfill(6) 
				+ str(int(time)).zfill(6), full=True)

	return epoch 


def dbg(msg, l=1, err=False):
	''' 
    if global debug is set to true, be more verbose 
    msg : str, Message to be printed
    l   : int, Debug level of message.  Set higher for lower level 
            messages.  As debug increases, noisiness should also.
	'''
	import inspect
	from libtools import to_string
	msg = to_string(msg)
	if hasattr(msg,'__iter__'): msg = ' '.join(msg)

	if DEBUG >= l:
		curf = inspect.currentframe()
		calf = inspect.getouterframes(curf,2)
		file, line, method = calf[1][1:4]
		file = '.'.join(file.split('/')[-1].split('.')[:-1])
		scream = '[%s.%s.%i] %s' % (file,method,line,msg)

		if not err:
			print scream
		else:
			raise RuntimeError, scream


if __name__ == '__main__':
	print '''
	FILE	hs3.py
	AUTHOR	Walter R Sessions
	USAGE	Don't run me.
	'''
	process_period()	#_SETUP FLIGHTS / COLLOCATE
##	simulate_flights()	#_LBLRTM/DIS
##	run_main()

