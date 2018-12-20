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


JOIN = os.path.join
DIR_PROD = os.environ['PRODUCTS']
DIR_HS3 = JOIN(os.environ['PRODUCTS'], 'hs3')
DIR_LBL = JOIN(os.environ['PRODUCTS'], 'LBL-RTM')
DIR_COL = JOIN(os.environ['PRODUCTS'], 'colloc')
DIR_LOG = os.environ['LOG']
DIR_OUT = JOIN(os.environ['PRODUCTS'], 'ssec_oe')
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

def get_caliop_times(fname):
	''' get start and end times out of caliop hdf file from SIPS '''
	from pyhdf.SD import SD, SDC
	from libtools import modis2epoch, dtg2epoch
	from numpy import floor

	#_open handle
	hdf = SD(fname, SDC.READ)
		
	#_get times and convert to dtg
	times = hdf.select('Profile_UTC_Time')[:]

	#_convert from float time to epoch
	tmin = times.min()
	tmax = times.max()

	#_pull out yymmdd
	dtg0 = '20{0:6.0f}00'.format(floor(tmin))
	dtg1 = '20{0:6.0f}00'.format(floor(tmax))

	#_get residual seconds as fraction of day in seconds
	sec0 = (tmin - floor(tmin)) * 86400.
	sec1 = (tmax - floor(tmax)) * 86400.

	#_add together
	e0 = dtg2epoch(dtg0) + sec0
	e1 = dtg2epoch(dtg1) + sec1
	
	return e0, e1


def get_airs_times(fname):
	''' yank out airs time range '''
	from pyhdf.SD import SD, SDC
	from libtools import modis2epoch, dtg2epoch

	#_open handle
	hdf = SD(fname, SDC.READ)

	#_yank out time variable
	times = hdf.select('Time')[:]
	times = modis2epoch(times)

	#_return start and end times
	return times.min(), times.max()
 

def create_collocation_files_airscal(dtg0=0, dtg1=1e9, dir_rad=DIR_HS3,
	pcfilt=True, dir_lid=JOIN(DIR_PROD, 'cpl'), dir_col=DIR_COL, **kwargs):
	''' 
	ON IRIS THIS CREATES A PICKLE TO BE USED BELOW TO GENERATE COLLOCATION FILES
	ON KEPLER

	From dir_rad, attempts to create collocationed segment files in
	dir_col (default $PRODUCTS/colloc) from the shiscpl script.
	The COLLOC files are then used to create segment files.

	THIS PART IS OUT OF DATE, FLIGHTS FILE DEPRECATED

	dtg0	str,	set if earliest date desired
	dtg1	str,	set if end date desired, default to last in dir_rad
	save_old bool,	true to skip already generated files	

	dir_col	str,	NOT CURRENTLY USED!  Check shislid script.
	'''
	import re, glob
	from libtools import dtg2iso, newdtg2
	from libtools import dtg2epoch as d2e
	from libtools import epoch2dtg as e2d 
	from numpy import array

	#_get first dtg for flights
	#_read in all available hours (what is last time?)
	re_dtg = re.compile('SHIS.CPL.COLLOC.(\d{12}).(\d{12}).hdf')
	re_rad = re.compile('AIRS.(\d{4}).(\d{2}).(\d{2}).(\d{3}).L1B.AIRS_Rad.v[\d.]+.G\d{11}.hdf')
	re_lid = re.compile('CAL_LID_L1-ValStage1-V3-30.(\d{4})-(\d{2})-(\d{2})T(\d{2})[\w\d.-]*hdf')
	
	#_collect start and end times for available SHIS files
	estart = []
	files_rad = glob.glob(dir_rad + '/AIRS*hdf')
	files_lid = array(glob.glob(dir_lid + '/CAL_LID_L1-ValStage1*hdf'))

	#_create list of lidar start and end times
	lid_bounds = [] # {key == filename of lidar : [list of start/end times]}
	for file_lid in files_lid:
		res = re_lid.search(file_lid)

		#_no match? look at next file
		if not res:
			continue

		#_pull out start/end times
		start, end = get_caliop_times(file_lid)

		lid_bounds.append({'bounds' : [start, end], 'file' : file_lid})

	#_MONDAY START HERE WALTER
	#_loop over radiance files to find overlapping periods
	collocation = {}
	for file_rad in files_rad:

		#_list that will contain any lidar files that overlap with this rad file
		lid_list = []
		
		#_pull out rad start and end times
		res	= re_rad.search(file_rad)

		#_file didn't match, move on
		if not res: 
			continue

		#_convert to epoch time
		str, end = get_airs_times(file_rad)
		from libtools import epoch2dtg as e2d

		#_loop over lid start and end times
		for values in lid_bounds:

			#_is ANY part of the file overlapping?
			lid_str, lid_end = values['bounds']
			if	(lid_str >= str and lid_str <= end) or	\
				(lid_end >= str and lid_end <= end) or	\
				(lid_str <= str and lid_end >= end):		#_lid encompasses
				lid_list.append(values['file']) 

		#_add rad file to collocation dictionary	
		if len(lid_list):
			collocation.update({file_rad : lid_list})

	#_put into pickle
	from pickle import dump
	pname = '{0}.pk'.format('airscal')
	dbg('dumping {0}'.format(pname))
	dump(collocation, open(pname, 'wb'))	


def create_collocation_files_airscal_kepler(dtg0=0, dtg1=1e9, dir_rad=DIR_HS3,
	pcfilt=True, dir_lid=JOIN(DIR_PROD, 'cpl'), dir_col=DIR_COL,
	pickle_file='airscal.pk', **kwargs):
	'''
	use pickle generated on IRIS to go through and create collocation files
	on kepler
	'''
	from pickle import load
	import re

	#_load in dictionary connecting rad files to lidar files
	collocation = load(open(pickle_file, 'rb'))

	#_loop over rad files, collocate and collect names of collocation files
	files_col = []
	err_files = []
	re_dtg = re.compile('colloc.airs_\d{8}T\d{4}.caliop_\d{8}T\d{4}.hdf')
	for file_rad, files_lid in collocation.iteritems():
		for file_lid in files_lid:
			#_swap out path
			file_rad = os.path.join(dir_rad, file_rad.split('/')[-1])
			file_lid = os.path.join(dir_lid, file_lid.split('/')[-1])

			try:
				#_run collocation code, get filename
				scp			= '/home/gregq/collopak/bin/airscal'
				cmd			= ' '.join((scp, file_rad, file_lid))
				dbg(cmd)
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
				err_files.append([file_rad, file_lid])
				continue
		
	#_write files that failed to error file
	from libtools import mkdir_p
	mkdir_p(dir_col)
	err_log = os.path.join(dir_col, 'failed_files')
	with open(err_log, 'w') as f:
		[f.write('{0} {1}\n'.format(s, c)) for s, c in err_files]


def create_segment_files_qsub(dir_col=DIR_COL, **kwargs):
	''' submit job to queue for generation '''
	from glob import glob
	from pickle import dump

    #_qsubmission prefix        
	env = ','.join(['='.join((var, os.environ[var])) \
		for var in ['PYTHONPATH', 'PRODUCTS', 'PATH', 'LOG']])
	qsub = ' '.join((   'qsub -v', env, '-o', DIR_LOG, '-e', DIR_LOG,
		'-cwd -S /opt/ShellB3/bin/python'))
	scpt = os.path.expanduser('~/lib/run_create_collocation.py')

	#_get list of collocation files
	files_col = glob(dir_col + '/SHIS.CPL.COLLOC.*hdf')
	for fname in files_col:
		#_create pickle file
		fn = fname.split('/')[-1]
		fkw = JOIN(dir_col, 'kwargs.{0:d}.{1}.pk'.format(os.getpid(), fn))
		dump(kwargs, open(fkw, 'wb'))

		#_put it all together and submit
		cmd = ' '.join((qsub, scpt, fname, fkw))
		dbg(cmd) 
		os.system(cmd)


def create_segment_files_sbatch(pname='/scratch/wsessions/updated_airscal.pk',
	segment_dt=30, dir_air=os.path.join(os.environ['PRODUCTS'],'sips/airs'),
	dtg0='20130824000000', dtg1='20130825000000',
	dir_cal=os.path.join(os.environ['PRODUCTS'], 'sips/caliop'),
	dir_col='/scratch/wsessions/colloc',
	dir_out='/scratch/wsessions/ssec_oe', **kwargs):
	'''
	Create segment files from updated pickle from kepler
	
	segment_dt	int,	length of segments in minutes
	'''
	
	from pickle import load	
	from libtools import epoch2dtg as e2d
	from libtools import dtg2epoch as d2e
	from pyhdf.SD import SD, SDC
	from libtools import mkdir_p
	from numpy import diff, vstack
	import re
	
	mkdir_p(dir_out)

	#_get times from collocation filenames
	def get_col_time(fname):
		from re import compile
		from libtools import dtg2epoch as d2e
		from libtools import epoch2dtg as e2d
		
		reg = compile('colloc.airs_(\d{8})T(\d{4}).caliop_(\d{8})T(\d{4}).hdf')
		res = reg.search(fname)		
		air = res.group(1) + res.group(2)
		cal = res.group(3) + res.group(4)

		return d2e(air, full=True), d2e(cal, full=True)

	#_open pickle (keys = collocation files, values = associated airs/caliop)
	col_dict = load(open(pname, 'rb'))

	#_put into float seconds
	segment_dt *= 60.

	#_pull out collocation files, organize temporally
	files_col = col_dict.keys()
	files_col.sort()

	#_loop over time segments
	e0 = d2e(dtg0, full=True)
	e1 = d2e(dtg1, full=True)
	epoch = e0
	while epoch < e1:
		''' A SEGMENT FILE SHOULD BE PRODUCED FOR EACH LOOP HERE '''

		#_set time range
		seg_str = epoch
		seg_end = epoch + segment_dt

		q = 0

		#_loop over and append all the crap together for the period
		for file_col in files_col:
			fname = os.path.join(dir_col, file_col)

			#_open file, get indices for satellite data	
			hdf = SD(fname, SDC.READ)
			idx_airX = hdf.select('Master_Index_1')[:].copy() - 1 	#_Xtrack (i think)
			idx_airL = hdf.select('Master_Index_2')[:].copy() - 1	#_Ltrack (i think)
			idx_cal = hdf.select('Follower_Index_1')[:].copy()	- 1			

			#_loop over the air/caliop files
			file_air = os.path.join(dir_air, 
				col_dict[file_col]['airs'].split('/')[-1])
			if re.search('08.23.240', file_air): continue
			air_str, air_end =  get_airs_times(file_air)

			files_cal = col_dict[file_col]['caliop']
			if len(files_cal) > 1: raise RuntimeError, 'nope'
			file_cal = os.path.join(dir_cal, files_cal[0].split('/')[-1])
			cal_str, cal_end = get_caliop_times(file_cal)

			#_see if these are within current segment range
			if ( air_end >= seg_str and air_end < seg_end) or \
				(air_str >= seg_str and air_str < seg_end) or \
				(air_str >= seg_str and air_end < seg_end) or \
				(air_str < seg_str and air_end > seg_end): 
				pass

			##	#_pull out appropriate times
			##	dbg(('IN SEG', e2d(seg_str, full=True),
			##		e2d(air_str, full=True), e2d(air_end, full=True),
			##		e2d(seg_end, full=True)))

					
			else:
				#_not within this time segment
				continue

			#_AIRS_####
			#_read in airs, limit to proper time range 
			tmp_airs, wave = read_airs(file_air, idx_airX, idx_airL, **kwargs)
			idx0 = tmp_airs.AIRS_epoch >= seg_str
			idx1 = tmp_airs.AIRS_epoch <  seg_end
			idx  = idx0*idx1
			if idx.sum() == 0:
				print 'none'
				continue
			tmp_airs = tmp_airs[idx0 * idx1].copy()
			
			if q > 0:
				#_append data to segment
				airs.resize(airs.size + tmp_airs.size)
				airs[-tmp_airs.size:] = tmp_airs[:].copy()
			else:
				#_init airs array
				airs = tmp_airs.copy()
				setattr(airs, 'AIRS_wavenumber', wave)

			#_CALIOP_####
			#_read in caliop, limit to proper time range 
			tmp_cal = read_caliop(file_cal, idx_cal, **kwargs)

			#_caliop is in a crazy shape, so flatten
			tim_cal = vstack(tmp_cal.CALIOP_epoch)[:,0]
			idx0 = tim_cal >= seg_str
			idx1 = tim_cal <  seg_end
			idx  = idx0 * idx1
			print 'BOUND ', tim_cal[idx==False]
			print 'BOUND0', e2d(tim_cal[idx == False], full=True)
			print 'BOUND1', e2d([seg_str, seg_end], full=True)
			tmp_cal = tmp_cal[idx]
			if q > 0:
				#_append data to segment
				caliop.resize(caliop.size + tmp_cal.size)
				caliop[-tmp_cal.size:] = tmp_cal[:].copy()
			else:
				#_init calip array
				caliop = tmp_cal.copy()
	
			#_GDAS_####
			#_read in GDAS profile
			tmp_gdas = read_profile_gdas(file_air, idx_airX, idx_airL, **kwargs)
			idx0 = tmp_gdas.GDAS_epoch >= seg_str
			idx1 = tmp_gdas.GDAS_epoch <  seg_end
			idx  = idx0 * idx1
			tmp_gdas = tmp_gdas[idx]
			if q > 0:
				#_append data to segment
				gdas.resize(gdas.size + tmp_gdas.size)
				gdas[-tmp_gdas.size:] = tmp_gdas[:].copy()
			else:
				#_init gdas array
				gdas = tmp_gdas.copy()

			#_append to arrays
			print airs.size, tmp_airs.size
			print caliop.size, tmp_cal.size
			print gdas.size, tmp_gdas.size
		
			#_testing location difference using great circ approx
			if 0:
				from libgeo import Point, great_circle
				from numpy import array

				#_convert locations to points
				a_coords = zip(airs.AIRS_latitude, airs.AIRS_longitude)
				c_coords = zip(caliop.CALIOP_latitude,caliop.CALIOP_longitude)
				g_coords = zip(gdas.GDAS_latitude, gdas.GDAS_longitude)
	
				a_pts = [Point(latitude=t, longitude=n) for t,n in a_coords] 
				c_pts = [Point(latitude=t, longitude=n) for t,n in c_coords] 
				g_pts = [Point(latitude=t, longitude=n) for t,n in g_coords]
				
				d_ac = [great_circle(a,b)/1e3 for a, b in zip(a_pts, c_pts)]
				d_ag = [great_circle(a,b)/1e3 for a, b in zip(a_pts, c_pts)]
				d_ac = array(d_ac)
				d_ag = array(d_ag)

				print 'MAX DISTANCE A/C:', max(d_ac)
				print 'MAX DISTANCE A/G:', max(d_ag)
				idx_ac = d_ac.argmax()	
				idx_ag = d_ag.argmax()	
				dbg(('MAX_DIST A/C/LON/LAT',
					caliop.CALIOP_longitude[idx_ac],
					airs.AIRS_longitude[idx_ac],
					caliop.CALIOP_latitude[idx_ac],
					airs.AIRS_latitude[idx_ac]))
				dbg(('MAX_DIST A/G/LON/LAT',
					gdas.GDAS_longitude[idx_ag],
					airs.AIRS_longitude[idx_ag],
					gdas.GDAS_latitude[idx_ag],
					airs.AIRS_latitude[idx_ag]))

			q += 1



		#_write segment file
		args = (e2d(seg_str, full=True)[2:], e2d(seg_end, full=True)[2:])
		out_name = 'AIRS.CALIOP.GDAS.COLLOC.{0}.{1}.nc'.format(*args)
		out_name = os.path.join(dir_out, out_name) 
		dbg(out_name)

		#_dump to netcdf
		write_collocated_airscal(airs, caliop, gdas, out_name, **kwargs)

		#_update time range	
		epoch += segment_dt


def precheck(pre_check='precheck.{0}.{1}.{2}.txt',
	dir_colloc=os.path.join(os.environ['WORK'], 'colloc'),
	latlon=[-90,90,-180,180], e0=0, e1=9e9, **kwargs):
	from libtools import epoch2dtg as e2d

	#_look if there is a prechecked list of files
	pre_check = os.path.join(dir_colloc, pre_check)
	latlon_str = '.'.join(['{0:.2f}'.format(ll) for ll in latlon])
	args = (e2d(e0, full=True), e2d(e1, full=True), latlon_str)
	pre_check = pre_check.format(*args)

	print pre_check, 'test'
	if os.path.exists(pre_check):
		files_colloc = [f.strip() for f in open(pre_check, 'r').readlines()]

		print 'MAKE SURE THIS WASNT A TEST FILE AND MAYBE ADD LAT/LON'
		pre_create = False

	else:
        #_get list of collocation files
		files_colloc = glob(dir_colloc + '/colloc*')

        #_create file
		pre_create = open(pre_check, 'w')

	#_return pre_create and files
	return files_colloc, pre_create


def get_collocation_indices(fname, **kwargs):
	''' depending on type of collocation, pass back indices '''
	import re
	from pyhdf.SD import SD, SDC

	if re.search('airs', fname) and re.search('caliop', fname):
		hdf = SD(fname, SDC.READ)
		i0 = hdf.select('Master_Index_1')[:] - 1 
		i1 = hdf.select('Master_Index_2')[:] - 1 
		i2 = hdf.select('Follower_Index_1')[:] - 1 
		#_missing caliop will be -1, up to 7# shots per AIRS footprint
		return i0, i1, i2


def get_collocation_files(fname, **kwargs):
	''' yank out names of collocated files '''
	import re
	from pyhdf.SD import SD, SDC

	fmt = re.compile('fname_(.*)')
	files = {}
	hdf = SD(fname, SDC.READ)
	for attr, value in hdf.attributes().iteritems():
		res = fmt.search(attr)
		if res:
			files.update({ res.group(1) : value })	

	return files


def read_caliop(fname, fidx, **kwargs):
	'''
	read in caliop data from SIPS

	Currently only yanking out time/location/backscatter

	To Add:
		Quality flags 
		Feature Mask
		Optical depth from L2 files
	'''
	from numpy import arange, array, recarray, ndarray, zeros
	from pyhdf.SD import SD, SDC
	from libsat import caliop2epoch as c2e
	from numpy.ma import masked_where
		
	xidx = 0 #_WHAT IS THE CROSS TRACK INDEX??

	#_open filehandle
	hdf = SD(fname, SDC.READ)
	n_airs, n_shot = fidx.shape
	n_fov, n_z = hdf.select('Total_Attenuated_Backscatter_532')[:].shape
	
	dtype = [	('CALIOP_total_attenuated_backscatter_532', ndarray),
				('CALIOP_epoch', ndarray),
				('CALIOP_longitude', ndarray), ('CALIOP_latitude', ndarray) ]
	caliop = recarray((n_airs), dtype)

##	fidx = fidx[:,0]
##	fidx = fidx.flatten()

	#_add in attributes
	tmp_back = hdf.select('Total_Attenuated_Backscatter_532')[:]
	tmp_time = hdf.select('Profile_UTC_Time')[:]
	tmp_lat = hdf.select('Latitude')[:]
	tmp_lon = hdf.select('Longitude')[:]
	
	#_build 3d array of tab532. (n_air, n_shot, n_z)
	for L in arange(n_airs):
		bak = zeros((n_shot, n_z)) - 9999.
		tim = zeros((n_shot)) - 9999.
		lat = zeros((n_shot)) - 9999.
		lon = zeros((n_shot)) - 9999.

		for S in arange(n_shot):
	##	for S in arange(1):
			I = S + L*n_shot
			I = int(fidx[L, S])
			bak[S, :] = tmp_back[I,:]
			tim[S] = c2e(tmp_time[I,0])
			lat[S] = tmp_lat[I,0]
			lon[S] = tmp_lon[I,0] 
		##	bak[S, :] = hdf.select('Total_Attenuated_Backscatter_532')[I,:]
		##	tim[S] = c2e(hdf.select('Profile_UTC_Time')[I,0])
		##	lat[S] = hdf.select('Latitude')[I,0]
		##	lon[S] = hdf.select('Longitude')[I,0]
		caliop.CALIOP_epoch[L]		= tim[:] 
		caliop.CALIOP_latitude[L]	= lat[:]
		caliop.CALIOP_longitude[L]	= lon[:]

		caliop.CALIOP_total_attenuated_backscatter_532[L] = bak[:]

##	caliop.CALIOP_total_attenuated_backscatter_532[:] = \
##		hdf.select('Total_Attenuated_Backscatter_532')[:,0][fidx]	
	return caliop


def read_oe(epochs, dir_lblrtm_fmt='.', experiment='hs3', out_label='.',
	lbldis_input='lbldis_input.{0}.final', **kwargs):
	''' read in out put from whatever '''
	from lblrtm_utils import Lbldis_input
	from libtools import epoch2dtg as e2d
	from numpy import recarray


##	dtype = [('tau', 'f4'), ('ref_wvnum', 'f4'), ('z_top', 'f4'), ('z', 'f4'),
##		('ref', 'f4'), ('epoch', 'f8'), ('habit', 'a12'), ('dbnum', 'i4')]
##	data = recarray((0,), dtype)

	#_init list 
	data = [None] * epochs.size

	#_loop over times, collect optical depth data
	for i, e in enumerate(epochs):

		#_convert time to dtg str for directory naming
		dtg = e2d(e, full=True)
		dir_lblrtm = dir_lblrtm_fmt.format(experiment, dtg)

		#_use outlabel to name the correct output file		
		fname = lbldis_input.format(out_label)
		fname = os.path.join(dir_lblrtm, fname)

		#_read in file, put into list
		try:
			input = Lbldis_input(fname)
			input.merge_layers()

			## for l in input.layers:
			##	habit = input.ssp_db[l['dbnum']].split('.')[-2]
			##	data.resize(data.size + 1)
			##	arg = (l['tau'][0], l['ref_wvnum'], l['z_top'], l['z'], l['ref'],
			##		e, habit, l['dbnum'])
					
			data[i] = input
		except IOError:
			print 'Missing', fname 

	#_return list
	return data
	

def read_airs(fname, idx_x=None, idx_l=None, **kwargs):
	'''
	read in airs data from SIPS

	Currently only yanking out time/location/radiances.

	idx_x	array,	cross track indices
	idx_l	array,	along track indices

	To Add:
	Quality flags 

	'''
	from numpy import arange, array, recarray, ndarray
	from pyhdf.SD import SD, SDC
	from pyhdf.HDF import *
	from pyhdf.VS import *

	from libtools import modis2epoch as m2e
	from libsat import AIRS_ancillary
	from numpy.ma import masked_where as mw 
	
	if not hasattr(idx_x, '__iter__'):
		idx_x, idx_l = array([idx_x]), array([idx_l])

	#_open filehandle
	hdf = SD(fname, SDC.READ)
	vdf = HDF(fname, HC.READ)
	vs = vdf.vstart()
	
	dtype = [ ('AIRS_radiances', ndarray), ('AIRS_epoch', 'f8'),
				('AIRS_longitude', 'f4'), ('AIRS_latitude', 'f4'),
				('AIRS_nen', ndarray) ]
	airs = recarray((len(idx_x),), dtype)

	#_load up bad channel mask
	mask = AIRS_ancillary().l2_ignore

	#_get location of error data
	edx = vs.find('NeN') #_WRS
	err = array(vs.attach(edx)[:]).flatten()
	err = mw(mask, err)
	
	#_pull out wavenumber from ancillary data and from hdf, compare
	wavenumber = AIRS_ancillary().srf_centroid

	for i in arange(idx_x.size):
		X, L = int(idx_x[i]), int(idx_l[i])
		airs.AIRS_epoch[i]		= m2e(hdf.select('Time')[L, X])
		airs.AIRS_latitude[i]	= hdf.select('Latitude')[L, X]
		airs.AIRS_longitude[i]	= hdf.select('Longitude')[L, X]
		airs.AIRS_radiances[i]  = mw(mask, hdf.select('radiances')[L, X, :]/1e3)
		airs.AIRS_nen[i]		= err[:] 

##	from numpy import vstack
##	radz = vstack(airs.AIRS_radiances)
##	print radz.shape
##	import matplotlib.pyplot as plt
##	plt.plot(airs.AIRS_radiances[0])
##	plt.savefig('test.png')

	return airs, wavenumber


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
	with Dataset(fname, 'r') as hdf:
		sidx = hdf.variables['SHIS_Index'][:] - 1#_the count seems to start at 1
		cidx = hdf.variables['CPL_Index'][:] - 1
		w = hdf.variables['Weights'][:]

		#_reformat as strictly strings
		f_shis	= '{0}'.format(hdf.SHIS_File)
		f_cpl	= '{0}'.format(hdf.CPL_File)	

	#_swap xdr file with nc file
	f_cpl = cpl_xdr2nc(f_cpl)

	#_read SHIS and CPL indices
	shis	= read_shis(f_shis, sidx, **kwargs)
	cpl		= read_cpl(f_cpl, cidx[:,0], **kwargs)
	prof	= read_profile_gdas(f_shis, sidx, **kwargs)

	#_check for jumps. Collocation code only looks at location,
	# does not seem to account for backtracking
	idx_jump = abs(shis.epoch - cpl.epoch) < 150 
##	idx_jump = diff(cpl.epoch) > 0
##	idx_jump = append(True, idx_jump)
	s_dict = shis.__dict__.copy()
	c_dict = cpl.__dict__.copy()
	g_dict = prof.__dict__.copy()
	shis = shis[idx_jump]
	cpl = cpl[idx_jump,:]
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


def write_collocated_airscal(airs, caliop, gdas, fname, dir_out=DIR_OUT,
	notes=None, **kwargs):
	''' write collocated flight segment file from read_collocated() '''
	from netCDF4 import Dataset
	from numpy import array, vstack
	import re
	from libtools import mkdir_p

	mkdir_p(dir_out)

	#_generate new flight segment file name
	fname = JOIN(dir_out, fname)	
	dbg(fname)

	#_pull out start and end times or default to arbitrary date after 1970
	res	= re.search('(\d{12}).(\d{12}).nc$', fname)
	dtg0, dtg1 = '20' + res.group(1), '20' + res.group(2)

	with Dataset(fname, 'w') as cdf:
		#_create dimensions
		nwvn	= airs.AIRS_wavenumber.size
		nfov	= airs.size
		narb, nz_caliop = \
			caliop[0].CALIOP_total_attenuated_backscatter_532.shape
	##	nz_airs = airs.AIRS_pressure.size
		nz_gdas = gdas.GDAS_pressure.size 

		cdf.createDimension('wavenumber', nwvn)
		cdf.createDimension('fov', nfov)
		cdf.createDimension('nz_caliop', nz_caliop) 
		cdf.createDimension('nz_gdas', nz_gdas) 
	##	cdf.createDimension('nz_airs', nz_airs)
		cdf.createDimension('dummy', narb) #_max shots from C for AIRS fov 

		#_write global attributes
	#	cdf.CALIOP_file	= caliop.CALIOP_fname	
	#	cdf.AIRS_file	= airs.AIRS_fname 
	#	cdf.GDAS_file	= gdas.GDAS_fname
		cdf.dtg0		= dtg0
		cdf.dtg1		= dtg1
 
		#_make sure everything is in temporal order
		idx	= airs.AIRS_epoch.argsort()
		rad_airs = vstack(airs.AIRS_radiances[idx])
	#	tem_airs = vstack(airs.AIRS_temperature[idx])
	#	ozo_airs = vstack(airs.AIRS_ozone_mixing_ratio[idx])
	#	rel_airs = vstack(airs.AIRS_relative_humidity[idx])
		tem_gdas = vstack(gdas.GDAS_temperature[idx])
		ozo_gdas = vstack(gdas.GDAS_ozone_mixing_ratio[idx])
		rel_gdas = vstack(gdas.GDAS_relative_humidity[idx])
		gpt_gdas = vstack(gdas.GDAS_geopotential_height[idx])
		pre_gdas = gdas.GDAS_pressure[0] 

		''' FIX HERE BECAUSE STACKING WON'T BE HAPPY '''	
		idx = vstack(caliop.CALIOP_epoch)[:,0].argsort()
		bak = vstack(caliop.CALIOP_total_attenuated_backscatter_532[idx])
		bak = bak.reshape(nfov, narb, nz_caliop)

		#_create variables
		#_write variables
		v = {'fill_value' : -9999.}

		c		= cdf.createVariable('CALIOP_total_attenuated_backscatter_532',
				'f4', ('fov', 'narb', 'nz_caliop'), **v)
		c[:]	= bak 
		c.units = 'km^-1 sr^-1'	

		c		= cdf.createVariable('CALIOP_latitude', 'f8', ('fov',), **v)
		c[:]	= caliop.CALIOP_latitude
		c.units = 'degrees_north'

		c		= cdf.createVariable('CALIOP_longitude', 'f8', ('fov',), **v)
		c[:]	= caliop.CALIOP_longitude
		c.units = 'degrees_east'

		c		= cdf.createVariable('CALIOP_epoch', 'f8', ('fov',), **v)
		c[:]	= caliop.CALIOP_epoch
		c.units = 'sec since Jan 1 1970'	

		c		= cdf.createVariable('AIRS_radiances','f8',
									('fov','wavenumber'),**v)
		c[:]	= rad_airs
		c.units = 'W/m2/sr/cm-1'
	
		c		= cdf.createVariable('AIRS_wavenumber','f4',('wavenumber'), **v)
		c[:]	= airs.AIRS_wavenumber
		c.units = 'cm-1'	

	##	c		= cdf.createVariable('AIRS_HBB_NESR','f4',('wavenumber'), **v)
	##	c[:]	= airs.hbb_nesr 

		c		= cdf.createVariable('AIRS_latitude', 'f8', ('fov',), **v)
		c[:]	= airs.AIRS_latitude
		c.units = 'degrees_north'

		c		= cdf.createVariable('AIRS_longitude', 'f8', ('fov',), **v)
		c[:]	= airs.AIRS_longitude
		c.units = 'degrees_east'
	
		c		= cdf.createVariable('AIRS_epoch', 'f8', ('fov',), **v)
		c[:]	= airs.AIRS_epoch
		c.units = 'sec since Jan 1 1970'	

		c		= cdf.createVariable('AIRS_altitude', 'f4', ('fov',), **v)
		c[:]	= [7077.75e3] * airs.size
		c.units = 'm'
		c.long_name = 'sensor altitude'

	##	c		= cdf.createVariable('AIRS_pressure', 'f4', ('nz_airs',), **v)
	##	c[:]	= airs.pressure 
	##	c.units = 'hPa'

		c		= cdf.createVariable('GDAS_latitude', 'f8', ('fov',), **v)
		c[:]	= gdas.GDAS_latitude
		c.units = 'degrees_north'

		c		= cdf.createVariable('GDAS_longitude', 'f8', ('fov',), **v)
		c[:]	= gdas.GDAS_longitude
		c.units = 'degrees_east'

		c		= cdf.createVariable('GDAS_pressure', 'f4', ('nz_gdas',), **v)
		c[:]	= pre_gdas #gdas.GDAS_pressure 
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
		c[:]	= gdas.GDAS_sfc_temperature 
		c.units = 'K'

		#_quick check
		print 'update this to use spherical code'
		if	abs(airs.AIRS_latitude - gdas.GDAS_latitude).max() > 0.5	\
			or abs(airs.AIRS_longitude - gdas.GDAS_longitude).max() > 0.5:
			raise ValueError, 'Collocation failed!'

		#_add arbitrary comments
		if notes is not None:
			cdf.notes = notes


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
		if	abs(shis.latitude - prof.latitude).max() > 0.5	\
			or abs(shis.longitude - prof.longitude).max() > 0.5:
			raise ValueError, 'Collocation failed!'

		#_add arbitrary comments
		if notes is not None:
			cdf.notes = notes


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



def read_cpl(fname, idx=None, dir_cpl=JOIN(DIR_PROD, 'cpl', 'nc'), **kwargs):
	'''
	fname	str,		filename of cpl data in netcdf format in dir_cpl
	dir_cpl	str,		path to cpl data files in netcdf format
	idx		ndarray,	collocated indices with SHIS for flight segment
	'''
	from netCDF4 import Dataset
	from numpy import ndarray, recarray, arange
	from libtools import julian2epoch as j2e 
	from numpy.ma import masked_where

	#_find filename for this dtg
	fname = JOIN(dir_cpl, fname)
	dbg(fname)
	h_cpl = Dataset(fname, 'r', format='NETCDF3_CLASSIC')

	#_read in cpl data
	dtype = [	('ext_532', ndarray),
				('tau_532', ndarray),
				('tau_type', ndarray),
				('layertop', ndarray),
				('layerbot', ndarray),
				('epoch', 'f8'),
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
	for t in xrange(idx.size):
		i = idx[t]
		cpl.ext_532[t]	= h_cpl.variables['Extprofile_sig_ch2'][i][::-1]
		cpl.tau_532[t]	= h_cpl.variables['Layertau_ch2'][i][:]
		cpl.tau_type[t]	= h_cpl.variables['Layertype'][i][:] 
	return cpl 


def read_gdas(epochs, lats, lons,
	path_gdas='/data/wsessions/GDAS_PROFILES', dummy=False, **kwargs):
	'''
	NOT ACTUALLY DONE.  copied from read_prof, then abandoned.

	This one DOES NOT use output from Greg Quinn's gdas2hdf script

	path_gdas	str,		path to gdas1 files output from gdas2hdf
	idx			ndarray,	index array for matched fovs with SHIS

	The time matching in this is pretty stupid.  Looks for roughly the
	nearest file within a 6 hour window.
	'''
	from numpy.ma import masked_where
	from libgeo import ppmv2kg, p2z, rh2wv
	from libclass import var
	from numpy import arange, array, recarray, zeros, append, tile
	import h5py, re
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
		with h5py.File(fname_pro, 'r') as hdf:
			sidx = hdf.get('Flight_Segment')[:] != 0
	else:
		with Dataset(fname_shis, 'r') as h_rad:
			#_if not specified, use all
			sidx = tile(True, (len(h_rad.dimensions['time']))) 

	#_if no idx specified for matching, load all
	idx = arange(sidx.sum()) if idx is None else idx

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
		handle = Dataset(fname_gdas, 'r')
		sfc.append(handle.variables['Surface_Temperature'][sidx])
		lat.append(handle.variables['Latitude'][sidx])
		lon.append(handle.variables['Longitude'][sidx])
		ozn.append(handle.variables['Ozone_Mixing_Ratio'][sidx,:])
		tmp.append(handle.variables['Temperature'][sidx,:])
		hum.append(handle.variables['Relative_Humidity'][sidx,:])
		pot.append(handle.variables['Geopotential_Height'][sidx,:])
		hdfs_gdas.append(handle)

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
	setattr(data, 'pressure', hdf.variables['Pressure'][:])

	#_close out filehandles
	[hdf.close() for hdf in hdfs_gdas]
	return data


def read_profile_gdas(fname_airs, adx_x=None, adx_l=None,
	pcfilt=True, path_gdas=os.path.join(os.environ['PRODUCTS'], 'gdas_prof'),
	dummy=False, **kwargs):
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
	from netCDF4 import Dataset
	from libtools import newdtg2, dtg2epoch, epoch2dtg, modis2epoch
	from pyhdf.SD import SD, SDC

	if not hasattr(adx_l, '__iter__'):
		adx_l = array([adx_l])
		adx_x = array([adx_x])
	print adx_l, 'WALTER', type(adx_l)
	res = re.search('AIRS.(\d{4}).(\d{2}).(\d{2}).(\d{3})', fname_airs)
	airs_dtg = '{0}{1}{2}00'.format(res.group(1), res.group(2), res.group(3))
	airs_n = res.group(4)

	#_dtg for today and tomorrow for gdas files
	dtg0, dtg1 = airs_dtg, newdtg2(airs_dtg, 24)
	dtgs_gdas = [ '{0}{1:02d}'.format(dtg_loop[:8], hr_loop)	\
				for dtg_loop in [dtg0, dtg1]					\
				for hr_loop in [0,6,12,18]]
	epochs_gdas = [ dtg2epoch(dtg_loop) for dtg_loop in dtgs_gdas ]

	#_initialize returned recarray
	dtype = [	('GDAS_latitude', 'f8'),
				('GDAS_longitude', 'f8'),
				('GDAS_epoch', 'f8'),
				('GDAS_sfc_temperature', 'f4'),
				('GDAS_ozone_mixing_ratio', ndarray),
				('GDAS_relative_humidity', ndarray),
				('GDAS_geopotential_height', ndarray),
				('GDAS_temperature', ndarray),
				('GDAS_pressure', ndarray)	]
	data = recarray((len(adx_l),), dtype)

	#_get time of profiles from radiance file
	dbg(fname_airs)
	h_rad = SD(fname_airs, SDC.READ)
	#_...
#	print '00', h_rad.select('Time')[:].copy().shape
#	print '01', adx_l.shape, adx_x.shape
#	print '02', adx_l
#	print '03', adx_x
	#_put into lists 

	epochs_airs = []
	n_l, n_x = h_rad.select('Time')[:].shape
	for L, X in zip(adx_l, adx_x):
		L, X = int(L), int(X)
		epochs_airs.append(modis2epoch(h_rad.select('Time')[L,X]))
	epochs_airs = array(epochs_airs)

	#_open all potential gdas files
	hdfs_gdas = []
	sfc, lat, lon, ozn, tmp, hum, pot, pre = [], [], [], [], [], [], [], []
	for gdas_dtg in dtgs_gdas:
		arg = (gdas_dtg, airs_n) #_WRS	
		fname_gdas = JOIN(path_gdas,'gdas_prof_airs.{0}.{1}.hdf'.format(*arg))

		#_open handle and add to list
		handle = SD(fname_gdas, SDC.READ)
		sfc.append(handle.select('Surface_Temperature')[:])
		lat.append(handle.select('Latitude')[:])
		lon.append(handle.select('Longitude')[:])
		ozn.append(handle.select('Ozone_Mixing_Ratio')[:])
		tmp.append(handle.select('Temperature')[:])
		hum.append(handle.select('Relative_Humidity')[:])
		pot.append(handle.select('Geopotential_Height')[:])
		pre.append(handle.select('Pressure')[:])
		hdfs_gdas.append(handle)

	#_loop over every AIRS time, use associated GDAS profile
	for j, epoch_airs in enumerate(epochs_airs):

		#_pull out appropriate collocated index
	##	i = idx[j]
		i = int(adx_x[j] + adx_l[j]*n_x)

		#_find closest timestamp between observation and GFS
		sec_diff = abs(epoch_airs - epochs_gdas)
		eidx = int(sec_diff.argmin())
		if sec_diff[eidx] > 3600 * 4. and not dummy:
			raise RuntimeError, 'Time between GFS and measurement > 4 hr'
		hdf = hdfs_gdas[eidx]

		#_add data into recarray
		data.GDAS_sfc_temperature[j]	= sfc[eidx][i]
		data.GDAS_latitude[j]			= lat[eidx][i]
		data.GDAS_longitude[j]			= lon[eidx][i] 
		data.GDAS_epoch[j]				= epochs_gdas[eidx]

		data.GDAS_ozone_mixing_ratio[j]	= ozn[eidx][i,:]
		data.GDAS_temperature[j]		= tmp[eidx][i,:]
		data.GDAS_relative_humidity[j]	= hum[eidx][i,:] 
		data.GDAS_geopotential_height[j]= pot[eidx][i,:]
		data.GDAS_pressure[j]			= pre[eidx][:]
	
	#_set pressure once since it is static
	setattr(data, 'GDAS_fname', fname_gdas[0])
##	setattr(data, 'GDAS_pressure', hdf.select('Pressure')[:])

	#_close out filehandles
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
	import h5py
	from time import sleep

	#_open file handle
	fname = JOIN(dir_shis, fname) if not os.path.exists(fname) else fname
	fname_pro = re.sub('rad{0}.nc'.format('_pcfilt'*pcfilt),
		'atm_prof_rtv_bc.h5', fname)
	#_principle component filtered and bias corrected

	dbg(fname)
	dbg(fname_pro)
	h_rad = Dataset(fname, 'r', format='NETCDF3_CLASSIC')

	if os.path.exists(fname_pro):
		n_pro = h5py.File(fname_pro, 'r')

		#_pull out profile segments from flight
		sidx = n_pro.get('Flight_Segment')[:] != 0
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
		setattr(shis, 'pressure', n_pro.get('Plevs')[:])
		shis.sfc_temperature[:]			= n_pro.get('TSurf')[sidx][idx]
		shis.cloud_top_temperature[:]	= n_pro.get('CTT')[sidx][idx]
		shis.cloud_top_pressure[:]		= n_pro.get('CTP')[sidx][idx]
		shis.cloud_optical_depth[:]		= n_pro.get('COT')[sidx][idx]
		tem = n_pro.get('TAir')[:,sidx]
		rel = n_pro.get('RelHum')[:,sidx]
		ozn = n_pro.get('O3VMR')[:,sidx]
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
			if not file_seg:
				return False 

	#	dbg(file_seg)
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
		cpl_fill = zeros((nc,nz))		

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
		out_label='hs3', **kwargs):
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

		#_field of view label string for prefixing directories 
		fov_str = dir_lblrtm_fmt.format(self.dtg0, self.dtg1, fov, experiment)

		#_current time
		date_fov = epoch2iso(self.SHIS_epoch[fov])

		#_gets overwritten if model run
		title = 'FOV={0:04d}, {1}'.format(fov, date_fov)
		path_fov = JOIN(DIR_PROD, 'LBL-RTM', fov_str) if dir_lblrtm is None \
				else dir_lblrtm	

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
		pname = 'dashboard_{2}_{0}_f{1:04d}.png'.format(self.dtg0,fov,out_label)
		pname = JOIN(dir_plot, pname)
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
