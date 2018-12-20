################################################################################
#_NAMES:	libaeronet.py						       #
#_PURPOSE:	Toolbox to read and manipulate aeronet data on NRL servers     #
#_USAGE:	Not a stand alone script.  Should only be imported, not run    #
#_AUTHOR:	Walter R Sessions, September 2011			       #
################################################################################
##_GENERAL_#####################################################################
################################################################################

import numpy as np
import libmeta as lm
import libtools as lt
import libclass as cl
import os
debug = lm.debug	#_Higher numbers result in noisier logs 

pnt_dict = lm.aeronet_points_icap()
mode_dict = lm.aeronet_modes()
mod_dict = lm.models()
dir_prod = lm.dir_prod

################################################################################
##_FUNCTIONS_###################################################################
################################################################################

def sort_rec( records, order=tuple ):
        ''' 
        Returns records sorted by order
        order:  A tuple of strings found in records (set by dtype, generally)
    
	This is duplicated, sorta, in libnva... check to see if there's a reason
	'''
        return records[ np.argsort( records, order=order ) ]

def read_aeronet( dtg_list, range=False, aeronet_level=15, path_aeroL2=dir_prod 
	+ '/AERO_L2', **kwargs ):
	''' 
	Read aeronet data for specified dates or period

	dtg_list	: list,	Contains dtgs(char*10) for period of interest
	range		: bool,	If True, all dtgs between earliest and latest
				in dtg_list are read. Otherwise, only ones in
				list are read.
				L2 always uses range.
	aeronet_level	: int,	15 or 2. Controls which aeronet level is read.
				L2 data should be recently downloaded and 
				converted to netcdf using la.write_aeronetL2_nc

	Returns cl.lidar() record array
	'''
	from multiprocessing import Process, Pipe
	from glob import glob

	#_initialize object
	aeronet = cl.lidar()

	#_make sure that dtg_list is, in fact, a list
	dtg_list = dtg_list if hasattr( dtg_list,'__iter__' ) else [dtg_list]
	dtg_list = sorted( dtg_list )
	dbg(( dtg_list[0], dtg_list[-1], aeronet_level ))  

	#_generate list of all possible dtgs needed
	dtg_read = []
	if range or aeronet_level == 2: #_treat dtg_list as range, no gaps
		dtg_loop = lt.newdtg( sorted( dtg_list )[0], -24 )
		dtg_end = lt.newdtg( sorted( dtg_list )[-1], 24 )
		while dtg_loop <= dtg_end:
			dtg_read.append( dtg_loop )
			dtg_loop = lt.newdtg( dtg_loop, 24 )

	else:	#_only dtgs in list
		for dtg in dtg_list:
			dtg_loop = dtg[:8] + '00' #_Files are daily
			dtg_read.append( dtg_loop )
			dtg_read.append( lt.newdtg( dtg_loop,  24) )
			dtg_read.append( lt.newdtg( dtg_loop, -24) )

	#_remove dupes	
	dtg_read = lt.unique( dtg_read ) 

	#_setup reading groups for level 1.5 data____________________________15_
	if aeronet_level == 15:
	    groups = lt.setup_groups( dtg_read, **kwargs )
 	    for group in groups:
		l = len( group )
		t = [None]*l
		p = [None]*l
		c = [None]*l
		
		data = [] 

		#_initialize forked processes
		for n in np.arange( l ):
			dtg = group[n]

			p[n], c[n] = Pipe()
			args = ( dtg, )
			kwargs.update({ 'pipe' : c[n] })
			t[n] = Process( target=parse_aeronet, 
					args=args,
					kwargs=kwargs )
			t[n].start()

		#_collect threads with file data
		for n in np.arange( l ):
			tmp = p[n].recv()

			if type(tmp) == np.recarray: data.append( tmp )

			t[n].join( 10 )

		#_append dtgs to master array	
		for dtg in data: aeronet = lt.merge(( aeronet, dtg ))

	#_setup reading groups for level 2 data_______________________________2_
	if aeronet_level == 2:
	    files = glob( path_aeroL2 + '/*_aeronet_l2.nc' )
	    groups = lt.setup_groups( files, **kwargs )
 	    for group in groups:
		l = len( group )
		t = [None]*l
		p = [None]*l
		c = [None]*l
		
		data = [] 

		#_initialize forked processes
		for n in np.arange( l ):
			file = group[n]

			p[n], c[n] = Pipe()
			args = ( file, dtg_read )
			kwargs.update({ 'pipe' : c[n] })
			t[n] = Process( target=parse_aeronetL2, 
					args=args,
					kwargs=kwargs )
			t[n].start()

		#_collect threads with file data
		for n in np.arange( l ):
			dtg = group[n]
			tmp = p[n].recv()

			if type(tmp) == np.recarray: data.append( tmp )

			t[n].join( 10 )

		#_append dtgs to master array	
		for aero_tmp in data: aeronet = lt.merge(( aeronet, aero_tmp ))

	#_return aeronet data
	return aeronet

def write_aeronetL2_nc( dir_aeroL2_raw=dir_prod+'AERO_L2/ALL_POINTS',**kwargs):
	from multiprocessing import Process, Pipe
	from glob import glob
	#_initialize object
	aeronet = cl.lidar()

	files = glob( dir_aeroL2_raw +'/*lev20' )
	dbg(( len( files ), 'file(s)' ))  

	#_setup reading groups
	groups = lt.setup_groups( files, **kwargs )
	for group in groups:
		l = len( group )
		t = [None]*l
		p = [None]*l
		c = [None]*l
		
		data = {} 

		#_initialize forked processes
		for n in np.arange( l ):
			file = group[n]

			p[n], c[n] = Pipe()
			args 	= ( file, )
			kwargs 	= { 'pipe' : c[n] }
			t[n] 	= Process( target=parse_aeronetL2_raw,
					args=args,
					kwargs=kwargs )
			t[n].start()

		#_collect threads with file data
		for n in np.arange( l ):
			dtg = group[n]
			tmp = p[n].recv()

			if type(tmp) == np.recarray: write_nc( tmp )
			t[n].join( 10 )
###
###		#_append dtgs to master array	
###		for dtg in data: aeronet = lt.merge(( aeronet, data[dtg] ))
###
###	#_return aeronet data
###	return aeronet

def readAeronet( *args, **kwargs ): read_aeronet_available( *args, **kwargs ) 

def read_aeronet_available( dtg_list, aeronet_level=None, finc=6, **kwargs ):
	'''
	read in a period of aeronet, whatever is handy. We ain't picky here.
	Use L2 when available, otherwise uses L15
	dtg_list	: list,	stand and end of period. always range.
	tolerance_dt	: int,	how far from valid time to accept as
				representative ob 
	'''
	import libnva as ln
	kwargs.update({ 'finc' : finc })

	#_read in aeronet data for period
	dtg_list = sorted( dtg_list )
	aeroL20 = read_aeronet( dtg_list, aeronet_level=2, **kwargs)
	aeroL15 = read_aeronet( dtg_list, range=True, aeronet_level=15,**kwargs)
##	aeroL15 = cl.lidar()
##	dbg('skipping 1.5')

	#_make list of dtgs we want obs for, then reduce 
	# aeronet data to only obs for those valid times
	dtg_filt = []
	dtg_str, dtg_end = dtg_list[0], dtg_list[-1]
	dtg_loop = dtg_str
	while dtg_loop <= dtg_end:
		dtg_filt.append( dtg_loop )
		dtg_loop = lt.newdtg( dtg_loop, finc )
	aeroL20 = filter_dtgs( dtg_filt, aeronet=aeroL20, **kwargs )
	aeroL15 = filter_dtgs( dtg_filt, aeronet=aeroL15, **kwargs )

	#_loop through aeronet L2
	dbg('check for overlap')
	aeronet = aeroL20 

	codes = lt.unique( aeroL15.code )
	for code in codes:
		aero_codeL15 = ln.subset( aeroL15, code=code )
		aero_codeL20 = ln.subset( aeroL20, code=code )
		for rec in aero_codeL15:
			#_see if we have it in level 2
			recL2 = ln.subset( aero_codeL20, dtg=rec.dtg,
				mode=rec.mode, wavelength=rec.wavelength )

			#_if not present in level 2, add to array
			if recL2.size == 0: aeronet = lt.merge(( aeronet, rec ))

	dbg(( aeroL15.size, aeroL20.size, aeronet.size ))
	return aeronet 

def write_nc( recs, dir_aeroL2=dir_prod + '/AERO_L2', **kwargs ):
	'''
	Takes aeronet level 2 from parse_aeronetL2_raw and puts it in ncdf file

	recs	: cl.lidar(),	Records to put in order then write to file
	'''
	from netCDF4 import Dataset
	import libnva as ln

	#_get site name for file
	code = lt.unique( recs.code, unique=True )
	nrl = lt.unique( recs.nrl_name, unique=True )
	lat = lt.unique( recs.lat, unique=True )
	lon = lt.unique( recs.lon, unique=True )

	#_initialize file object
	lt.mkdir_p( dir_aeroL2 )
	file = '/'.join(( dir_aeroL2, code + '_aeronet_l2.nc' )) 
	ncdf = Dataset( file, mode='w', format='NETCDF3_CLASSIC' )
	dbg( file )

	ncdf.lon = lon
	ncdf.lat = lat
	ncdf.code = code
	ncdf.nrl_name = nrl 
	ncdf.wavelength = lt.unique( recs.wavelength, unique=True )

	#_pull out mode and put in temporal order
	epochs 	= lt.unique ( recs.epoch )
	nt	= len( epochs )
	recs 	= ln.sort_rec( recs, order=('epoch',) )

	#_create dimensions
	ncdf.createDimension( 'epoch', nt ) 
		
	#_create dimension variable
	cdf = ncdf.createVariable( 'epoch', 'f8', 'epoch' )
	cdf.units = 'second(s) since Jan 1 1970'
	cdf[:] = epochs

	for mode in ['total','fine','coarse']:
		#_pull out mode and put in temporal order
		data = ln.subset( recs, mode=mode )
		data = ln.sort_rec( data, order=('epoch',) )

		#_create variable object and write data
		cdf = ncdf.createVariable( mode, 'f4', 'epoch',
			fill_value=-9999 )
		cdf[:] = data.tau

	ncdf.close()

def parse_aeronet( dtg, dir_aero='/aerosol_ops/native/aeronet', only_key=True,
	wavelength=550, pipe=None, **kwargs ):
   	'''
	Passed a list of dtgs, returns record array (class.lidar())
   	'''
	import libnva as ln
	import re, calendar, time

	#_initialize object
	aeronet = cl.lidar()

	#_initialize list to carry lidar objects.  Converted at return
	#_Do a song and dance and find which file is available
	GZIP = False	
	ZIP = False
	file_in = '_'.join(( dtg[:4], dtg[4:6], dtg[6:8], 'level1.5' )) 
	dir0 = '/'.join(( dir_aero, dtg[:6] ))
	dir1 = '/'.join(( dir_aero, dtg[:4] ))
	if os.path.exists( dir0 ):
		file_in = '/'.join(( dir0, file_in ))
	elif os.path.exists( dir1 ):
		file_in = '/'.join(( dir1, file_in )) 

	if os.path.isfile( file_in ): 
		pass
        elif os.path.isfile( file_in + '.gz' ):
		res = os.system( 'gunzip ' + file_in + '.gz' )
		if res != 0: raise RuntimeError, 'failed to gunzip'
		GZIP = True
	elif os.path.isfile( file_in + '.Z' ):
		res = os.system( 'zcat ' + file_in + '.Z > ' + file_in )
		if res != 0: 		
			dbg(( 'unable to extract data from', file_in ), l=3 )
			if pipe == None: 
				return -1
			else: 
				pipe.send( None )
				pipe.close()
				os._exit(-1)
		ZIP = True
        else:
          	dbg(( file_in, 'is missing' ), l=3 )
		if pipe == None: 
			return -1
		else: 
			pipe.send( None )
			pipe.close()
			os._exit(-1)	

	lamb_key	= str( int( wavelength ))
	index_site	= None	
	ncol		= None
	tau		= {}

###	site2code	= lm.site2code() 

	dbg(file_in, l=1)
	handle	= open(file_in,'r')
	lines	= handle.readlines()
	handle.close()

	if GZIP == True: 
		dbg(( 'compressing', file_in ), l=3 )
		os.system( 'gzip ' + file_in )
		os.system( 'chgrp aerosol_users ' + file_in + '.gz')

	#_Read AERONET file into temporary vars
	for line in lines:
	  	try:
	  		columns = line.split(',')
	  		if ncol != None and len(columns) < ncol: 
				break	#_Skip broke data
	  	except ValueError:
	  		dbg(( 'Incomplete aeronet line',file_in,'\n',line ),l=3)
	  		continue

		#_get indices of columns based upon header
		columns = line.split(',')
		if 'Site' in columns:	#_Define keys
			header = [ col.strip() for col in columns ]
			ncol = len(columns)
    			index_time = columns.index('Time') 
    			index_date = columns.index('Date') 
    			index_site = columns.index('Site')
    			index_lats = columns.index('Latitude')
    			index_lons = columns.index('Longitude')
    			#_Get list of AODS
    			for n in xrange( len(columns) ):
      				if columns[n][:4] == 'AOT_': 
					tmp = columns[n][4:].strip()
					tau.update( { tmp : '' } ) 
    			continue #_Skip to next line
		
		#_skip if site index not present
		elif index_site == None:
			dbg(( 'Missing header in file', file_in ), l=3 )
			dbg(( str( len(lines) )), l=3 )
			break

		#_get name from preset site dictionary
		site = columns[index_site]
		code = name2code(site)
		long_name = site
###		try:
###			code = site2code[site]
###			long_name = pnt_dict[code][1]
###		except KeyError:
###			code = 'UNKNOWN'
###			long_name = 'UNKNOWN'
 
		lat = float( columns[index_lats] )
		lon = float( columns[index_lons] )
 
  		#_use subsequent lines to fill in data 
  		date = re.search( '(\d{2}):(\d{2}):(\d{4})',
			columns[index_date])
		time = re.search( '(\d{2}):(\d{2}):(\d{2})',
			columns[index_time])
		try: 
			epoch = calendar.timegm((
				int( date.group(3) ), int( date.group(2) ),
				int( date.group(1) ), int( time.group(1) ),
				int( time.group(2) ), int( time.group(3) )))
			dtg_aero = lt.epoch2dtg( epoch )
		except ValueError: 
			dbg(( 'Invalid date/time', str( date.group(0) ), 
				str( time.group(0) )), l=3 )
			continue
  		hh = time.group(1) if time.group(2) < 30 \
			else str( int( time.group(1) )+1 ).zfill(2) 

		#_read in tau at various wavelentghs.  Ignoring 24
		for l in tau:
			index 	= header.index('AOT_' + str(l).zfill(3))
			tmp 	= columns[index]
			tmp_aod	= float( tmp ) \
				if re.search( 'N/A',tmp ) == None \
				else -9999.
	  	 	tau[l] 	= tmp_aod if tmp_aod >= 0 else -9999. 

		#_check to make sure there are at 
		# least 4 good wavelengths to pass 
		lambs = [] 	#_Fill with available wavelengths
		taus = [] 	#_Just go with me on this
		for lamb in tau:
			#_skip the wavelength we're baking out
			if lamb == lamb_key: continue
 
			if tau[lamb] > 0:
		
				#-Put into Micrometers
				tmp = float( tau[lamb] )
	      			lambs.append( int(lamb)/1000. )
	      			taus.append( tmp )

				#_Add total tau for lambda
				aeronet.resize( aeronet.size + 1 )
				aeronet[-1] = ( tmp, epoch, dtg_aero, code,
					site, long_name, int(lamb),
					'total', lat, lon )

	  	if len(lambs) < 4:
			dbg(( 'not enough data to derive 550nm tau at', site ))
   			continue
  		elif max(lambs) > 20: raise RuntimeError, 'not in micrometers'

		#_Compute the spectral polynomial (ln-ln fit)
		# e.g., ln(taus) = 	p[3]*ln(lambs)**3 +   
		#			p[2]*ln(lambs)**2 +  
		#			p[1]*ln(lambs) +  
		#			p[0] for deg 3
		if 0 in taus: continue #-don't take log of zero, boyo 

		#_calculate total, fine and coarse aod modes
		kwargs.update({ 'site' : site, 'wavelength' : wavelength })
		tau_ref, tau_f, tau_c = calc_tfc( taus, lambs, **kwargs )

		#_add total at 550
		aeronet.resize( aeronet.size + 1 )
		aeronet[-1] = ( tau_ref, epoch, dtg_aero, code, site, 
			long_name, lamb_key, 'total', lat, lon )

		#_add coarse mode at 550
		aeronet.resize( aeronet.size + 1 )
		aeronet[-1] = ( tau_c, epoch, dtg_aero, code, site, 
			long_name, lamb_key, 'coarse', lat, lon )

		#_add fine mode at 550
		aeronet.resize( aeronet.size + 1 )
		aeronet[-1] = ( tau_f, epoch, dtg_aero, code, site, 
			long_name, lamb_key, 'fine', lat, lon )
	
	#_If desired, only return results at wavelength
	if only_key: aeronet = ln.subset( aeronet, wavelength=wavelength )

	if pipe == None:
		return aeronet
	else:
		pipe.send( aeronet )
		pipe.close()

def name2code( name ):
	''' convert full name to six string code '''
	import re
	#_special rules for sites with repeated prefixes
	re_repeat = re.compile('^(ceilap|dragon)-(.*)', re.IGNORECASE )

	res = re_repeat.search( name )
	if res:
		prefix = res.group(1)
		site = res.group(2)
		code = prefix[:3] + site[:3] 
	else:
		code = re.sub('[\s_,]', '', name )[:6]

	strL = len( code )
	if strL < 6: code = code + code[-1]*(6-strL)

	return code.lower()

def parse_aeronetL2( file, dtg_list, pipe=None, **kwargs ):
	'''
	Reads aeronet level2 data from netcdf files produced by 
	la.write_aeronetL2_nc()
	
	file	: str,	Full path to aeronet L2 netcdf file
	dtg_list: list,	Two item list containing date range to read
	pipe	: multiprocessing object
	'''
	from netCDF4 import Dataset
	dtg_list = sorted( dtg_list )

	#_initialize object
	aeronet = cl.lidar()

	#_open netcdf file and get global attributes
	ncdf = Dataset( file, mode='r', format='NETCDF3_CLASSIC' )
	site = ncdf.nrl_name
	code = ncdf.code
	lamb = ncdf.wavelength
	lat = ncdf.lat
	lon = ncdf.lon

	#_read in time variable and get indices of time period of interest
	epoch0 = lt.dtg2epoch( dtg_list[0] )
	epoch1 = lt.dtg2epoch( dtg_list[-1] )
	epochs = ncdf.variables['epoch'][:]

	nidx = np.where( epochs >= epoch0 )[0]	
	xidx = np.where( epochs <= epoch1 )[0]
	idx  = np.array( lt.intersection([ nidx, xidx ]) )
	if idx.size == 0: return lt.exit_pipe( aeronet, pipe=pipe )
	epochs = epochs[idx]
	for mode in ncdf.variables:
		if mode == 'epoch': continue

		#_pull times of interest from file
		data = ncdf.variables[mode][idx]
		nrec = data.size

		#_initialize temporary recarray for data
		aero_mode = cl.lidar()
		aero_mode.resize( nrec )

		#_fill recarray
		aero_mode.tau 	= data
		aero_mode.epoch = epochs
		aero_mode.dtg	= lt.epoch2dtg( epochs )
		aero_mode.mode 	= [mode] * nrec
		aero_mode.code 	= [code] * nrec
		aero_mode.lat 	= [lat] * nrec
		aero_mode.lon 	= [lon] * nrec
		aero_mode.wavelength 	= [lamb] * nrec
		aero_mode.nrl_name 	= [site] * nrec
		aero_mode.long_name 	= [site] * nrec
		aeronet = lt.merge(( aeronet, aero_mode ))
	if pipe == None:
		return aeronet
	else:
		pipe.send( aeronet )
		pipe.close()

def parse_aeronetL2_raw( file, wavelength=550, pipe=None, **kwargs ):
	import libnva as ln
	import re, calendar, time

	#_initialize object
	aeronet = cl.lidar()

	lamb_key	= str( int( wavelength ))
	index_site	= None	
	ncol		= None

	site2code	= lm.site2code() 

	dbg( file )
	handle	= open( file, 'r' )
	lines	= handle.readlines()
	handle.close()

	#_get metadata from header
	re_header = re.compile('Location=(.*?),.*long=(.*?),.*lat=(.*?),')
	site, lon, lat = heading = re_header.findall( lines[2] )[0]

	code = name2code( site )
	long_name = site 

	#_get column index of date/time
	header = lines[4].split(',')
    	index_time = header.index('Time(hh:mm:ss)') 
    	index_date = header.index('Date(dd-mm-yy)')

    	#_get AOD_wavelength column indices
	tau_idx = {}
    	for n in xrange( len(header) ):
      		if header[n][:4] == 'AOT_':
			tmp_lamb = int( header[n][4:].strip() )
			tau_idx[ tmp_lamb ] = n 

	#_read AERONET file into temporary vars
	for line in lines[5:]:
		#_get name from preset site dictionary
		columns = line.split(',')
 
  		#_use subsequent lines to fill in data
		date =  columns[index_date]
		time =  columns[index_time]
		dd, mm, yy = date[:2], date[3:5], date[6:10] 
		hh, MM, ss = time[:2], time[3:5], time[6:8]  
		try:
			dbg(( yy, mm, dd, hh, MM, ss, date ), l=9) 
			epoch = calendar.timegm((
				int( yy ), int( mm ), int( dd ), 
				int( hh ), int( MM ), int( ss ) ))
			dtg_aero = lt.epoch2dtg( epoch )
		except ValueError: 
			dbg(( 'invalid date/time', date, time ))
			continue
		
		#_read in tau at various wavelentghs.
		#_check to make sure there are at 
		# least 4 good wavelengths to pass 
		lambs = [] 	#_Fill with available wavelengths
		taus = [] 	#_Just go with me on this
		for lamb, idx in tau_idx.iteritems():
			#_skip the wavelength we're baking out
			if lamb == lamb_key: continue
			if lamb == 1640 or lamb == 340: continue
 
			tmp 	= columns[idx]
			#tmp_aod	= float( tmp )
			try:
				tau = float( tmp )
			except ValueError:
	  	 		tau = -9999.

			if tau != -9999.:
				#_convert to micrometers
		      		lambs.append( int(lamb)/1000. )
		      		taus.append( tau )

		#_calculate total, fine and coarse aod modes
		kwargs.update({ 'site' : site, 'wavelength' : wavelength })
		vals = calc_tfc( taus, lambs, **kwargs )
		if vals.count(-9999) > 0: continue
		tau_ref, tau_f, tau_c = vals 

		#_add total at 550
		aeronet.resize( aeronet.size + 1 )
		aeronet[-1] = ( tau_ref, epoch, dtg_aero, code, site, 
			long_name, lamb_key, 'total', lat, lon )

		#_add coarse mode at 550
		aeronet.resize( aeronet.size + 1 )
		aeronet[-1] = ( tau_c, epoch, dtg_aero, code, site, 
			long_name, lamb_key, 'coarse', lat, lon )

		#_add fine mode at 550
		aeronet.resize( aeronet.size + 1 )
		aeronet[-1] = ( tau_f, epoch, dtg_aero, code, site, 
			long_name, lamb_key, 'fine', lat, lon )

	#_remove records that share the same epoch 
	keep = cl.lidar()	
	for mode in ['total','coarse','fine']:
	 	data = ln.subset( aeronet, mode=mode )
		data = ln.sort_rec( data, order=('epoch',) )
		e = data.epoch
		dt = e[1:] - e[:-1]

		if dt.tolist().count(0) == 0: 
			keep = lt.merge(( keep, data )) 
		else:
			dbg(( 'trimming', code, data.size ), l=3 )
			idx = np.where( dt != 0 )[0]
			keep = lt.merge(( keep, data[idx] ))
			keep = lt.merge(( keep, data[-1] ))
			dbg(( 'trimmedd', code, data[idx].size ), l=3 )
	aeronet = keep

	if pipe == None:
		return aeronet
	else:
		pipe.send( aeronet )
		pipe.close()

def parse_aeronet_icap( file, **kwargs ):
	from netCDF4 import Dataset
	import numpy as np
	import libtools as lt
	dbg(file)

	ncdf    = Dataset( file, mode='r', format='NETCDF3_CLASSIC' )
	nt      = len( ncdf.dimensions['time'] )
	modes   = ncdf.modes.split(',')

	#_construct recarray datatyp
	values = {}
	for mode in modes:
	        tmp = ncdf.variables[mode][:]
	        tmp = np.ma.masked_where( tmp == -9999., tmp )
	        mask = tmp.mask == False
	        values.update({ mode : tmp[mask] })

	nt = mask.tolist().count(True)

	dtype = [('epoch','i8'),('code', 'a6'),('lat','f4'),('lon','f4')]
	[ dtype.append((str(mode),'f4')) for mode in modes ]
	records         = np.recarray((nt,), dtype=dtype)
	records.dtg_start = ncdf.dtg_start
	records.dtg_end = ncdf.dtg_end
	records.code    = [ncdf.code] * nt
	records.modes   = modes
	records.lat     = [ncdf.lat] * nt
	records.lon     = [ncdf.lon] * nt

	#_some of the files didn't save the time correctly
	# this corrects that.
	epochs = ncdf.variables['time'][:].copy()
	if (epochs[1:] - epochs[:-1]).max() == 1:
	        epoch_start = np.min( epochs )
	        epochs -= epoch_start
	        epochs *= 21600
	        epochs += epoch_start
	epochs = epochs[mask]
	records.epoch = epochs[:]

	records.total   = values['total']
	records.coarse  = values['coarse']
	records.fine    = values['fine']
	return records


def calc_tfc( taus, lambs, site='UNKNOWN', wavelength=550,
	physical_forcing=True, compute_errors=True, **kwargs ):
	'''
	tau	: dict,	Contains aod values at various wavelengths

	REFERENCES
	1.    O'Neill, N. T., Eck, T. F., Smirnov, A., Holben B. N., S. 
		Thulasiraman, (2003), Spectral discrimination of coarse and 
		fine mode optical depth, JGR, Vol. 108, No. D17, 4559.
	2.    O'Neill, N. T., Dubovik, O., Eck, T. F., (2001), A modified 
		Angstrom coefficient for the characterization 
		of sub-micron aerosols, App. Opt., Vol. 40, No. 15, pp. 
		2368-2374.
	3.    O'Neill, N. T., Eck, T. F., Smirnov, A., Holben B. N., S. 
		Thulasiraman, (2005), Spectral Deconvolution 
		algorithm Technical memo.
	
	** WARNING **
	The algorithm has only been physically validated at a reference 
	wavlength of 0.5 um. The input AOD spectrum must have at least 
	4 elements. The wavelengths employed in the 0.5 um reference 
	wavelength validation exercise were [1.02 0.870 0.670 0.5 0.44 0.38] 
	um (in other words 0.34 um was systematically excluded)

	'''
	if len(lambs) < 4:
		dbg(('not enough data to derive 550nm tau at',site),l=3)
		return -9999, -9999, -9999
  	elif max(lambs) > 20: raise RuntimeError, 'not in micrometers'
	physical_forcing = True if compute_errors else physical_forcing	

	#_aeronet L2 wavelengths, ignoring 1640 and 340nm
	lamb_ref	= wavelength / 1000.
	lamb_key	= str( int( wavelength ) )
	deg_fit		= 2.
	Dtau		= 0.05

	#_Not used, but could loop this 
	alph_c		= -0.15
	alphp_c		= 0.0

	#_fine mode
	a_top		= -0.22
	b_top		= (10**(-0.2388)) * ((lamb_ref)**( 1.0275))
	c_top		= (10**( 0.2633)) * ((lamb_ref)**(-0.4683))
	a_bot		= -0.3
	b_bot		= 0.8
	c_bot		= 0.63
	a		= (a_bot + a_top) / 2.
	b		= (b_bot + b_top) / 2.
	c		= (c_bot + c_top) / 2.
	b_star		= b + 2. * alph_c * a
	c_star		= c - alphp_c + b * alph_c + a * alph_c**2 
	p = np.polyfit( np.log(lambs), np.log(taus), deg_fit )[::-1] 	
						#_rev. to keep 
						# in IDL order 
						# (for sanity)
	#_calculate AOD at desired lambda
	tau_ply	= np.sum( p*np.log(lamb_ref)**np.arange(deg_fit+1) ) 
	tau_ref	= np.exp( tau_ply ) 
	dtau_reg = 0

	#_compute alpha at reference lambda
	alph_ref = 0.
	for i in np.arange(deg_fit)+1:
		exponent = deg_fit - i
		alph_ref = alph_ref - (exponent+1) * p[ deg_fit + 1 - i] \
			* np.log(lamb_ref)**exponent

	#_compute alpha prime at reference lambda
	alphp = 0
	if deg_fit >= 2:
		for i in np.arange(deg_fit-1)+2:
			exponent = deg_fit - i
			alphp = alphp - (exponent+1) * (exponent+2) 	\
				* p[(deg_fit+1)-(i-1)] 			\
				* np.log(lamb_ref)**exponent

	#_compute derived fine and coarse mode parameters
	alph_offset	= alph_ref - alph_c
	alphp_offset	= alphp - alphp_c
	t		= alph_offset - alphp_offset/alph_offset
	alph_f		= ((t+b_star) 			\
			+ np.sqrt(( t+b_star )**2 + 4.	\
			* (1-a) * c_star)) 		\
			/ (2.*(1-a)) + alph_c
	alphap_f	= a*alph_f**2 + b*alph_f + c
	alph_f_offset	= alph_f - alph_c
	eta		= alph_offset/alph_f_offset

	tau_f		= eta * tau_ref	#_fine + coarse == total
	tau_c		= tau_ref - tau_f

	#_bias correction.  Correct the alphp bias and 
	# propagaate this correctio through all derived 
	# parameters ( as per appendix of ref 1)
	# Bias logic is x = x_true + x_bias and therfore 
	# x + -x_bias = (x_true + x_bias) + -x_bias = x_true
	#_wrs - can't this be done the first go around?
	alphp_bias_corr	= 0.65 * np.exp( -((eta-0.78)**2) / (2.*0.18**2) ) \
			if deg_fit == 2	else 0.0
	alphp_bias_corred	= alphp + alphp_bias_corr
	t_bias_corred		= alph_offset - (alphp_bias_corred - alphp_c) \
				/ alph_offset
	alphf_bias_corred	= ((t_bias_corred + b_star) \
				+ np.sqrt((t_bias_corred + b_star)**2 + 4. \
				*(1-a)*c_star))/(2.*(1-a)) + alph_c
	eta_bias_corred		= alph_offset / (alphf_bias_corred - alph_c)
	tau_f_bias_corred	= eta_bias_corred * tau_ref
	tau_c_bias_corred	= tau_ref - tau_f_bias_corred

	alphp 			= alphp_bias_corred
	alphp_offset		= alphp - alphp_c 
	eta			= eta_bias_corred
	alph_f			= alphf_bias_corred
	alphp_f			= a * alph_f**2 + b*alph_f + c
	
	tau_f = tau_f_bias_corred
	tau_c = tau_c_bias_corred
	alph_f_offset = alph_f - alph_c

	if compute_errors:
		k2 = -2.5  
		factor = alph_ref < 2   
		k1 = factor*(10 + 10/np.exp( ((alph_ref-2)/(np.sqrt(2)*.75))**2\
			 )) + 20*(1 - factor)

		# rms errors over an ensemble of measurements (see ref. 1)
		Dalpha_c = .15 
		Dalphap_c = Dalpha_c 

		# error in alphap_f due to fine-mode model 
		# uncertainty (see visible.xls for details)
		Dalphap_f = (a_top - a_bot)/2.*alph_f**2 + (b_top - b_bot) \
			/2.*alph_f + (c_top - c_bot)/2.
		D = np.sqrt((t + b_star)**2 + 4*(1 - a)*c_star)
		alphap_f_est = a*alph_f**2 + b*alph_f + c 
		temp1 = alphp_offset/alph_offset**2 + 2*a 
    		temp2 = 2*( 1 - a ) 
		dD_dalpha_c = ((t + b_star)*temp1 + 4*a*alph_c*(1 - a))/D
    		dalpha_f_dalpha_c = 1 + (temp1 + dD_dalpha_c)/temp2 
		dalpha_f_dalphap_c = (1 + (t + b_star - temp2*alph_offset)/D)\
			/(temp2*alph_offset)
		dalpha_f_dalphap_f = 1/np.sqrt(t**2 \
			+ 4*(alphap_f_est - alphp_c)) 
		dalpha_f_dalphap = -alph_f_offset/alph_offset/D
		dalpha_f_dalpha = alph_f_offset*(1 \
			+ alphp_offset/alph_offset**2)/D
		deta_dalpha_c = -1/alph_f_offset*(eta*dalpha_f_dalpha_c \
			+ (1 - eta))
		deta_dalphap_c = -1/alph_f_offset*(eta*dalpha_f_dalphap_c)
		deta_dalphap_f = -1/alph_f_offset*(eta*dalpha_f_dalphap_f) 
		deta_dalphap = -1/alph_f_offset*(eta*dalpha_f_dalphap) 
		deta_dalpha = 1/alph_f_offset*(1 - eta*dalpha_f_dalpha) 
		Dtau_rel = Dtau/tau_ref  
		Dalphap = k1*Dtau_rel 
		Dalpha = k2*Dtau_rel  

		# compute the individual contributions to the total rms errors. 
                # alpha_f contribution to the rms error
		Dalpha_f_alphap = dalpha_f_dalphap*Dalphap 
		Dalpha_f_alpha = dalpha_f_dalpha*Dalpha 
		Dalpha_f_alphap_f = dalpha_f_dalphap_f*Dalphap_f 
		Dalpha_f_alphap_c = dalpha_f_dalphap_c*Dalphap_c 
		Dalpha_f_alpha_c = dalpha_f_dalpha_c*Dalpha_c 

		# eta contribution to the rms error
		Deta_alphap = deta_dalphap*Dalphap
		Deta_alpha = deta_dalpha*Dalpha
		Deta_alpha_alphap = Deta_alphap + Deta_alpha
		Deta_alphap_f = deta_dalphap_f*Dalphap_f
		Deta_alphap_c = deta_dalphap_c*Dalphap_c
		Deta_alpha_c = deta_dalpha_c*Dalpha_c
                                
		# tau_f contribution to the rms error
		Dtau_f_alphap = Deta_alphap*tau_ref
		Dtau_f_alpha = Deta_alpha*tau_ref
		Dtau_f_alpha_alphap = Dtau_f_alphap + Dtau_f_alpha
		Dtau_f_tau = eta*Dtau 
		Dtau_f_alphap_f = Deta_alphap_f*tau_ref
		Dtau_f_alphap_c = Deta_alphap_c*tau_ref
		Dtau_f_alpha_c = Deta_alpha_c*tau_ref 
                                
		# tau_c contribution to the rms error
		Dtau_c_alphap = -Dtau_f_alphap
		Dtau_c_alpha = -Dtau_f_alpha
		Dtau_c_alpha_alphap = -Dtau_f_alpha_alphap
		Dtau_c_tau = Dtau - Dtau_f_tau
		Dtau_c_alphap_f = -Dtau_f_alphap_f
		Dtau_c_alphap_c = -Dtau_f_alphap_c
		Dtau_c_alpha_c = -Dtau_f_alpha_c 

		# compute the total rms stochastic errors
		Deta = np.sqrt(Deta_alpha_alphap**2 + Deta_alphap_f**2 \
			+ Deta_alphap_c**2 + Deta_alpha_c**2)
		Dalpha_f = np.sqrt((Dalpha_f_alphap + Dalpha_f_alpha)**2 \
			+ Dalpha_f_alphap_f**2 + Dalpha_f_alphap_c**2 \
			+ Dalpha_f_alpha_c**2)
		Dtau_f = np.sqrt((Dtau_f_alpha_alphap + Dtau_f_tau)**2 \
			+ Dtau_f_alphap_f**2 + Dtau_f_alphap_c**2 \
			+ Dtau_f_alpha_c**2)
		Dtau_c = np.sqrt((Dtau_c_alpha_alphap + Dtau_c_tau)**2 \
			+ Dtau_c_alphap_f**2 + Dtau_c_alphap_c**2 \
			+ Dtau_c_alpha_c**2) 
	else:
		Dalpha = 0 
		Dalphap = 0 
		Dtau_f = 0             
		Dalpha_f = 0          
		Dtau_c = 0           
		Deta = 0

	if physical_forcing:
		alpha_f_max = alph_f + Dalpha_f 
		alpha_f_min = alph_f - Dalpha_f 
		alpha_c_max = alph_c + Dalpha_c 
		alpha_c_min = alph_c - Dalpha_c 
	
		# from visible.xls (combined with the Rayleigh limit)
		alpha_f_max_theoretical = np.min([4,10**(.18*\
			np.log10(wavelength)+.57)])

		# bubble of unphysical problems
		if (alpha_f_min < alph_ref) or (alpha_c_max > alph_ref):
			# near eta = 1 modify the nominal alpha_f 
			# starting at the point where its lower 
			# error bar fouls alpha

			# alpha in between alpha_f_min and alpha_f_max
			if alpha_f_min < alph_ref and alpha_f_max > alph_ref:
				alph_f = (alph_ref + alpha_f_max)/2. 
				alpha_f_min = alph_ref 
			else:  
				# alpha is above alpha_f_max
				if (alpha_f_max <= alph_ref):
					alph_f = alph_ref
					alpha_f_max = alph_ref
			# alpha in between alpha_c_min and alpha_c_max
			if alpha_c_max > alph_ref and alpha_c_min < alph_ref:
				alph_c = (alph_ref + alpha_c_min)/2.
			else:
				# alpha is below alpha_c_min
				if alpha_c_min >= alph_ref: alph_c = alph_ref

		# upper limit and possibly lower limit 
		# of alpha_f is beyond theoretical max.
		if alpha_f_max > alpha_f_max_theoretical:
			# alpha_f_min and alpha_f_max straddle theoretical max.
			if alpha_f_max > alpha_f_max_theoretical and \
			alpha_f_min < alpha_f_max_theoretical:
				alph_f = (alpha_f_max_theoretical + \
					alpha_f_min)/2.
				alpha_f_max = alpha_f_max_theoretical ;
			else:
				# all of the error bar is above 
				# alpha_f_max_theoretical
				if alpha_f_min > alpha_f_max_theoretical: 
					alph_f = alpha_f_max_theoretical
					alpha_f_min = alpha_f_max_theoretical

		# this could be zero if, before the corrections 
		# above, alpha_f < alpha < alpha_c ... a very unlikely event
		alph_f_offset 	= alph_f - alph_c 
		alph_offset	= alph_ref - alph_c 
		eta 		= alph_offset/alph_f_offset
    		tau_f 		= eta*tau_ref  
		tau_c 		= tau_ref - tau_f  

		# force alphap_f to be coherent with changes in alpha_f
		alphap_f = a*alph_f**2 + b*alph_f + c

	fmt = '%10.4f'
	ftt = '%-10s'
 	tt, ff = [site], [site]
	for n in xrange( len(taus) ): 
		tt.append( taus[n] )
		ff.append( lambs[n] )
		ftt += fmt
	dbg( fmt*3 % (tau_ref, tau_c, tau_f), l=7 )
	dbg( ftt % tuple( tt ), l=7 )
	dbg( ftt % tuple( ff ), l=7 )
	if tau_ref < tau_f or tau_ref < tau_c:
		dbg( 'negative tau value', l=2)
		dbg(( 'site', site, eta, tau_ref, tau_c, tau_f ), l=2)
		dbg(( 'taus', site,  taus ), l=2)		
		dbg(( 'lambs', site, lambs ), l=2)		
		return -9999, -9999, -9999 
	else:
		return tau_ref, tau_f, tau_c

def plot_aeronet_period( records, path_out='./plots' ):
	''' plots timeseries for a aeronet data '''
	import libnva as ln
	from numpy import ma
	import matplotlib.pyplot as plt
	#_Make output directory
	lt.mkdir_p( path_out )

	dtg_start = yyyymm + '0100'

	#_Find last day of month. I don't use a lookup table because
	# I'd rather it all leap year nonsense be handled by the calendar
	# module.  It's a matter of preference.
	dtg_loop = dtg_start
	ym_loop = dtg_loop[:6]
	dtg_list = []
	while ym_loop == yyyymm:
		dtg_list.append( dtg_loop )
		dtg_loop = lt.newdtg( dtg_loop, 24 )
		ym_loop = dtg_loop[:6]
	else:
		#_Once the month changes, go back a day
		dtg_end = lt.newdtg( dtg_loop, -24 )	

	#_Read aeronet data if not passed
	if aeronet == None:
		aeronet = read_aeronet( dtg_list )

	#_Check for datapoints outside of month of interest
	ind_end = np.where( aeronet.dtg < dtg_end )[0]
	ind_str = np.where( aeronet.dtg >= dtg_start )[0]
	indices = lt.intersection([ ind_end, ind_str ])
	aeronet = aeronet[indices]

	#_Loop over available sites, plot
	sites_avail = lt.unique( aeronet.code )
	for code in sites_avail:
		if code == 'UNKNOW': continue

		#_Subset only those sites
		aero_sub = ln.subset( aeronet, code=code )

		#_Sort by epoch
		aero_srt = sort_rec( aero_sub, order=('epoch') )

		#_Make abscissa seconds for each month
		epoch_start = lt.dtg2epoch( dtg_start )
		epoch_end = lt.dtg2epoch( dtg_end ) + 86399
		nt = epoch_end - epoch_start 
		x = np.arange( nt )

		#_Label ids
		cs = {}

		#_Setup aspect ratio
		h, w = plt.figaspect( 20 )
		fig = plt.figure( figsize=(w,h))
		ax = plt.subplot( 111 )

		#_Loop over modes
		for mode in mode_dict:
			#_Pull out total/fine/coarse
			aero_mode = ln.subset( aero_srt, mode=mode )	

			#_Initialize plotting array
			y = np.empty((nt,))
			y[:] = -9999. 

			#_Use epochs mine epoch_start as indexing values
			y[ aero_mode.epoch - epoch_start ] = aero_mode.tau

			#_Plot data
			col = mode_dict[mode]['color']
			y = ma.masked_where( y == -9999., y )
			cs[mode] = ax.scatter( x, y , color=col, s=1.,
				marker='x', edgecolors=col )

		#_Max plot y value
		ymax_tau = 4.

		#_Limit size
		ax.set_xlim( 0, nt )
		ax.set_ylim( 0, ymax_tau )

		#_Setup axis label and header
		title = 'AERONET TAU for ' + code.upper()
		ax.set_title( title, size='x-small' )
		ax.set_ylabel('TAU', size='xx-small' )

		#_X frequency in days and seconds
		frq_x = 2.5 
		frq_xs = frq_x * 86400
		frq_y = 0.5
		mjl = plt.MultipleLocator( frq_xs )	#_Major Label
		ax.xaxis.set_major_locator( mjl )
		mjl = plt.MultipleLocator( frq_y )	#_Major Label
		ax.yaxis.set_major_locator( mjl )
		ax.grid( True )#, zorder=-1 )

		#_Setup grid and tick labels
		tick_lab = dtg_list[1::5]
		tick_loc = x[86400::5*86400]
		plt.xticks( tick_loc, tick_lab, size='xx-small' )
		tick_lab = np.linspace( 0, ymax_tau, ymax_tau / frq_y + 1 )
		plt.yticks( tick_lab, size='xx-small' )

		#_Create output file name
		file_out = path_out + '/' + yyyymm + '_' + code + '.png'

		#_Save figure
		dbg( file_out )
		plt.savefig( file_out )
		plt.close()
	
def plot_aeronet_month( yyyymm, aeronet=None, path_out='./plots' ):
	''' plots timeseries for a month of aeronet data '''
	import libnva as ln
	from numpy import ma
	import matplotlib.pyplot as plt
	#_Make output directory
	lt.mkdir_p( path_out )

	dtg_start = yyyymm + '0100'


	#_Find last day of month. I don't use a lookup table because
	# I'd rather it all leap year nonsense be handled by the calendar
	# module.  It's a matter of preference.
	dtg_loop = dtg_start
	ym_loop = dtg_loop[:6]
	dtg_list = []
	while ym_loop == yyyymm:
		dtg_list.append( dtg_loop )
		dtg_loop = lt.newdtg( dtg_loop, 24 )
		ym_loop = dtg_loop[:6]
	else:
		#_Once the month changes, go back a day
		dtg_end = lt.newdtg( dtg_loop, -24 )	

	#_Read aeronet data if not passed
	if aeronet == None:
		aeronet = read_aeronet( dtg_list )

	#_Check for datapoints outside of month of interest
	ind_end = np.where( aeronet.dtg < dtg_end )[0]
	ind_str = np.where( aeronet.dtg >= dtg_start )[0]
	indices = lt.intersection([ ind_end, ind_str ])
	aeronet = aeronet[indices]

	#_Loop over available sites, plot
	sites_avail = lt.unique( aeronet.code )
	for code in sites_avail:
		if code == 'UNKNOW': continue

		#_Subset only those sites
		aero_sub = ln.subset( aeronet, code=code )

		#_Sort by epoch
		aero_srt = sort_rec( aero_sub, order=('epoch') )

		#_Make abscissa seconds for each month
		epoch_start = lt.dtg2epoch( dtg_start )
		epoch_end = lt.dtg2epoch( dtg_end ) + 86399
		nt = epoch_end - epoch_start 
		x = np.arange( nt )

		#_Label ids
		cs = {}

		#_Setup aspect ratio
		h, w = plt.figaspect( 20 )
		fig = plt.figure( figsize=(w,h))
		ax = plt.subplot( 111 )

		#_Loop over modes
		for mode in mode_dict:
			#_Pull out total/fine/coarse
			aero_mode = ln.subset( aero_srt, mode=mode )	

			#_Initialize plotting array
			y = np.empty((nt,))
			y[:] = -9999. 

			#_Use epochs mine epoch_start as indexing values
			y[ aero_mode.epoch - epoch_start ] = aero_mode.tau

			#_Plot data
			col = mode_dict[mode]['color']
			y = ma.masked_where( y == -9999., y )
			cs[mode] = ax.scatter( x, y , color=col, s=1.,
				marker='x', edgecolors=col )

		#_Max plot y value
		ymax_tau = 4.

		#_Limit size
		ax.set_xlim( 0, nt )
		ax.set_ylim( 0, ymax_tau )

		#_Setup axis label and header
		title = 'AERONET TAU for ' + code.upper()
		ax.set_title( title, size='x-small' )
		ax.set_ylabel('TAU', size='xx-small' )

		#_X frequency in days and seconds
		frq_x = 2.5 
		frq_xs = frq_x * 86400
		frq_y = 0.5
		mjl = plt.MultipleLocator( frq_xs )	#_Major Label
		ax.xaxis.set_major_locator( mjl )
		mjl = plt.MultipleLocator( frq_y )	#_Major Label
		ax.yaxis.set_major_locator( mjl )
		ax.grid( True )#, zorder=-1 )

		#_Setup grid and tick labels
		tick_lab = dtg_list[1::5]
		tick_loc = x[86400::5*86400]
		plt.xticks( tick_loc, tick_lab, size='xx-small' )
		tick_lab = np.linspace( 0, ymax_tau, ymax_tau / frq_y + 1 )
		plt.yticks( tick_lab, size='xx-small' )

		#_Create output file name
		file_out = path_out + '/' + yyyymm + '_' + code + '.png'

		#_Save figure
		dbg( file_out )
		plt.savefig( file_out )
		plt.close()

def plot_model_timeseries( records, aeronet=None, plot_members=True, path='.',
	codes=None, aspect=40, finc=1, **kwargs ): 
	'''
	Plot aeronet values with model timeseries at aeronet points
	records		: model_object(), det or ensemble
	aeronet		: lidar(), if called multiple times, pass aeronet data
	plot_members	: bool, determines if ensemble members are plotted
	path		: str,	Output directory for plots
	code		: list,	Aeronet codes to plot, None for all 
	dtg_start	: str*10
	dtg_end		: str*10
			: These two could be derived from various versions
				of dtg_init and dtg_vald, but the naming 	
				requirements for certain output resulted
				in it being easier to include these to 
				explicitly give these values.  Primary
				influence on filenames.
	fcst		: bool,	Special case that affects file naming;
				when True, the last dtg_init available
				is used as dtg instead of dtg_start
	finc		: int,	hours between datapoints (xvalues)
	If multiple .models are passed, they will be put on the same plot.

	Ensembles will also have all members plotted with mean. 

	ToDo: 	Cleanup the kludiness with naming between verif, forecasts,
		and straight up analyses 
	'''
	import matplotlib.pyplot as plt
	from numpy import ma
	import libicap as li
	import libnva as ln
	from scipy.interpolate import interp1d
	import pylab, re, time

	dbg(( 'records', str( records.size ) ), l=3 )

	codes = pnt_dict if codes == None else codes

	#_Create list of all dtg we might need aeronet data for
	dtg_list = lt.unique( records.dtg_vald )
	dtg_inits = lt.unique( records.dtg_init )

	#_Define first and last valid times and convert to epochs
	dtg_start = dtg_list[0]
	dtg_end	= dtg_list[-1] 
	epoch_start = lt.dtg2epoch( dtg_start )
	epoch_end = lt.dtg2epoch( dtg_end )

	fincs = finc * 3600
	nt = (epoch_end - epoch_start) / fincs + 1 
	#_Get list of species to plot
	species = lt.unique( records.variable )

	#_make list of dtgs
	epoch_x = np.arange(nt) * fincs + epoch_start
###	epoch_x = np.linspace( epoch_start, epoch_end, 

	#_read in and bin aeronet
	if aeronet == None: aeronet = read_aeronet( dtg_list )
	aeronet = filter_dtgs( lt.epoch2dtg( epoch_x ), aeronet=aeronet, 
		tolerance_dt=finc/2., **kwargs )

	#_Fcst is used to define naming schemes of files for ICAP,
	# if the std of them is greater than 12, plot for general ICAP.
	#_The mean is not used to allow for similar plots using, for example,
	# 72-96 hour forecast data instead of 0-24
	#_anal == true if dtg_inits > 1
	#_Get list of dtgs, should be 1 for daily fcsts, and the last will
	# be the naming key for the verif plots
	fcst = True if np.std( lt.unique( records.fhr )) > 12 else False
	anal = True if len( dtg_inits ) > 1 else False

	dtg = dtg_inits[-1] if fcst else dtg_inits[0] 

	#_Get string of forcast length (max)
	max_fhr = str( np.max( records.fhr )).zfill(3)

	#_Limit to timer period of interest (could also masked_outside())
	ind_end = np.where( aeronet.epoch < epoch_end )[0]
	ind_str = np.where( aeronet.epoch >= epoch_start )[0]
	indices = lt.intersection([ ind_end, ind_str ])
	aeronet = aeronet[indices]

	#_Initialize array used for plotting
	x = np.arange( nt ) 
	y = np.empty((nt,))
	y[:] = -9999. 

	#_Get list of available AERONET size modes
	modes_avail = lt.unique( aeronet.mode )

	#_Create re object to strip member names of _m##
	re_mem = re.compile( '_m\d{2,3}' )

	#_loop over species
	for s in species:
		#_only do these for total_aod right now
		if anal and s != 'total_aod': continue

		#_subset model data by species
		fcst_spec = ln.subset( records, variable=s )

		#_if more than one model is passed, 
		# give it generic name for output 
		models = [ tmp.upper() \
			for tmp in lt.unique( fcst_spec.model ) ]
		if len( models ) != 1: 
			model = 'CUSTOM'
			dir_m = 'custom_01'
		else:
			model = models[0]
			dir_m = mod_dict[model]['web']

		#_plot MODEL data
		for m in models:
			#_subset by model, order by valid time
			fcst_mod = ln.subset( fcst_spec, model=m )
			fcst_mod = ln.sort_rec( fcst_mod, order=('dtg_vald',) )
 
			#_check if ensemble, setup special vars
			if lt.unique( fcst_mod.ensemble ) != [0]:
				nens = fcst_mod.ensemble[0]
				id = fcst_mod.dimname[0].index('member')
				nens_max = fcst_mod.dimsize[0][id]

				is_ens = True
			else:	
				is_ens = False

			#_get coordinates for model at point
			lats = fcst_mod.values[0].lat
			lons = fcst_mod.values[0].lon
			mbs = fcst_mod.values[0].member.tolist()

			#_get list of epochs to act as indices
			epochs = lt.dtg2epoch( fcst_mod.dtg_vald )

			#_find where to mask
			gaps = epochs[1:] - epochs[:-1]
			dt = np.min( gaps )
			gaps_to_mask = np.where( gaps > (dt+1e-4) )[0]

			#_warn if multiple values for single dtg
			if len( lt.unique( fcst_mod.dtg_vald )) != \
				len( lt.unique( epochs) ):
				#_Did you have to search for this gem?
				dbg( 'multiple values', l=3 )

			#_flatten recarray into +1D array
			#_this takes a long ass time. don't
			fcst_flt = ln.join_values( fcst_mod )

			#_take average if ensemble
			if is_ens:
				fcst_plt = fcst_flt.mean( axis=id+1 ) 
			else:
				fcst_plt = fcst_flt	

			#_initialize plot array
			fcst_idx = ( epochs - epoch_start ) / fincs  

			#_Generate mask array values
			mask = y.copy()
			mask[:] = False 
			n_ext = 0
			for gap_idx in gaps_to_mask:
				#_add compounding offset
				gap_idx += n_ext

				#_where does the gap start?
				idx0 = gap_idx * dt

				#_how many indices is gap?
				n_idx = int( gaps[gap_idx] / dt )

				idx1 = ( gap_idx+n_idx ) * dt 

				#_set mask between those to true
				mask[idx0:idx1] = True

				#_increment offset			
				n_ext += n_idx - 1 

			#_Loop over points, plot timeseries
			for p in pnt_dict:
				#_skip undesired
				if p not in codes: continue

				#_Setup aspect ratio
				if anal and not fcst: 
					h, w = plt.figaspect( aspect )
					fig = plt.figure( figsize=(w,h))
				ax = plt.subplot( 111 )

				#_Subset aeronet by point
				aero_pnt = ln.subset( aeronet, code=p )
	
				#_Get long_name
				long_name = pnt_dict[p][1]
	
				#_Get coordinates
				lon, lat = pnt_dict[p][0] 
				i, j = lt.ll2ij( lat, lon, lats, lons )

				#_Create
				cs = {}

				#_output file name
				dir_out = '/'.join((path, dir_m, 'spag_pnt',p, 
					dtg[:6]))
				lt.mkdir_p( dir_out )
				if fcst: #_Standard forecast Timseries
					file_out = '_'.join(( dir_out+'/'+dtg,
						'f' +max_fhr, s, p, 
						model.lower()+'.png' ))
				elif anal:	#_Timeseries with funky names
					file_out = '_'.join(( dir_out+'/'+dtg
						+'-' +dtg_end, s, p,
						model.lower()+'.png' ))

				#_Plot AERONET data
				for mode in modes_avail:
					#_subset by mode
					aero_mode = ln.subset( aero_pnt,
						mode=mode )	

					#_loop if no data
					if aero_mode.size == 0: continue

					#_filter set dtg to mode bin time
					aero_vald = lt.dtg2epoch(aero_mode.dtg)

					#_Initialize plot array
					y_a = y.copy()
					aidx = (aero_vald - epoch_start) / fincs
					y_a[ aidx ] = aero_mode.tau

					#_Mask missing values and plot
					col = mode_dict[mode]['color']
					y_a = ma.masked_where( y_a==-9999.,y_a )
					cs[mode] = ax.scatter( x, y_a, 
						color=col, s=10., marker='x',
						edgecolors=col, lw=0.6 )

				#_Plot ensemble members
				if is_ens and plot_members:
					#_For non-icap ensembles
					ens_col = lt.colors_spag( nens_max )

					#_Loop over ensemble members, plot
					for mem in mbs:
						e = mbs.index( mem )
						k = re_mem.sub( '', mem )

						#_Interpolate ndtg to nepoch
						f = interp1d( fcst_idx, 
							fcst_flt[:,e,j,i] )
						y_e = ma.masked_where(mask,f(x))

						#_Select plot color
						if model == 'ICAP':	
							c = mod_dict[k]['color']
						else:
							c = ens_col[e]

						#_Mask missing values and plot
						e_str = str(e).zfill(3)
						cs[mem] = ax.plot( x, 
							y_e, color=c, ls='--')

				#_interpolate ndtg to nepoch
				f = interp1d( fcst_idx, fcst_plt[:,j,i] )

				#_mask data
				y_m = ma.masked_where( mask, f(x) )

				#_plot Ensemble Mean 
				col = mod_dict[m]['color']
				cs[m] = ax.plot( x, y_m, color=col, ls='-' )

				#_setup legend
				pylab.rcParams.update( { 'legend.fontsize' : 8})
				leg_labels = [ re_mem.sub( '', k ) for k in cs ]
				leg_labelz = leg_labels[:]

				leg_labels = lt.unique( leg_labels ) 
				a, b = [], []
				for name in leg_labels:
					nlabel = name if name not in mod_dict\
						else mod_dict[name]['label'] 
##					nlabel = mod_dict[name]['label']\
##					if name in mod_dict else name
					a.append( nlabel )
					b.append( cs[name]  )
	
				a, b = tuple(a), tuple(b) 
				l1 = ax.legend( b, a, loc=1, ncol=2,
					borderpad=0.5, fancybox=True,
					labelspacing=0.1 )
	
				#_setup axes size
				y_max = np.max( [ np.max( fcst_plt[:,j,i])+0.2,
					3. ] )
				y_max = 3 if y_max == -9999. else y_max 
				x_max = nt-1
				ax.set_ylim([ 0, y_max ])
				ax.set_xlim([ 0, x_max ])	
	
				#_If doing verification plot, plot vertical line
				# at last dtg_init
				if fcst:
					xv = ( lt.dtg2epoch( dtg_inits[-1] ) \
						- epoch_start ) / fincs
					ax.plot( [xv, xv], [0, y_max], 'black',
						ls=('-') )
					frq_x = 86400 * 2 / fincs
	
					#_change filename
					if anal: file_out = file_out.replace( s,
							'verif')
				else:
					frq_x = 86400 * 1 / fincs 
	
				#_setup x-grid
				mjl = plt.MultipleLocator( frq_x )
				mnl = plt.MultipleLocator( frq_x / 6 )
				ax.xaxis.set_major_locator( mjl )
				ax.xaxis.set_minor_locator( mnl )
	
				#_setup x tickmarks
				tick_loc = np.arange( frq_x, x_max, frq_x )	
				tick_lab = lt.epoch2dtg( tick_loc*fincs 
					+ epoch_start )
				if (x_max % frq_x) == (frq_x-1): 
					tick_loc = np.append( tick_loc, x_max )
					tick_lab = np.append( tick_lab, dtg_end)
				tick_lab = lt.date_label( tick_lab )
				ax.set_xticks( tick_loc ) #, size='xx-small' )
				ax.set_xticklabels( tick_lab, size='xx-small' )

				#_setup y-grid
				frq_y = 0.25
				mjl = plt.MultipleLocator( frq_y )
				mnl = plt.MultipleLocator( frq_y / 5 )
				ax.yaxis.set_major_locator( mjl )
				ax.yaxis.set_minor_locator( mnl )
	
				#_setup y tickmarks (default is actually ok)
				tick_lab = np.round( np.arange( 0, y_max, 
					frq_y ), 2 )
				ax.set_yticklabels( tick_lab, size='x-small' )
	
				#_turn on Grid
				ax.grid( True )
	
				#_setup header
				title = 'Modeled ' + s.upper() + ' fields at ' \
					+ long_name + '\n' +lt.human_date( dtg )
				ax.set_title( title, va='bottom', ha='left',
					position=( 0.0, 1.0 ), size='small' )
		
				#_setup footer
				n       = time.gmtime( time.time() )    
       				cdtg    = str( n.tm_year ) \
					+ str( n.tm_mon ).zfill(2) \
					+ str( n.tm_mday ).zfill(2) \
					+ str( n.tm_hour ).zfill(2)
                	        footer  = 'Plot generated ' \
					+ lt.human_date( cdtg ) \
                	                + ' NRL/Monterey Aerosol Modeling\n'
                	        plt.figtext( 0.97, 0.85, footer, size='x-small',
                	                rotation='vertical' )

				title_y = s.upper() + ' [550nm]'
				ax.set_ylabel( title_y, size='small' )

				#_Save figure
				dbg( file_out )
				try:
					plt.savefig( file_out )
					lt.make_readable( file_out )
				except:
					dbg(( 'failed', file_out ))
				plt.close()

def plot_forker( records, aeronet=None, nproc=4, **kwargs ):
	'''
	When passed records

	aeronet: recarray containing aeronet data
	records: forecast recarray 

	To do forecast timeseries, just past a single forecast
	to this module

	To do verif style plots, you know, pass the whole thing.  I don't know
	'''
	import libnva as ln

	kwargs.update({ 'nproc' : nproc })
	codes = pnt_dict.keys()
	groups = []
	for i in xrange( nproc ): groups.append( codes[i::nproc] )

	#_read in aeronet data if not passed
	dtg_list = lt.unique( records.dtg_vald )
	if aeronet == None: aeronet = read_aeronet( dtg_list )

	dbg(( str( len(groups) ), 'group(s)' ))

	#_loop over processing groups
	for group in groups:
		dbg(( str( len( group )), 'processes'), l=7 )
	
		#_initialize list to keep pids of children procs
		children = []
		pid = os.fork()
		if pid != 0:
			dbg(( group, str(pid)), l=7 )
			children.append( pid )
		elif pid == 0:
			aero = ln.subset( aeronet, code=group )
			plot_model_timeseries( records, aeronet=aero,
				codes=group, **kwargs )
			os._exit(0)
		
	#_pause until all children processes finish
	for kid in children: os.waitpid( kid, 0 )

####	#_Split all potential points into forkable groups
####	variables = lt.unique( records.variable ) 
####	
####	groups = lt.setup_groups( variables, **kwargs )
####	dbg(( str( len(groups) ), 'group(s)' ))
####
####	#_Loop over processing groups
####	for group in groups:
####		dbg(( str( len( group )), 'processes'), l=7 )
####	
####		#_Initialize list to keep pids of children procs
####		children = []
####		for variable in group:
####			pid = os.fork()
####			if pid != 0:
####				dbg(( variable, str(pid )), l=7 )
####				children.append( pid )
####			elif pid == 0:
####				recs = ln.subset( records, variable=variable )
####				plot_model_timeseries( recs, aeronet=aeronet,
####					**kwargs )
####				os._exit(0)
####		
####		#_Pause until all children processes finish
####		for kid in children: os.waitpid( kid, 0 )
####
###	#_Split all potential points into forkable groups
###	codes = pnt_dict.keys()
###	codes = codes if code == None else lt.intersection([code, codes])
###	groups = lt.setup_groups( codes, **kwargs )
###	
###	dbg(( str( len(groups) ), 'group(s)' ))
###
###	#_Loop over processing groups
###	for group in groups:
###		dbg(( str( len( group )), 'processes'), l=7 )
###	
###		#_Initialize list to keep pids of children procs
###		children = []
###		for code in group:
###			pid = os.fork()
###			if pid != 0:
###				dbg(( code, str(pid )), l=7 )
###				children.append( pid )
###			elif pid == 0:
###				plot_model_timeseries( records, aeronet=aeronet,
###					code=code, **kwargs )
###				os._exit(0)
###		
###		#_Pause until all children processes finish
###		for kid in children: os.waitpid( kid, 0 )

def synchronous_aeronet_data( dtg_vald, aeronet=None, tolerance_dt=2, **kwargs):
	'''
	Returns boolean based upon if there is aeronet data within 
	tolerance_dt hours of valid time

	Ok, returns index of aeronet record, which POTENTIALLY could be
	0, which creates a bit of a problem.  Now return a list, which
	won't trigger False
	'''
	#_Find smallest difference between aero.epoch and fcst time
	epoch = lt.dtg2epoch( dtg_vald )
	index = np.argmin( np.fabs( aeronet.epoch - epoch ) )
	min_diff = np.fabs( aeronet.epoch[index] - epoch )

	#_False if outside of window
	if min_diff / 3600 > tolerance_dt:
		return False
	#_Add dtgs to list to use if within window
	else:
		return [ index, min_diff ]


def closest_aeronet_site( lat, lon ):
	"""
	Given coordinates, find nearest aeronet site (that is in our dictionary)
	site_code_string, distance(km) 
		= closest_aeronet_site(lat_float, lon_float)
	
	Site codes are six character identifiers we use throughout the website.
	e.g., capove => Capo Verde

	This is a really bad way of doing this, but it's Good Enough For Now(tm)
	"""
        min_dist = [999999999999999,'VOID']
        for site in pnt_dict:	#_Checks ALL sites
          	lon_site, lat_site = pnt_dict[site][0]
          	lon_diff = lon_site - lon
          	lat_diff = lat_site - lat
          	if abs(lon_diff) > 180:	#_Work across date line
            		lon_shift = lon_site + 360 if lon_site < 0 \
				else lon_site - 360
            		lon_diff = lon_shift - lon
	  	km_x = lon_diff * 111 * np.cos(np.pi*lat/180)
	  	km_y = lat_diff * 111
          	dist = km_x**2 + km_y**2 	#_no need to sqrt all 
						# since this is a mag chk
          	if dist < min_dist[0]: min_dist = [ dist, site ]
	distance = np.round( np.sqrt(min_dist[0]), 2 )	#_only do sqrt once 
	site_code = min_dist[1]
	coords = pnt_dict[site][0] 
	return site_code, distance, coords 

def representativeRecords( aeronet ):
	''' return recarray with each site represented only once '''
	aerorep = cl.lidar()
	for rec in aeronet:
		if rec.code not in aerorep.code:
			aerorep = lt.merge(( aerorep, rec ))
	return aerorep

def closestAeronet( lat, lon, aeronet, degree=1 ):
	'''
	finds nearest two aeronet sites... sorta...
	lat	: flt,	latitude of point
	lon	: flt, 	longitude of point
	aeronet	: cl.lidar(),	recarray containing aeronet data after being
				processed by representativeRecords()
	degree	: int,	how many of the closest sites to return (not inc. min)
	'''
	if lon > 180: lon -= 360.
	code = aeronet.code
	latc = aeronet.lat
	lonc = aeronet.lon

	#_calculate the distance in km in the crudest way possible
	km_x = ( lat - latc ) * 111.
	km_y = ( lon - lonc ) * 111. * np.cos( np.pi*latc/180. )
	dist = np.sqrt( km_x**2 + km_y**2 )

	#_sort the codes and distances
	sort = dist.argsort()
	code = code[sort]
	dist = dist[sort]

	latc = latc[sort]
	lonc = lonc[sort]
	dbg(( lat, latc[:degree] ), l=5 )
	dbg(( lon, lonc[:degree] ), l=5 )

	return code[:degree], dist[:degree]

def generate_components( data, dtg_start='0', dtg_end='2050013100',aeronet=None,
        pipe=None, error_models=None, score_members=True, bin=None, 
	**kwargs ):
        '''
        data		: libnrl.model_object()	Contains record array forecasts

        dtg_start	: str*10,	Used to limit which dtgs are scored
        dtg_end		: str*10,	These are generally NOT used because
					this will only score for the DTG_VALDS 
					in data
        aeronet		: class.aeronet() 	object, If calling multiple 
				times, 
				best to read aeronet 
                        externally and pass it using this keyword       

	Returns residuals, gross, and normalized residuals for scoring against
	aeronet as aeronet_comps object

	If you want to limit what gets scored, use ln.subset() to reduce 
	the number of records passed
        '''
	import time
	import libnva as ln

	dbg(( data.size, 'record(s)' ))

	#_I/O should be top loop, and right now it's sitting in ln.write_stats
	# which should be changed later, but not today...
	#_When that is done, thresholds can be placed at a less awkward 
	# part of the code

	#_Get 
	thresholds = lm.thresh_bins( bin=bin )

	if error_models == None:
		error_models = ['diagnostic','prognostic']
	else:
		error_models = [error_models] \
			if type(error_models) == str else error_models

        #_Get list of dtgs we have data for_
        dtg_list = lt.unique( data.dtg_vald )	#_For initial testing
        dtg_valds = []				#_For score iteration
        for dtg in dtg_list:
                if dtg >= dtg_end or dtg <= dtg_start: continue
                dtg_valds.append(dtg)
	dbg(( dtg_list[0], str( len( data )), 'records' ), l=3 )

        #_Read aeronet data if not passed
        if aeronet == None: aeronet = read_aeronet( dtg_valds )

	#_Filter out observations too far from vald times
	aeronet = filter_dtgs( dtg_valds, aeronet=aeronet )

	#_Create list of points both available and in ICAP list
	points_avail = lt.unique( aeronet.code )
	points = lt.intersection( [ pnt_dict, points_avail ] )

	#_Get aeronet modes present
	modes = lt.unique( aeronet.mode )

	#_Ensure we're only using 550nm
	number_of_wavelengths = len( set( aeronet.wavelength ) )
	if number_of_wavelengths > 1:
                raise ValueError, 'More than 550nm in aeronet data'

	name_dbug = '_'.join(( dtg_list[0], dtg_list[-1], 'debug.txt' ))
	f_dbug = open( name_dbug, 'w' )

        #__________________________________________________CALCULATE_COMPONENTS_
	#_loop over fcst records
	comps = cl.aeronet_comps()
	for rec in data:
		dtg_vald 	= rec.dtg_vald
		variable 	= rec.variable
		fhr 	= rec.fhr
		model 	= rec.model
		epoch 	= lt.dtg2epoch( dtg_vald )
		full	= rec.values	

		#_get model dimensions
        	lats = full.lat
        	lons = full.lon

		#_crude way of checking ensemble
        	is_ens = False if rec.ensemble == 0 else True
        	if is_ens:
			ens_idx = rec.dimname.index( 'member' )
			nens 	= rec.ensemble #_dimsize[ens_index]
			members = tuple( full.member ) #.tolist()
     		        fcst 	= full.mean( axis=ens_idx )
       		else:
       			nens = 0        #_Skips loops later
            		fcst = full 
		
		#_For sanity
		dbg(( model, dtg_vald, variable, str(fhr) ), l=3 )

		#_limit to aeronet points at this time
		# filter_dtgs() sets dtg of lidar object to 
		# DTG of associated fcst
		aero_tmp = ln.subset( aeronet, dtg=dtg_vald )
		#_MAKE THIS A 2D MAP OF AERONET VALUES, THIS LOOP IS EXPENSIVE 
		# THEN MASK THAT OBJECT BASED UPON ERROR MODEL TYPE

		#_loop over observations
		for obs_rec in aero_tmp:
		   p		= obs_rec.code
		   mode 	= obs_rec.mode
		   tau_obs 	= obs_rec.tau.copy() 
		   epoch_obs 	= obs_rec.epoch
	
		   #_Find nearest lat/lon and get fcst value.  
		   # Take mean of ensemble, if necesary??? 
		   # Or just use closest
		   lon, lat = obs_rec.lon, obs_rec.lat
		   i, j = lt.ll2ij( lat, lon, lats, lons )

		   #_Pull out forecast value to score	
		   tau_fcst = fcst[j,i]

		   #_Loop over scoring models
		   for err in error_models:
		     for ( label, thresh ) in thresholds.iteritems():
			
			#_Loop over thresholds bins (high/med/low/etc)

			#_For testing events within certain magnitudes
			#_The bins are based upon the model forecast for prog
			# against the aeronet for diagnostic
			if err == 'prognostic'	:	thresh_fcst = tau_fcst
			elif err == 'diagnostic':	thresh_fcst = tau_obs 
			else: 
				error = 'error model not defined '
				raise RuntimeError, error + str(err)

			#_check to see if within threshold
			if (thresh_fcst < thresh[0] or thresh_fcst > thresh[1]):
					continue 
			dbg(( 'LLIJ:', str(fcst.shape), str(lat), str(lats[j]),
				str(lon), str(lons[i])), l=3 )
			dbg(( 'FCST:', rec.dtg_vald, p, str(tau_fcst)), l=3 )
			dbg(( 'OBSV:', rec.dtg_vald, p, str(tau_obs)), l=3 )
			dbg(( 'THSH:', rec.dtg_vald, p, str(thresh_fcst)), l=3 )
			dbg(( 'DTGV:', rec.dtg_vald, lt.epoch2dtg(epoch_obs)),
				l=3 )
			dbg( '____', l=3 )

	                #_calculate residual, squared
	                residual 	= tau_fcst - tau_obs
			gross 		= tau_fcst + tau_obs
                	res_nrml 	= residual / tau_obs 
	
			dbug = ( rec.dtg_vald, epoch_obs, model, tau_fcst, 
				p, tau_obs, residual )
			f_dbug.write( str(dbug) + '\n' )

	                #_calculate individual ensember member error
			if score_members:
			    for mem in members:
				e = members.index( mem ) 
				
				#_skip missing data
				if hasattr( full, 'mask' ):
					#_skip masked datapoints. 
					# this catches spots where things can
					# be expected to be masked, such as
					# NGAC not having smoke, or just loose
					# error reports in general
					try:
						if full.mask[e,j,i]:
							dbg(( 'masked at', 
							variable, mem, 
							dtg_vald, e, j, i,
							full.mask[e,j,i] ),l=3) 
							continue
					except: 
						dbg(( 'skipped at', variable, 
						mem, dtg_vald, e, j, i ), l=3) 
						continue 

        	                tau_fcst = full[e,j,i].copy()
				
				#_Calculate member comps	
        	                res_e = tau_fcst - tau_obs 
				grs_e = tau_fcst + tau_obs
				nrl_e = res_e / grs_e	

				#_rank member
				midx = full[:,j,i].argsort().tolist().index(e)
				rnks = np.zeros( nens )
				rnks[midx] += 1

				#_Add compoents to array
				comps.resize( comps.size + 1 )
				comps[-1] = ( p, lat, lon, epoch_obs, mode, 
					variable, mem, fhr, res_e, grs_e, 
					nrl_e, -9999., rnks, err, False, -9999.,
					label)
			else: dbg( 'not scoring ensemble members', l=3 )

			#_calculate ensemble only scoring (briers, ranks)
			if is_ens:
				#_Calculate binary score for HALF-BRIERS 
				#_Find how many members are within threshold
				abv = np.where( full[:,j,i].data > thresh[0])[0]
				bel = np.where( full[:,j,i].data < thresh[1])[0]
				n_hits = len( lt.intersection( [abv, bel] ) )

				#_Calculate probability of event
				pE = float( n_hits ) / nens
			
				#_Create observation binary
				if tau_obs > thresh[0] and tau_obs < thresh[1]:
					oE = 1	#_event occured 
				else: 
					oE = 0	#_Event did not occur

				#_Calculate residual for half-Brier's score
				brier_res =  pE - oE 

				#_RANK IN ENSEMBLE
				sample = np.append( tau_obs, full[:,j,i].copy())
				idx = sample.argsort().tolist().index(0)

				#_Incriment appropate rank
				ranks = np.zeros( sample.size )
				ranks[idx] += 1

			#_Add compoents to array
			comps.resize( comps.size + 1 )
			comps[-1] = ( p, lat, lon, epoch_obs, mode, variable,
				model, fhr, residual, gross, res_nrml, 
				brier_res, ranks, err, True, members, label )

	f_dbug.close()
	
	#_Use pipe if threaded
	if pipe == None:
		return comps 
	else:
		pipe.send( comps )
		pipe.close()

def aeronet_score( components ):
        ''' 
	components	: class.aeronet_comps from stat_comps() 

	to score individual points or regions, subset aeronet_comps()
	externally
	'''
	dbg(( components.size, 'component(s)'), l=6 )

	#_Number of observations
	res = components.residual
	grs = components.gross
	nrl = components.nrml_res
	bri = components.briers_res 
	rnk = components.rank.sum()
        n = len( res )

	res = np.ma.masked_where( res == -9999., res )
	grs = np.ma.masked_where( grs == -9999., grs )
	nrl = np.ma.masked_where( nrl == -9999., nrl )
	bri = np.ma.masked_where( bri == -9999., bri )
	
	#_if completely missing, blargh
###	res = -9999. if res.mask.all() else res
###	grs = -9999. if grs.mask.all() else grs
###	nrl = -9999. if nrl.mask.all() else nrl
###	bri = -9999. if bri.mask.all() else bri
	#_______________________________________________________________________
	bias = np.sum( res ) / n #_This just mean absolute error... sort it 

	#_______________________________________________________________________
        #_ROOT MEAN SQUARE ERROR
        total = np.sum( np.array( res )**2 )
        rmse = np.sqrt( total / n )

	#_______________________________________________________________________
	#_RMSD
	total = np.sum ( ( np.array( res ) - bias )**2 )
	rmsd = np.sqrt( total / n )
	
	#_______________________________________________________________________
	#_CALC BRIERS HALF-SCORE
	tmp = np.array( bri )**2
	briers = np.sum( tmp ) / n

	#_______________________________________________________________________
	#_MEAN ABSOLUTE ERROR
	total = np.sum( np.fabs( res ) )	
	mae = total / n

	#_______________________________________________________________________
	#_NORMALIZED ROOT MEAN SQUARE ERROR
        total = np.sum( np.array( nrl )**2 )
        nrmse = np.sqrt( total / n )

	#_______________________________________________________________________
	#_FRACTIONAL GROSS ERROR
	total = np.sum( np.fabs( res / grs ) )
	fgs = 2*total/n

	#_Put stats into object to pass back to calling script
	stats = cl.aeronet_stats((n, bias, rmse, rmsd, nrmse, mae, fgs, briers, rnk))

	return stats 

def filter_dtgs( dtgs, aeronet=None, interpolation='neighbor', tolerance_dt=3,
	**kwargs ):
	'''
	dtgs		: list, 	List of date-time-groups
	aeronet		: cl.lidar()	recarray of aeronet observations
	tolerance_dt 	: int		distance from vald time of dtgs to allow
	interpolation	: string,	neighbor or median	

	Given the dtgs passed, filter will return records of nearest
	aeronet observations within tolerance_dt hours.

	No checking is done that dt of desired obs times is greater than 
	the window, so tolerance == 24, dt = 3 will result in the same obs
	being used for median interpolation.	
	'''
	import libnva as ln
	from scipy import stats
	median = stats.mstats.mquantiles
	dbg( dtgs, l=3 )

	#_If not passed, read in aeronet
	if aeronet == None: aeronet = read_aeronet( dtgs, **kwargs )

	#_Get list of codes present
	codes = lt.unique( aeronet.code )
	names = lt.unique( aeronet.long_name )
	modes = lt.unique( aeronet.mode ) 
	aero_out = cl.lidar()

	#_Loop over dtgs
	#_Loop over aeronet sites
###	for code in codes:
	for name in names:
	    for mode in modes:
		aero_sub = ln.subset( aeronet, long_name=name, mode=mode )
	
		#_Loop over valid times
		for dtg in dtgs:
		    epoch = lt.dtg2epoch( dtg )

		    #_return nearest observation within window__________________
		    if interpolation == 'neighbor':

			#_Find smallest difference and index of that record
			idx = np.argmin( np.fabs( aero_sub.epoch - epoch ) )
			min_diff = np.fabs( aero_sub.epoch[ idx ] - epoch )

			#_Add to output if within window
			if min_diff/3600. < tolerance_dt:
				dtg_fail = lt.epoch2dtg( aero_sub.epoch[idx] )
				dbg(( dtg_fail, dtg ), l=3 )

				aero_sub[idx].dtg = dtg
				aero_out = lt.merge(( aero_out, aero_sub[idx] ))
			
		    #_return median of values within window_____________________
		    elif interpolation == 'median':
			#_get start and end times of window
			epoch0 = epoch - tolerance_dt*3600. 	
			epoch1 = epoch + tolerance_dt*3600.	

			#_create mask of points outside window
			idx0 = np.where( aero_sub.epoch >= epoch0 )[0]
			idx1 = np.where( aero_sub.epoch < epoch1 )[0]
			idx = lt.intersection([idx0,idx1])

			#_take median within window
			taus = aero_sub.tau.copy()[idx]
			tau = median( taus, 0.5 )

			#_if there are no data, it will be a masked array	
			if len(idx) > 0:
				name = lt.unique( aero_sub.long_name,
					unique=True )
				lat = lt.unique( aero_sub.lat )[0]
				lon = lt.unique( aero_sub.lon )[0]
				wvl = lt.unique( aero_sub.wavelength,
					unique=True )
				code = name2code(name)

				#_expand output array, add new values
				aero_out.resize( aero_out.size + 1 )
				aero_out[-1] = ( tau, epoch, dtg, code, 
					name, name, wvl, mode, lat, lon )
		    else:
			raise ValueError, 'invalid option for interpolation'
				
	#_Return filtered data
	return aero_out

def test_map(  dir_map='./', show=True ):
	import matplotlib.pyplot as plt
	import libnva as ln
	import time
        lt.mkdir_p( dir_map )
        file_out= dir_map + '/map.png'
        m = ln.draw_map( [-80,81,-180,181], [-90,90,-180,180], 'cyl', 20 )

        #_Setup plot header and footer text
        plt.title('test',va='bottom',ha='left',position=(0.0,1.0),size='medium')
        n       = time.gmtime(time.time())
        cdtg    = str(n.tm_year)                \
                + str(n.tm_mon).zfill(2)        \
                + str(n.tm_mday).zfill(2)       \
                + str(n.tm_hour).zfill(2)

	x,y     = m(0,20)
	m.scatter( x, y, color='r', s=350*111, edgecolors='none' )
	if show:
		plt.show()
	else:
		dbg(file_out)
		plt.savefig( file_out, dpi=(120) )
		lt.make_readable( file_out )
        plt.close()
	
def create_meanmed( dtg_range, finc=6, write_period=False, **kwargs ):
	'''
	write_period	: bool,		writes raw aeronet data to intermediary
					files
	'''
	import libnva as ln
	from netCDF4 import Dataset
	from scipy import stats

	path = kwargs.get('path')
	dtg_range = sorted( dtg_range )
	aeronet = read_aeronet_available( dtg_range, **kwargs )

	#_reduce to every six hours
	dtg_filt = []
	dtg_str, dtg_end = dtg_range
	dtg = dtg_str
	while dtg < dtg_end:
		dtg_filt.append( dtg ) 
		dtg = lt.newdtg( dtg, finc )

	#_limit the dtgs to closest to model valid times
	aeronet = filter_dtgs( dtg_filt, aeronet=aeronet, **kwargs )
	codes = lt.unique( aeronet.code )
	modes = lt.unique( aeronet.mode )

	#_make recarray with one record per present site
	# only used to find the nearest aeronet location
	aeroRep = representativeRecords( aeronet )
	
	#_for each site, get mean value and size of circle
	dtype = [	('code','a6'),
			('mode','a6'),
			('lat','f4'),
			('lon','f4'),
			('mu','f4'),
			('median','f4'),
			('geomu','f4'),
			('circle','f4'),
			('dtg_str', 'a10'),
			('dtg_end', 'a10'),
			('nobs','i4') ]
	circles = np.recarray((0,),dtype=dtype)

	tmp = np.array( lt.dtg2epoch( lt.unique( aeronet.dtg ) ))
	min_dt = int( np.min( tmp[1:] - tmp[:-1] ))
	dbg(('min', min_dt))

	epoch_start, epoch_end = lt.dtg2epoch(dtg_str), lt.dtg2epoch(dtg_end)
	nt = (epoch_end - epoch_start)/min_dt + 1
	for code in codes:
	    aero_code = ln.subset( aeronet, code=code, wavelength=550 )
	    if write_period:
		lt.mkdir_p( path )
		file = '-'.join(( code, dtg_str, dtg_end+'.nc' )) 
		file = '/'.join(( path, file ))
		dbg(file)	
		ncdf = Dataset( file, mode='w', format='NETCDF3_CLASSIC' )

		ncdf.createDimension('time', nt)
		cdf = ncdf.createVariable('time','f8',('time'))
		cdf[:] = np.arange( nt )*min_dt + epoch_start

		ncdf.epoch_start 	= epoch_start
		ncdf.epoch_end 		= epoch_end
		ncdf.dtg_start 		= dtg_str
		ncdf.dtg_end 		= dtg_end
		ncdf.code 		= code
		ncdf.long_name		= aero_code[0].long_name 
		ncdf.modes		= ','.join(modes)
		ncdf.lat		= aero_code[0].lat
		ncdf.lon		= aero_code[0].lon

	    for mode in modes:
		aero = ln.subset( aero_code, mode=mode, wavelength=550)
		aero = sort_rec( aero, order=('epoch',) )
		lat, lon = aero[0].lat, aero[0].lon

		#_find closest aeronet site
		sites, dists = closestAeronet( lat, lon, aeroRep, degree=2 )

		#_find mean tau value for period
		mu 	= np.mean( aero.tau )
		median 	= np.median( aero.tau )
		geomu 	= stats.mstats.gmean( aero.tau )
		obs 	= aero.size

		#_add to circles recarray
		circles.resize( circles.size + 1 )
		circles[-1] = ( code, mode, lat, lon, mu, median, geomu,
			dists[1], dtg_str, dtg_end, obs )

		#_write the raw data, including times when not available
		if write_period:
			#_put values into array the length of the period
			values = np.ones(nt) * -9999.
			idx = (lt.dtg2epoch(aero.dtg)-epoch_start)/min_dt
			values[idx] = aero.tau

			#_create ncdf variable and put data into file
			cdf = ncdf.createVariable(mode,'f4',('time'), 
				fill_value=-9999.)
			cdf[:] = values

	    if write_period: ncdf.close()
	
	return circles

def write_aeronet_meanmed( records, path='.', label=None, **kwargs ):
	'''
	write aeronet mean as read by readAeronetMean()
	records = np.recarray(), dtype = 
	'''
	from netCDF4 import Dataset
	import libnva as ln

	#_get record dimensions
	codes = lt.unique( records.code )
	modes = lt.unique( records.mode )
	dtg_str = lt.unique( records.dtg_str )[0]
	dtg_end = lt.unique( records.dtg_end )[-1]

	nsites = len( codes )

	#_generate file name and open output netcdf
	file = '-'.join(('aeronet_meanmed',dtg_str,dtg_end+'.nc'))
	file = '/'.join(( path, file ))
	dbg(file)
	ncdf = Dataset( file, mode='w', format='NETCDF3_CLASSIC' )
	ncdf.createDimension( 'code', nsites )
	for mode in modes:
		recs = ln.subset( records, mode=mode )
		cdf = ncdf.createVariable( mode+'_mean', 'f4', ('code',) )
		cdf[:] = recs.mu
		cdf = ncdf.createVariable( mode+'_median', 'f4', ('code',) )
		cdf[:] = recs.median
		cdf = ncdf.createVariable( mode+'_geomu', 'f4', ('code',) )
		cdf[:] = recs.geomu
	cdf = ncdf.createVariable( 'lat', 'f4', ('code',) )
	cdf[:] = recs.lat
	cdf = ncdf.createVariable( 'lon', 'f4', ('code',) )
	cdf[:] = recs.lon
	cdf = ncdf.createVariable( 'distance', 'f4', ('code',) )
	cdf[:] = recs.circle
	cdf.long_name = 'distance to nearest aeronet site'
	cdf.units = 'km'
	cdf = ncdf.createVariable( 'nobs', 'f4', ('code',) )
	cdf[:] = recs.nobs
	cdf.long_name = 'number of observations'
	ncdf.codes = ','.join( codes )
	ncdf.modes = ','.join( modes )
	ncdf.dtg_str = dtg_str
	ncdf.dtg_end = dtg_end
	ncdf.close()

def read_aeronet_meanmed( file, **kwargs ):
	''' read mean aeronet values into recarray, previously calculated '''
	from netCDF4 import Dataset
	dbg(file)
	ncdf = Dataset( file, mode='r', format='NETCDF3_CLASSIC' )

	codes = ncdf.codes.split(',')
	modes = ncdf.modes.split(',')

	#_final sanity check that the code dimension is correct
	if len(ncdf.dimensions['code']) != len(codes):
		raise RuntimeError, 'error in write_aeronet_mean, ncode wrong'

	try:
		dtg_str = ncdf.dtg_str
		dtg_end = ncdf.dtg_end
	except AttributeError:
		import re
		res = re.search('(\d{10})-(\d{10}).nc', file)
		dtg_str = res.group(1)
		dtg_end = res.group(0)

	#_initialize recarray
	dtype = [	('code','a6'),
			('mode','a6'),
			('lat','f4'),
			('lon','f4'),
			('mu','f4'),
			('median','f4'),
			('geomu','f4'),
			('circle','f4'),
			('dtg_str', 'a10'),
			('dtg_end', 'a10'),
			('nobs','i4') ]
	circles = np.recarray((0,),dtype=dtype)

	#_iterate over sites and read data into recarray
	for code in codes:
		#_fine site's array index
		idx = codes.index(code)

		#_pull general metadata
		lat = ncdf.variables['lat'][idx]
		lon = ncdf.variables['lon'][idx]
		nobs = ncdf.variables['nobs'][idx]
		km = ncdf.variables['distance'][idx]

		for mode in modes: 
			#_put tau values in variable
			mean = ncdf.variables[mode+'_mean'][idx]
			median = ncdf.variables[mode+'_median'][idx]
			geomu = ncdf.variables[mode+'_geomu'][idx]

			#_put data into recarray
			circles.resize( circles.size + 1 )
			circles[-1] = ( code, mode, lat, lon, mean, median,
				geomu, km, dtg_str, dtg_end, nobs )	
	ncdf.close()

	return circles

def aeronet_meanmed_map( circles, path_map='.', file='map.png', show=False,
	title=None, metric='median', **kwargs ):
	'''
	plot map over avarage aeronet values over prescribed period
	circles	: np.recarray,	information on average tau, location and size
				of scatter circle	
	'''
	import matplotlib.pyplot as plt
	import libnva as ln
	import time

        dellat = 20 # 90
        dellon = 20 #45
        labx = [0,0,0,0]
        laby = [1,0,0,0]
        grid = [-90, 91, -180, 181]
        grid = [-80, 81, -180, 181]

	modes = lt.unique(circles.mode)
###	fig = plt.figure()

        lt.mkdir_p( path_map )
###	file_out = '/'.join(( path_map, file ))

	ax = {}
	nx, ny = 360, 180
###	pidx = 0
	for mode in modes:
		file_out = '/'.join(( path_map, mode+ '-' + file ))
###		ax[pidx] = fig.add_subplot( 1, 3, pidx+1 )
        	m = ln.draw_map( grid, [-90,90,-180,180], 
			'cyl', 1, dellat=dellat, dellon=dellon, labx=labx,
			laby=laby )

	        #_Setup plot header and footer text
		head = ' '.join(( title, mode ))
   		plt.title( head ,va='bottom', ha='left',
			position=(0.0,1.0), size='xx-small' )
###   		ax[pidx].set_title( title ,va='bottom', ha='left',
###			position=(0.0,1.0), size='xx-small' )

        	n       = time.gmtime(time.time())
        	cdtg    = str(n.tm_year)                \
        	        + str(n.tm_mon).zfill(2)        \
        	        + str(n.tm_mday).zfill(2)       \
        	        + str(n.tm_hour).zfill(2)
	
		rgb = lt.rgbgen('aod')
		lev = [0.01, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 9]

		circ_mode = ln.subset( circles, mode=mode )	
		for rec in circ_mode: #circles:
			x,y = m(rec.lon,rec.lat)
			#_size = rec.circle / 2. / 5.
###			size = 2.		
			size = 20.		
	
			for i in xrange( len(rgb) ):
				if rec.tau > lev[i]: color = rgb[i]
	 
###			ax[pidx].scatter( x, y, color=color, s=size, 
###				edgecolors='none')
			m.scatter( x, y, color=color, s=size, edgecolors='none')

		#_make fake plot for colorbar. yes, this is silly.	
		y2d=np.ones((ny,nx))*np.linspace( -89.5, 89.5,ny).reshape(ny,1)
		x2d=np.ones((ny,nx))*np.linspace(-179.5,179.5,nx).reshape(1,nx)
		d = np.zeros((ny,nx)) 
		cb = m.contourf( x2d, y2d, d, levels=lev, colors=rgb ) 
		bar = plt.colorbar( cb, orientation='horizontal', aspect=40,
			pad=0.05, shrink=0.9 ) 
###		[ t.set_fontsize(4) for t in bar.ax.get_xticklabels() ]

###		pidx += 1

		if not show:
			dbg(file_out)
			plt.savefig( file_out, dpi=(360) )
			lt.make_readable( file_out )
		else:
			plt.show()
	        plt.close()

def aero_geostd( values ):
	''' calculates geometric standard deviation of array given values '''
	from scipy import stats
	import numpy as np
	#_calculate geometric mean
	geomu = stats.mstats.gmean( values )
	dbg( geomu )

	#_mask missing and fill with geomu value so that it won't affect sum
	data = np.ma.masked_where( values == 0, values.copy() )

	#_calculate number of data points available
	miss = np.ones( data.shape )
	if data.mask.any() != False: miss[data.mask] = 0
	n = miss.sum()

	#_fill in masked values to not affect summation
#	data = np.ma.MaskedArray.filled( data, geomu )

##	dbg(( data.max(), data.min(), n, miss.shape ))
	tmp_log = np.log(data/geomu)**2
	tmp_sum = np.sum( tmp_log/n )
	tmp_sqr = np.sqrt( tmp_sum )
	tmp_gsd = np.exp( tmp_sqr )
	
	geosd = np.exp( np.sqrt( np.sum( np.log( data/geomu )**2 )/n ))
	dbg(( tmp_log.max(), tmp_sum, tmp_sqr, n ))
	dbg(( geosd, tmp_gsd ))
	return geosd

def dbg( msg, l=1 ):
	import inspect
	msg = lt.to_string( msg )
        if hasattr( msg, '__iter__'): msg = ' '.join( msg )

        if debug >= l:
                curf = inspect.currentframe()
                calf = inspect.getouterframes( curf, 2 )
                file, line, method = calf[1][1:4]
		file = file.split('/')[-1]
                print '[%s.%s.%i] %s' % ( file, method, line, msg )
