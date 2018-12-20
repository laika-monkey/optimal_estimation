#############################################################################_80
# NAME		libtools.py														   #
# PURPOSE	This is the first stab at a Python general purpose toolbox		   #
#			for NVA and use within the aerosol group.						   #
# USAGE		Not a stand alone script.  Should only be imported, not run		   #
# AUTHOR	Walter R. Sessions, September 2011								   #
#############################################################################_80

debug = 1 #_Turns on messages sent to dbg 


def int_with_commas(n):
	''' returns integer as a string formatted with commas '''
	if type(n) not in [int]:
		raise TypeError, 'Only uses integers or strings to avoid roundoff err'

	import locale
	locale.setlocale(locale.LC_ALL, 'en_US')
	return locale.format('%d', n, grouping=True)


def newdtg(dtg_old, fhour):
	'''
	USAGE: 	OUTPUT 		= newdtg(<DTG_STRING>,<INT>)
		'2011040200' 	= newdtg('2011040100',24)
	'''
	dtg_old	= str(dtg_old)
	fhour = int(fhour)
	a = dtg2epoch(dtg_old) 
	a += 3600 * fhour 
	dtg_new	= epoch2dtg(a)
	dbg((dtg_old, fhour, dtg_new), l=9 )
	return dtg_new 


def epoch2iso(e):
	from datetime import datetime
	if hasattr(e, '__iter__'):
		return [datetime.utcfromtimestamp(i).isoformat() for i in e]
	else:
		return datetime.utcfromtimestamp(e).isoformat()


def epoch2local(e):
	from time import localtime, mktime 
	from datetime import datetime

	#_create initial time tuple
	str_time = localtime(e)

	#_convert to seconds
	time_sec = mktime(str_time)
	
	#_convert to datetime
	dt = datetime.fromtimestamp(time_sec)

	return dt.isoformat()


def newdtg2(dtg_old, dt, inc=3600., **kwargs):
	''' 
	increment time forward

	dtg_old		str,	start date-time-group
	dt			int,	number of increments forward
	inc			int,	size of increments in seconds, defaults to hours 

	'''
	epoch = dtg2epoch(dtg_old, **kwargs)
	epoch += dt * inc
	return epoch2dtg(epoch, **kwargs) 


def julian2dtg(year, jday, full=False):
	''' convert julian to date-time-group '''
	import datetime as dt
	from numpy import array
	fmt = '%Y%m%d%H%M%S' if full else '%Y%m%d%H'

	if hasattr(jday, '__iter__'):
		d = []
		for j in jday:
			dtg = dt.datetime(int(year), 1, 1)
			dtg += dt.timedelta(float(j)-1)
			d.append(dtg.strftime(fmt))
		d = array(d)
	else:
		d = dt.datetime(int(year), 1, 1) + dt.timedelta(float(jday)-1)
		d = d.strftime(fmt) 
	return d 


def julian2epoch(year, jday):
	''' convert julian to date-time-group '''
	from calendar import timegm
	import datetime as dt
	from numpy import array
	if hasattr(jday, '__iter__'):
		d = []
		for j in jday:
			dtg = dt.datetime(int(year), 1, 1)
			dtg += dt.timedelta(float(j)-1)
			d.append(timegm(dtg.timetuple()))
		d = array(d)
	else:
		d = dt.datetime(int(year), 1, 1) + dt.timedelta(float(jday)-1)
		d = timegm(d.timetuple()) 
	return d


def dtg2iso(dtg, full=False):
	''' convert YYYYMMDDHH to YYYY-mm-ddTHH:MM:SS '''
	import datetime
	fmt = '%Y%m%d%H' if not full else '%Y%m%d%H%M%S'
	dt = datetime.datetime.strptime(dtg, fmt)
	return dt.isoformat()


def dtg2julian(dtg, full=False):
	''' Convert datetimegroup to julian day '''
	import datetime
	fmt = '%Y%m%d%H' if not full else '%Y%m%d%H%M%S'
	dt = datetime.datetime.strptime(dtg, fmt)
	tt = dt.timetuple()

	return tt.tm_yday+(tt.tm_hour*3600 + tt.tm_min*60 + tt.tm_sec) / 86400.


def find_runlength( dtgs, dtge ):
	'''Find time between dtgs in seconds'''
	import calendar
	import numpy as np
	
	ys,ms,ds,hs	= dtgsplit( dtgs, type=int )
	ye,me,de,he	= dtgsplit( dtge, type=int )
	start 		= calendar.timegm((ys,ms,ds,hs,0,0))
	end 		= calendar.timegm((ye,me,de,he,0,0))
	dbg(( dtgs, dtge, start-end ), l=9 )
	if start > end	: dbg( 'Start time after end.' )
	return 		np.abs(start-end)


def dtgsplit( dtg, type=str ):
	'''just a wasted module to return yyyy,mm,dd,hh'''
	dbg(( dtg, str(type)), l=9 )
	if type == str	: return dtg[:4], dtg[4:6], dtg[6:8], dtg[8:10]
	else		: return int( dtg[:4] ), int( dtg[4:6] ),	\
				int( dtg[6:8] ), int( dtg[8:10] )


def human_date( dtg ):
	'''
	produce human readable date string
	result	= human_date(yyyymmddhh)
	'''
	import calendar
	dbg( dtg, l=9 )
	weekdays=['Monday ','Tuesday ','Wednesday ','Thursday ',\
		'Friday ','Saturday ','Sunday ']
	months	= ['Null ','January ','February ','March ',	\
		'April ','May ','June ','July ','August ',	\
		'September ','October ','November ','December ']
	if hasattr( dtg, '__iter__' ):
		dates = []
		for day in dtg:
			day	= str( day )
			y 	= int( day[:4] )
			m 	= int( day[4:6] )
			d 	= int( day[6:8] )
			h 	= day[8:10]
    
			weekday	= weekdays[ calendar.weekday( y, m, d ) ]
			month	= months[m]
			date	= weekday+str(d)+' '+month+str(y)+' '+h+'UTC' 
			dates.append(date)
		return dates
	else:
		dtg_old	= str( dtg )
		y 	= int( dtg_old[:4] )
		m 	= int( dtg_old[4:6] )
		d 	= int( dtg_old[6:8] )
		h 	= dtg_old[8:10]
    
		weekday	= weekdays[ calendar.weekday( y, m, d ) ]
		month	= months[m]
		date	= weekday+str(d)+' '+month + str(y) + ' ' + h + 'UTC' 
		return date 


def to_string( a ):
	if hasattr( a, '__iter__' )	: return [ str(b) for b in a ]
	else						: return str( a )


def date_label( dtgs, icap=False ):
	'''
	produce human readable date string, abbreviated for tickmarks
	REPLACE THIS WITH CALENDAR/DATETIME TOOLS PLEASE
	'''
	months = [	'Null','Jan','Feb','Mar','Apr','May','Jun', \
			'Jul','Aug','Sep','Oct','Nov','Dec']
	dates = []
	for dtg in dtgs:
	  	dtg_old	= str(dtg)
	  	y = int(dtg_old[:4])
	  	m = int(dtg_old[4:6])
	  	d = int(dtg_old[6:8])
	  	h = dtg_old[8:10]
    
	  	month = months[m]
		if icap == True:
	  		date = str(y) +' '+ month 
		else:
	  		date = h + 'Z ' + str(d) + ' ' + month 
	  	dates.append(date)
	dbg((dtgs, dates), l=9 )
	return dates 


def epoch2dtg(s, full=False):
	'''
	convert from seconds since 1 jan 1970 to dtg.  loss of res 
	s	int	seconds
	full	bool	if true, includes minutes, seconds

	returns
	dtg	str	date-time-group of length 10|14 depending on full
	'''
	import time
	import numpy as np
	if not hasattr(s, '__iter__'):
		n = time.gmtime(float(s))
		dtg = '{0:04d}{1:02d}{2:02d}{3:02d}'	\
			.format(n.tm_year, n.tm_mon, n.tm_mday, n.tm_hour)

		if full:
			dtg += '{0:02d}{1:02d}'.format(n.tm_min, n.tm_sec)
	else:
		dtype = 'a10' if not full else 'a14'
		dtg = np.array((), dtype=dtype)
		for seconds in s:
			n = time.gmtime(float(seconds))
			d = '{0:04d}{1:02d}{2:02d}{3:02d}'	\
			.format(n.tm_year, n.tm_mon, n.tm_mday, n.tm_hour)
			if full:
				d += '{0:02d}{1:02d}'.format(n.tm_min, n.tm_sec)
			dtg = np.append(dtg, d)
	dbg(( s, dtg ), l=9)
	return dtg


def remote_copy(src, passarg=''):
	''' copy over via scp '''
	import re
	import os

	#_just assume no colons in file.  We aren't using WRF output.
	result = re.search('^(.*?):(.*)$', src)

	if not result:
		return False

	host = result.group(1)
	fnam = result.group(2)
	cmd = 'scp {1} {0} .'.format(src, passarg)
	print cmd
	print '\n\n'
	os.system(cmd)
	return (fnam.split('/')[-1], host)


def modis2epoch(s):
	''' atomic time, offset below. THIS IS APPROXIMATE '''
	s += dtg2epoch('1993010100')
	return s


def epoch2date( s ):
	''' convert from seconds since 1 jan 1970 to dtg.  loss of res '''
	import time
	import numpy as np
	if hasattr( s, '__iter__' ) == False:
		n = time.gmtime(float(s))
		date = '.'.join(( str(n.tm_year), str(n.tm_mon).zfill(2),
				str(n.tm_mday).zfill(2) ))
		date += '_' 
		date +=	':'.join(( str(n.tm_hour).zfill(2),
			str(n.tm_min).zfill(2), str(n.tm_sec).zfill(2) ))
	else:
		date = np.array((), dtype='a10' )
		for seconds in s:
			n = time.gmtime( float( seconds ))
			d = '.'.join(( str(n.tm_year), str(n.tm_mon).zfill(2),
				str(n.tm_mday).zfill(2) ))
			d += '_' 
			d += ':'.join(( str(n.tm_hour).zfill(2),
				str(n.tm_min).zfill(2), str(n.tm_sec).zfill(2)))
			date = np.append( date, d )
	return date 


def dtg2epoch(dtg, full=False):
	''' 
	convert from date-time-group to seconds since 1 jan 1970 
	dtg	str(10|14)	date time group in yyyymmddHH(MMSS) fmt
	full	bool		if true, include minutes and seconds in dtg
	'''
	import calendar, datetime
	import numpy as np
	fmt = '%Y%m%d%H' if not full else '%Y%m%d%H%M%S'

	if not hasattr(dtg, '__iter__'):
		dt = datetime.datetime.strptime(dtg, fmt)
		y = dt.year 
		m = dt.month
		d = dt.day
		h = dt.hour
		timegm = [dt.year, dt.month, dt.day, dt.hour, 0, 0]
		if full:
			timegm[-2:] = dt.minute, dt.second
		a = calendar.timegm(tuple(timegm))
	else:
		a = np.array((), dtype='i4')
		for d in dtg:
			dt = datetime.datetime.strptime( d, fmt )
			timegm = [dt.year, dt.month, dt.day, dt.hour, 0, 0]
			if full:
				timegm[-2:] = dt.minute, dt.second
			a = np.append(a, calendar.timegm(tuple(timegm)))
	return  a


def ecmwf_day2dtg( hours ):
	'''	
	USAGE:	<DTG_STRING>	= ecmwf_day2dtg(<INT>)
		Convert from ECMWF timestamp to DTG.  
		ECMWF uses hours since 1900
	'''
	import numpy as np
        seventy	= 613608
	if hasattr( hours, '__iter__' ) == False:
		epoch = ( hours-seventy )*3600 
		a = epoch2dtg( epoch )
	else:
		a = np.array((), dtype='a10')
		for h in hours:
			epoch = ( hours-seventy )*3600
			dtg = epoch2dtg( epoch )
			a = np.append( a, dtg )
	dbg(( hours, a ), l=9 )
        return a 

def ncep_day2dtg( hours ):
        seventy	= 876576 + 613608 #_(1800-1900, 1900-1970)
	epoch	= (hours-seventy)*3600 
	dtg_new	= epoch2dtg(epoch)
	dbg(( hours, dtg_new ), l=9 )
        return dtg_new

def gsfc_day2dtg( days ):
	'''
	USAGE:	<DTG_STRING>	= gsfc_day2dtg(<INT>)
		Convert from GSFC timestamp to DTG.  They use days since 1970.	
	'''
	import numpy as np
	seventy	= 719164.
	if hasattr( days, '__iter__' ) == False:
		epoch = round( (days-seventy)*86400, -2 )
		a = epoch2dtg( epoch )
	else:
		a = np.array((), dtype='a10')
		for d in days:
			epoch = round( ( d-seventy ) * 86400, -2 )
			dtg = epoch2dtg( epoch )
			a = np.append( a, dtg )
	dbg(( days, a ), l=9 )
	return a 

def gsfc_dtg2day( dtg ):
	'''
	Usage: <float, days since 1970> = gsfc_day2dtg(<str, dtg>)
	'''
	seventy = 719164.
	epoch = dtg2epoch( dtg )
	days = (epoch/86400.) + seventy
	dbg(( dtg, days ), l=9 )
	return days

def get_firstlast( yyyymm ):
	''' returns first and last dtgs for month '''
	yyyymm = str( yyyymm )
	dtg0 = yyyymm + '0100'
	dtg1 = newdtg( newdtg( dtg0, 32*24 )[:6] + '0100', -24 )
	dbg(( yyyymm, dtg0, dtg1 ), l=9 )
	return dtg0, dtg1

def strip_aod( string ):
        ''' removes "_aod" from strings... kludge for transition '''
        import re
        re_aod = re.compile( '_aod', re.IGNORECASE )
        out = re_aod.sub( '', string )
	dbg(( string, out ), l=9 )
	return out

def setup_groups(objects, nproc=1, **kwargs):
	'''
	USAGE:	GROUPS 			= setup_groups( <ITER_OBJ>, <INT> )
		[[1,2,3],[4,5,6]]	= setup_groups( range(6), 3 )

	Break iterable object into groups of NP size to be forked.
	objects	: iterable object 
	nproc	: number of processors, int
	'''
	n_obs	= len(objects)
	ng		= int(n_obs / nproc)			#_Number of full groups
	extra	= 1 if (n_obs % nproc) else 0	#_Check for residual 
	groups 	= ['none'] * (ng + extra)		#_Setup 2D list 
	for i in range(ng):
		groups[i] = objects[(nproc*i):nproc*(i+1)]	
	if extra != 0: 
		groups[ng] = objects[nproc*ng:]
	dbg((objects, groups), l=9)
	return groups

def intersection( a ):
	''' Pass a 1d list and it will return common elements '''
	#_ set(a).intersection(b)
	x = a.pop()
	y = a.pop()

	x = x if hasattr( x, '__iter__' ) else [x]
	y = y if hasattr( y, '__iter__' ) else [y]

	inter = list( set(x) & set(y) )
	if len(a) > 0:
		for i in xrange( len(a) ):
			inter = list( set(inter) & set(a[i]) )
	dbg(( a, inter ), l=9 )
	return inter

def mkdir_p( path ):
	'''
	USAGE: 	mkdir_p(<FULL_PATH_STRING>)
		Similar to using `mkdir -p /directory/to/make/`
	'''
	if path == '.': return 0
	import errno, os
	dbg( path, l=9 )
	def mkdir( p ):
		try: 
			os.makedirs( p )
			os.chmod( p, 0775 )
		except OSError as exc: # Python >2.5
	        	if exc.errno == errno.EEXIST: pass
	          	else: raise

	if hasattr( '__iter__', path ):
		for p in path: mkdir( p )
	else:
		mkdir( path )

		
def set_log_axis(ax, xloc=None, xlab=None):
     ''' setup spectral axis for heating and flux plots '''
     import matplotlib
	 
     #_explicit x ticks
     #for b in xloc:
    ##     if b >= 1:
  ##           xlab.append( str(int(b)) )
##        else: 
##             xlab.append(str(b)[1:])
     
     #_setup x scale and limits
     ax.set_xscale('log')
     ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
     ax.set_xticks(xloc)   
     ax.set_xticklabels(xlab)


def latlon_to_xy( map, lat, lon ):
	'''
	Not used, but is here for reference

	returns x/y coords for plotting
	'''
	import numpy as np
	dbg(( lat, lon, type(map)), l=9 )
	if len(lat.shape) == 1: return map(*np.meshgrid(lon,lat))
	else: return map(lon,lat)


def ll2ij( lat, lon, lats, lons ):
	'''
	lat	: float
	lon	: float
	lats	: 1d array of latitudes
	lons	: 1d array of longitudes
	Ghetto and often wrong way to convert between lat/lon and 
	i/j indices
	meshgrid methods also functional
	'''
	import numpy as np
	lats = np.array( lats )
	lons = np.array( lons )

	if hasattr( lat, '__iter__' ) and hasattr( lon, '__iter__' ):
		j, i = [], []
		ncoords = len(lat)
		if len(lon) != ncoords: raise RuntimeError, 'come on'

		for n in xrange(ncoords):
			jt = int(np.argmin( np.fabs( lats - lat[n] )))
			it = int(np.argmin( np.fabs( lons - lon[n] )))
			j.append(jt)
			i.append(it)
	else:
		j = int(np.argmin( np.fabs( lats - lat )))
		i = int(np.argmin( np.fabs( lons - lon )))
	dbg(( lat, lon, i, j ), l=9 )
	return i, j

def masked_sum0( a, b, axis=None ):
	''' there really must be a better option '''
	#_Make sure they are masekd arrays and replace masked with 0s
	import numpy as np
	a = np.ma.MaskedArray.filled( np.ma.masked_array( a ), 0 )
	b = np.ma.MaskedArray.filled( np.ma.masked_array( b ), 0 )
	dbg(( type(a), type(b) ), l=9 )
	return a + b

def masked_sum( a, b ):
	''' there really must be a better option '''
	#_Make sure they are masekd arrays and replace masked with 0s
	import numpy as np
	shape = a.shape
	if shape != b.shape: raise SizeError, 'incompatable '+shape+' '+b.shape
	shape = list( shape )
	shape.insert(0,1)
	data = np.append( a.reshape(shape), b.reshape(shape), axis=0 )
	data = np.ma.masked_array( data )

	mask0 = np.zeros((shape))
	if not hasattr(a,'mask'):
		a = np.ma.masked_array( a )
		mask0[:] = False
	else:
		if type( a.mask ) == np.bool_: mask0[:] = a.mask
		else: mask0[:] = a.mask	

	mask1 = np.zeros((shape))
	if not hasattr(b,'mask'):
		b = np.ma.masked_array( b )
		mask1[:] = False
	else:
		if type( b.mask ) == np.bool_: mask1[:] = b.mask
		else: mask1[:] = b.mask	

	mask = np.append( mask0, mask1,axis=0 )
	data.mask = mask
	return data.sum(axis=0) 

def geostd( a, gmu, axis=None ):
	import numpy as np
	a = np.ma.MaskedArray.filled( np.ma.masked_array( a ), gmu )
	return np.log( a/gmu )**2


def masked_prod( a, b, axis=None ):
	''' there really must be a better option '''
	#_make sure they are masekd arrays and replace masked with 1s
	import numpy as np
	a = np.ma.MaskedArray.filled( np.ma.masked_array( a ), 1 )
	b = np.ma.MaskedArray.filled( np.ma.masked_array( b ), 1 )
	dbg(( type(a), type(b) ), l=9 )
	return a * b


def round_thresh( a, thresh ):
	''' round to a specified threshold '''
	from numpy import empty, round
	if hasattr(a, '__iter__'):
		out = empty(0)
	    	for x in a:
			out = append(out, round(float(x) / thresh) * thresh )
	else:
		out = round(float(a) / thresh) * thresh  
	dbg((a, thresh, out), l=9)
	return out


def join_records():
	''' placeholder for eventually making time a coordinate '''
	#_Purpose is to go from ([nens], ny, nx) => (nt, [nens], ny, nx)

	#_Initialize model object
	#_Join values, setting time var
	#_Update dimname/dimsize, pass everything to new record
	#_Return updated record
	pass


def colors_spag( ncol ):
	'''produces rainbow colors for spaghetti plots'''
	from numpy import array, arange, append
##	codes = ['0000FF','00FF00','00FFFF','FF0000','FF00FF','FFFF00']
	codes = ['FF0000','FF9900','FFFF00','00FF00','0000FF','000066','CC00CC']
	spect = array(())
##	for i in arange(len(codes)-1):
	#_convert codes to hex and fill in integer steps between
	for i, c in enumerate(codes):
		try:
			steps = arange(int(c, 16), int(codes[i+1], 16))
			spect = append(spect, steps)
		except IndexError:
			pass

	spect   = spect[0:-1:spect.size/500]
	spect   = spect[0:spect.size:spect.size/ncol]
	colors  = []
	for n in arange(spect.size):
		colors.append('#'+str(hex(int(spect[n]))[2:].zfill(6)))
	dbg(( ncol, colors ), l=9 )
	return colors


def rgbgen( field, nens=0 ):
	''' Produces rgb triplets for contour plots '''
	import numpy as np
        if field == 'spag':
                rgb = colors_spag( nens )
                rgb.append('#FF0000')   #-Mean color
                rgb.append('#000000')   #-Deterministic color
        else:
                red =   {
                'aod' : [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
                'warn': [1.,1.],
                'corr': [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
                'cmyk': [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0],
                'musig':[1,1,1,1,1]     }
                blue =  {
                'aod' : [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                'warn': [0,0],
                'corr': [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                'cmyk': [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                'musig':[1,1,1,1,1]     }
                green = {
                'aod' : [0.0, 0.5, 1.0, 1.0, 0.5, 0.0, 0.0], 
                'warn': [0,0],
                'corr': [0.0, 0.5, 1.0, 1.0, 0.5, 0.0, 0.0],
                'cmyk': [1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0], 
                'musig':[1,1,1,1,1]     }
                r       = np.array( red[field] )
                b       = np.array( blue[field] )
                g       = np.array( green[field] )
                rgb     = np.ones(( r.size, 3 ))

                for i in np.arange( r.size ):
                        rgb[i] = np.ones(3) * [r[i],g[i],b[i]]
	dbg(( field, nens, rgb ), l=9 )
        return rgb


def cb_fontsize( size, figure ):
	''' 
	size	: int, size of desired font
	figure	: matplotlib.pyplot.colorbar() object
	'''
	dbg(( size, type(figure)), l=9 )
	[ t.set_fontsize( size ) for t in figure.ax.get_xticklabels() ]
 
def unique( list_of_items, unique=False ):
	''' creates list of unique items '''
	a = sorted( set ( list_of_items ))
	dbg(( list_of_items, unique ), l=9 )
	if not unique:			return a 
	elif unique and len(a) == 1:	return a[0]
	elif unique and len(a) > 1:	
		raise ValueError, 'Too many values present when single exp'
	else:
		dbg(list_of_items)
		raise ValueError, 'Well this shouldn\'t happen'

def regrid_field( x_nat, y_nat, data, x_new, y_new ):
	#from scipy.interpolate import interp2d
	from scipy.interpolate import RectBivariateSpline as interp 
	"""
	Function returns data on x_nat/y_nat grid at x,y locations.
	len(x) * len(y) == len(data)
	"""
	#_Some of the data is flipped and flopped, so allow for transposition
	try		: f = interp( x_nat, y_nat, data )
	except TypeError: f = interp( x_nat, y_nat, data.transpose() )
	help(f)
	#_Actual Regridding
	regrid	= f( x_new, y_new )
	dbg(( x_nat.shape, y_nat.shape, data.shape, 
		x_new.shape, y_new.shape, regrid.shape ), l=9 )

	return regrid

#############################################################################_80
######_FILE-IO_#################################################################
################################################################################

def binary_append( field, file, bs=False ):
	"""
	Appends fortran-like unformatted records to binary
	field	== dict { 
				'data' :  ,
				'size' : <bytes>,
				'type' : <i,f,S#>,
			}
	file	== filename to append
	bs == byteswap
	USAGE:	binary_append(<DICT_FIELD>,<FILE_STR>)
	"""
	import numpy as np
	dbg( file, l=5 )
	f	= open(file,'ab')
	type	= field['type']
	size	= field['size']
	if bs == False:
		data		= np.array( size, dtype='i' )
		data.tofile(f)	#_Leading unformated FORTRAN record desriptor
		data		= np.array( field['data'], dtype=type )
		data.tofile(f)	#_Data
		data		= np.array( size, dtype='i' )
		data.tofile(f)	#_Closing unformatted FORTRAN record descriptor
	elif bs == True:
		data		= np.array( size, dtype='i' ).byteswap()
		data.tofile(f)	#_Leading unformated FORTRAN record desriptor
		data		= np.array( field['data'],dtype=type).byteswap()
		data.tofile(f)	#_Data
		data		= np.array( size, dtype='i' ).byteswap()
		data.tofile(f)	#_Closing unformatted FORTRAN record descriptor
	f.close()
	dbg(( field, file ), l=9 ) 

def make_readable( file ):
	"""
	For some reason, on heck, output images are set to 600 instead of 644
	If this is the case, set to 644
	"""
	import os
	dbg( file, l=9 )
	os.chmod( file, 0664 )

def rec_read( f, dtype ):
	''' Read in FORTRAN unformatted write '''
	from numpy import fromfile
	sz0 = byte_order( fromfile(f, dtype='i4', count=1))[0]
	dat = byte_order( fromfile(f, dtype=dtype, count=1))[0]
	sz1 = byte_order( fromfile(f, dtype='i4', count=1))[0]

	''' SOME VERSIONS HAVE DTYPE MINUTES AT END, so are not i8 for first '''
	#_get byte size at start of record
##	sz0 = fromfile(f, dtype='i4', count=1)

	#_the nature of this test is flawed.  Please delete
	#_read in unformated FORTRAN data 
##	dat = fromfile(f, dtype=dtype, count=1)[0]

	#_get byte size at end of record, to verify
##	sz1 = fromfile(f, dtype='i4', count=1)
	if sz0 != sz1:
		raise MemoryError('Malformed record {0} {1}'.format(sz0, sz1))

	dbg((f, dtype), l=9)
	return dat

def byte_order( x ):
	''' Checks system endian-ness. Swaps order of binary data if little '''
	import sys
	if (sys.byteorder == 'little'): x = x.byteswap()
	return x

def merge( arrays ):
	''' Takes two recarrays, merges them together and keeps attributes '''
	from numpy import hstack	
	return arrays[0].__array_wrap__(hstack(arrays))

def merge2( arrays ):
	ar0, ar1 = arrays[0].copy(), arrays[1].copy()
	sz0, sz1 = ar0.size, ar1.size
	ar0.resize( sz0 + sz1 )
	ar0[-ar1.size:] = ar1.copy()
	return ar0

def dir_size( path='.' ):
	''' returns the size of path on disk '''
	import os
	total = 0
	for dirpath, dirnames, filenames in os.walk( path ):
		for f in filenames:
			fp = os.path.join( dirpath, f )
			if os.path.exists( fp ): total += os.path.getsize(fp)
	return total

def exit_pipe( return_object, pipe=None, attrv=[None] ):
	''' function to cleanly exit optional pipes '''
	if pipe == None: 
		return return_object 
	else:
		pipe.send(( return_object, list(attrv) ))
		pipe.close()    

def print_opts( **flags ):
	'''
	flags	: dict
	
	used in a lot of my standalone scripts to print 
	dictionaries standing in as namelists.
	'''
        fmt = '%25s | %-25s'
        for k in sorted( flags.keys() ):
                v = flags[k]
                if not hasattr( v, '__iter__' ): #type(v) != list:
                        print fmt % ( str(k), str(v) )
                else:
                        print fmt % ( str(k), '' )
                        if type( v ) == dict:
                                for item in v:
                                        print fmt % ( str(item), str(v[item]) )

                        elif type( v ) == list:
                                for item in v:
                                        print fmt % ( '', str(item) )
                        print fmt % ( '', '' )
        print ''


def shrink_ticks(ax, size=6):
	''' reduce font size of ticks '''   
##	[ t.label.set_fontsize(size) for t in ax.yaxis.get_major_ticks() ]
##	[ t.label.set_fontsize(size) for t in ax.xaxis.get_major_ticks() ]
	for t in ax.yaxis.get_major_ticks():
		t.label1.set_fontsize(size)
		if hasattr(t, 'label2'):
			t.label2.set_fontsize(size)

	for t in ax.xaxis.get_major_ticks():
		t.label1.set_fontsize(size)
		if hasattr(t, 'label2'):
			t.label2.set_fontsize(size)
     

def capture_stdout(cmd):
	#_WTF leave in for now, but where did this garbage come from?
##	from os import open
##	p = open(cmd, 'r')
##	p.readline()
	import subprocess as s
	return s.Popen(cmd.split(), stdout=s.PIPE).communicate()[0].split('\n')


def subset( d, unique=False, **kwargs ):
	'''
	Database selection tool from np.recarray()

	If multiple keywords are pass, it is taken as AND operation
	If list for single option is passed, OR operation (fhr=[0,12,15])
	'''
	options=locals()
	import numpy as np
	for descriptor, desired in kwargs.iteritems(): #options.iteritems():
			#_if desired is an array, convert
			if type(desired) == np.ndarray:
				desired = desired.tolist()

			#_skip special
			if descriptor in ['d',None,'unique','kwargs']: continue

			#_read in values of attributes from data
			values=d.__getattribute__(descriptor)

			#_convert values for comparison
			values=np.array(values)

			#_makes it easier to keep consistent
			indices=np.array((0))
			if type(desired) != list:
				indices=np.where(values==desired)[0]
			else:
				indices=[]
				for d_loop in desired:
					[indices.append(i) for i \
						in np.where(values==d_loop)[0]]
			d = d[indices]

	if unique and d.size==1: return d[0]
	elif unique and d.size==0:
		raise ValueError, 'No values present when single rec expected'
	elif unique and d.size>1:
		print unique, d.size
		raise ValueError, 'Too many values present when one rec exp'
	else: return d


def tighten_plot(axes, rounded=False, onlyx=False, onlyy=False):
	'''
	make plot area as small as possible
	axes,			matplotlib.pyplot.axes
	rounded,		bool,	round axes up
	'''
	xmin, xmax = axes.xaxis.get_data_interval()
	ymin, ymax = axes.yaxis.get_data_interval()
	if onlyx	: axes.set_xlim(xmin,xmax)
	elif onlyy	: axes.set_ylim(ymin,ymax)
	else		: axes.axis([xmin,xmax,ymin,ymax])
	if rounded: raise ValueError, 'rounded not implemented'

############################################################################_80
#_stats_#######################################################################

def eof( d0, d1, centered=True, idx=0, **kwargs ):
	'''
	compute regression?  eofs?
	d0, d1	np.matrix,	Datasets, properly aligned d0.shape[0] == d1.shape[0]
	centered,	bool,	Has mean already been removed from data?
	idx,		int,	state space index	
	'''
	import numpy as np
	import libnva as ln
	import  matplotlib.pyplot as plt
	
	''' I haven't removed the variance at each gridpoint so this doesn't make sense '''
	#_that is, the var in time is 0
	
	lat, lon = np.linspace(-90,90,181), np.linspace(-180,180,361)
	lat2d, lon2d = np.meshgrid( lat, lon )
	
	#_remove mean if not centered datasets
	if not centered:
		d0 -= d0.mean()
		d1 -= d1.mean()
	
	if idx:
		raise RuntimeError, 'Not implemented yet. Order them yourself for now.'
		#_make array of indices 0-ndim
		#_remove idx from array, put at front
		#_.transpose((order))
	
	#_reshape datasets and convert to matrixes? Matrices?
	shape = d0.shape[1:] #_???
	d0 = np.matrix( d0.reshape( d0.shape[0], np.prod(d0.shape[1:]) ))
	d1 = np.matrix( d1.reshape( d1.shape[0], np.prod(d1.shape[1:]) ))
	n, ndummy = d0.shape
	if d0.shape[0] != d1.shape[0]:
		raise RuntimeError, 'Matrices not alligned '+str((d0.shape, d1.shape))

	covar = d0.T * d1 / (n-1.)
	covar = covar.reshape(shape)
	
	m = ln.draw_map()
	m.contourf( lon2d, lat2d, covar, latlon=True )
	plt.colorbar()
	plt.show()
	
def vector_len(x,**kwargs):
	''' 
	return length of vector 
	np.linalg.norm(x) does the same operation
	'''
	
	import numpy as np
	if type(x) != np.matrixlib.defmatrix.matrix:
		raise RuntimeError, 'only use on matrix'
	x = x.reshape(x.size,1)
	return np.sqrt(x.T*x)
	
def vector_angle(x,y,radian=True,**kwargs):
	''' find angle between to vectors '''
	import numpy as np
	if type(x) != np.matrixlib.defmatrix.matrix:
		raise RuntimeError, 'only use on matrix'
		
	x = x.reshape(x.size,1)
	y = y.reshape(y.size,1)
	
	#_wilkes 9.14
	theta = np.arccos( x.T*y / np.linalg.norm(x) / np.linalg.norm(y) )
	if not radian: theta *= 180/np.pi
	
	return theta[0,0]
	
def vector_shadow(x,y):
	''' lenth of x in direction of y '''
	import numpy as np
	if type(x) != np.matrixlib.defmatrix.matrix:
		raise RuntimeError, 'only use on matrix'
		
	x = x.reshape(x.size,1)
	y = y.reshape(y.size,1)
	
	#_wilkes 9.16
	L = x.T*y / np.linalg.norm(y)
	
	return L[0,0]
	
def covar_matrix(x,y,**kwargs):
	'''
	np.cov() does this, do not use.
	Wilkes 9.30
	
	Computer the covariance matrix of datasets
	x, y	np.matrix,	data to covar, must have same sample space dimension
	centered,	bool,	was mean removed across sample space
	sample,		int,	dimension of sample space 
	
	Sample space is currently assume to be the first axis
	'''
	x0, x1 = x.shape
	y0, y1 = y.shape
	
	if x0 != y0:
		raise RuntimeError, 'sample space not aligned'+str((x.shape,y.shape))
		
	x = center_data(x, **kwargs)
	y = center_data(y, **kwargs)
		
	#_calculate covariance matrix
	return x.T*y / (x0-1.)
	
def correlation_matrix(x,**kwargs):
	'''
	calc correlation matrix 
	x	np.array/matrix,	axis=0, sample space; axis=1, state space
	
	not sure what this actually is
	
	Wilkes 9.31
	'''
	import numpy as np
	
##	x = center_data(x, **kwargs)
	
	#_get covariance matrix and a diagonal of only the state variances
	s = np.sqrt( np.cov(x) )
	d = np.matrix( np.identity( np.sqrt(s.size) ) * s )

	#_calculate corr matrix
	return d.I*s*d.I
	
def center_data(x, centered=True, sample_index=0, **kwargs):
	''' remove mean from data if not centered, assumes sample spaceaxis=0 '''
	if not centered: x -= x.mean(axis=sample_index)
	return x
	
def remove_climo( data, dt=12 ):
	''' 
	returns data with dt signal removed
	
	data	np.ndarray,	Contains data with time dimension in 0 axis
	dt		int,		number to skip in cycle (e.g., 12 for monthly)
	'''
	import numpy as np
	
	#_take mean along first axis at dt rate
	climo = np.array([ data[i::dt].mean(axis=0) for i in xrange(12) ])
	
	#_build array in which to repeat climatology for removal
	shape = tuple( [data.shape[0]/dt]+[1] * (data.ndim-1) )

	#_remove signal form data
	data -= np.tile( climo, shape )

	return data, climo
	
def detrend( data, fit='linear', order=1, breakpoint=False ):
	''' 
	remove trends from data vector.  If data.ndim > 1, detrend columns
	
	fit,	string,		'linear' or 'constant'   Linear removes a best fit
						line while constant removes the mean
	breakpoint,	int,	Not implemented.  Allow for trending regions to be
						removed instead of over whole set
	order		int,	power of fit
	
	'''
	import numpy as np
	
	if fit == 'linear':
		size = data.shape			#_store to return data to size
		
		#_if only one vector, add dummy dimension
		if data.ndim == 1	: column = np.array( data.reshape( data.size, 1 ))
		else				: column = np.array( data.copy() )
		
		#_loop over each column and remove trend
		for n in range( column.shape[1] ):
			y_o = column[:,n]			#_pull y values
			x_o = range( y_o.size )		#_generate x values

			fit = np.polyfit( x_o, y_o, order )	#_get first order fit params
			f_x = np.poly1d( fit )				#_use them to build a function
			column[:,n] = y_o - f_x(x_o)		#_remove from data

		column = column.reshape(size)	#_return to original shape
		
	elif fit == 'constant':
		pass
	else:
		raise RuntimeError, error, 'Undefined type of trend'
	
	return column

	
def autocorrelation(x, y=None, lag=10, **kwargs):
	'''
	x	np.ndarray	autocorrelate at lag n
	y	np.ndarray	crosscorrelate at lag n
	lag	int			interval between corr
	
	autocorrelations return one tail
	'''
	from numpy import array, corrcoef, arange
	from scipy.stats import pearsonr
	
	if y == None:
		return array( [1] + [pearsonr(x[:-i], x[i:])[0] 
			for i in arange(1,lag+1) ] )
			
	else:	#_perform cross correlation
		tmp = []
		for n in arange(lag*2+1)-lag:
			if n<0		: 
				tmp.append( pearsonr(x[-n:],y[:n])[0] )
			elif n==0	: 
				tmp.append( pearsonr(x,y)[0] )
			elif n>0	: 
				tmp.append( pearsonr(x[:-n],y[n:])[0] )
			else		: raise RuntimeError, 'What?'
		return array(tmp)		


def combination(arrays, out=None, dtype=None):
	'''
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> combination(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

	'''
	import itertools
	return  list(itertools.product(*arrays))
	'''
	import numpy as np
	arrays = [np.asarray(x) for x in arrays]

	n = np.prod([x.size for x in arrays])
	if out is None and dtype is None:
		out = np.tile(None, (n, len(arrays)))
	elif dtype is not None:
		out = np.zeros((n, len(arrays)), dtype=dtype)

	m = n / arrays[0].size
	out[:,0] = np.repeat(arrays[0], m)
	if arrays[1:]:
		combination(arrays[1:], out=out[0:m,1:])
		for j in xrange(1, arrays[0].size):
			out[j*m:(j+1)*m,1:] = out[0:m,1:]

	return out
	'''


def dummy_sbatch(cmd,
	job_name='dummy',
	sname='sbatch_{0}_{1}',
	unique_id='id-fill',
	time_sbatch='2:00:00',
	nice=False,
	**kwargs):
	''' write a basic sbatch to get something going '''
	from os import getpid, environ, path
	from time import gmtime, strftime

	time = strftime('%Y-%m-%dT%H-%M-%S', gmtime())
	pid = getpid()
	unique_id = '{0}_{1}_p{2}'.format(unique_id, time, pid)
	sname = sname.format(job_name, unique_id)
	
	out =  '#!/bin/bash\n'
	out += '#SBATCH --job-name={0}_{1}\n'.format(job_name, unique_id)
	out += '#SBATCH --partition=all\n'
	out += '#SBATCH --share\n'
	out += '#SBATCH --time={0}\n'.format(time_sbatch)
	out += '#SBATCH --ntasks=1\n'
	out += '#SBATCH --cpus-per-task=1\n'
	out += '#SBATCH --output=/odyssey/scratch/%u/logs/dummy_%A.txt\n'
	if nice:
		out += '#SBATCH --nice={0:d}\n'.format(nice)
	out += 'module purge\n'
	out += 'module load license_intel\n'
	out += 'module load impi\n'
	out += 'module load intel/15.0-2\n'
	out += 'module load hdf5/1.8.14\n'
	out += 'module load netcdf4/4.3.3\n'
	out += 'module load anaconda27/base_2015a_sharppy13\n'
##	out += 'export TMPDIR=${{SLURM_JOB_NAME}}.${{SLURM_JOB_ID}}\n'
##	out += 'mkdir -p $TMPDIR\n'
	out += '\n'.join(['export {0}={1}'.format(var, environ[var]) \
			for var in ['PYTHONPATH', 'PRODUCTS', 'WORK', 'PATH', 'LOG']])
	out += '\n'
	out += 'source activate my_root\n'
	out += 'echo `which python`\n'
	out += '{0}\n'.format(cmd)

	#_get unique identifiers
	sname = path.join(environ['WORK'], 'batch_scripts', sname)
	with open(sname, 'w') as f:
		f.write(out)

	print sname
	return sname


def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))


def strictly_decreasing(L):
    return all(x>y for x, y in zip(L, L[1:]))


def non_increasing(L):
    return all(x>=y for x, y in zip(L, L[1:]))


def non_decreasing(L):
    return all(x<=y for x, y in zip(L, L[1:]))


def open_url(url):
	'''remmmiiiinndder'''
	import urllib2
	return urllib2.urlopen(url).read()


def dbg( msg, l=1 ):
	''' if global debug is set to true, be more verbose '''	
	import inspect
	msg = to_string( msg )
	if hasattr( msg, '__iter__'): msg = ' '.join( msg )

	if debug >= l:
		curf = inspect.currentframe()
		calf = inspect.getouterframes( curf, 2 )
	#	file, method = calf[1][1], calf[1][3]
		file, line, method = calf[1][1:4]
		print '[%s.%s.%i] %s' % ( file, method, line, msg )

