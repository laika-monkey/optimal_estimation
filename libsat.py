#############################################################################_80
#_NAME: 	libsat.py						       #
#_PURPOSE:	Split this off of former libnrl.py to try to keep model        #
#		specific code for plotting and recarray manipulation here      #
#									       #
#_AUTHOR:	Walter R. Sessions, September 2011			       #
#############################################################################_80
##########_GENERAL_#############################################################
################################################################################
import libtools as lt
import libmeta 	as lm
import libclass as cl 
import numpy as np
import os, matplotlib
if 'DISPLAY' not in os.environ:
	matplotlib.use('Agg')	#_This is needed for plotting 
				# on systems with no active X service 
				# or through cron
import matplotlib.pyplot as plt

debug 		= 1 #_Turns on messages sent to dbg 
sen_dict	= lm.sensors()

################################################################################
######_FILE-IO_#################################################################
################################################################################


def read_airs_hdfeos(fname, fidx=None, **kwargs):
	'''
	2016.04.11, Walter Sessions
	Moving to IRIS, need to read in radiances for retrieval.

	fname	str,	full path to hdf-eos file

	Wait... I should be generating files and making new collocation deals
	'''
	from pyhdf import SD, SDC
	from numpy import recarray, ndarray

	#_
	hdf = SD(fname, SDC.READ)

	#_define record array
	dtype = [	('AIRS_radiance', ndarray),
				('AIRS_epoch', ndarray),
				('AIRS_latitude', ndarray),
				('AIRS_longitude', ndarray) ]
	airs = recarray((0,), dtype)
	setattr(airs, 'fname', fname)

	#_	
	idx = arange(hdf.select('Time')[:].size) if fidx is None else fidx	


def caliop2epoch(t):
	''' convert from yyyymmdd.FLOAT to epoch seconds '''
	
	def convert(n):
		from numpy import floor
		from libtools import dtg2epoch as d2e
		dtg = '20{0:6.0f}00'.format(floor(n))
		sec = (n - floor(n)) * 86400.
		return d2e(dtg) + sec
		
	if hasattr(t, '__iter__'):
		from numpy import array
		return  array([convert(n) for n in t])

	else:
		return convert(t)	



class AIRS_ancillary(object):
	'''
	return AIRS ancillary data
	
	Paolo provided L2.chan_prop.2003.11.19.v8.1.0.tobin.anc which may
	be reduncanty with AIRS_srftables_51118v4.hdf, but I have not checked
	as of May 2014.  Primarily loading this for RAD QUAL/L2_Ingore columns
	'''
	def __init__(self, nchan=2378, AIRS_ancfile=
		'/home/wsessions/data/L2.chan_prop.2003.11.19.v8.1.0.tobin.anc',
		**kwargs):
		#_from numpy import genfromtxt (not fixed ncol, so can't use)
		from numpy import tile, nan

		self.chanid			= tile(nan, nchan)
		self.l2_ignore		= tile(nan, nchan)
		self.srf_centroid	= tile(nan, nchan)

		#_read in data
		with open(AIRS_ancfile) as f:
			for n, line in enumerate(f.readlines()[123:]):
				cols = line.split()
				self.chanid[n]		= int(cols[0])
				self.srf_centroid[n]= float(cols[1])
				self.l2_ignore[n]   = int(cols[12])


def find_satfiles( dtg, sensor='modis', path=sen_dict['modis']['dir'] ):
	'''
	NOT UPDATED, DO NOT USE
	

	Find the most appropriate satellite file for a given dtg and sensor
	For missing files, loop backward 24 hours until one is found
	For forecasts, start the day prior to the analysis timestep

	NOTE: As of the lastest version of flambe, no longer go back.
	'''
	satellites = sen_dict[sensor]['sats']
	filelist = []

	#_loop over satellites, build file name, add to return list 
	for satellite in satellites:
		filename = '_'.join(( 'smoke', satellite, dtg+'00.dat' ))
		file = '/'.join(( path, satellite, dtg[:6], filename ))

		#_if the file exists, add to return list
		if os.path.exists( file ): 
			pass
		elif os.path.exists( file + '.Z' ):
			file += '.Z'
		elif os.path.exists( file + '.gz' ):
			file += '.gz'
		else:
			continue

		filelist.append( file )
		dbg( file, l=10 )

	return filelist

def goes_codes( dtg ):
	'''
	Returns goes east and west codes based upon dtg
	dtg	: str*10,	date time group of interest

	Just having it loop instead.... if they're there, they're there
	'''
	epoch = lt.dtg2epoch( dtg )

	if epoch < lt.dtg2epoch('2003021300'):	#_2000.09.01 to 2003.02.13
                code_east = 'g8'
                code_west = 'g10'
	elif epoch < lt.dtg2epoch('2006062600'):#_2003.02.13 to 2006.06.26
		code_east = 'g12'
		code_west = 'g10'
	elif epoch < lt.dtg2epoch('2010040100'):#_2006.06.26 to 2010.04.01
		code_east = 'g12'
                code_west = 'g11'
        elif epoch < lt.dtg2epoch('2011120600'):#_2010.04.01 to 2011.12.06
                code_east = 'g13'
                code_west = 'g11'
        elif epoch < lt.dtg2epoch('2012092400'):#_2011.12.06 to 2012.09.24
		code_east = 'g13'
		code_west = 'g15'
        elif epoch < lt.dtg2epoch("2012102300"):#_2012.09.24 to 2012.10.23
		code_east = 'g14'
		code_west = 'g15'
        else: #epoch < lt.dtg2epoch('2030120600') #_2012.10.23 to ???
		code_east = 'g13'
		code_west = 'g15'

	return code_west, code_east

def read_obsnew( dtg, 
	dir='/shared/aerosol_ops2/modis_obsnew_1d6h/landocean/terraaqua/',
	skip=1, name='MODIS', **kwargs ):
	"""
	UPDATE FOR RECORDS
	Read obs as of 2011.06.08
	Are on half degree marks (0.5, 1.5, etc)

	This module returns a big-ole-grid full of missing values
	If you want a single point, use get_modis_value

	MODIS column information:
	> dtg (i-14) YYYYMMDDHH
	> lon (f10.3) Longitude (-180 to 180)
	> lat (f10.3) Latitude (-90 to 90)
	> tau (f10.4) AOD
	> tausd (f10.4) AOD error standard deviation
	The over-land data files include the following additional items:
	> taustdL2 (f10.4) Standard deviation of L2 AOD within grid cell
	> obs_err (f10.4) Observation error of MODIS AOD
	> n (i4) Number of MODIS L2 retrievals in grid cell

	This was written when first starting with Python, so there
	are a number of things done in 'non-pythonic' ways

	"""
	import numpy as np

	#_ DON"T USE GENFROM TEXT! Ed's multiple ways of storing this crap
	# makes it too hard to generalize

	if name == 'MODIS':
		pre = '/'.join(( dir, dtg[:6], dtg + '_obsnew' ))
	elif name == 'NPP VIIRS':
		pre = '/'.join(( dir, dtg[:4], dtg[:6], dtg + '_obsnew' ))

	if os.path.exists(pre + '.txt'):
		file = pre + '.txt'
	elif os.path.exists(pre):
		file = pre 
	else:
		raise OSError, str(pre) + ' does not exist.'	
	dbg(file)

	#_read data into array
	hand = np.genfromtxt( file, dtype=None, skiprows=skip, comments=name )
	nobs = hand.size
	ncol = len(hand.dtype)

###	#_subset by file structure
###	fields_all = ['dtg','lon','lat','tau','tausd',
###			'tausdL2','obs_err','grid_n',
###			'1','2','3','4','5','6']
###	fields = fields_all[:ncol]

	#_put into recarray
	sat = cl.sat_object()
	sat.resize(nobs)
	sat.dtg[:]	= hand['f0']
	sat.lon 	= hand['f1']
	sat.lat 	= hand['f2']
	sat.tau[:]	= hand['f3']
	sat.tausd[:]	= hand['f4']
	sat.sensor[:]	= [name] * nobs
	sat.wavelength[:] = [550] * nobs
	sat.level[:] 	= [3] * nobs
	sat.units[:]	= [''] * nobs
	sat.long_name[:]= ['Aerosol Optical Depth at 550nm'] * nobs
	return sat 

def read_obsnew_loop( dtgs, dtge, finc=6, **kwargs ):
	taus = cl.sat_object()
	dtg = dtgs
	while dtg <= dtge:
		tau = read_obsnew( dtg, **kwargs )
		taus = lt.merge([taus, tau])
		dtg = lt.newdtg( dtg, finc )

	return taus

def plot_forker( recs, **kwargs ):
	'''	forks plots based on dtg	'''
	dtgs = lt.unique( recs.dtg )
	groups = lt.setup_groups( dtgs, **kwargs )

	#_loop over proc groups
	for group in groups:
		children = []
		for dtg in group:
			pid = os.fork()
			if pid != 0:
				children.append(pid)
			elif pid == 0:
				#_pull out day
				day = ln.subset( recs, dtg=dtg )

				#_pass day to plotter
				plot_sat( day, **kwargs )
				os._exit(0)

		for kid in children: os.waitpid(kid,0)

def plot_sat( recs, field='aod', path='.', **kwargs ):
	'''
	rec	: np.recarray, 

	plots all sensor data and differenc plots and whatever else ed wanted
	'''
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	from matplotlib.colors import LogNorm
	import libnva 	as ln
	import libcmap 	as lc
	reg_dict = lm.regions()
	plt_dict = lm.fields()


	#_pull out metadata for title
	dtg = lt.unique( recs.dtg, unique=True )
	sensor = lt.unique( recs.sensor, unique=True )
	long_name = lt.unique( recs.long_name, unique=True )
	wavelength = lt.unique( recs.wavelength, unique=True )
	level = lt.unique( recs.level, unique=True )
	title = ' '.join(( lt.human_date(dtg), sensor, long_name ))
	plt.title( title, size='small', ha='left', va='bottom',position=(0.,1.))

	#_plot specifics
	corn = reg_dict['global']['corn']
	grid = reg_dict['global']['grid']
	delt = reg_dict['global']['delta']

	#_setup colors
	icap = lc.rgbcmap( 'icap', N=256 )

	#_get specifics for field
	field = field.lower()
	levs = [0.1,0.2,0.4,0.8]
	norm = LogNorm( levs[0], levs[-1] )

	#_setup output name
	dir_out = '/'.join(( path, 'sat_01', 'global', dtg ))
	file = '_'.join(( dtg, 'total_aod', str(wavelength), sensor.lower(), 
		'l'+str(level)+'.png'))
	file = '/'.join(( dir_out, file ))
	lt.mkdir_p( dir_out )

	#_generate map to plot over
	m = ln.draw_map( grid, corn, 'cyl', delt, fill='white', **kwargs )

	#_put into 2d array
	tau2d = obsnew_field( recs, **kwargs )

	#_put data onto grid
	d, x, y = ln.map_data( tau2d, m )

	#_plot data
	CS = m.pcolormesh( x, y, d, cmap=cm.jet, norm=norm, rasterized=True )

	#_generate colorbar
	bar = plt.colorbar( CS, orientation='horizontal', shrink=0.8,
		aspect=40, pad=0.05, ticks=levs )
	bar.set_ticklabels( levs )

	dbg(file)
	plt.savefig(file)
	plt.close()
	
def obsnew_field( recs, nx=360, ny=180, **kwargs ):
	'''
	Takes recarray of obsnewdata and puts it into a 2d field for plot

	recs	: cl.sat_object
	
	Will crash if more than 1 dtg pass
	'''
	dtg = lt.unique( recs.dtg, unique=True )

	lat = np.linspace( -89.5, 89.5, ny )
	lon = np.linspace( -179.5, 179.5, nx )
	
	#_create dummy array to hold values
	tau = np.zeros((ny,nx)) - 9999.
	
	#_convert lats and lons to indices
	i, j = lt.ll2ij( recs.lat, recs.lon, lat, lon )

	tau[j,i] = recs.tau
	idx = np.where( tau == -9999. )
	tau[idx] = 0.
##	tau = np.ma.masked_where( tau == -9999., tau )
##	tau = np.ma.MaskedArray.filled( tau, 0. )
##	tau = np.linspace(0.1,0.8,nx*ny).reshape(ny,nx)
	tau = cl.var( tau, attrv=(lat,lon) )
	return tau

def error_msg( msg, pipe=None, **kwargs ):
        ''' 
        if global debug is set to true, be more verbose 
        msg     : str, Message to be printed
        l       : int, Debug level of message.  Set higher for lower level 
                        messages.  As debug increases, noisiness should also.
        '''
        import inspect
        msg = lt.to_string( msg )
        if hasattr( msg, '__iter__'): msg = ' '.join( msg )

        if pipe != None:
                pipe.send(False)
                pipe.close()

        curf = inspect.currentframe()
        calf = inspect.getouterframes( curf, 2 )
        file, line, method = calf[1][1:4]
        raise RuntimeError, '[%s.%s.%i] %s' % ( file, method, line, msg )

def dbg( msg, l=1 ):
        ''' 
        if global debug is set to true, be more verbose 
        msg     : str, Message to be printed
        l       : int, Debug level of message.  Set higher for lower level 
                        messages.  As debug increases, noisiness should also.
        '''
        msg = lt.to_string( msg )
        if hasattr( msg, '__iter__'): msg = ' '.join( msg )

        import inspect
        if debug >= l:
                curf = inspect.currentframe()
                calf = inspect.getouterframes( curf, 2 )
                file, line, method = calf[1][1:4]
                print '[%s.%s.%i] %s' % ( file, method, line, msg )

def plot_modis(data,corners=[-90,90,-180,180], grid=[-80,81,-180,181],
			delta=20, polar=0, foot='', title=None,
			levs=[.1,.2,.4,.8,1.6,3.2,6.4,10.], show=False, 
			rgb=None, fout='default.png', field='aod', cb=1):
	''' NOT UPDATED FOR RECORDS '''
	import matplotlib.pyplot as plt
	import time
	rgb	= rgbgen(field)
  	proj 	= 'cyl' if polar != 1 else 'npstere'
	m	= draw_map(grid,corners,proj,delta)
	n	= time.gmtime(time.time())
	d,x,y	= shift_data(data,m) 
	head_tail(title,foot,vert=polar)
	CS	= m.contourf(x,y,d,levels=levs,colors=rgb)
        if cb == 1:
 	  plt.colorbar(CS,orientation='horizontal',shrink=0.9,aspect=40,pad=0.05)
	if fout == False:
	  plt.show()
	else:
 	  print 'Writing '+fout
	  plt.savefig(fout,dpi=(120))
	  make_readable(fout)
	plt.close()
	
def satcode2name(code):
	code	= int(code)
	table	= { 	251	: 'GOES-11',
			252	: 'GOES-12',
			253	: 'GOES-13',
			254	: 'GOES-14',
			255	: 'GOES-15',
			499	: 'Meteosat-8',
			500	: 'Meteosat-9',
			600	: 'MTSAT-1R',
			601	: 'MTSAT-2',
			783	: 'MODIS-Terra',
			784	: 'MODIS-Aqua'		}
	try:
		satname	= table[code]
	except KeyError:
		satname	= 'unknown'
	return satname

def dbg( msg, l=1 ):
        ''' 
        if global debug is set to true, be more verbose 
        msg     : str, Message to be printed
        l       : int, Debug level of message.  Set higher for lower level 
                        messages.  As debug increases, noisiness should also.
        '''
        msg = lt.to_string( msg )
        if hasattr( msg, '__iter__'): msg = ' '.join( msg )

        import inspect
        if debug >= l:
                curf = inspect.currentframe()
                calf = inspect.getouterframes( curf, 2 )
                file, line, method = calf[1][1:4]
                file = file.split('/')[-1]
                print '[%s.%s.%i] %s' % ( file, method, line, msg )
