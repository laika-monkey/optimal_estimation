#############################################################################_80
#_NAME: 	libnva.py						       #
#_PURPOSE:	Split this off of former libnrl.py to try to keep model        #
#		specific code for plotting and recarray manipulation here      #
#									       #
#		Current focused on AOD, so the var class limits to an very     #
#		lenient set of potential physical values.  However, it does    #
#		limit.... so know your data and know what's going on in here.  #
#_USAGE:        Not a stand alone script.  Should only be imported, not run    #
#_AUTHOR:	Walter R. Sessions, September 2011			       #
#############################################################################_80
##########_GENERAL_#############################################################
################################################################################
''' 
TODO: 
	- Fix the read scripts to be smarter. Use dimensions and indexing to
	only read in the necessary values.
	- Remove all the extraneous 'write|read_period_whatevers that were
	one offs for paper

'''
import libtools		as lt
import libaeronet 	as la
import libmeta 		as lm
import libclass		as cl
import numpy as np
from numpy import ma
import os, matplotlib
if 'DISPLAY' not in os.environ: matplotlib.use('Agg')	
				#_This is needed for plotting 
				# on systems with no active X service 
				# or through cron
import matplotlib.pyplot as plt
import matplotlib.cm as cm

debug 		= 1 		#_Turns on messages sent to dbg 
dir_prod 	= lm.dir_prod 
dir_nva		= lm.dir_ops 
ifmt		= 'png' #_'gif'
sec_per_hour	= 3600
sec_per_10_year	= 86400 * 3650

#_Metadata used for models and plotting
mod_dict = lm.models()
reg_dict = lm.regions()
plt_dict = lm.fields()
pnt_dict = lm.aeronet_points()
sen_dict = lm.sensors()

#############################################################################_80
########_GENERAL_###############################################################
################################################################################

def subset( d, 	model=None, 
		code=None,		#_AERONET 6 char site name
		day=None,
		dtg_init=None, 
		dtg_vald=None, 
		dtg_str=None, 
		dtg=None, 
		dtgs=None, 
		dtge=None, 
		error_model=None,
		ensemble=None,
		epoch=None, 
		fhr=None, 
		fhr0=None, 
		fhr1=None, 
		label=None,
		location=None,
		long_name=None, 
		mean=None,
		mode=None, 
		nrl_name=None, 
		region=None,
		stdev=None,
		tau=None, 
		units=None, 
		variable=None, 
		wavelength=None,
		unique=False ):
	'''
	Creates a subset of data depending on desired keywords

	If multiple keywoards are passed, it is taken as an "AND" operation
	If a list of options for a single option is passed, it is an "OR" 
	operation

	Could make the argument a dictionary so I don't have to keep manually
	expanding what is filterable, but I like the failability.

	Potentially, could merge with sub_region.

	NEVER enable this with **kwargs.  Ever ever ever. It will not do what
	you think it is.
	'''
	options = locals()
	for descriptor, desired in options.iteritems():
		if descriptor=='d' or desired == None or descriptor == 'unique':
			continue
		#_Read in values of attribute from data
		values = d.__getattribute__( descriptor )

		#_Convert values for comparison [REMOVED - too damn slow ]
		values = np.array( values )

		#_Makes this easier to keep consistent
		indices = np.array((0))
		if type(desired) != list:
			indices = np.where( values == desired )[0]
		else:
			indices = [] 
			for d_loop in desired:
				[ indices.append(i) for i 	\
				in np.where( values == d_loop )[0] ]
		d = d[indices]

	if unique and d.size == 1:
		return d[0]
	elif unique and d.size == 0:
		raise ValueError, 'No values present when single rec expected'
	elif unique and d.size > 1:
		raise ValueError, 'Too many values present when single rec exp'
	else:
		return d 

def sub_region( data, corn=[-90,90,-180,180], region=None, label='custom' ):
	'''
	Make a new model_object(), subset by region
	Will update region name
	corners = list( llat, ulat, llon, ulon )
	'''
	if region != None:
		corn = reg_dict[ region ]['corn']
		label = region

	llat, ulat, llon, ulon = corn

	#_Initialize object to return
	data_sub = cl.model_object()	

	#_Loop over each record
	for record in data:
		if data.size == 1: record = data
		dtg_init = record.dtg_init
		dtg_vald = record.dtg_vald
		variable = record.variable
		ensemble = record.ensemble
		long_name = record.long_name
		units 	= record.units
		attrn 	= record.dimname
		region	= label
		vhr	= record.fhr
		model 	= record.model
		values 	= record.values
		dimname = record.dimname

		#_Shift data if it crosses the International Date Line
		if ulon > 180 or region == 'npolar':
			tmp_val, attrv = slide_grid( record, llon, reg=region )	
			values = cl.var( tmp_val, attrv=attrv )

		lats = values.lat
		lons = values.lon

		#_Get i, j values for corners
		li, lj = lt.ll2ij( llat, llon, lats, lons )
		ui, uj = lt.ll2ij( ulat, ulon, lats, lons ) 

		#_Sanity check, then store new lats/lons
		ui = np.min( [ui+1, len(lons)] )
		uj = np.min( [uj+1, len(lats)] )

		#_Setup attribute dimension values
		attrv = []
		for name in dimname:
			dim = values.__getattribute__( name )
			if name == 'lat': dim = lats[lj:uj]
			if name == 'lon': dim = lons[li:ui]
			attrv.append( dim )	

		#_Subset actual data
		if ensemble:
			values_reg = cl.var( values[:,lj:uj,li:ui], attrn=attrn,
					attrv=attrv )
		else:
			values_reg = cl.var( values[lj:uj, li:ui], attrv=attrv,
					attrn=attrn )

		#_Add record to returnable array		
		ndim = values_reg.ndim
		dimsize = values_reg.shape
		data_sub.resize( data_sub.size + 1 )
		data_sub[-1] = ( values_reg, model, dtg_init, dtg_vald, 
			vhr, label, variable, ensemble, dimname, dimsize,
			units, long_name )
 
	return data_sub 

def sub_point( aod, lon, lat ):
	'''
	Make a new model_class object, subset by point 
	aod = model_class() 
	'''
	raise RuntimeError, 'Do not use'

	lats	= aod.lat
	lons 	= aod.lon
	dtg 	= aod.dtg
	model 	= aod.model
	variable = aod.var
	
	i, j = lt.ll2ij( lat, lon, lats, lons )

	if aod.ndim == 2:
		meta = [ dtg, variable, [lat], [lon], model ]
		subset = cl.var( aod[i,j], meta=meta )
	if aod.ndim == 3:
		nens = aod.shape[2]
		meta = [ dtg, variable, [lat], [lon], model, 'members' ]
		subset = cl.var( aod[i,j,:], meta=meta )
	
	return subset

def error_msg( msg, pipe=None, **kwargs ):
	''' 
	if global debug is set to true, be more verbose 
	msg	: str, Message to be printed
	l	: int, Debug level of message.  Set higher for lower level 
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
	msg	: str, Message to be printed
	l	: int, Debug level of message.  Set higher for lower level 
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

def join_values( records, newattr='dummy', newattrv=None, order=None ):
        ''' 
        When passed recarrays, flatten array 
        mainly doing this for VALUES attr, work on rest later 

        This could be used to shove multiple forecast times into one
        array to output for ncdf.  If so, probably want to change order of
        dummy.  And allow for passing of attrv (for coordinates)

        newattr:        New dimensional coordinate name 
        newattrv:       New dimenstional coordinates
                        If not passed, integer array put in
        order:          Set the order of the axes to return.  Defaults
                        to ( (shape_in), new_axis ), eg.( [nens], ny, nx, ndum )

        If using to merge time dimenions, make sure input records are in
        temporal order and have same timestep between them.
        Don't pass dtgs reversed when one model is three hourly and the other
        hourly.
        '''
        import numpy as np
        return_array = None
        attrv, attrn, dummy = [], [], []

        #_To join these, value arrays need to be of the same shape
        variables = lt.unique( records.variable )
        shape = lt.unique( records.dimsize )
        dimname = lt.unique( records.dimname )[0]

        #_Return array (without dummy) if a single record
        if records.size == 1:
                attrn, attrv = get_attr( records[0] )
                return cl.var( records[0].values, attrv=attrv, attrn=attrn )

        if len( shape ) != 1:
                err = 'cannot join arrays of different size! ' + str(shape)
                raise ValueError( err )
        shape_in = list( shape[0] )

        #_Also don't allow the merging of different variables yet
        if len( variables ) != 1:
                err = 'cannot join arrays of different vars! ' + str(variables)
                raise ValueError( err )

        #_Add dummy index 
        ndum = records.size             #_size of new dimension
        shape_in.insert( 0, ndum )      #_shape in prepended with dummy

        #_Fill attribute lists for var class
        values = records.values
        order_out = []
        index = 1
        for name in records.dimname[0]:
                dim = values[0].__getattribute__( name )
                attrn.append( name )
                attrv.append( dim )

                order_out.append( index )
                index += 1

        #_Add new dim attribute name and values to attrn/attrv arrays
        attrn.insert( 0, newattr )
        if newattrv == None:
                attrv.insert( 0, np.arange( ndum ))
        else:
                attrv.insert( 0, np.array( newattrv ))

        #_Explicit dictation of output dimensinoal order
        order_out.insert( 0, 0 )#_Default: put new index on end for transpose
        if order != None:       #_Non-default
                ndim = len( attrv )
                attrn_tmp = [None] * ndim
                attrv_tmp = [None] * ndim
                for n in np.arange( ndim ):
                        i = order.index( n )
                        attrv_tmp[i] = attrv[n]
                        attrn_tmp[i] = attrn[n]
                order_out = order
                attrv = attrv_tmp
                attrn = attrn_tmp

        #_Stack values, change shape, put dummy on end
        #_NETCDF, we normally put dummy/time first.... change?
        return_array = np.vstack( values )
        return_array.shape = shape_in
        return_array = return_array.transpose( order_out )

        return cl.var( return_array, attrv=attrv, attrn=attrn )

def get_attr( rec ):
        ''' Returns list of rec.values dimenions and dimensional names '''
        import numpy as np
        if type( rec ) == np.recarray:
                raise ValueError, 'get_attr is for use on single records'

        attrn, attrv = [], []
        for name in rec.dimname:
                dim = rec.values.__getattribute__( name )
                attrn.append( name )
                attrv.append( np.array(dim) )
        return list( attrn ), list( attrv )


def sort_rec( records, order=('dtg_vald',) ):
        ''' 
        Returns records sorted by order
        
        order   : tuple, of strings describing record attributes
                        (set by dtype, generally)
                        IF A SINGLE ITEM, MAKE SURE TO HAVE A COMMA INCLUDED 
                        (1,)
        Passing more than one does not appear to work well right now.
        '''
        import numpy as np
        for column in order[::-1]:
                idx = np.argsort( records.__getattribute__( column ) )
                records = records[ idx ] 
        return records


#############################################################################_80
########_NVA-CORE_##############################################################
################################################################################

	
#############################################################################_80
########_PLOT_##################################################################
################################################################################

def plot_rec( rec, field='aod', water='white', res='c', path='.', show=False,
	file=None, dpi=120, **kwargs ):
	''' 
	USAGE: plot_rec( model_object() )

	rec	: model_object() [.size == 1]
	field	: string, field to be plotted.  'aod','spag','musig','warn'
	water	: string, color of water in plots
	res	: string, resolution of plots ('c'oarse,'l'ow, etc )
	path	: string, where to plot

	For a single record, plots prescribed field.
	When passed a model_object with .size > 1, will plot only the first
	Please refer to plot_forker() when wanting to plot multiple records

	'''
	import matplotlib.pyplot as plt
	import time, re
	#_Ensure we're receiving single record, or if not, use first
	if type( rec ) == np.core.records.recarray: 
		dbg( 'passed array - using first record.', l=2 )
		rec = rec[0]

	#_If data is ensemble, average 
  	variable	= rec.variable
	vhr 		= rec.fhr
	model 		= rec.model.upper()
	dtg_init	= rec.dtg_init
	dtg_vald 	= rec.dtg_vald
	values 		= rec.values
	dimnames	= rec.dimname
	dimsizes	= rec.dimsize
	label		= mod_dict[model]['label']

	#_Get regional information
	reg 	= rec.region
	corners = reg_dict[reg]['corn']
	grid 	= reg_dict[reg]['grid']
	delta 	= reg_dict[reg]['delta']

	#_Check if record is for ensemble data
	if rec.ensemble != 0:
		is_ens = True 
		ens_index = dimnames.index( 'member' )
		lat_index = dimnames.index( 'lat' )
		lon_index = dimnames.index( 'lon' )
		nens = rec.ensemble 
		ny = values.shape[ lat_index ]
		nx = values.shape[ lon_index ]
	else:
		is_ens = False
		nens = None
	
		#_alrighty, until we add more
		if field != 'aod': return 0

	#_Get plot specifics for field
	field	= field.lower()
	f_title = plt_dict[field]['title'] 
	levs	= plt_dict[field]['levs']
	rgb 	= lt.rgbgen( field, nens=nens )

	#_Strip out '_AOD' from variable name
	respec = re.compile( '_aod', re.IGNORECASE )
	spec = respec.sub( '', variable ) 

	#_Setup output file and directory
	dir_m = mod_dict[model]['web']
	if file == None:
		dir_out	= '/'.join(( path, dir_m, reg, dtg_init )) 
		file = '_'.join(( dtg_init, dtg_vald, 'f'+str(vhr).zfill(3), 
			spec.lower(), field, '550', reg.lower(),
			model.lower()+'.'+ifmt ))
		file = '/'.join(( dir_out, file ))
	else:
		dir_out = '/'.join( file.split('/')[:-1] )

	#_Make output directory
	if not show: lt.mkdir_p( dir_out )

	#_Setup output file name
	file = file.replace('_global','')
      	i_s = lt.human_date( dtg_init )
	v_s = lt.human_date( dtg_vald )

	#_Setup plot header
###    	title = i_s + ' ' + rec.model + ' Forecast t+' + str(vhr).zfill(3) \
       	title = i_s + ' ' + label + ' Forecast t+' + str(vhr).zfill(3) \
		+ '\n' + v_s + ' Valid Time\n' + spec.upper() + ' ' + f_title
	title = title.replace('NMMB-','NMMB/')
	if is_ens: title = ' '.join(( title, '(', 'nMEM', '=', str(nens), ')' ))

	#_Get time of plot generation
	n	= time.gmtime(time.time())
	cdtg    = str(n.tm_year) + str(n.tm_mon).zfill(2)  \
		+ str(n.tm_mday).zfill(2) + str(n.tm_hour).zfill(2)

	foot    = 'Plots Generated ' + lt.human_date(cdtg) \
		+ ' NRL/Monterey Aerosol Modeling\n' \
		+ mod_dict[model]['foot']
	if reg == 'byzantium': foot = ''

	#_Setup map
	if rec.region != 'npolar':
		proj = 'cyl'
		polar = 0 
	else:
		proj = 'npstere'
		polar = 1 

	#_Create MAP object
	m = draw_map( grid, corners, proj, delta, water=water, res=res )
	n = time.gmtime( time.time() )

	#_Write header and footer
	head_tail( title, foot, vert=polar )

	#_When taking individual slices, meta attrs need to be passed
	attrn, attrv = [], []
	for name in dimnames:
		if name == 'member': continue
		dim = values.__getattribute__( name )
		attrn.append( name )
		attrv.append( dim )
	#_______________________________________________________________________
	#_AOD___________________________________________________________________
	if field == 'aod':
		#_Get mean if data is an ensemble
		if is_ens:
			vmean = values.mean( axis=ens_index )
			values = cl.var( vmean, attrv=attrv, attrn=attrn )

		#_Plot
		d, x, y	= map_data( values, m )
		CS = m.contourf( x, y, d, levels=levs, colors=rgb )

		#_Plot Colorbar
 		plt.colorbar( CS, orientation='horizontal', shrink=0.9,
			aspect=40, pad=0.05 )

	#_SPAGHETTI PLOTS_______________________________________________________
	elif field == 'spag' and is_ens:
		#_Mask nans; this object now lacks metdat
		ndum, ny, nx = values.shape

		#_Setup colors
		rgb = lt.colors_spag( nens )
		rgb.append( '#FF0000' )
		rgb.append( '#000000' )

		#_Plot individual members
		for e in np.arange( nens ):
			col = rgb[e]
			ens_values = cl.var( values[e,:,:], attrv=attrv )
 
			d, x, y = map_data( ens_values, m )
			CS = m.contour( x, y, d, levels=levs, colors=col, 
				linewidths=0.5 ) 

		#_Plot ensemble mean
		ens_mean = cl.var( values.mean( axis=ens_index ), attrv=attrv )
		d, x, y = map_data( ens_mean, m )
		CS = m.contour( x, y, d, levels=levs, colors='#FF0000', 
			linewidths=1.2 )


	#_MEAN / SIGMA PLOTS____________________________________________________
	elif field == 'musig' and is_ens:
		#_Mask nans; this object now lacks metdat
		ndum, ny, nx = values.shape

		#_Find ensemble mean and map
		ens_mean = cl.var( values.mean( axis=ens_index ), attrv=attrv )
		mu, x, y = map_data( ens_mean, m )
		CS = m.contour( x, y, mu, levels=levs, lindwidths=1, colors=rgb)
		plt.clabel( CS, fontsize=8, inline=1, inline_spacing=1, 
			fmt='%1.1f' )

		#_Find ensemble standard deviation and normalize it, then map
		ens_std = values.std( axis=ens_index )
		std_2 = ens_std.std()
		if std_2 == 0.:
			ens_tmp = np.zeros(( ny, nx ))
		else:
			ens_tmp = ( ens_std - ens_std.mean() ) / std_2

		ens_nml = cl.var( ens_std, attrv=attrv )
		sig, x, y = map_data( ens_nml, m )
		CS = m.contourf( x, y, sig ) #, cmap=cm.jet )

	#_THRESHOLD
	elif field == 'warn' and is_ens:
		#_Above pulling of shape values not working for some reason
		ndum, ny, nx = values.shape

		#_Initialize threshold value.  warn_array is final plotting arr
		thresh = levs[0] 
		mask_array = np.zeros( values.shape )
		warn_array = cl.var( np.zeros(( ny, nx )), attrv=attrv  )

		#_Find where models are above threshold value
		indices = np.where( values.data >= thresh )

		#_Set those locations to 1
		try:
			mask_array[indices] = 1
		except:	
			dbg( 'No values above threshold ', l=2 )

		#_Could stop here, but need at least half to be above threash
		sum_array = np.sum( mask_array, axis=ens_index )
		indices = np.where( sum_array >= ( nens/2. ) )
		warn_array[indices] = 1
	
		#_Map and plot
		d, x, y = map_data( warn_array, m ) 
		CS = m.contourf( x, y, d, levels=levs, colors=rgb )

	#_FIELD NOT IMPLEMENTED
	else:
		dbg(( field, 'not implemented' ))
		plt.close()
		return

	#_Save image to file 
	dbg( file )
	if show:
		plt.show()
	else:
		plt.savefig( file, dpi=(dpi) )
		lt.make_readable( file )
	plt.close()

def plot_forker( *args, **kwargs ):
	'''
	Takes model_object() and splits the records in groups
		to be processed in parallel
	records	: np.recarray
	'''
	#_Put each record into threading groups
	groups = lt.setup_groups( *args, **kwargs )
        dbg( str( len(groups) ) +' group(s)' )

	#_Loop over processing groups
	for group in groups:
        	dbg( str( len(group) ) +' processes', l=3 )
		
		#_Initialize list to keep pids of child processes
		children = []
	  	for rec in group:
	   		pid = os.fork()
	   		if pid != 0: children.append(pid)
			elif pid == 0:
				try:
      					plot_rec( rec, **kwargs )
				except:
					dbg(( 'Failed to plot', rec.model,   
						rec.variable, rec.fhr,
						rec.dtg_init, rec.dtg_vald,
						rec.region ))

				#_Send exit signal
				os._exit(0)

		#_Pause until all children processes finish
	  	for kid in children: os.waitpid( kid, 0 )

def plot_epsgrams( d, path_out='.', label='epsgram' ):
	'''
	NOT UPDATED FOR STAT RECARRAY 
	Take in data in the form 
	d[point] = ndarray(nt, 7) 
	7 = max 90 75 med 25 10 min
	'''
	import matplotlib.pyplot as plt
	import libaeronet as la
	lt.mkdir_p( path_out )

	dtg_init = d.dtg_init
	fhr = str(d.fhr).zfill(3)
	finc = d.finc
	frq = 4
	frqy = 1.

	dtgs = []
	t=0
	while t <= int(fhr):
		dtgs.append( lt.newdtg( dtg_init, t ) )
		t += finc

	nt = len(dtgs)
	x = np.arange( nt ) + 0.5 

	points = d.keys()

	for p in points:
	   species = d[p].keys()
	   for s in species:
		coords, name = pnt_dict[p]
		coords = np.round(coords, 2)
		file_out = path_out + '/' + '_'.join(( label, dtg_init, 'f'+fhr,
			s, p+'.'+ifmt ))
		title = '_'.join(( label.upper(), s.upper(), 'AOD Meteogram at',
			name, str(coords), 'for', lt.human_date( dtg_init ) ))

		p100 = d[p][s][:,0] 
		p090 = d[p][s][:,1] 
		p075 = d[p][s][:,2] 
		p050 = d[p][s][:,3] 
		p025 = d[p][s][:,4] 
		p010 = d[p][s][:,5] 
		p000 = d[p][s][:,6] 

		h, w = plt.figaspect( 10 )
		fig = plt.figure( figsize=( w, h ) )
		ymax = np.ceil( np.max([3., p100.max()]) )

		tick_lab = dtgs[frq::frq]
		tick_loc = np.linspace( frq, nt, nt/frq ) 
		plt.xticks( tick_loc, tick_lab, size='small' )

		tick_lab = np.round( np.linspace( 0, ymax, 11 ), 2 )
		plt.yticks( tick_lab, size='small' )

		if x_coord != 'label':	plt.xlim( 0, nt )
		plt.ylim( 0, ymax )

		ax = plt.subplot( 111 )
		mjl = plt.MultipleLocator( frqy )
		mnl = plt.MultipleLocator( frqy/2 )
		ax.yaxis.set_major_locator( mjl )
		ax.yaxis.set_minor_locator( mnl )
		ax.grid( True, zorder=-1 )

		plt.title( title )
		fig.set_xlabel = ( 'DTG' )
		fig.set_ylabel = ( 'AOD' )

		plt.bar( x, p100-p000, width=.00001, bottom=p000, facecolor='#000000' ) 
		plt.bar( x-.2, p090-p010, width=.4, bottom=p010, facecolor='#6bffde' ) 
		plt.bar( x-.3, p075-p050, width=.6, bottom=p050, facecolor='#22dbb2' ) 
		plt.bar( x-.3, p050-p025, width=.6, bottom=p025, facecolor='#22dbb2' ) 

		dbg( file_out )
		plt.savefig( file_out )
		plt.close()

def head_tail( title, footer, vert=0 ):
	'''Add standard header and footer to plots'''
	import matplotlib.pyplot as plt
	plt.title( title, va='bottom', ha='left', position=(0.0,1.0), \
		size='small')
        if vert == 0:
		plt.figtext( 0.16, 0.1, footer, size='x-small' )
 	elif vert == 1:
		plt.figtext( 0.97, 0.85, footer, size='x-small', \
			rotation='vertical' )

def numeric_rank( ens, mem, metrics=['rmse'], path='.', **kwargs ):
	'''
	based on scores
	'''
	import numpy as np
	import time

	fhrs = lt.unique( ens.fhr )
	vars = lt.unique( ens.variable )
	dtgs = lt.unique( ens.dtgs )
	dtge = lt.unique( ens.dtge )
	mods = lt.unique( ens.mode )
	errs = lt.unique( ens.error_model )
	locs = lt.unique( ens.location )
	label = lt.unique( ens.label, unique=True )

	mems = lt.unique( ens.members, unique=True )
	model = lt.unique( ens.model, unique=True )
	all = [model]
	all.extend(mems)

	#_ensemble+ members

	path_out = '/'.join(( path, 'numeric_ranks' ))
	lt.mkdir_p(path_out)

	#_pull out appropriate matching metrics
	file = '-'.join(( dtgs[0], dtge[0], model.upper(), label+'.txt' ))
	file = '/'.join(( path_out, file ))
	f = open( file, 'w' )
	header = '%10s'*len(all) % tuple(all)
	f.write(header+'\n')
	dbg(vars)
	dbg(fhrs)
	dbg(dtgs)
	dbg(dtge)
	dbg(mods)
	dbg(locs)
	dbg(metrics)
#	time.sleep(3)
	for metr in metrics:
	  for var in vars:
	    for fhr in fhrs:
	      rank = np.zeros( len(all) )
	      for dtg0 in dtgs:
	        for dtg1 in dtge:
		  for mode in mods:
		    for loc in locs:
			try:
			  e = subset( ens, dtgs=dtg0, dtge=dtg1, 
				variable=var, fhr=fhr, mode=mode,
				location=loc, unique=True )
			  m = subset( mem, dtgs=dtg0, dtge=dtg1, 
				variable=var, fhr=fhr, mode=mode,
				location=loc )
			except ValueError:
				#print 'NONE:', dtg0, dtg1, var, fhr, mode, loc
				continue
			dbg((fhr,dtg0,dtg1,mode,var,loc,metr), l=3)	
			#_check to see if members are ordered the same as
			# we expect in mems
			mlist = e.members
			if mlist != mems: 
				raise RuntimeError, 'meh, make sure theyre in the same order'

			#_initialize array to temporarily put metrics into
			rank_tmp = np.zeros( len(all) ) - 9999.

			#_put ensemble score in first
			rank_tmp[0] = e.__getattribute__( metr )

			#_loop over members to get scores
			for member in mems:
				#_put into expected order in array
				midx = all.index(member)
				tmp2 = subset(m,model=member, 
					unique=True).__getattribute__(metr)
				
				if np.isnan(tmp2): continue

				#-put into array for ranking
				rank_tmp[midx] = tmp2

			#_mask missing values
			rank_tmp = np.ma.masked_where( rank_tmp == -9999., rank_tmp )

			order = rank_tmp.argsort()
			rank[order[0]] += 2
			rank[order[1]] += 1 
			
	      scores = '%10i'*len(all) % tuple(rank)
  	      f.write(scores + '\t' + var + '\t' + str(fhr) + '\t' + metr + '\n')

			#Rank models, score and add to rank
	f.close()
		



def plot_rankhistogram_point( stats, path='.', show=False, **kwargs):
	'''
	Attempts to plot rank histograms for ensemble data

	stats	: class.regional/point_stats
	show	: print to screen
	individual : make plots standalone

	Will make grid based on labels.
	Fhrs will be on same plot
	|___________________________________...|
	|________|__pnt0__|__pnt1__|__pnt2__...|
	|_label0_|___hr___|________|________...|
	|_label1_|________|________|________...|
	|_label2_|________|________|________...|
	|_label3_|________|________|________...|
	|...
	Different species on different pages
	'''
	import matplotlib.pyplot as plt
	import libaeronet as la
	var2mode = lm.var2mode

	#_Only the ensemble has rank data, so if passed full stat object
	# drop it down to something more useful
	if type( stats ) != np.recarray: stats = stats.ensemble.copy()

	#_Get all Iterations of Labels, fhrs, and regions
	error_models 	= lt.unique( stats.error_model ) 
	variables 	= lt.unique( stats.variable )
	labels 		= lt.unique( stats.label )
	model 		= lt.unique( stats.model, unique=True )
	modes 		= lt.unique( stats.mode ) 
	locs 		= lt.unique( stats.location )
	fhrs 		= lt.unique( stats.fhr )
	dtg 		= lt.unique( stats.dtgs )[0]

	points = lt.intersection([ pnt_dict.keys(), locs ]) 
	 
	#_Dimensions of plotting grid
	ncol 	= len( points )
	nrow 	= len( labels )
	nbar 	= len( fhrs )
  	colors 	= lt.colors_spag( nbar )
	nrow 	= len( labels )
	lw 	= 9. / ncol / (nbar/2.)

	#_For small nbar, some of the default colors are hard to see
	scol = ['red','black','orange','green','blue']
	scol = ['red','orange','green','blue','violet']
	nt = len( fhrs )
	if nt <= 4: colors = scol[:nt] 
	dbg("TEMPORARY WORKAROUND")
	colors = scol[:nt] 

	matplotlib.rcParams.update( { 'font.size' : 5 } )
	matplotlib.rcParams.update( { 'legend.fontsize' : 6 } )

	#_prog vs. diag shouldn't matter for rank vs obs
	stats = subset( stats.copy(), error_model=error_models[0] )

	dbg( points )
	dbg( variables )
	dbg( modes )
	dbg( fhrs )
	dbg( locs )
	path_out = '/'.join(( path, 'rank' ))
	lt.mkdir_p( path_out )

	#_loop over variable pages and create new plot
	for variable in variables:
		mode = var2mode[variable]

		#_Output file name
		file = path_out + '_'.join(( '/rhist_reg', model, variable,
			mode+'.'+ifmt )) 

		#_Create artist object
		h, w = plt.figaspect( 6 )
		fig = plt.figure( figsize=( w, h ) )
		#fig = plt.figure()
		idx = 1 
		a, b, cs = [], [], {}
 
		#_Loop over regions
		nrank = len( stats[0].ranks )
		x = np.arange( nrank )
		xmax = nrank-.5
		dbg('TEMPFIX');
		if variable != 'dust_aod': xmax-=1
		dbg((variable,xmax,nrank))
###		ymax = 0.
		axT = {}
		for label in labels:
			ymax = 0.
			ax = {}

			#_Loop over labels
			for point in points:
        		#h, w = plt.figaspect( 10 )
        		#fig = plt.figure( figsize=( w, h ) )
				ax[idx] = fig.add_subplot( nrow, ncol, idx )

				#_Create bar offset
				x_m = np.linspace( -0.4, 0.4, nbar )

				#_Loop over hours, plot vertical lines
				for fhr in fhrs:
					i = fhrs.index( fhr )

					#_Subset data records
					try:
						stats_sub = subset( stats, 
						fhr=fhr, variable=variable,
  						location=point, label=label, 
						unique=True, mode=mode )
					except ValueError:
						dbg('Should have more data',l=2)
						continue
				
					#_Get the ens+1 size and setup abscissa
 					y = stats_sub.ranks

					#_Plot vertical bars
					cs[fhr] = ax[idx].vlines( x+x_m[i], 
						[0], y, color=colors[i], lw=lw )
					tmp = np.max( np.array(y) )
					ymax = np.max( [tmp, ymax] )

					#_First time, make legend items
					if str(fhr) not in a:
						a.append( str(fhr) )
						b.append( cs[fhr] )

				#_Put ylabels on far left plots

				if not (idx-1) % ncol == 0: 
					plt.setp( ax[idx].get_yticklabels(),
						visible=False )
				else:
					ax[idx].set_ylabel( label.upper(), 
						rotation=90, size='small' )
						
				#_Set xlabel on top 
				if not idx > ncol*(nrow-1): 
					plt.setp( ax[idx].get_xticklabels(),
						visible=False )
				else:
					ax[idx].set_xlabel( point, size='small')

				#_Increment Index of Plot
				idx += 1

			#_Expand x/y area
			for n in ax:
				ax[n].set_xlim([ -.5, xmax ])
				ax[n].set_ylim([ 0, ymax+1  ])
				ax[n].tick_params( axis='y', which='both',
					left='off',right='off',labelbottom='off')
				ax[n].tick_params( axis='x', which='both',
					top='off',bottom='off')

			axT.update(ax)
		ax = axT

		#_Expand x/y area(A)
#		for n in ax:
#			ax[n].set_xlim([ -.5, nrank ])
#			ax[n].set_ylim([ 0, ymax+1  ]);ax[n].tick_params(axis='y',which='both',left='off',right='off',labelbottom='off')

###			for idx in ax:
###	                        #_add N values
###       	                x_srt = sorted( x.keys() )
###				ymax = 5 if ymax == 0 else ymax
###	                        #_setup grid
###				dbg( (nrank, ymax, ymax/5.) )
###	                        mjl = plt.MultipleLocator( nrank )
###	                        ax[idx].xaxis.set_major_locator( mjl )
###	                        mjl = plt.MultipleLocator( ymax / 5 )
###	                        ax[idx].yaxis.set_major_locator( mjl )
###	                        ax[idx].grid( False )
###	                        lab, loc = [], []
###                       		ax[idx].set_xticks( loc  )
###                        	ax[idx].set_xticklabels( lab )
###
###                        	if (idx-1) % ncol:
###                        	        plt.setp( ax[idx].get_yticklabels(),
###                        	                visible=False )
###	                        if not idx > ncol*(nrow-1):
###	                                plt.setp( ax[idx].get_xticklabels(),
###	                                        visible=False )
		
		#_reorder legend entries by value
		order = []
		for fhr in sorted( fhrs ): order.append( a.index(str(fhr)) )
		a = np.array(a)[order]
		b = np.array(b)[order]
	
		#_plot legend on last img, often clipped
		plt.legend( b, a, loc=2, bbox_to_anchor=(1.05,1), frameon=False) 

		#_write image and close
		if show == False:
			dbg( file )
			plt.savefig( file, dpi=(200) )
			lt.make_readable( file )
		else:
			plt.show()
		plt.close()	

def plot_rankhistogram( stats, path='.', show=False, **kwargs ):
	'''
	Attempts to plot rank histograms for ensemble data

	stats	: class.regional/point_stats
	show	: print to screen
	individual : make plots standalone

	Will make grid based on labels.
	Fhrs will be on same plot
	|___________________________________...|
	|________|__reg0__|__reg1__|__reg2__...|
	|_label0_|___hr___|________|________...|
	|_label1_|________|________|________...|
	|_label2_|________|________|________...|
	|_label3_|________|________|________...|
	|...
	Different species on different pages
	'''
	import matplotlib.pyplot as plt
	import libaeronet as la

	#_Only the ensemble has rank data, so if passed full stat object
	# drop it down to something more useful
	if type( stats ) != np.recarray: stats = stats.ensemble.copy()

	#_Get all Iterations of Labels, fhrs, and regions
	error_models 	= lt.unique( stats.error_model ) 
	variables 	= lt.unique( stats.variable )
	labels 		= lt.unique( stats.label )
	model 		= lt.unique( stats.model, unique=True )
	modes 		= lt.unique( stats.mode ) 
	locs 		= lt.unique( stats.location )
	fhrs 		= lt.unique( stats.fhr )
	dtg 		= lt.unique( stats.dtgs )[0]

	regs = lt.intersection([ reg_dict.keys(), locs ])
	points = lt.intersection([ pnt_dict.keys(), locs ]) 
	 
	#_Dimensions of plotting grid
	ncol = len( regs )
	nrow = len( labels )
	nbar = len( fhrs )
  	colors = lt.colors_spag( nbar )
	nrow = len( labels )
	dbg(( ncol, nbar ))
	lw = 9. / ncol / (nbar/2.)

	#_For small nbar, some of the default colors are hard to see
	scol = ['red','black','orange','green']
	nt = len( fhrs )
	if nt <= 4: colors = scol[:nt] 

	matplotlib.rcParams.update( { 'font.size' : 5 } )
	matplotlib.rcParams.update( { 'legend.fontsize' : 6 } )

	#_prog vs. diag shouldn't matter for rank vs obs
	stats = subset( stats.copy(), error_model=error_models[0] )

	dbg( regs )
	dbg( variables )
	dbg( modes )
	dbg( fhrs )
	dbg( locs )
	dbg( lt.intersection( [pnt_dict.keys(), locs]) )

	#_loop over variable pages and create new plot
	for variable in variables:
	    for mode in modes:
		#_Make output directory
		path_out = '/'.join(( path, mode ))
		lt.mkdir_p( path_out )

		#_Output file name
		file = path_out + '_'.join(( '/rhist_reg', model, variable,
			mode+'.'+ifmt )) 

		#_Create artist object
		fig = plt.figure()
		idx = 1 
		a, b, cs = [], [], {}
		first = True
 
		#_Loop over regions
		nrank = len( stats[0].ranks )
		x = np.arange( nrank )
		for label in labels:
			ymax = 0.
			ax = {}

			#_Loop over labels
			for reg in regs:
				ax[idx] = fig.add_subplot( nrow, ncol, idx )

				#_Create bar offset
				x_m = np.linspace( -0.4, 0.4, nbar )

				#_Loop over hours, plot vertical lines
				for fhr in fhrs:
					i = fhrs.index( fhr )

					#_Subset data records
					try:
						stats_sub = subset( stats, 
						fhr=fhr, variable=variable,
  						location=reg, label=label, 
						unique=True, mode=mode )
					except ValueError:
						dbg('Should have more data',l=2)
						continue
				
					#_Get the ens+1 size and setup abscissa
 					y = stats_sub.ranks

					#_Plot vertical bars
					cs[fhr] = ax[idx].vlines( x+x_m[i], 
						[0], y, color=colors[i], lw=lw )
					tmp = np.max( np.array(y) )
					ymax = np.max( [tmp, ymax] )

					#_First time, make legend items
					if first:
						a.append( str(fhr) )
						b.append( cs[fhr] )

				first = False
				#_Put ylabels on far left plots
				if (idx-1) % ncol == 0: 
					ax[idx].set_ylabel( label.upper(), 
							rotation=90, 
							size='small' )
				else: 
					ax[idx].yaxis.set_visible( False )

				#_Set xlabel on top 
				if idx > ncol*(nrow-1): #_could just use reg[0]
					ax[idx].set_xlabel( reg, size='small' )
				else:
					ax[idx].xaxis.set_visible( False )

				#_Increment Index of Plot
				idx += 1

			#_Expand x/y area
			for n in ax:
				ax[n].set_xlim([ -1.5, nrank ])
				ax[n].set_ylim([ 0, ymax+1  ])
		
		#_plot legend on last img, often clipped
		plt.legend( b, a, loc=2, bbox_to_anchor=(1.05,1), frameon=False ) 

		#_write image and close
		if show == False:
			dbg( file )
			plt.savefig( file, dpi=(200) )
			lt.make_readable( file )
		else:
			plt.show()
		plt.close()	

####	#_plot individual aeronet sites?
####	dbg( points )
####	ncol = len( variables )
####	nrow = len( labels )
####	lw = 9. / ncol / (nbar/2)
####	groups = lt.setup_groups( points, **kwargs )
####	dbg(( len(points), 'points', len(groups)))
####	for group in groups:
####	    #_initialize list to keep ids of subprocs
####	    children = []
####	    for pnt in group:
####		pid = os.fork()
####		if pid != 0: children.append( pid )
####		elif pid == 0:
####		   for mode in modes:
####			#_output file name
####			path_out = '/'.join(( path, mode ))
####			lt.mkdir_p( path_out )
####			file = path_out + '_'.join(( '/rhist_pnt', model, 
####				pnt, mode+'.'+ifmt )) 
####
####			#_Create artist object
####			fig = plt.figure()
####			a, b, cs = [], [], {}
####			first = True
#### 
####			ymax = 0.
####			ax = {}
####			idx = 1 
####
####			#_Setup x axis size
####			nrank = len( stats[0].ranks )
####			x = np.arange( nrank )
####			for label in labels:
####	  		    for variable in variables:
####				#_Loop over labels
####				ax[idx] = fig.add_subplot( nrow, ncol, idx )
####
####				#_Create bar offset
####				x_m = np.linspace( -0.4, 0.4, nbar )
####
####				#_Loop over hours, plot vertical lines
####				for fhr in fhrs:
####					i = fhrs.index( fhr )
####
####					#_Subset data records
####					try:
####						stats_sub = subset( stats, 
####						fhr=fhr, variable=variable,
####  						location=pnt, label=label, 
####						unique=True, mode=mode )
####					except ValueError:
####						dbg('Should have more data',l=2)
####						continue
####				
####					#_Get the ens+1 size and setup abscissa
#### 					y = stats_sub.ranks
####
####					#_Plot vertical bars
####					cs[fhr] = ax[idx].vlines( x+x_m[i], 
####						[0], y, color=colors[i], lw=lw )
####					tmp = np.max( np.array(y) )
####					ymax = np.max( [tmp, ymax] )
####
####					#_First time, make legend items
####					if first:
####						a.append( str(fhr) )
####						b.append( cs[fhr] )
####
####				first = False
####				#_Put ylabels on far left plots
####				if (idx-1) % ncol == 0: 
####					ax[idx].set_ylabel( label.upper(), 
####							rotation=90, 
####							size='small' )
####				else: 
####					ax[idx].yaxis.set_visible( False )
####
####				#_Set xlabel on top 
####				if idx > ncol*(nrow-1): #_could just use reg[0]
####					ax[idx].set_xlabel( variable, 
####						size='small' )
####				else:
####					ax[idx].xaxis.set_visible( False )
####
####				#_Increment Index of Plot
####				idx += 1
####
####			#_Expand x/y area
####			for n in ax:
####				ax[n].set_xlim([ -1.5, nrank ])
####				ax[n].set_ylim([ 0, ymax+1  ])
####
####			#_Plot legend on last img, often clipped
####			plt.legend( b, a ) 
####
####			#_Write image and close
####			if show == False:
####				dbg( file )
####				plt.savefig( file, dpi=(1000) )
####				lt.make_readable( file )
####			else:
####				plt.show()
####			plt.close()	
####
####		   os._exit(0)
####
####	    #_Wait for plots to wrap up and return exit code 0
####	    for kid in children: os.waitpid( kid, 0 )

def plot_errors( stats_ens, stats_mem, path='./error_plots', show=False,
	metric='rmse', row_name='label', col_name='region', line='member',
	x_coord='fhr', loc_type='region', size=3.2, lw=.4, anon=False,**kwargs):
	'''
	Attempts to plot rank histograms for ensemble data

	stats	: class.regional/point_stats
	metric	: str,	metric used ('rmse','fge', 'mae')

	loc_type: str,	whatever you're using for location, set here.
			do not mix and match...

	row	: str,	attribute iterated by row
	col	: str,	attribute iterated by column
	line	: str,	attribute iterated by line
	x_coord	: str,	attribute to define xcoordinate

	loc_type: str,	region | point, and also superfluous. Do through type
 
	Will make grid based on labels.
	|_________________________________...|
	|______|__col0__|__col1__|__col2__...|
	|_row0_|__mods__|________|________...|
	|_row1_|________|________|________...|
	|...
	Vars not in row_name, col_name, line or x_coord are on different pages
	'''
	import matplotlib.pyplot as plt
	import itertools
	import libaeronet as la

	#_make output directory
	path = '/'.join(( path, 'error' ))
	lt.mkdir_p( path )


	#_upper plot bounding by metric
	ymax_metric = { 'rmse'	: 0.5,
			'mae'	: 1.0,
			'fge'	: 2.0,
			'nrmse'	: 2.0,
			'briers': 1.0, }
	n_max = 10000

	#_get all Iterations of Labels, fhrs, and regions
	model 			= lt.unique( stats_ens.model, unique=True )
	dtgs 			= lt.unique( stats_ens.dtgs )
	disc 			= {}
	disc['error_model'] 	= lt.unique( stats_ens.error_model ) 
	disc['variable'] 	= lt.unique( stats_ens.variable )
	disc['label']		= lt.unique( stats_ens.label )
	disc['fhr']		= lt.unique( stats_ens.fhr )
	disc['member']		= lt.unique( stats_mem.model )
	locations		= lt.unique( stats_ens.location )

	#_allow region or poitns to be used as location
	loc 			= {}
	loc['region']		= lt.intersection([ reg_dict.keys(), locations])
	loc['point']		= lt.intersection([ pnt_dict.keys(), locations])

	if metric == 'briers': disc['member'] = []

	if col_name in ['region','point']: 
		disc['location'] = loc[col_name]
		col_name = 'location'
	elif row_name in ['region','point']:
		disc['location'] = loc[col_name]
		row_name = 'location'
	elif line in ['region','point']:
		disc['location'] = loc[col_name] 
		line = 'location'

	#_get ensemble size from names of members 
	nens = len( disc['member'] )

	#_add model to members if that is what defines the lines
	if line == 'member': 
		disc['member'].insert( 0, model )

	#_build x coords
	stats_ens = sort_rec( stats_ens.copy(), order=( x_coord, ) )
	x = {}
	if x_coord == 'label': 
		thresholds = lm.thresh_bins() 
		for label in disc['label']: 
			x[label] = np.min([ np.mean(thresholds[ label ] ), 1.5])
	else:
		for key in disc[ x_coord ]:
			x[key] = disc[ x_coord ].index( key )

	#_dimensions of plotting grid
	rows = lt.unique( disc[ row_name ] )
	cols = lt.unique( disc[ col_name ] ) 
	ncol = len( cols )
	nrow = len( rows )
	nlin = len( disc[ line ] )
	nx = len( x.keys() )

	#_setup line colors
	colors = lt.colors_spag( nlin )
	if model == 'ICAP' and line == 'member':
		for mem in disc['member']:	
			i = disc['member'].index( mem ) 
			colors[i] = mod_dict[mem]['color']	

	matplotlib.rcParams.update( { 'font.size' : 8 } )
	matplotlib.rcParams.update( { 'legend.fontsize' : 8 } )

	#_c-c-c-c-c-combo
	#_full list of what could be looped over
	loops_all = ['label','member','location','variable','fhr','error_model']

	#_these fields will each be their own image file 
	loop = []
	for opt in loops_all:
		if opt not in [ row_name, col_name, line, x_coord ]:
			loop.append( opt )
	loops = []
	tmp = itertools.product( dtgs, disc[loop[0]], disc[loop[1]] )
#	[ loops.append( grp ) for grp in tmp ]
	loops.extend( tmp )
	groups = lt.setup_groups( loops, **kwargs )

	dbg(groups)
	dbg(( 'x_coord', x_coord, 'row_name', row_name, 'col_name', col_name )) 

	#_fork over each image type 
	for group in groups:
	  children = []
	  for dtg, l0, l1 in group:
	    pid = os.fork()
	    if pid != 0: children.append( pid )
	    elif pid == 0:
		#_output file name
		file = path + '/' + '_'.join(( metric, dtg, model,
			str(l0).zfill(3), str(l1).zfill(3)+'.'+ifmt )) 
		#_create artist object
		fig = plt.figure()
		idx = 1 
		a, b, cs = [], [], {}
		first = True

		x_max = 0.0
		y_max = ymax_metric[ metric ] #0.750
		ax = {}
		ay = {}

		#_subset data
		kwargs	= {	'dtgs' : dtg,
				loop[0] : l0,
				loop[1] : l1,	}
		if 'member' in kwargs: del kwargs['member'] 
		sub_pge = subset( stats_ens, **kwargs )
		sub_pgm = subset( stats_mem, **kwargs )

		#_loop over rows 
		for row in rows:
			#_loop over columns 
			for col in cols:
				ax[idx] = fig.add_subplot( nrow, ncol, idx )
				ax[idx].__setattr__('nobs', {})
				ax[idx].set_ylim( 0.0, y_max )
				ay[idx] = ax[idx].twinx()
				ay[idx].set_zorder(0)
				ay[idx].set_ylim( 0, n_max )
		
				#_y values for plottin gvertical lines for nobs
				nobs_y = []

				kwargs	= {	'dtgs' : dtg,
						row_name : row,
						col_name : col }
				if 'member' in kwargs: del kwargs['member'] 

				#_initialize coordinates to plot
				sub_ens = subset( sub_pge, **kwargs )
				sub_mem = subset( sub_pgm, **kwargs )
	
				#_for each line
				for le in disc[ line ]:
				  #_linestyle shoudl be dots unless is mean
				  lss = ':'

				  #_loop over bins, get x/y values
				  i = disc[line].index( le )
				  xe, ye = [], []

				  #_if ensemble, pass that; if member, that
				  if line == 'member':
					if le == model:
						sub_tmp = sub_ens.copy()
						lss = '-'
					else: 
						tmp = sub_mem.copy()
						sub_tmp = subset( tmp,model=le )
				  #_if members not line, don't bother with them	
				  else:
					sub_tmp = sub_ens.copy()

				  #_build x/y values
				  x_srt = sorted( x.keys() )
				  for k in x_srt:
					v = x[k]
					try:
						kwargs = {	x_coord : k,
								'unique' : True}
						st_bin = subset( sub_tmp,
							**kwargs )
					except:
						dbg( 'not enough data', l=2 )
						ax[idx].nobs[k] = 0 
						continue
				
					#_Put error in y, and x in x	
					value = st_bin.__getattribute__(metric)
					xe.append( v )
					ye.append( value )
			
					#_Store nobs (of first line)
					if i == 0:
						ax[idx].nobs[k] = st_bin.nobs
						nobs_y.append( st_bin.nobs )
						
				  #_no data, skip 
				  if len(xe) == 0: continue	

				  xe = np.ma.masked_where( xe == -9999., xe ) 
				  ye = np.ma.masked_where( ye == -9999., ye ) 

				  #_plot ensemble error
				  c = colors[i]
				  tru_arr = [True]*len(ye)
				  if type( ye.mask ) != np.bool_ and \
					ye.mask.tolist() == tru_arr:
					dbg( 'NO DATA!', l=2 )
					continue

				  #_create independent yaxis with shared x
				  if i == 0:
					  ay[idx].vlines( xe, [0], nobs_y, 
						lw=3, color='#b1a177' )
					  tick_loc = np.linspace( 100, n_max, 5)
					  ay[idx].set_yticks( tick_loc )

				  #_setup fit line
##				  fit = np.polyfit( xe, ye, 2 )
##				  fit_fn = np.poly1d( fit )
				  cs[le] = ax[idx].plot( xe, ye, 
					color=c, ls=lss, lw=.4 )
				  try:
				    ax[idx].scatter( xe, ye, marker='x', 
					color=c, s=size, lw=lw,
					facecolors='none' )
				  except ValueError:
				    dbg( 'ngac probbaly nan', l=4 )
				  x_max = np.max( [np.max(xe), x_max] )
###				  y_max = np.max( [np.max(ye), y_max] )

				  #_for Legend	
				  if first:
					print 'AWEE', le
					if anon and le != 'ICAP':
						continue
					a.append( le )
					b.append( cs[le] )

				#_done adding legend labels	
				first = False

				#_put ylabels on far left plots
				if (idx-1) % ncol == 0: 
					ax[idx].set_ylabel( row, rotation=90, 
							size='small' )

				#_put ylabels on far right plots
				if idx % ncol == 0:
					ay[idx].set_ylabel('nobs', 
							rotation=270, 
							size='small' )

				#_set xlabel on top  
				if idx > ncol*(nrow-1): 
					ax[idx].set_xlabel( str(col).zfill(2), 
						size='large' )
						#size='x-small' )

				#_increment Index of Plot
				idx += 1

		#_set all plots to same grid
		for idx in ax:
			#_add N values
			x_srt = sorted( x.keys() )
			if hasattr( ax[idx], 'nobs' ):
				nobs = ax[idx].nobs
				x_n, y_n = 0.1, y_max + 0.05
				obs_str = ','.join( str(nobs[k]) for k in x_srt)
				ax[idx].text( x_n, y_n, 'N='+obs_str, 
					size='xx-small')
			#_setup grid
			mjl = plt.MultipleLocator( x_max / 12 )
			ax[idx].xaxis.set_major_locator( mjl )
			mnl = plt.MultipleLocator( x_max / 24 )
			ax[idx].xaxis.set_minor_locator( mnl )
			mjl = plt.MultipleLocator( y_max / 5 )
			ax[idx].yaxis.set_major_locator( mjl )
			mnl = plt.MultipleLocator( y_max / 50 )
			ax[idx].yaxis.set_minor_locator( mnl )
			ax[idx].grid( False )
			lab, loc = [], []
			for k in x_srt:
				lab.append( str( k ) ) 
				loc.append( x[k] )
			ax[idx].set_xticks( loc  )
			ax[idx].set_xticklabels( lab )
			 
			if (idx-1) % ncol: 
				plt.setp( ax[idx].get_yticklabels(), 
					visible=False )
			if idx % ncol or rows.index( row ) != 0:
				plt.setp( ay[idx].get_yticklabels(),
					visible=False )

			if not idx > ncol*(nrow-1):
				plt.setp( ax[idx].get_xticklabels(), 
					visible=False )
			#_Set to log scale
			if x_coord == 'label': 
				ax[idx].set_xscale( 'log' )
				ax[idx].set_xlim( 0.0, x_max + 0.01 )
			else:
				ax[idx].set_xlim( -0.01, x_max+.01 )

		#_plot legend on last img, often clipped
		plt.legend( b, a, borderaxespad=0., frameon=False, 
			bbox_to_anchor=(1.05,1.), loc=2 ) 
		
		#_write image and close
		if show == False:
			dbg( file )
			plt.savefig( file, dpi=(200) )
			lt.make_readable( file )
		else:
			plt.show()
		plt.close()

		os._exit(0)

	  for kid in children: os.waitpid( kid, 0 )		

#############################################################################_80
######_PLOT-MAPPING_############################################################
################################################################################

def addcyclic( record ):
    	'''
	``arrout, lonsout = addcyclic(arrin, lonsin)``
	adds cyclic (wraparound) point in longitude to ``arrin`` and ``lonsin``.
	'''
	arrin = record.values
	nlons = arrin.lon.size #arrin.shape[1]
	nlats = arrin.lat.size #arrin.shape[0]
	lonsin = arrin.lon

	#_Find lon dimension position
	lonidx = record.dimname.index( 'lon' )

	#_Get in shape, then add one to lon dimension
	arr_shape = list( record.dimsize )
	arr_shape[ lonidx ] = arr_shape[ lonidx ] + 1
	arr_shape = tuple( arr_shape )

	#_Actually add extra band of data
	arrout = np.zeros( arr_shape, arrin.dtype )
	if arrin.ndim == 2:
		arrout[:,0:nlons] 	= arrin[:,:]
		arrout[:,nlons] 	= arrin[:,0]
	elif arrin.ndim == 3:
		arrout[:,:,0:nlons] 	= arrin[:,:,:]
		arrout[:,:,nlons] 	= arrin[:,:,0]

	#_Add extra longitude value
	lonsout = np.zeros( arr_shape[ lonidx ], lonsin.dtype )
	lonsout[0:nlons] = lonsin[:]
	lonsout[nlons] = lonsin[-1] + lonsin[1]-lonsin[0]

	return arrout, lonsout

def slide_grid( record, lon0, reg=None ):
    	'''
	Shift global lat/lon grid east or west.
	assumes wraparound (or cyclic point) is included.
	'''
	#_Add cyclic point
	values_in, lons_in = addcyclic( record )
	
	polar = True if reg[-5:] == 'polar' else False
	if np.fabs( lons_in[-1]-lons_in[0]-360. ) > 1.e-4: # and not polar:
		print lons_in[-1],  lons_in[0], lons_in[-1]-lons_in[0]-360.
        	raise ValueError, 'cyclic point not included'
	if (lon0 < lons_in[0] or lon0 > lons_in[-1]) and not polar:
        	print lon0, lons_in[0], lons_in[-1]
        	raise ValueError, 'lon0 outside of range of lonsin'

	#_Setup output arrays
	values_out = np.zeros( values_in.shape, values_in.dtype )
	lons_out = np.zeros( lons_in.shape, lons_in.dtype )

	#_Find where the llcorner will be of shift
	i0 = np.argmin( np.fabs( lons_in - lon0 ) )

	lons_out[ 0:lons_in.size-i0 ] = lons_in[ i0: ]
	lons_out[ lons_in.size-i0: ] = lons_in[ 1:i0+1 ] + 360.

	if values_out.ndim == 2:
		values_out[ :, 0:lons_in.size-i0 ] 	= values_in[ :, i0: ]
		values_out[ :, lons_in.size -i0: ] 	= values_in[ :, 1:i0+1 ]
	elif values_out.ndim == 3:
		values_out[:, :, 0:lons_in.size-i0 ] 	= values_in[:,:,i0: ]
		values_out[:, :, lons_in.size -i0: ] 	= values_in[:,:,1:i0+1 ]

	#_Settup attribute values to pass back
	attrv = []
	for name in record.dimname:
		dim = record.values.__getattribute__( name )
		if name == 'lon': dim = lons_out 
		attrv.append( dim )	

	return values_out, attrv 

def map_data( values, m, latlon=None ):
	''' produces x/y coords and adjusted indices for particular map	'''
	if latlon == None:
		lats = values.lat
		lons = values.lon
		ny = lats.size 
		nx = lons.size
	else:
		lats, lons = latlon
		ny, nx = [ len(n) for n in latlon ]
	
	#_make meshgrid of lat/lons
	lats2d = np.ones(( ny, nx )) * lats.reshape( ny, 1 )
	lons2d = np.ones(( ny, nx )) * lons.reshape( 1, nx )

	#_map lats/lons to x/y values on the map
	x, y = m( lons2d, lats2d )
	return values, x, y


def draw_map(grid=[-80,81,180,181], corners=[-90,90,-180,180], proj='cyl', delta=40, 
	water='white', fill='#dcbc80', res='c', dellat=None, dellon=None, 
	labx=None, laby=None, **kwargs):
	from mpl_toolkits.basemap import Basemap

	if (proj == 'npstere'):		#-Currently only supporting NPolar
		m = Basemap(projection=proj, boundinglat=grid[0], lon_0=0, 
			resolution=res, **kwargs)
	  	labelsx	= [1,1,0,0]
	  	labelsy	= [0,0,0,0]
	elif (proj == 'cyl'):		#-Or Cylindrical Equidistant
	  	m = Basemap(projection=proj, 
			llcrnrlat=corners[0], urcrnrlat=corners[1],
			llcrnrlon=corners[2], urcrnrlon=corners[3],
			resolution=res, **kwargs)
	  	labelsx	= [0,0,0,1]
	  	labelsy	= [1,1,0,0]
	labelsx = labelsx if labx == None else labx
	labelsy = labelsy if laby == None else laby
	
	m.drawcoastlines(linewidth=0.5)
	m.fillcontinents( color=fill, lake_color=water, zorder=0 )
	if ( corners[0] != -90 and corners[1] != 90 ): m.drawcountries()
	dellat = delta if dellat == None else dellat
	m.drawparallels( np.arange(grid[0], grid[1], dellat ), linewidth=0.5,
		dashes=[4,1], labels=labelsy, size='xx-small' )

	dellon = delta if dellon == None else dellon
	m.drawmeridians( np.arange(grid[2], grid[3],dellon), linewidth=0.5,
		dashes=[4,1], labels=labelsx, size='xx-small' )

	return m


#############################################################################_80
######_FILE-IO_#################################################################
################################################################################

def parse_env( path_env='path.env', parameter_env='scripts/parameters.env',
	**kw ):
	'''
	Read in NVA path.env or parameters.env file and add them to 
	return dictionary. Relies on the variable name and the semicolon
	ending the line to be on the same line.

	file	: str,	path to file to be parsed
	'''
	import re
	#_initialize file options
	opts = {}
	for file in [ path_env, parameter_env ] :
  	    with open( file, 'r' ) as f:
		for line in f:
			#_skip comments
			if re.match('#', line): continue
	
			#_search search line for a key and value
			dbg( line, l=9 )	
			res = re.search( '\$([\w\d_]+)\s*=\s*(.*?)\s*;', line )
			if res:
				key = res.group(1)
				value = res.group(2)
	
				#_check for perl style environment variables
				env = re.search('(\$ENV{["\'](.*)["\']})',value)
				if env:
					pl_env = env.group(1)
					py_env = os.environ[ env.group(2) ]
	
					#_put environment variable into string
					value = value.replace( pl_env, py_env ) 
					
				#_strip white space and quotes	
				value = re.sub( '["\'\s]', '', value )
	
				#_put search results into return dict
				opts[key] = value
		
	#_loop over until all nested variables fleshed out
	local_variables = True
	while local_variables:
		#_if any present, switches to True, allow +1 iteration
		# necessary to find any nesting
		local_variables = False
		for k, v in opts.iteritems():
			env = re.search( '(\$(.*?))[\s/;]', v )

			#_if value contains a variable, find explicit value
			if env:
				var_path = env.group(1)
				var_repl = opts[ env.group(2) ]
		
				#_peel one layer of onion
				dbg(( var_path, var_repl ), l=9 )
				v = v.replace( var_path, var_repl )

				#_put back into return dictionary
				opts.update({ k : v })

				local_variables = True 

	#_make everything lowercase
	for k, v in opts.iteritems():
		del opts[k]
		opts[k.lower()] = v

	return opts

def write_smoke( dtg, inc_emis=1, dir_emis='smoke',
	dtg_start='1970010100', arc_modis_hourly=None, arc_goes_hourly=None,
	smoke_sensors=['modis','goes'], forecast=False, **flags ):
	'''
	Generate smoke input files for NAAPS
	dtg	: str*10, 	date time group of smoke valid time
	dtg_strt: str*10, 	initialization time of model, used for forecasts
	inc_emis: int,		timestep between emission files
	dir_emis: str,		path to naaps input directory
	'''
	import libsat as ls
	import re
	
	#_if arc_modis_hourly or arc_goes_hourly not defined,
	# look for path.env file
	if arc_modis_hourly == None or arc_goes_hourly == None:
		flags.update( parse_env( **flags ) )
		arc_modis_hourly = flags.get('arc_modis_hourly')
		arc_goes_hourly = flags.get('arc_goes_hourly')

	#_name output smoke file, and if it already exists, don't remake
	file_smoke = '/'.join(( dir_emis, dtg+'_smoke.dat' ))
	if os.path.exists( file_smoke ): return file_smoke
	lt.mkdir_p( dir_emis )	

	#_calculate first dtg to include
	dtg_emis = lt.newdtg( dtg, -inc_emis )
	inc_hour = 1

	#_initialize list of files to include in file_smoke
	sensor_files = []

	#_loop over sensor data
	for sensor, metadata in sen_dict.iteritems():
		#_skip if not in list
		if sensor not in smoke_sensors: continue

		#_go through each hour and make list and local copies of
		# sensor files to concatenate into a single input
		dtg_loop = dtg_emis

		#_again, this is a dumb way of doing this
		if sensor == 'modis': dir_sens = arc_modis_hourly
		elif sensor == 'goes': dir_sens = arc_goes_hourly
		else: raise RuntimeError, 'no directory defined for sensor'

		#_inputs are in all sorts of compression forms,
		# so make a local copy to do with as we please
		# before concatenation. There will be one local copy
		# for each sensor class
		local_copy = '/'.join(( dir_emis, dtg+'_'+sensor+'.txt' ))

		while dtg_loop < dtg:
			#_forecasts should use the previous day of fires,
			# maintaining its diurnal cycle	

			#_package keyword args for finding satfiles
			kwgs = { 'sensor' : sensor, 'path' : dir_sens }

			#_if forecast, use the smoke files from the 
			# day prior to init
			if not forecast:
				satfiles = ls.find_satfiles( dtg_loop, **kwgs )
			elif forecast:
				dtg_tmp = dtg_start[:8] + dtg_loop[8:]
				dtg_tmp = lt.newdtg( dtg_tmp, -24 )
				satfiles = ls.find_satfiles( dtg_tmp, **kwgs )

			#_dump all satefiles into local_copy
			for satfile in satfiles:
				dbg( satfile, l=5 )
				if satfile[-2:] == '.Z':
					cmd = 'zcat '+satfile+' >> '+local_copy
				else:
					cmd = 'cat '+satfile+' >> ' +local_copy

				#_could use gzip modules, but get 'er done
				os.system( cmd )

			#_increment foreward by one hour
			dtg_loop = lt.newdtg( dtg_loop, inc_hour )	

		#_add this dtg_sensor file to files to concatenate
		sensor_files.append( local_copy )

	#_dump smoke files into file_smoke
	dbg( file_smoke )
	with open( file_smoke, 'w' ) as fw:
		for sensor_file in sensor_files:
			with open( sensor_file, 'r' ) as fr:

				for line in fr:
					#_skip empty lines 
					if re.search( '^\s*$', line ): continue

					#_skip modis west of 30W
					if re.search( 'modis.txt$',sensor_file):
						cols = line.split()
						if cols[1] < -30: 
							dbg(cols)
							continue	

					fw.write( line )

	#_return smoke file name
	return file_smoke

def read_ieee( file='2007013100_ieee_press' ):
	''' UPDATE OR JUST DELETE '''
	f = open(file,'rb')

	nx, ny, nz, nz_nog = 360, 180, 25, 31
	fields, order = lm.ieee_vars(nz, ny, nx)
	out = {}

	for v in order:
		size 	= np.multiply.reduce(fields[v]['shape'])
		dt	= np.dtype([( 'a', str(size)+'f4' )])
		tmp 	= lt.byte_order( np.fromfile( f, dt, count=1 ) )
		meta	= [ '1900010100', v, 'ieee' ]
		out[v]= nrl_var(tmp['a'].reshape(fields[v]['shape']),meta=meta)	
		out[v].setdims(fields[v]['dims'])
		out[v].setunits(fields[v]['units'])
		out[v].setlong_name(fields[v]['long_name'])
		out[v].setfillvalue(fields[v]['_FillValue'])
	f.close()
	return out

def write_conc( records, dir_out='.', faux=False ):
	'''
	fields must be object with same format as read in by read_conc()
	'''
	#_Get metadata
	conc = subset( records, variable='conc' )
	ns, nz, ny, nx = conc[0].dimsize
	dtg = lt.unique( records.dtg_vald )[0]
	fhr = lt.unique( records.fhr )[0]
	
	#_Setup output file
	file = dir_out + '/' + dtg + '_conc'
	if os.path.exists(file): 
		os.unlink( file )
	dbg( file )

	#_Add temporal metadata to file
	f = open( file,'wb' )
	d = np.array( 8, dtype='i' ).byteswap()		; d.tofile(f)
	d = np.array( dtg, dtype='i' ).byteswap()	; d.tofile(f)
	d = np.array( fhr, dtype='f' ).byteswap()	; d.tofile(f)
	d = np.array( 8, dtype='i' ).byteswap()		; d.tofile(f)
	f.close()

	#_Add spatial metadata to file
	w = { 'data' : [nx, ny, nz, ns], 'type' : 'i', 'size' : 16 }
	lt.binary_append( w, file, bs=True )	

	#_Get order and metadata for individual vars, then write to file
	vars, order = lm.conc_vars( ns, nz, ny, nx )
	for v in order:
		if v == 'sigma':
			sig_a = subset( records, variable='sig_a' )[0]
			sig_a = sig_a.values.flatten()
			sig_b = subset( records, variable='sig_b' )[0]
			sig_b = sig_b.values.flatten()
			d = np.append( sig_a, sig_b )
		else:
			d = subset( records, variable=v )[0]
			d = d.values.flatten()
		s = len(d)
		w = { 'data' : d, 'type' : 'f', 'size' : s }
		lt.binary_append( w, file, bs=True )

''' converts file to netcdf in same directory '''
def conc2ncdf( file ): to_ncdf( read_conc( file ) )

def read_conc(file, pipe=None, ncdf=False, ensemble=0, **kwargs):
	'''
	USAGE: <recarray> = read_conc(<FILENAME_STR>,[Pipe obj])
		Reads standard NAAPS files, assumes little endian 

	ensemble should always be false. Reading of an ensemble of conc
	files will repeatedly call this and set ensemble True on the approp
	master array, not the individual members	
	'''
	dbg( file )
	model = 'NAAPS'
	if os.path.exists(file + '.gz'): os.system('gunzip ' + file + '.gz')

	if ncdf: 
	  return read_conc_ncdf( file, pipe=pipe, ensemble=ensemble, **kwargs )
 
	f = open(file,'rb')

	#_Initialize concentration object
	conc = cl.model_object()

	#_Get temporal metadata
	dtype = [('dtg','i4'), ('fhr','f4')]#, ('minutes','i4')]
	rec = lt.rec_read(f, dtype)
	dtg_vald = rec[0]
	fhr = rec[1]
	dtg_init = lt.newdtg( dtg_vald, -fhr )

	#_Get spatial metadata 
	dtype = [('dims','4i4') ]
	rec = lt.rec_read( f, dtype ) 
	nx, ny, nz, nspec = rec['dims']
	cbyte = str( nx*ny*nz*nspec )

	#_Add species names
	species = np.array(['so2','sulfate','dust','smoke','seasalt'])
	conc.resize( conc.size + 1 )
	conc[-1] = ( species[:nspec], model, dtg_init, dtg_vald, fhr, 
		'global', 'species', ensemble, ('species',), (nspec,), '', 
		'name of model species' )

	#_Get table of conc file metadata and begin reading
	var_meta, var_order = lm.conc_vars( nspec, nz, ny, nx ) 
	for var_name in var_order:
		dimsize = np.array( var_meta[var_name]['shape'] ) 
		dtype = var_meta[var_name]['dtype'] 

		#_Read record length, data, then closing record length
		size_exp = dimsize if dimsize.size == 1 else \
			 np.multiply.reduce( np.array(dimsize) )

		#_Read record
		data = lt.rec_read( f, dtype )

		#_Add data to return object
		for name in data.dtype.names:
			dimname = var_meta[name]['dims']
			dimsize = var_meta[name]['shape']
			units = var_meta[name]['units']
			long_name = var_meta[name]['long_name']
			ndim = len(dimsize)

			#_Increase recarray size and reshape data
			tmp = data[name].reshape( dimsize )
	
			#_Make list of dimension metadata
			attrv = []
			for dim in dimname:
				if ndim < 2: break
				attr = subset( conc, variable=dim ).values[0]
				attrv.append( attr )

			#_Add sigma midpoint values
			if name == 'sig_b':
				conc.resize( conc.size + 1 )
				sig_b = tmp[:] 
				sigma = sig_b[:nz]-(sig_b[:nz] - sig_b[1:])/2
				conc[-1] = ( sigma, model, dtg_init, \
					dtg_vald, fhr, 'global', 'sigma', 
					ensemble, ('sigma',), (nz,), '', 
					'hybrid level midpoints' )

			#_Add record and fill
			if ndim > 1:
				tmp = cl.var( tmp, attrv=attrv, attrn=dimname )
	
			conc.resize( conc.size + 1 )
			conc[-1] = ( tmp, model, dtg_init, dtg_vald,
				fhr, 'global', name, ensemble, dimname, dimsize,
				units, long_name )

	#_POST_METADATA________________________________________________________
	try:	#_Newer versions of conc have extra records 
		# at the end of the file. It will produce a 
		# MemoryError if it is an old version and 
		# simply move on
		dtype =[('z0','i4'), ('version','a20'), ('z1','i4'),
			('z2','i4'), ('fhour','i4'), ('z3','i4'),
			('z4','i4'), ('n_obs','i4'), ('z5','i4'),
			('z6','i4'), ('modis','a430'), ('z7','i4')	]
	  	data = np.fromfile( f, dtype, count=1 )
		conc.__setattr__( 'naaps_version', data['version'][0] )
		conc.__setattr__( 'fhr', data['fhour'][0] )
		conc.__setattr__( 'n_obs', data['n_obs'][0] )
		conc.__setattr__( 'modis_version', data['modis'][0] )

		#_This is set up to accept an ever expanding flag length 
	  	#_Currently, 2 == 'n/a', 1 == 'True', 0 == 'False'
	  	# 111
	  	# |||_sedimentaion
	  	# ||__ginoux dust source
	  	# |___cmorph precipitation
          	dtype = [('bflen','i4')]		#_get size of 
          	tmp = np.fromfile( f, dtype, count=1 )	# flags rec
		bflen = str( tmp['bflen'][0] )
          	dtype = [ ('bitflag', 'a'+bflen), ('z','i4') ]
          	data = np.fromfile( f, dtype, count=1 )	#_Read that rec
		conc.__setattr__( 'bitflag', data['bitflag'][0] )

	except MemoryError:
		dbg('No metadata for ' + file, l=2)

	except RuntimeError:
		pass

	except IndexError:
		pass

	try:    #_More stapled on metadat - this time for SRCSPEC output
          	dt = [	('z6','i4'),('srcspec','a300'),		('z7','i4'),
			('z0','i4'),('flambe_version','a7'), 	('z1','i4')]
          	data = np.fromfile( f, dt, count=1 )
          	#data = lt.rec_read( f, dt )
		conc.__setattr__( 'srcspec', data['srcspec'][0] )
		conc.__setattr__( 'flambe_version', data['flambe_version'][0] )

	except MemoryError:
		dbg( 'No srcspec metadata for ' + file, l=2 )

	except RuntimeError:
		dbg( 'numpy version not fully supported' )

	except IndexError:
		dbg('Not exactly sure, but failed on metadata')


	if pipe == None	: return conc 
	else		: pipe.send(conc)

def read_conc_ncdf( file, pipe=None, ensemble=0, **kwargs ):
	'''
	USAGE: <recarray> = read_conc(<FILENAME_STR>,[Pipe obj])
		Reads standard NAAPS files, assumes little endian 

	ensemble should always be false. Reading of an ensemble of conc
	files will repeatedly call this and set ensemble True on the approp
	master array, not the individual members	
	'''
	from netCDF4 import Dataset
	dbg( file )
	model = 'NAAPS'

	#_Initialize concentration object
	conc = cl.model_object()

	if os.path.exists(file + '.gz'): os.system('gunzip ' + file + '.gz') 
	f = Dataset( file, mode='r', format='NETCDF3_CLASSIC' ) 

	dtg_vald = f.dtg_vald
	fhr = f.runhours
	dtg_init = lt.newdtg( dtg_vald, -fhr )
###	minutes = f.minutes

	#_get dimension metadata 
	nx = len( f.dimensions['lon'] )
	ny = len( f.dimensions['lat'] )
	nz = len( f.dimensions['sigma'] )
	nspec = len( f.dimensions['specie'] )

	#_add species names
##	try:
##		species = np.array( f.species.split(',') )
##		''' check to see what these... actually... are... '''
##	except:
	species = np.array(['so2','sulfate','dust','smoke','seasalt'])

	conc.resize( conc.size + 1 )
	conc[-1] = ( species[:nspec], model, dtg_init, dtg_vald, fhr, 
		'global', 'species', ensemble, ('species',), (nspec,), '', 
		'name of model species' )

	#_make sigma levels
	sig_b = f.variables['sig_b'][:]
	sigma = sig_b[:nz]-(sig_b[:nz] - sig_b[1:])/2
	conc.resize( conc.size + 1 )
	conc[-1] = ( sigma, model, dtg_init, dtg_vald, fhr, 'global', 'sigma', 
		ensemble, ('sigma',), (nz,), '', 'hybrid level midpoints' )

	#_get table of conc file metadata and begin reading
	var_meta = lm.conc_vars_ncdf( nspec, nz, ny, nx ) 
	for var_name in f.variables:
		dimsize = f.variables[ var_name ].shape
		dimname = f.variables[ var_name ].dimensions
		ndim = f.variables[ var_name ].ndim

		#_still need to add these to FORTRAN 
		units = var_meta[ var_name ]['units']
		long_name = var_meta[ var_name ]['long_name']

		#_increase recarray size and reshape data
		tmp = f.variables[ var_name ][:]
	
		#_make list of dimension metadata
		attrv = []
		if ndim > 1:
			for dim in dimname:
				if dim == 'specie':
					attr = subset( conc, 
						variable='species' ).values[0]
				elif dim == 'sigma':
					attr = subset( conc, 
						variable='sigma' ).values[0]
				else:
					attr = f.variables[dim][:] 
				attrv.append( attr )

			#_add record and fill
			tmp = cl.var( tmp, attrv=attrv, attrn=dimname )
	
		conc.resize( conc.size + 1 )
		conc[-1] = ( tmp, model, dtg_init, dtg_vald,
			fhr, 'global', var_name, ensemble, dimname, dimsize,
			units, long_name )

	if pipe == None	: return conc 
	else		: pipe.send(conc)

def convert_conc2kg( data, units='gridbox', rho=None ):
	'''
	UPDATE FOR RECARRAYS
	ALSO DO NOT TRUST.  NOT TESTED.
	units = {'gridbox' | 'm2'}

	CMAP == Area of each grid cell (nx,ny) [m2 / gridbox]

	Taken from cmass.f
	c calculate mass kg for species n1 through n2
c
c species 1 - 2:
c model units are  ppbvs / density, at stp, or
c          [s atoms / m**3 air] [ 1 / billion air atoms / m**3 @ stp ] / density
c multiply by grid volume, billions air atoms / m**3 @ stp, and density, and 
c convert s atmos to kg s,
c or
c mass (kg) = model_conc * delt-x * delta-y * delta-z * bn_air_stp
c             * density * m_sulfur / avagadro's number
c
c species 3 - 5:
c model units are kg / m**3 / density
c multiply by grid volume and density
c or
c mass (kg) = model_conc * delt-x * delta-y * delta-z * density
c
	real avagadro, m_sulfur, bn_air_stp, const
	avagadro = 6.023e23 ! # molecules / mole
	m_sulfur = 32.0e-3  ! kg / mole sulfur
	bn_air_stp = 2.65e16 ! # air molecules / cubic meter at stp
c
	cmass=0.0
c
	do k=n1,n2
	  cmassk(k) = 0.0
	  if ( k .le. 2 ) then
	   const = bn_air_stp * m_sulfur / avagadro
	  else if ( k .ge. 3 .and. k .le. nequat ) then
	   const =  1.0
	  else
	   write(6,*) 'cmas: nbin > ',nequat,' not allowed'
	   stop 4
	  endif
	  do  l=1,nz
	    do  i=1,ny*nx
	     rhodz = sp(i)*(bsig(l)-bsig(l+1))/ 9.8
	     cmassk(k) = cmassk(k) +c(i,l,k)*rhodz*cmap(i) * const
	    enddo
	  enddo
	 cmass = cmass + cmassk(k)
	 if (debug) write(6,*) 'cmas: k, cmass(kg) = ',k,cmassk(k)
	enddo
c
	return
	end
	'''
	import re
	conc = data['conc'] #_Species concentrations
	sp = data['sfc_pres'] #_Surface Pressure
	bsig = data['sig_b'] #_b-sigma coords
	nspec, nz, ny, nx = conc.shape 
	conc_kg = np.zeros((nspec, nz, ny, nx))
	rhodz = np.zeros((nz, ny, nx))
	cmap = np.zeros((ny)) #_Acrea / gridbox
	K = np.zeros((nspec))
	pi = np.pi
	r_e = 6370000.	#_Radius of the earth, m
	lats = data['lats']
	lons = data['lons']

	bn_air_stp = 2.65e16 #_air molecules / m3 (@stp)
	avagadro = 6.023e23 #_Avogadro's Number
	m_sulf = 32.0e-3 #_ kg / mole, sulfur

	for j in xrange(ny):
		#_Calculate area per gridbox (cmap)
		rprime = r_e * np.cos( lats[j] * pi / 180. )
		dxx = 2. * pi * rprime / float(nx) #_pi*d=circ/nx=arclen
		dyy = pi * r_e / float(ny)
		cmap[j] = dxx * dyy
	cmap = cmap.reshape(1,1,ny,1)

	for s in xrange(nspec):	
		if s <= 1:
			K[s] = 100. * bn_air_stp * m_sulf / avagadro
		elif s > 1 and s < nspec:
			K[s] = 100. * 1.0
		else:
			exit("ERROR: Invalid species index -"+str(s))
	K = K.reshape(nspec,1,1,1)
	bsig1 = []
	[bsig1.append(bsig[n+1]) for n in xrange(nz)] 
	bsig1 = np.array( bsig1 )
	bsig0 = bsig[0:nz].reshape(nz,1,1)
	bsig1 = bsig1.reshape(nz,1,1)
	sp = sp.reshape(1,ny,nx)
	rhodz = sp * (bsig0-bsig1) / 9.8
	
	if units == 'gridbox':	
		conc_kg = conc * rhodz * cmap * K
	elif units == 'm3':
		if rho == None: exit('ERROR: Density required for kg/m3')
		if rho.shape != (nz, ny, nx): rho = rho[:nz,:,:]
		conc_kg = conc * rho * K
	elif units == 'm2':
		conc_kg = conc * rhodz * K
	else:
		dbg( 'Invalid "units" entry ' + units, l=2 )

	data['conc'] = conc_kg
	return data

def query_conc( file_in ):
	''' 
	Reads in concentration metadata and prints to stdout
	file_in	: str, full path to _conc binary file
	'''
	import re
	conc = read_conc(file_in)
	dtg = lt.unique( conc.dtg_vald, unique=True )	

	#_Print main metadata
	print '\n'
	print '%25s : %-25s (%s)' % ('VALID_TIME',lt.human_date( dtg ), dtg )
	if hasattr( conc, 'naaps_version' ):
		#_The 'fhr' in standard conc metadata is not useful
		version = conc.naaps_version
		fhr = lt.unique( conc.fhr, unique=True )
		print '%25s : %-25s' % ('FORECAST_HOUR', fhr )

	shape = subset( conc, variable='conc', unique=True ).values.shape
	print '%25s : %-25s' % ( 'NS,NZ,NY,NX', shape )

	#_Loop over variables
	print '%25s : %-25s' % ( 'VARIABLES', '' )
	for v in lt.unique( conc.variable ):
		rec = subset( conc, variable=v )
		if not hasattr( rec, 'units' ): continue
		print '{0:28}'.format('') + \
			'{0:{fill}{align}49}'.format( v, fill='_', align='<' )
		print '%30s %-9s : %-25s' % ( '', 'units', rec.units )
		print '%30s %-9s : %-25s' % ( '', 'long_name', rec.long_name )

	#_METADATA FIRST SET (2011.02)
	try:
 	    if hasattr( conc, 'naaps_version' ):
		print '%25s : %-25s' % ( 'NAAPS_VERSION', version )
		print '%25s : %-25s' % ( 'N_OBS_TOTAL', conc.n_obs )
		tmp 	= conc.modis_version
		re_type = re.compile('(MODIS.*?v[\d.]+)')
		re_sat 	= re.compile('SAT\s+=\s+(\w+)')
		re_nobs = re.compile('N\s+=\s+(\d+)')
		re_vlab = re.compile('VLABEL\s+=\s+([\w\d_.-]+)')
		tmp_modis = re_type.findall( tmp )

		#_No modis metadat
	  	if len( tmp_modis ) == 0: 
			print '%25s : %-25s' % ( 'MODIS_TYPE', 'n/a' )
			print '%25s : %-25s' % ( 'SATELLITE', 'n/a' )
			print '%25s : %-25s' % ( 'N_OBS', 'n/a' )
			print '%25s : %-25s' % ( 'VLABEL', 'n/a' )
		else:
			print '%25s : %-25s' % ( 'MODIS_TYPE',
						',\n\t\t\t    '.join(tmp_modis))
			print '%25s : %-25s' % ( 'SATELLITE',
						', '.join(re_sat.findall(tmp)))
			print '%25s : %-25s' % ( 'N_OBS',
						', '.join(re_nobs.findall(tmp)))
			print '%25s : %-25s' % ( 'VLABEL',
						', '.join(re_vlab.findall(tmp)))
		tmp = conc.bitflag
		tf = { '1' : 'True' , '0' : 'False' , '2' : 'n/a'}
		print '%25s : %-25s' % ( 'CMORPH', tf[tmp[0:1]] )
		print '%25s : %-25s' % ( 'GINOUX_DUST', tf[tmp[1:2]] )
		print '%25s : %-25s' % ( 'SEDIMENTATION', tf[tmp[2:3]] )
		print '%25s : %-25s' % ( 'BITFLAG', tmp )

	    #_METADATA SECOND SET (2011.06)
	    if hasattr( conc, 'srcspec' ):
		tmp = conc.srcspec
		re_sat = re.compile( 'SAT\s+:\s+(\d+)' )
		re_inj = re.compile( 'INJ\s+:\s+([\d.+E]+)' )
		print '%25s : %-25s' % ( 'SMOKE_SATCODES',
					','.join(re_sat.findall(tmp)))
		print '%25s : %-25s' % ( 'SMOKE_INJECTION',
					','.join(re_inj.findall(tmp)))
		print '%25s : %-25s' % ( 'FLAMBE_VERSION',
					conc.flambe_version )
	except:
		dbg('numpy version not fully supported')

def diff( file0, file1 ):
	'''
	Opens two concentration files, finds the differences, outputs it 
	netcdf file
	'''
	import re
	#_Read data
	con = re.compile( 'conc$' )
	aod = re.compile( 'aod$' )
	if con.search( file0 ) and con.search( file1 ):
		f0 = read_conc( file0 )
		f1 = read_conc( file1 )
	elif aod.search( file0 ) and aod.search( file1 ):
		f0 = read_naapsaod( file0 )
		f1 = read_naapsaod( file1 ) 
	else:
		dbg(( 'Invalid files', file0, file1 ))
		return -1

	#_Get list of variables
	vars = lt.unique( f0.variable )
	#_Initialize output recarray
	diff = cl.model_object()

	#_Loop over vars, assume they have the same ones for now
	for v in vars:
		v0 = subset( f0, variable=v )
		v1 = subset( f1, variable=v )

		dimname = v0.dimname
		dimsize = v0.dimsize

		ensemble = v0.ensemble
		long_name = v0.long_name
		units = v0.units
		try:	#_Many things simply cannot be subtracted
			tmp = v0.values - v1.values
		except:
			dbg('cannot diff ' + v )
			#_Skip to next in cases such as species names
			continue

		attrn = v0.dimname
		attrv = []

		#_Add to recarray object with some dummy metadata
		if hasattr( v0.values, 'lat') \
		and hasattr( v1.values, 'lon'):
			for attrname in attrn:
				attr = v0.values.__getattribute__( attrname )
				attrv.append( attr ) 
			tmp = cl.var( tmp, attrv=attrv )

		#_Add difference into master array
		diff.resize( diff.size + 1 )  
		diff[-1] = ( tmp, 'DIFFS', 1, 1, 1, 'reg', v, ensemble, 
			dimname, dimsize, units, long_name )

	plot_default( diff )
	to_ncdf( diff, filename='diff.nc' )


def to_ncdf( data, filename='default.nc' ):
	'''
	Creates netcdf file when passed a recarray with these fields:
		values		(nd.array, data)
		long_name	str
		units		str
		dtg_init	str
		dtg_vald	str
		variable	str	
		model		str	
		dimname		tuple( str )
		dimsize		tuble( int )
	'''
	from netCDF4 import Dataset

	#_Get list of variables to write
	vars = lt.unique( data.variable )
	if len(vars) < data.size:
		raise ValueError, 'Currently not setup to handle multiple ' \
					+ 'instances of the same variable'
	dimname_list 	= data.dimname 
	dimsize_list	= data.dimsize 
	ntypes		= len( dimname_list )

	#_Loop over dimnames, which will have repeats, but I don't know
	# a clean way to match these up otherwise
	dimensions = {}
	for n in np.arange( ntypes ):
		dimnames = dimname_list[n]
		dimsizes = dimsize_list[n]
		ndims = len( dimnames )

		for i in np.arange( ndims ):
			name = dimnames[i]
			size = dimsizes[i]
			if name not in dimensions:
				dimensions[ name ] = size
			else:
				if size != dimensions[name]:
					raise TypeError( 'Size mismatch' )

	#_Create dimensions
	ncdf = Dataset( filename, 'w', format='NETCDF3_CLASSIC' )
	dbg( filename )
	for dim in dimensions:
		ncdf.createDimension( dim, dimensions[dim] )

	#_Create variables and write var attributes 
	for rec in data:
		cdf = ncdf.createVariable( rec.variable, 'f4', rec.dimname )
		cdf.model	= rec.model
		cdf.dtg_init	= rec.dtg_init
		cdf.dtg_vald	= rec.dtg_vald
		cdf.long_name	= rec.long_name
		cdf.units	= rec.units
		try:
			cdf[:] = rec.values
		except ValueError:
			#_This was specifically put in to add species as ncattr
			cdf[:] = np.arange( len( rec.values ) )
			try:
				val_str = ', '.join( rec.values.tolist() )
			except:
				val_str = str( rec.values )
			ncdf.setncattr( rec.variable, val_str )
	ncdf.close()
	lt.make_readable( filename )

def to_ncdf2( data, filename='default.nc' ):
	'''
	THIS FUNCTION IS A STOP GAP MEASURE.  NEEDS TO BE REPLACED!

	Creates netcdf file when passed a recarray with these fields:
		values		(nd.array, data)
		long_name	str
		units		str
		dtg_init	str
		dtg_vald	str
		variable	str	
		model		str	
		dimname		tuple( str )
		dimsize		tuble( int )
	This was made as a quick kludge to get multiple times into single file
	Should be used to dump ICAP netcdf, but isn't.

	TODO: Make different function This needs to be updated to handle i
	concentration metadata
	'''
	from netCDF4 import Dataset
	from time import gmtime, strftime

	#_Get list of variables to write
	vars = lt.unique( data.variable )
	dimname_list 	= data.dimname 
	dimsize_list	= data.dimsize 
	ntypes		= len( dimname_list )

	#_Loop over dimnames, which will have repeats, but I don't know
	# a clean way to match these up otherwise
	dimensions = {}
	for n in np.arange( ntypes ):
		dimnames = dimname_list[n]
		dimsizes = dimsize_list[n]
		ndims = len( dimnames )

		#_Loop over each possible dim
		for i in np.arange( ndims ):
			name = dimnames[i]
			size = dimsizes[i]
	
			#_Create dictionary of dim name => dimvalues 
			if name not in dimensions:
				dimensions[ name ] \
				= data.values[n].__getattribute__(name)

	epoch = lt.dtg2epoch( lt.unique( data.dtg_vald ))
	nt = len( epoch )

	#_Create dimensions
	dbg( filename )
	ncdf = Dataset( filename, 'w', format='NETCDF3_CLASSIC' )
	for dim in dimensions:
		#_Add dimension
		size = dimensions[ dim ].size
		ncdf.createDimension( dim, size )

		#_Add corresponding variable
		cdf = ncdf.createVariable( dim, 'f4', (dim,) )
		cdf[:] = dimensions[ dim ]

	#_Add dummy dimension 
	ncdf.createDimension( 'epoch', nt )
	cdf = ncdf.createVariable( 'epoch', 'f8', ('epoch') )
	if hasattr( data[0].values, 'member'):
		 ncdf.members = ','.join( data[0].values.member.tolist() )

	#_Add dimensions as variables

	#_Create variables and write var attributes 
	for variable in vars:
		sub_var = subset( data, variable=variable )
		sub_srt = sort_rec( sub_var, order=('dtg_vald',) )

		atn, atv = get_attr( sub_srt[0] )	
		sub_full = join_values( sub_srt, newattr='epoch',newattrv=epoch)

		atn.insert( 0, 'epoch' )
		atv.insert( 0, epoch )	#_used?

		cdf = ncdf.createVariable( variable, 'f4', atn )
		cdf.model	= lt.unique( sub_srt.model, unique=True ) 
		cdf.dtg_start	= lt.epoch2dtg( epoch[0] ) 
		cdf.dtg_end	= lt.epoch2dtg( epoch[-1] ) 
		cdf.long_name	= lt.unique( sub_srt.long_name, unique=True ) 
		cdf.units	= lt.unique( sub_srt.units, unique=True ) 
		cdf[:] = sub_full 

	ncdf.generated = strftime( "%Y-%m-%d %H:%M:%S UTC", gmtime() ) 

	ncdf.close()
	lt.make_readable( filename )

def read_metmod( file, plevs=True, rho=True, rh=False, pipe=None, ensemble=0,
	**kwargs ):
	'''	
	UPDATE TO WORK WITH RECARRAYS 
	USAGE: <recarray> = read_nogaps(<filename_str>)

	plevs	: bool,	Set to calculate full pressure levels
	rho	: bool, Set to calculate full density on all levels
	rh	: bool, Just don't set.
	ensemble: int,	Number of total ensemble members, set to 0 for det

	Ensemble should always be false.  When reading an ensemble forecast,
	this will be looped over and the master function will set True
	to the appropriate array.
	'''
	import re	
	dbg( file )

	#_Crude test based on filename
	model = 'NAVGEM' if re.search( 'nvg', file ) else 'NOGAPS'

	#_Initial data object
	metmod = cl.model_object()

	#_Open binary file
	f = open( file,'rb' )

	#_Read in metadata
	dtype	= np.dtype( [ 	('a0', 'i4'), 
					('dtg', 'i4'),
					('unk0', 'i4'),
				('a1', 'i4'),
				('b0', 'i4'),
					('dims', '3i4'),
				('b1', 'i4'),
				('c0', 'i4'),
					('unk1', 'i4'),
					('unk2', 'i4'),
				('c1', 'i4'),	] )
	data = lt.byte_order( np.fromfile( f, dtype, count=1 ))
	dtg_vald = str( data['dtg'][0] )
	dtg_init = dtg_vald #_Don't see this information in metadata
	fhr = lt.find_runlength( dtg_init, dtg_vald ) #_and this is also bunk
	ny, nx, nz = data['dims'][0]

	#_Read in coordinate data
	var_meta, var_order = lm.nog_vars( nz, ny, nx )
	mtype = [('scalar','f4'), ('offset','f4')] #_Used in ndarrays
	for v in var_order:
		dimsize = np.array( var_meta[v]['shape'] )	
		dtype = np.dtype( var_meta[v]['dtype'] )


		#_Make list of variable names.  Important when more than
		# one variable is in a single write
		names = dtype.names

		#_Read in record length, data, then closing record length
		size_exp = dimsize if dimsize.size == 1 else	\
			np.multiply.reduce( np.array( dimsize ) )

		ndim = dimsize.size
		#_Multi-dimensional variables are prefaced by scalars 
		# and offsets in an record of length 8
		if ndim == 1:
			#_Read in data 
			data = lt.rec_read( f, dtype )

		elif ndim == 2:	
			#_Read in metadata
			meta = lt.rec_read( f, mtype )
			scalar = meta['scalar']
			offset = meta['offset']

			#_Read in data and apply offset/scalar
			temp = lt.rec_read( f, dtype )
			data = {}
			data[v] = temp[v] * scalar + offset 
	
		elif ndim == 3:
			#_Each level has metadata/offset/scalar
			# initialize a 3d array and populate it
			d = np.zeros(( nz, ny, nx ))

			for k in np.arange( nz ):
				#_Read in metadata
				meta = lt.rec_read( f, mtype )
				scalar = meta['scalar']
				offset = meta['offset']

				#_Read in data and apply offset/scalar
				tmp = lt.rec_read( f, dtype )
				tmp = tmp[v] * scalar + offset
				d[k,:,:] = tmp.reshape( ny, nx )

			#_Put in dict to be added by name below 
			data[v] = d

		#_Add data to return object
		for name in names:
			dimname = var_meta[name]['dims']
			dimsize = var_meta[name]['shape']
			units = var_meta[name]['units']
			long_name = var_meta[name]['long_name']

			#_Add record and fill
			metmod.resize( metmod.size + 1 )
			tmp = data[name].reshape( dimsize )
			metmod[-1] = ( tmp, 'NOGAPS', dtg_init, dtg_vald,
				fhr, 'global', name, ensemble, dimname, dimsize,
				units, long_name )
	f.close()

	sigma_b		= subset( metmod, variable='sig_b' ).values[0]
	p_sfc		= subset( metmod, variable='sfc_pressure' ).values[0] \
			* 100	#_Convert hPa -> Pa
	p_top		= subset( metmod, variable='sig_a' ).values[0].min() 
	Rd 		= 287.058	#_Gas constant, dry air (J/kgK)
	Rv 		= 461.5		#_Gas constant, moist air (J/kgK)
	Lv		= 2.454e6	#_Latent heat of vaporization (J/kg)

	#_Calculate pressure levels
	if plevs or rho or rh:
		name = 'pressure'
		#_sigma = (p - p_sfc) / (p_top - p_sfc)
		#_p = sigma * (p_top - p_sfc) + p_sfc

		#_Get midpoint
		sigma_m	= sigma_b[:nz]-(sigma_b[:nz] - sigma_b[1:])/2	
		sigma_m3d = np.ones((nz,ny,nx)) * sigma_m.reshape(nz,1,1)
		pressure = sigma_m3d * (p_top - p_sfc) + p_sfc

		dimname = var_meta[name]['dims']
		dimsize = var_meta[name]['shape']
		units = var_meta[name]['units']
		long_name = var_meta[name]['long_name']

		#_Add record and fill
		metmod.resize( metmod.size + 1 )
		P = pressure.reshape( dimsize )
		PhPa = P / 100
		metmod[-1] = ( PhPa, model, dtg_init, dtg_vald,
			fhr, 'global', name, ensemble, dimname, dimsize, units,
			long_name )

	#_Calculate 3D air density
	if rho or rh:
		#_The actual calculations
		T = subset( metmod, variable='temperature', 
			unique=True ).values
		rho = ( P / T ) / Rd 	#_True value overwritten

		#_Get variable metadata
		name = 'rho'
		dimname = var_meta[name]['dims']
		dimsize = var_meta[name]['shape']
		units = var_meta[name]['units']
		long_name = var_meta[name]['long_name']

		#_Put into records arrawy
		metmod.resize( metmod.size + 1 )
		metmod[-1] = ( rho, model, dtg_init, dtg_vald,
			fhr, 'global', name, ensemble, dimname, dimsize, units,
			long_name )

	#_Calculate surface relative humidity
	if rh:	#_I don't trust this calculation, by the way.
		# It mixes SFC values with LOWEST SIGMA LEVEL values...so....
		Es = 6.11 * np.e**(Lv/Rv*(1/273.15 - 1/T))
		#_Claussius-Clapyron (sat vp)

		#_Surface specific humidity
		sfc_sh = subset( metmod, variable='spec_humidity', 
			unique=True ).values[0]
		sfc_rh = sfc_sh * p_sfc / (0.622 + sfc_sh) / Es 
		rhov = -sfc_sh * rho / (sfc_sh-1)
		vapor_p	= T * 4.612 * rhov 	#_vapor pressure
		vapor_ps = 6.112 * np.e**(6816 * ((1./273.15)-(1./T)) \
			+ 5.1309*np.log(273.15/T))

		#_Get variable metadata
		name = 'sfc_rh'
		dimname = var_meta[name]['dims']
		dimsize = var_meta[name]['shape']
		units = var_meta[name]['units']
		long_name = var_meta[name]['long_name']

		#_Put into records arrawy
		sfc_rh = (vapor_p / vapor_ps) * 100
		metmod.resize( metmod.size + 1 )
		metmod[-1] = ( sfc_rh, model, dtg_init, dtg_vald,
			fhr, 'global', name, ensemble, dimname, dimsize, units,
			long_name )
	#_Return data
	if pipe == None:
		return metmod 
	else:
		pipe.send( metmod )

def read_nogaps_thread( *args, **kwargs ):
	from multiprocessing import Process,Pipe
	groups	= lt.setup_groups( *args, **kwargs )
	data	= [None]*len(list) 
	for group in groups:
		length	= len(group)
 		thread	= [None]*length
 		p		= [None]*length
  		c		= [None]*length
		for i in range(length):
			file	= group[i]
			p[i],c[i]	= Pipe()
			thread[i]	= Process( 	target=read_nogaps,
							args=( file, c[i] ) )
			thread[i].start()

	for i in range(length):
        	data[i]	= p[i].recv() 
		thread[i].join()  

        return data

def read_ensfcst( dtg, path_ens=dir_prod, model='ENAAPS-NAV', fhr=120,
	pipe=None, modes=False, modes_only=False, **kwargs ):
	'''
	Reads netcdf file of AODs for ensemble forecasts.  
	If expected file is missing, attempts to generate ncdf file
	using ASCII (or individual model files for ICAP) input.
	
	dtg	: str, 	start dtg
	path	: str,	directory of MASINGAR files
	model	: str, 	name of ensemble model
	fcst_l	: int,	used to limit number of times read
	'''
###	def exitf( aod, pipe, attrv=None):
###		'''
###		function to cleanly exit optional pipes
###		'''
###		if pipe == None:
###			return aod
###		else:
###			pipe.send(( aod, list(attrv) ))
###			pipe.close()	

	from netCDF4 import Dataset,Variable

	modes = True if modes_only else modes
 
	file_netcdf = '_'.join(( model.lower(), dtg, 'aod.nc' ))
	file = '/'.join(( path_ens, model, dtg[:6], file_netcdf ))

	#_I do not like this kludge.  ONE WAY, MAN,  ONE WAY
	if model == 'ICAP':
		import libicap as li
		current_icap = lm.current_icap()
	else:
		current_icap = [] 

	#_If file has not been created, make a local ICAP file
	if not os.path.exists( file ):
		dbg(( file, 'is missing' ))
		dbg(( 'Creating', model, 'file', dtg ))
		try:
 		    if model == 'ENAAPS-NAV':
			raw = read_ensfcst_raw( dtg )
		    elif model == 'ICAP':
			raw = li.read_icapfcst_raw( dtg, remove_members=True )
		    elif model == 'ENAAPS-DART':
			raw = read_ensfcst_raw( dtg, model=model, 
				path_ens=path_ens )
		    else:
			raise ValueError, 'Unknown ensemble type ' + model
		except:
		    raise IOError, 'failed to read raw data' 
		#return lt.exit_pipe( -1, pipe )

		#_Write a single netcdf file for forecast for faster reader
		write_ensfcst( raw )
	
		if not os.path.exists( file ):
			raise IOError, 'failed to create file '+ file 
			return lt.exit_pipe( -1, pipe )
###			dbg( 'Failed to write ' + file )
###			raise RuntimeError, 'Failed to write' #return -1

	dbg((model, dtg, path_ens)) 
	dbg( file, l=7 )

	aod = cl.model_object()
	vars = lm.aod_vars()

	handle	= Dataset( file, mode='r', format='NETCDF3_CLASSIC' )

	null,nens,ny,nx = handle.variables['total_aod'].shape
	lons	= handle.variables['lon'][:]
	lats	= handle.variables['lat'][:]
	times	= handle.variables['time'][:]
	dtg_init= lt.epoch2dtg(times[0]) 		
	dtg_fhr	= lt.newdtg( dtg, fhr ) 	#_Don't always want all
	nt	= len(times) 	
	specs 	= [ s.lower() for s in mod_dict[model]['specs'] ]
	members = handle.members.split(',')

	#_Loop over each variable, store in dictionary
	for t in np.arange(nt):
		s = times[t]
	        dtg_loop = lt.epoch2dtg(s)

		modes_tmp = {}
		for spec in specs: 
			if spec == 'lat' or spec == 'lon' or spec == 'time': 
				continue
			long_name = vars[spec]['long_name']
			tmp = handle.variables[spec]
			if dtg_loop > dtg_fhr: break

			#_mask missing variables
			tmp = handle.variables[spec][t,:,:,:].copy()

			#_for ICAP, older files will not have more recent
			# additions to the consensus. If reading an old file,
			# append masked arrays to keep them the same size
			nens_expected 	= len( current_icap )
			nens_infile 	= len( members )
			mem_loc 	= members[:]
			if model == 'ICAP' and ( nens_expected != nens_infile ):
				nan_arr = np.empty(( 1, ny, nx ))
				nan_arr[:] = -9999. 
				for mem in current_icap:
					#_don't append if present
					if mem in members: continue
					dbg( 'appending nan array', l=7 )

					#_This assumes mem axis is 0, as below
					tmp = np.append( tmp, nan_arr, axis=0 )	
					mem_loc.append( mem )

			#_make sure order of models is correct
			if model == 'ICAP':	
				#_check order of members in file array vs exp
				order = []
				##for mem in mem_loc:
				for mem in current_icap:
				##	idx = current_icap.index( mem )
					idx = mem_loc.index( mem )
					order.append( idx )

				#_Reorder array to match members dim, if nec
				tmp = tmp[order]

				#_Set members to current_icap()
				mem_loc = current_icap[:]
			else:
				mem_loc = members[:]

			#_if members have same name as model, append
			# differentiator. this will likely break everything.
			for member in mem_loc:
				if member != model: continue

				en = mem_loc.index( member ) 
				e2 = str(en).zfill(2)
				mem_loc[en] = member + '_m' + e2
				
			#_put in variable class and mask 
			attrv = (mem_loc, lats, lons,)
  	        	data = cl.var( tmp, attrv=attrv )
			nens_p = check_members( data ) 

			dimname = list( vars[spec]['dims'] )
			dimname.insert( 0, 'member' )
			dimname = tuple( dimname )
			dimsize = data.shape

			#_get current forecast hour	
			vhr = lt.find_runlength( dtg_init, dtg_loop ) / 3600

			aod.resize( aod.size + 1 )
			aod[-1] = ( data, model, dtg_init, dtg_loop, vhr, 
				'global', spec, nens_p, dimname, dimsize,
				'', long_name )
			if modes: modes_tmp[spec] = aod[-1].copy()
		else:
		    if modes:
			sulf = modes_tmp['sulfate_aod'].values
			smok = modes_tmp['smoke_aod'].values
			dust = modes_tmp['dust_aod'].values
			salt = modes_tmp['seasalt_aod'].values
			
			'''need to skip days where one or the other missing'''
			aod.resize( aod.size + 1 )
###			fine = lt.masked_sum( sulf, smok ) #_allows one to be missing
			fine = sulf + smok #_requires both
			fine = cl.var( fine, attrv=attrv )
			units = vars['fine_aod']['units']
			long_name = vars['fine_aod']['long_name']
			nens = check_members( fine ) 
			aod[-1] = ( fine, model, dtg_init, dtg_loop, vhr, 
				'global', 'fine_aod', nens, dimname, dimsize,
				'', long_name )

			aod.resize( aod.size + 1 )
			coarse = lt.masked_sum( dust, salt ) #_allows salt to be missing	
###			coarse = dust + salt 
			coarse = cl.var( coarse, attrv=attrv )
			units = vars['coarse_aod']['units']
			long_name = vars['coarse_aod']['long_name']
			nens = check_members( coarse ) 
			aod[-1] = ( coarse, model, dtg_init, dtg_loop, vhr, 
				'global', 'coarse_aod', nens, dimname, dimsize,
				'', long_name )
	handle.close()

	if modes_only: 
		aod = subset(aod,variable=['fine_aod','coarse_aod','total_aod'])


	#_mask where not expected
	mask_members( aod )

	#_if ICAP and strict, filter out records with 
	if model == 'ICAP': aod = li.filter( aod, **kwargs )
	if aod.size > 0: attrn, attrv = get_attr( aod[0] )
	return lt.exit_pipe( aod, pipe, attrv=attrv )

def read_ensanal( dtg_start, dtg_end, fcst_finc=24, finc=6, **kwargs ):
        '''
        TODO:   Write Help
                Generalize to all ensembles
                Allow for reading of other fcst hours?
        '''
	dbg((dtg_start, dtg_end ))
        #_Which forecast hours do we want 
        fhrs = range( finc, fcst_finc+finc, finc ) #_[6, 12, 18, 24]

	if 'fhr' not in kwargs:	kwargs.update({ 'fhr' : fcst_finc })

        #_Initialize recarray
        aod = cl.model_object()

        dtg = dtg_start
        while dtg < dtg_end:
		try:
                        anal = read_ensfcst( dtg, **kwargs )
                        anal = subset( anal, fhr=fhrs )
                        aod = lt.merge(( aod, anal ))
		except:
                        dbg(( 'error reading/merging ensemble analysis', dtg ))

                #_Incriment date
	  	dtg = lt.newdtg( dtg, fcst_finc )

        return aod


def read_ensfcst_raw( dtg, nens=20, path_ens=dir_prod, model='ENAAPS-NAV',
	**kwargs):
	'''
	Read ensemble forecast data from netcdf or _aod files
	ensembles var deprecated	
	'''
	from netCDF4 import Dataset,Variable
	dbg( model + dtg )
	pre = cl.model_object()
	vars = lm.aod_vars()

	#_Get list of expected species
	species = mod_dict[model]['specs']

	#_Generate input directory
	path_in = path_ens + '/' + model

	#_Loop over ensemble members and read in forecasts
	for e in np.arange(nens) + 1:
		if model == 'ENAAPS-NAV':
			tmp = read_naapsfcst( dtg, path_nva=path_in, e=e,
				model=model, **kwargs )
		elif model == 'ENAAPS-DART':
			''' UPDATE THIS TO WORK'''
			pass
###			tmp = read_dartfcst( dtg, e=e, model=name )
		else:
			raise ValueError, 'Ensemble type not setup ' + model 
		pre = lt.merge(( pre, tmp ))

	#_Loop over variables and valid times and join along member axis
	aod = cl.model_object()
	for variable in lt.unique( pre.variable ):
		tmp_var = subset( pre, variable=variable )

		#_Get list of dtgs, loop over and add to records
		dtgs = lt.unique( tmp_var.dtg_vald )
		for dtg_loop in dtgs:
			tmp = subset( tmp_var, dtg_vald=dtg_loop )
			vhr = tmp[0].fhr	
			long_name = tmp[0].long_name
			
			attr = 'member'
			value = [model]*nens
			attrn, attrv = get_attr( tmp[0] )
			attrn.insert( 0, attr )
			attrn = tuple( attrn )
			attrv.insert( 0, value )
			data = join_values( tmp, newattr=attr, newattrv=value )
			data = cl.var( data, attrn=attrn, attrv=attrv )
			dimsize = data.shape
			aod.resize( aod.size + 1 )
			aod[-1] = ( data, model, dtg, dtg_loop, vhr, 'global',
				variable, nens, attrn, dimsize, '', long_name)	
	return aod
	
def read_fakefcst( nens=5, nx=360, ny=180, nt=3, model='ICAP', fhrs=[6,12],
	dtg='2012010100', variables=['fine_aod','coarse_aod','total_aod'] ):
	''' Static 4d fields to test geomu/geostd '''
	records = cl.model_object()

	lats = np.linspace(-89.5, 89.5, ny)
	lons = np.linspace(-179.5, 179.5, nx)
	mems = ['NAAPS','GEOS5','MACC','MASINGAR','NGAC']

	attrn = ('member','lat','lon')
	attrv = (mems,lats,lons)

	n2 = ny*nx

###	flip = 1 
	for variable in variables:
	    for t in xrange(nt):
	    	dtg_init = lt.newdtg( dtg, t*24 )
		for fhr in fhrs:
			dtg_vald = lt.newdtg( dtg_init, fhr )

			values = np.ones((nens,ny,nx))
			for e in xrange( nens ):
###				values[e] = np.arange(n2)[::flip].reshape(ny,nx)
				values[e] += t

			#_values = np.ones((nens,ny,nx))
			values = cl.var( values, attrv=attrv )

			records.resize( records.size + 1 )
			records[-1] = ( values, model, dtg_init, dtg_vald,
				fhr, 'global', variable, nens, attrn,
				values.shape, '', variable )
###	    	flip *= -1
	return records

def write_ensfcst( data, path_out=dir_prod, **kwargs ):
	'''
	Creates netcdf file when passed a recarray with these fields:
		values		(nd.array, data)
		long_name	str
		units		str
		dtg_init	str
		dtg_vald	str
		variable	str	
		model		str	
		dimname		tuple( str )
		dimsize		tuble( int )

	The ensemble forecasts take up ridiculous amounts of room,
	so this puts AOD data into single netcdf files
	'''
	from netCDF4 import Dataset
	from time import gmtime, strftime

	#_Make sure only one model passed
	model = lt.unique( data.model )
	if len( model ) != 1:
		raise RuntimeError, 'Can only handle one model per file'	
	model = model[0]

	#_Get list of variables to write
	vars 		= lt.unique( data.variable )
	dimname_list 	= data.dimname 
	dimsize_list	= data.dimsize 
	ntypes		= len( dimname_list )

	#_Get 
	dtg_init = lt.unique( data.dtg_init )
	if len( dtg_init ) != 1 or data.size == 0: 
		raise RuntimeError, 'Only one forecast per file'
	dtg_init = dtg_init[0]

	#_Loop over dimnames, which will have repeats, but I don't know
	# a clean way to match these up otherwise
	dimensions = {}
	for n in np.arange( ntypes ):
		dimnames = dimname_list[n]
		dimsizes = dimsize_list[n]
		ndims = len( dimnames )

		for i in np.arange( ndims ):
			dimname = dimnames[i]
			size = dimsizes[i]
			if dimname not in dimensions:
				dimensions[ dimname ] = size
			else:
				if size != dimensions[dimname]:
					raise TypeError( 'Size mismatch' )

	#_Filename
	dir_out = path_out + '/' + model.upper() + '/' + dtg_init[:6]
	lt.mkdir_p( dir_out )
	filename = dir_out + '/' + model.lower() + '_' + dtg_init + '_aod.nc'

	#_Create dimensions
	ncdf = Dataset( filename, 'w', format='NETCDF3_CLASSIC' )
	dbg( filename )
	for dim in dimensions:
		ncdf.createDimension( dim, dimensions[dim] )

	#_Create time dimension
	ref_var = subset( data, variable='total_aod' )
	nt = ref_var.size
	ncdf.createDimension( 'time', nt )

	#_Write dimensional variables
	cdf = ncdf.createVariable( 'lat', 'f4', ('lat',) )
	cdf[:] = ref_var.values[0].lat 
	cdf = ncdf.createVariable( 'lon', 'f4', ('lon',) )
	cdf[:] = ref_var.values[0].lon 
	cdf = ncdf.createVariable( 'time', 'f8', ('time',) )
	dtgs = lt.unique ( data.dtg_vald )
	cdf[:] = lt.dtg2epoch( dtgs )

	#_Create variables and write var attributes 
	first = True
	for variable in vars:
		#_Subset by variable
		var_data = subset( data, variable=variable ) 

		#_Get attributes
		model = lt.unique( var_data.model )
		dtg_init = lt.unique( var_data.dtg_init )
		long_name = lt.unique( var_data.long_name )
		units = lt.unique( var_data.units )
		members = var_data.values[0].member

		#_Stack into single array and sort by dtg_vald if multiple vals
		# and not a dimension
		attrn, attrv = get_attr( var_data[0] )
		var_srt = sort_rec( var_data, order=('dtg_vald',) )
		epoch = lt.dtg2epoch( var_srt.dtg_vald )
		attrn = list(attrn)
		attrv = list(attrv)
		attrn.insert( 0, 'time' )
		attrv.insert( 0, epoch )
		var_stk = join_values( var_srt, newattr='epoch', newattrv=epoch)

		#_Write variables
		cdf = ncdf.createVariable( variable, 'f4', attrn,
			fill_value=-9999.)
		cdf.model	= model
		cdf.dtg_init	= dtg_init
		cdf.long_name	= long_name
		cdf.units	= units
		cdf[:] 		= var_stk

		first = False

	#_Add global variables. Use last aod vars meta data... eerrrooorrr
	ncdf.members = ','.join( members )
	ncdf.dtg_init = dtg_init
	ncdf.generated = strftime( "%Y-%m-%d %H:%M:%S", gmtime() ) 
	ncdf.close()
	lt.make_readable( filename )

def read_naapsanal( dtg_s, dtg_e, finc=6, path=dir_prod + '/NVA', **kwargs ):
	'''
	Reads in NAAPS analysis files for prescribed period

	dtg_s	: str, DTG start of period
	dtg_e	: str, DTG end of period
	finc	: int, time between analysis files
	path	: str, equivalent to DIR_WORK in NVA
	'''
	from multiprocessing import Process,Pipe
	dtgs = []
	dtg_loop = dtg_s
	while dtg_loop <= dtg_e:
		dtgs.append(dtg_loop)
		dtg_loop = lt.newdtg( dtg_loop, finc )

	model		= 'NAAPS'
	groups		= lt.setup_groups( dtgs, **kwargs )
	species		= mod_dict[model]['specs']
	aod		= cl.model_object()
	for group in groups:
          	l           = len(group)
          	thread      = [None]*l
          	p           = [None]*l 
          	c           = [None]*l
          	for i in range(l):  #-Pipe init, open connections 
            		dtg_loop	= group[i] 
            		file 		= path + '/NAAPSAOD/' + dtg_loop[:6] \
                			+ '/' + dtg_loop + '_aod' 
            		p[i],c[i] 	= Pipe()
            		args      	= ( file, c[i] )
            		thread[i] 	= Process( target=read_naapsaod, 
					args=args )
            		thread[i].start()
         
		for i in range(l):  #-Pipe recv, close connections
      			tmp, [lat, lon] = p[i].recv()

			for rec in tmp:
            			if (len(tmp) == 2 and tmp[0] == 1): 
					err = 'Failed to read ' + tmp[1] 
              				raise RuntimeError, err 

				dtg_tmp = rec.dtg_vald

				#_Get metadata
				spec = rec.variable 
				data = cl.var( rec.values, attrv=[lat,lon] )
				dimname = rec.dimname
				dimsize = rec.dimsize 
				long_name = rec.long_name
				dtg_init = rec.dtg_vald
				dtg_vald = rec.dtg_vald
				vhr = 0
	
				#_Add to recarray
				aod.resize( aod.size + 1 )
				aod[-1] = ( data, model, dtg_init, 
					dtg_vald, vhr, 'global', 
					spec, 0, dimname, dimsize,
					'', long_name )

       		     	thread[i].join()

	if lt.unique( aod.fhr )[0] != 0:
		raise ValueError, 'No analysis data read', lt.unique( aod.fhr )

	return aod

def read_naapsfcst( dtg, fhr=120, finc=6, fstrt=0, nt=None, 
	path_nva=dir_prod + '/NVA', e=None, model='NAAPS', **kwargs ):
	"""
	Defaults are those being used for ICAP procedures.
	aod[dtg][spec]	= read_naapsfcst(dtg,**keywords)
		dtg	= string, 	'2011050100'	Init date
		fhr	= integer,	120		How far past dtg to look
		finc	= integer,	6		Timestep up to fhr
		fstrt	= integer,	0		Number of hrs from dtg 
							to start. Almost always
							leave at zero
		nt	= integer,	None		If specific number of 
							timesteps
							required instead of 
							proceeding to fhr
		path	= string,	path_to_files	Where them files at?
	"""
	dbg(( dtg, model, 'ENS', str(e) ))
	from multiprocessing import Process,Pipe
	dtg_init	= lt.newdtg(dtg,-fstrt)
	nt		= (fhr - fstrt) / finc + 1 if nt == None else nt
	groups		= lt.setup_groups( np.arange(nt), **kwargs )
	aod		= cl.model_object()
	vars		= lm.aod_vars()
	species		= mod_dict['NAAPS']['specs']
	#_If ensemble, create member number string*2
	if e != None:
		e2 = str(e).zfill(2)
		estr = 'E' + e2 + '/'
		model += '_m' + e2
	else:
		estr = ''

	for group in groups:
          	l           = len(group)
          	thread      = [None]*l
          	p           = [None]*l 
          	c           = [None]*l
          	for i in range(l):  #-Pipe init, open connections 
            		t 		= group[i]
            		vhr		= finc*t + fstrt
            		dtg_loop	= lt.newdtg(dtg_init,vhr)
            		file 		= '/'.join(( path_nva, 'NAAPSAOD',
					'FCAST_' + dtg_init + '_f' 
					+ str(fhr).zfill(3), estr, dtg_loop[:6],
					dtg_loop + '_aod' ))
            		p[i],c[i] 	= Pipe()
            		args      	= ( file, )
			kwargs.update( { 'pipe' : c[i] } )
            		thread[i] 	= Process( target=read_naapsaod, \
					args=args, kwargs=kwargs)
            		thread[i].start()
         
		for i in range(l):  #-Pipe recv, close connections
       			t         = group[i]
            		vhr       = finc*t + fstrt
           		dtg_loop  = lt.newdtg( dtg_init, vhr )
			try:
	      			tmp, [lat, lon] = p[i].recv()
			except:
              			raise RuntimeError, 'Failed to read '+ dtg_loop 

			#_Loop over each record (should be npspec)
			for rec in tmp:
				dtg_tmp = rec.dtg_vald
				if dtg_tmp != dtg_loop: 
					dbg(( 'DTG mismatch', dtg_tmp, 
						dtg_loop ), l=3 )

				#_Get metadata
				spec = rec.variable 
				data = cl.var(rec.values, attrv=(lat, lon,))
				dimname = rec.dimname
				dimsize = rec.dimsize 
				long_name = rec.long_name

				aod.resize( aod.size + 1 )
				aod[-1] = ( data, model, dtg_init, 
					dtg_loop, vhr, 'global', 
					spec, 0, dimname, dimsize,
					'', long_name )
       		     	thread[i].join()
	return aod

def read_naapsaod( file, pipe=None, ncdf=False, **kwargs ):
	'''
	Read NAAPS AOD ASCII file
	aod = read_naapsaod( filename, pipe=<multiprocessing obj> )

	As long as the order of species is consistent, this works.  If total
	gets pushed back further, we run into problems.
	'''
	import re
	if ncdf or file[-3:] == '.nc': 
		return read_naapsaod_ncdf( file, pipe=pipe, **kwargs )

	model = 'NAAPS'
	species_all = [ s.lower() for s in mod_dict[model]['specs'] ]
	dbg( file, l=2 )
        if not os.path.exists( file ) and       \
        not os.path.exists(file + '.gz'):               #-File check
                dbg(( file, 'is missing '))
                return lt.exit_pipe( False, pipe=pipe )
        elif os.path.exists(file + '.gz'):      #-File check
                dbg( 'GUNZIPPING', l=3 )
                test = os.system('gunzip ' + file + '.gz')
                if test != 0:
                        dbg(( file, '.gz is missing' ))
                        return lt.exit_pipe( False, pipe=pipe )

	#_Open file, get coordinates and metadata
	hand = open(file,'r')

	#_Read in lines
	lines = hand.readlines()
	hand.close()
	lines = [ line.split() for line in lines ]
	nrecs = len(lines)

        #_read in all values for lat
        lats = [ float(n[1]) for n in lines ]

        #_lats repeat in the first column, the number of repeated times
        # will equal the number of meridians
        nx = lats.count(lats[0])

        #_the number of parallels can then be derived   
        ny = nrecs / nx

        #_reduce lats to non-repetition and read in lons
        lats = lats[:ny]
        lons = [ float(n[0]) for n in lines[::ny] ]

	#_Convert to more sensible array structures
	lats = np.array( lats )
	lons = np.array( lons )
 
	#_oddly...convert data to float.
	for line in lines:
		for i in xrange(len(line)): line[i] = float(line[i])

	nspec = len(lines[0]) - 2 
	species	= species_all[:nspec]
	tmp = np.array(lines).reshape( nx, ny,nspec+2 )
	dtg_vald = re.search( '(\d{10})_aod$', file ).group(1)
	dtg_init = dtg_vald[:] # Gets set appropriately if known 

	#_initialize object to return
	aod = cl.model_object()
	vars = lm.aod_vars()
	vhr = 0

	#_loop over species, put data into aod[spec] = np.ndarray form
	for s in xrange( nspec ):
		spec = species[s].lower()
		long_name = vars[spec]['long_name']

		data = cl.var( tmp[:,:,s+2].transpose(), attrv=(lats,lons,) )
		dimsize = data.shape  
		dimname = vars[spec]['dims'] 

		#_add to record array and substitute dummy dtg_init (unknown)
		aod.resize( aod.size + 1 )
		aod[-1] = ( data, model, dtg_init, dtg_vald, vhr, 
			'global', spec, 0, dimname, dimsize, '', long_name )

	return lt.exit_pipe( aod, pipe=pipe, attrv=[lats,lons] )

def read_naapsaod_ncdf( file, pipe=None, lamb=0.550, **kwargs ):
	'''
	Read NAAPS AOD NCDF file
	aod = read_naapsaod( filename, pipe=<multiprocessing obj> )

	As long as the order of species is consistent, this works.  If total
	gets pushed back further, we run into problems.
	'''
	from netCDF4 import Dataset
	model = 'NAAPS'
	species_all = [ s.lower() for s in mod_dict[model]['specs'] ]
	dbg( file, l=2 )
        if not os.path.exists( file ) and       \
        not os.path.exists(file + '.gz'):               #-File check
                dbg(( file, 'is missing '))
                return lt.exit_pipe( False, pipe=pipe )
        elif os.path.exists(file + '.gz'):      #-File check
                dbg( 'GUNZIPPING', l=3 )
                test = os.system('gunzip ' + file + '.gz')
                if test != 0:
                        dbg(( file, '.gz is missing' ))
                        return lt.exit_pipe( False, pipe=pipe )

	#_open file, get coordinates and metadata
	hand = Dataset( file, mode='r', format='NETCDF3_CLASSIC' )

	#_Convert to more sensible array structures
	lats = hand.variables['lat'][:] 
	lons = hand.variables['lon'][:] 
	ny = len( lats ) 
	nx = len( lons ) 

	#_create species list
	species = [ s.lower() for s in mod_dict[model]['specs'] ]
	nspec = len( species )

	#_get dtg_vald from metadata
	dtg_vald = hand.dtg_vald 
	dtg_init = dtg_vald[:] 	# eventually set this up to actually 
				# be in there. read_naapsfcst attempts to set

	#_initialize object to return
	aod = cl.model_object()
	vars = lm.aod_vars()
	vhr = 0

	#_loop over species, put data into aod[spec] = np.ndarray form
	for s in xrange( nspec ):
		spec = species[s].lower()
	
		#_read in values, put in var class
###		if spec == 'total_aod':
###			idx = np.where( hand.variables['lambda'][:] == lamb )
###			idx = idx[0][0]
###			values = hand.variables[ spec ][idx]	
###		else:
		values = hand.variables[ spec ][:]
		data = cl.var( values, attrv=(lats,lons,) )

		#_create metadata
		dimsize = data.shape  
		dimname = vars[spec]['dims'] 
		long_name = vars[spec]['long_name']

		#_add to record array and substitute dummy dtg_init (unknown)
		aod.resize( aod.size + 1 )
		aod[-1] = ( data, model, dtg_init, dtg_vald, vhr, 
			'global', spec, 0, dimname, dimsize, '', long_name )

	return lt.exit_pipe( aod, pipe=pipe, attrv=[lats,lons] )

##	if pipe == None:
##		return aod
##	else:
##		pipe.send( (aod, [lats,lons]) )
##		pipe.close()

#############################################################################_80
######_STATS_###################################################################
################################################################################

def generate_period_comps( dtg_start, dtg_end, fhrs=None, fcst_finc=24, 
	summation=True, product=True, components=True, aeronet=None, **kwargs ):
	'''
	Uses model and aeronet data to produce components to be
	converted into gross statistics by other modules.

	Namely, write_period_statistics()
	Aeronet mode can be controlled by prefiltering aeronet = values
	To reduce what is scored, use subset() to filter records

	dtg_start	: str*10, 	Start of Scoring Period
	dtg_end		: str*10, 	End of Scoring Period
	fcst_inc	: int,		Hours between model initializations
	fhrs		: list||str,	Fhrs to score

	summation	: bool,		Return summation of all present data
	product		: bool,		Return product of all present data
	components	: bool,		Return statistical residuals, etc for
					all present data
	'''
	from multiprocessing import Process,Pipe
	import libaeronet as la
	import libicap as li

	dbg(( dtg_start, dtg_end ))

	#_don't read aeronet if only doing summations
	if components: 
		#_if aeronet not passed, read it in
		# this variable is the FULL period aeronet.  Subsets are read
		# for individual scoring until I implement a temporal range 
		# sub-set for aeronet data taht's smart
		if aeronet == None:
			aeronet = la.read_aeronet( [dtg_start, dtg_end], 
					range=True )

		#_Filter out points not in data or ICAP dictionary
	        points_avail = lt.unique( aeronet.code )
		points = lt.intersection( [ pnt_dict, points_avail ] )
		aeronet = subset( aeronet, code=points )
		kwargs.update({ 'aeronet' : aeronet })

	#_Get plist of dtgs to read in
	dtgs = []
	dtg_loop = dtg_start
	while dtg_loop <= dtg_end:	
		dtgs.append( dtg_loop )
		dtg_loop = lt.newdtg( dtg_loop, fcst_finc )

	#_By point 
	components = cl.aeronet_comps() 
	prod = cl.cum_mean()
	summ = cl.cum_mean()	#_A funky recarray 
				#_Running mean only when there is AERONET data
				#_Seems useless as aeronet data is not always
				# off and on synchronously between stations
	#_READ_IN_DATA__________________________________________________________
	#_Scoring is now done on the fly, which requires a whole lot more
	# disk I/O.  Reading is now threaded, as it happens so damn many times
	groups = lt.setup_groups( dtgs, **kwargs )
	for group in groups:
		l 	= len( group )	#_Can't use nproc if grp%np!=0	
		t 	= [None]*l	#_List of processing threads
		pi	= [None]*l	#_List of pipes in
		po	= [None]*l	#_List of pipes out

		comps = {}
		fcst_recs = {}

		#_reading of ensemble data
		for n in np.arange( l ):
			dtg = group[n]	#_Name of point

			pi[n], po[n] = Pipe()
			args = ( dtg, )
			kwargs.update( { 'pipe'	: po[n] } )
			t[n] = Process( target=read_ensfcst,
					args=args, 
					kwargs=kwargs )
			try:
				t[n].start()
			except:
				print 'KLUDGE'
		
		#_collect multi-processing threads of forecast data 
		for n in np.arange( l ):
			dtg = group[n]

			#_store in temporary dictionary. Missing data loops.
			tmp_fcst, attrv = pi[n].recv()
			if type( tmp_fcst ) != np.recarray: 
				dbg(('error reading model', dtg), l=3)
				continue

			#_i am not good at python... my dim attributes aren't
			# being transferred through the pipe (since it pickles
			# them), so I reset them here. Can be fixed with 
			# Manager?
			for rec in tmp_fcst:
				rec.values = cl.var( rec.values, attrv=attrv )

			#_if ICAP and strict, filter out records with 
			tmp_fcst = li.filter( tmp_fcst, **kwargs )

			#_filter by FHR
			tmp_fcst = subset( tmp_fcst, fhr=fhrs )

			#_If all records weren't filtered out, continue 
			if tmp_fcst.size > 0: fcst_recs[n] = tmp_fcst	
			else: dbg(( 'no model records', dtg ), l=3 )
			t[n].join( 60 )

		#_this loop is done separately so that the rest can be
		# skipped if only summation/means desired
		for n in fcst_recs:
			#_keep running summation of aod data
			summ = running_summation( fcst_recs[n], summ, **kwargs )
			prod = running_product( fcst_recs[n], prod, **kwargs )

		#_skip aeronet scoring
		if not components: continue 

		#_loop over gathered data and thread to scoring function
		for n in fcst_recs:
			dtg = group[n]
			if type(fcst_recs[n]) != np.recarray: 
				dbg(('non-record array',type(fcst_recs[n])),l=3)
				continue
	
			#_if data was read in, score it and add to component
			pi[n], po[n] = Pipe()
			args = ( fcst_recs[n], )
			kwargs.update({	'pipe' : po[n] })
			t[n] = Process( target=la.generate_components,
					args=args, 
					kwargs=kwargs )

			t[n].start()

		#_merge aeronet_comp() records after all are computed 
		for n in fcst_recs:
			#_get residuals through pipe
			comps[n] = pi[n].recv()
			t[n].join( 3600 )

			#_merge components 
			components = lt.merge(( components, comps[n] ))

	if len( components ) == 0: dbg( 'No comps calculated, check data', l=3 )

	#_Add dtg_init dtg_vald for later
	summ.__setattr__('dtg_init', dtg_start)	
	summ.__setattr__('dtg_vald', dtg_start)	

	return components, summ, prod

def generate_period_components( dtg_start, dtg_end, fhrs=None, fcst_finc=24, 
	summation=True, product=True, components=True, aeronet=None, **kwargs ):
	'''
	Uses model and aeronet data to produce components to be
	converted into gross statistics by other modules.

	Namely, write_period_statistics()
	Aeronet mode can be controlled by prefiltering aeronet = values
	To reduce what is scored, use subset() to filter records

	dtg_start	: str*10, 	Start of Scoring Period
	dtg_end		: str*10, 	End of Scoring Period
	fcst_inc	: int,		Hours between model initializations
	fhrs		: list||str,	Fhrs to score

	summation	: bool,		Return summation of all present data
	product		: bool,		Return product of all present data
	components	: bool,		Return statistical residuals, etc for
					all present data
	'''
	from multiprocessing import Process,Pipe
	import libaeronet as la
	import libicap as li

	dbg(( dtg_start, dtg_end ))

	#_don't read aeronet if only doing summations
	if components: 
		#_if aeronet not passed, read it in
		# this variable is the FULL period aeronet.  Subsets are read
		# for individual scoring until I implement a temporal range 
		# sub-set for aeronet data taht's smart
		if aeronet == None:
			aeronet = la.read_aeronet( [dtg_start, dtg_end], 
					range=True )

		#_Filter out points not in data or ICAP dictionary
	        points_avail = lt.unique( aeronet.code )
		points = lt.intersection( [ pnt_dict, points_avail ] )
		aeronet = subset( aeronet, code=points )
		kwargs.update({ 'aeronet' : aeronet })

	#_Get plist of dtgs to read in
	dtgs = []
	dtg_loop = dtg_start
	while dtg_loop <= dtg_end:	
		dtgs.append( dtg_loop )
		dtg_loop = lt.newdtg( dtg_loop, fcst_finc )

	#_By point 
	components = cl.aeronet_comps() 
	prod = cl.cum_mean()
	summ = cl.cum_mean()	#_A funky recarray 
				#_Running mean only when there is AERONET data
				#_Seems useless as aeronet data is not always
				# off and on synchronously between stations
	#_READ_IN_DATA__________________________________________________________
	#_Scoring is now done on the fly, which requires a whole lot more
	# disk I/O.  Reading is now threaded, as it happens so damn many times
	groups = lt.setup_groups( dtgs, **kwargs )
	for group in groups:
		l 	= len( group )	#_Can't use nproc if grp%np!=0	
		t 	= [None]*l	#_List of processing threads
		pi	= [None]*l	#_List of pipes in
		po	= [None]*l	#_List of pipes out

		comps = {}
		fcst_recs = {}

		#_reading of ensemble data
		for n in np.arange( l ):
			dtg = group[n]	#_Name of point

			pi[n], po[n] = Pipe()
			args = ( dtg, )
			kwargs.update( { 'pipe'	: po[n] } )
			t[n] = Process( target=read_ensfcst,
					args=args, 
					kwargs=kwargs )
			try:
				t[n].start()
			except:
				print 'KLUDGE'
		
		#_collect multi-processing threads of forecast data 
		for n in np.arange( l ):
			dtg = group[n]

			#_store in temporary dictionary. Missing data loops.
			tmp_fcst, attrv = pi[n].recv()
			if type( tmp_fcst ) != np.recarray: 
				dbg(('error reading model', dtg), l=3)
				continue

			#_i am not good at python... my dim attributes aren't
			# being transferred through the pipe (since it pickles
			# them), so I reset them here. Can be fixed with 
			# Manager?
			for rec in tmp_fcst:
				rec.values = cl.var( rec.values, attrv=attrv )

			#_if ICAP and strict, filter out records with 
			tmp_fcst = li.filter( tmp_fcst, **kwargs )

			#_filter by FHR
			tmp_fcst = subset( tmp_fcst, fhr=fhrs )

			#_If all records weren't filtered out, continue 
			if tmp_fcst.size > 0: fcst_recs[n] = tmp_fcst	
			else: dbg(( 'no model records', dtg ), l=3 )
			t[n].join( 60 )

		#_this loop is done separately so that the rest can be
		# skipped if only summation/means desired
		for n in fcst_recs:
			#_keep running summation of aod data
			summ = running_summation( fcst_recs[n], summ, **kwargs )
			prod = running_product( fcst_recs[n], prod, **kwargs )

		#_skip aeronet scoring
		if not components: continue 

		#_loop over gathered data and thread to scoring function
		for n in fcst_recs:
			dtg = group[n]
			if type(fcst_recs[n]) != np.recarray: 
				dbg(('non-record array',type(fcst_recs[n])),l=3)
				continue
	
			#_if data was read in, score it and add to component
			pi[n], po[n] = Pipe()
			args = ( fcst_recs[n], )
			kwargs.update({	'pipe' : po[n] })
			t[n] = Process( target=la.generate_components,
					args=args, 
					kwargs=kwargs )

			t[n].start()

		#_merge aeronet_comp() records after all are computed 
		for n in fcst_recs:
			#_get residuals through pipe
			comps[n] = pi[n].recv()
			t[n].join( 3600 )

			#_merge components 
			components = lt.merge(( components, comps[n] ))

	if len( components ) == 0: dbg( 'No comps calculated, check data', l=3 )

	#_Add dtg_init dtg_vald for later
	summ.__setattr__('dtg_init', dtg_start)	
	summ.__setattr__('dtg_vald', dtg_start)	

	return components, summ, prod

def write_period_points( records, points, label=None, path='.', **kwargs ):
	'''
	given a period of model data, writes netcdf files
	with only the model data for that point.
	
	records	: np.recarray()
	points	: dict, 	contains point name and coordinates
	label	: sting,	normally aeronet site code
	'''
	from netCDF4 import Dataset
	path = '/'.join(( path, 'points' ))
	lt.mkdir_p( path )

	#_get model data
	model = lt.unique( records.model, unique=True ).lower()
	lats = records[0].values.lat
	lons = records[0].values.lon
	yidx = records[0].dimname.index('lat')
	xidx = records[0].dimname.index('lon')
	members = records[0].values.member
	variables = lt.unique( records.variable )

	#_put records in temporal order
	records = sort_rec( records, order=('dtg_vald','fhr',) )
	epochs	= lt.dtg2epoch( lt.unique( records.dtg_init ))
	fhrs	= lt.unique( records.fhr )
	nt	= len( epochs )
	nens	= len( members )
	nhrs	= len( fhrs )

	#_initialize netcdf outputs
	ncdf = {}
	vcdf = {}	
	for code, coords in points.iteritems():
		lat, lon = coords
		i, j = lt.ll2ij( lat, lon, lats, lons )

		#_build filename
		if label == None:
			dtgs = lt.unique( lt.epoch2dtg( epochs ))
			slat = str( int(lat) )
			slon = str( int(lon) )
			file = '-'.join(( model, code, dtgs[0], dtgs[1],
				slat, slon+'.nc' ))
		else:
			file = '-'.join(( model, code, label+'.nc' ))
		file = '/'.join(( path, file ))
		dbg(file)

		#_open output file
		ncdf[code] = Dataset( file, mode='w', format='NETCDF3_CLASSIC' )

		#_create time dimension and coordinate variable, epoch
		ncdf[code].createDimension( 'time', nt )
		ncdf[code].createDimension( 'fhr', nhrs )
		ncdf[code].createDimension( 'member', nens )
		ncdf[code].member = ','.join( members )
		ncdf[code].code = code
		ncdf[code].lat = lat 
		ncdf[code].lon = lon 
		ncdf[code].model = model.upper()
 
		cdf = ncdf[code].createVariable( 'epoch', 'f8', ('time') )
		cdf[:] = epochs
		cdf.long_name = 'unix time for forecast initialization'
		cdf = ncdf[code].createVariable( 'fhr', 'i4', ('fhr') )
		cdf[:] = fhrs 

		vcdf[code] = {}

	#_loop over variables and write to file
	for variable in variables:
	    for fhr in fhrs:	
		fidx = fhrs.index(fhr)
	    	rec_var = subset( records, variable=variable, fhr=fhr )
	    	stack = join_values( rec_var )
	    	ndum, nens, ny, nx = stack.shape

	    	for code, coords in points.iteritems():
			if variable not in vcdf[code]:
				vcdf[code][variable] =\
					ncdf[code].createVariable(variable,'f4',
					('time','fhr','member'),
					fill_value=-9999. )

			lat, lon = coords
			i, j = lt.ll2ij( lat, lon, lats, lons )
	
			#_get values from model at point
			slice = stack[:,:,j,i]

			#_build index list
			tidx = []
			for e_init in lt.dtg2epoch( rec_var.dtg_init ):
				tidx.append( epochs.tolist().index(e_init) )

			#_make initial array of missing values, then fill
			tmp = np.zeros(( nt, nens )) - 9999.
			tmp[tidx,:] = slice[:]

			#_full out aod at
			vcdf[code][variable][:,fidx,:] = tmp[:] #slice[:]

	#_close all files
	for code in points: ncdf[code].close()

def read_period_points( file, pipe=None, **kwargs ):
	''' 
	reads points file as written by write_period_points() 
	this is a throwaway. not meant for the larger system

	'''
	from netCDF4 import Dataset
	import libclass as cl
	dbg(file)
	ncdf = Dataset( file, mode='r', format='NETCDF3_CLASSIC' )
	dtype = [       ('values', cl.var),
        		('model', 'a16'),
			('epoch', 'i8'),
			('fhr', 'i4'),
			('code', 'a10'),
			('lat', 'f4'),
			('lon', 'f4'),
			('variable', 'a20'),
			('ensemble', 'i4'),
			('dimname', tuple),
			('dimsize', tuple) ]
	records = np.recarray( (0,), dtype=dtype )

	model 	= ncdf.model.upper()
	code 	= ncdf.code
	lat	= ncdf.lat
	lon	= ncdf.lon
	members = ncdf.member.split(',')
	fhrs 	= ncdf.variables['fhr'][:]
	epoch_inits = ncdf.variables['epoch'][:]

	nens	= len( members )
	attrn	= ('member',)
	attrv	= (members,)
	dimsize = ( nens, )

	try:
		if len(lat) > 1: lat = lat[0]
	except TypeError:
		pass
	try:
		if len(lon) > 1: lon = lon[0]
	except TypeError:
		pass
	
	for variable in ncdf.variables:
	    for fhr in fhrs:
		fidx = fhrs.tolist().index(fhr)
		if variable in ['fhr','epoch']: continue

		#_get data from netcdf file and put into var class
		values = ncdf.variables[variable][:,fidx,:]		

		#_calculate dtg vald list
		epochs = epoch_inits + fhr * 3600
		epochs = [ int(i) for i in epochs ]

		#_loop over valid times, put in recarray
		for epoch in epochs:
			didx = epochs.index( epoch )			

			vals = values[didx,:]
			vals = np.ma.masked_where( vals == -9999., vals )
			vals = cl.var( vals, attrv=attrv, attrn=attrn )
			
			records.resize( records.size + 1 )
	#		print ( vals, model, epoch, fhr, code, lat,
	#			lon, variable, nens, attrn, dimsize )
			records[-1] = ( vals, model, epoch, fhr, code, lat,
				lon, variable, nens, attrn, dimsize )
	ncdf.close()

###	return lt.exit_pipe( aod, pipe=pipe, attrv=[lats,lons] )

	if pipe == None:
		return records
	else:
		pipe.put(( records, list(attrv) ))

def write_period_statistics( components, summ, path='.', **kwargs ):
	'''
	Writes ascii output of statistical properties of forecast model
	when compared to AERONET. Uses data from generate_period_components() 

	dtg_start	: str*10, 	Start of Scoring Period
	dtg_end		: str*10, 	End of Scoring Period
	finc		: int, 		Hours between scoring periods
	fcst_inc	: int,		Hours between model initializations
	model		: str,		Which ensemble to score
	'''
	from scipy import stats
	from multiprocessing import Process, Pipe
	import libaeronet as la
	percentile = stats.scoreatpercentile
			

	#_Only pass ONE ENSEMBLE at a time
	model = lt.unique( subset( components, ensemble=True ).model, 
					unique=True )

	#_Get first and last model valid times
	dtgs_all = []
	dtg_inits = []
	[[ dtgs_all.append( dtg ) for dtg in dtgs ] for dtgs in summ.dtg_valds ]
	[[ dtg_inits.append( dtg ) for dtg in dtgs ] for dtgs in summ.dtg_inits]

	dtg_start = dtg_inits[0]
	dtg_end = dtg_inits[-1]
	dbg(( model, dtg_start, dtg_end ))

	#_Setup output directory
	dir_out = '/'.join(( path, dtg_start[:6], 'stat' )) 
	lt.mkdir_p( dir_out )

	#_Get lat/lon arrays
	lats, lons = summ[0].values.lat, summ[0].values.lon

	#_Loop over scoring distinctions
	fhrs 		= lt.unique( components.fhr )
	points 		= lt.unique( components.code )
	modes 		= lt.unique( components.mode )
	vars 		= lt.unique( components.variable )
	thresholds 	= lt.unique( components.label )
	error_models 	= lt.unique( components.error_model )

	#_Create dictionary of which points are in which region
	pnt_reg = {}
	for region in reg_dict:
		#_initialize regional list of points
		pnt_reg[ region ] = []
		
		#_Get edge coordinates for region
		llat, ulat, llon, ulon = reg_dict[ region ]['corn']
	
		#_Loop over points we have component scores for	
		for p in points:
			pnt = subset( components, code=p )[0]
			lat, lon = pnt.lat, pnt.lon
			if 	lat <= ulat and lat >= llat and \
				lon <= ulon and lon >= llon:
				pnt_reg[ region ].append( p )	

	dbg(( 'FHRS', fhrs ))
	dbg(( 'PNTS', points ))
	dbg(( 'MODS', modes ))
	dbg(( 'VARS', vars ))
	dbg(( 'THRS', thresholds ))
	dbg(( 'ERRS', error_models ))

	#_loop over Fhrs and Fork to output writing
	groups = lt.setup_groups( fhrs, **kwargs )
	for group in groups:
	  #_initialize list to hold child process PIDs
	  children = []
	
	  #_Loop over fhr
	  for fhr in group:
	    #_fork child process and get PID
	    pid = os.fork()
	
	    #_if parent, save pid
	    if pid != 0: children.append( pid )
	
	    #_if child, loop over variable, mode and error model, writing stats
	    elif pid == 0:

	      #_loop over AOD variables to score
	      for variable in vars:
		#_subset summation
		sub_mn = subset( summ, variable=variable, fhr=fhr, 
			unique=True ).copy()

		#_Get attributes of summation record and apply to mean
		attrn, attrv = get_attr( sub_mn )  
		tmp_mu = sub_mn.values / sub_mn.values.n 
		sub_mn.values = cl.var( tmp_mu, attrn=attrn, attrv=attrv )
		ens_idx = attrn.index( 'member' )

		#_loop over and score AERONET modes { total, coarse, fine }
		for mode in modes:

		  #_loop over error model types
		  for err in error_models:

		    #_labels for thresholds (as defined in la.aerone_comps() )
		    for label in thresholds:

			#_subset components by everything but model
			sub_comp = subset( components, fhr=fhr, 
				variable=variable, error_model=err, 
				label=label, mode=mode )

			#_subset by model
			sub_mod = subset( sub_comp, model=model, ensemble=True )

			#_if not enough data were available to score
			if sub_mod.size == 0: 
				dbg(( 'cannot score', fhr ), l=2 )
				continue

			#_get ensemble size
			members = lt.unique( sub_mod.members )[0] 
			mem_str = ','.join( members  )
			mfmt = str( len( mem_str ) + 3 )
			nens = len( members )
		
			#_setup output file name
			file_out = dir_out + '/' + '_'.join(( label, dtg_start, 
				dtg_end, 'f'+str(fhr).zfill(3), variable,
				 mode, 'pnt', err, model.lower()+'.txt' ))
			dbg( file_out )
			f = open( file_out, 'w' )

			head = '%-8s ' 			#_pnt 
			head += '%10s ' * 16 		#_Stats
			head += '%10s %12s %-'+mfmt+'s'#_nens model members
			head += '%7s' * (nens+1)
			labels = [ 'AERONET', 'MEAN', 'STDEV', 'MAX', '90TH',
			'75TH', '50TH', '25TH', '10TH', 'MIN', 'NOBS', 'RMSE', 
			'BRIERS', 'BIAS', 'MAE', 'NMRSE', 'FGE', 'NMEMBERS', 
			'MODEL', 'MEMBERS' ]

			#_Add Headers for Ranks
			[ labels.append(str(i).zfill(2)+'_RNK')  \
				for i in np.arange(nens+1) ]

			header = head % tuple(labels)
			f.write( header + '\n' )

			#_Loop over components and write to file
			for p in points:
				'''
				pnt__8_mean__8
				capove 01.1111 
				'''
				#_pull out point
				sub_pnt = subset( sub_mod, code=p )
				if sub_pnt.size == 0:
					dbg(( 'cannot score', p ), l=5 )
					continue

				#_pass components to calculate scores
				stats = la.aeronet_score( sub_pnt ) 

				dbg(( model, p, sub_pnt.size, fhr, variable, 
					err, label, mode ), l=5 )

				#_Get model indices for point lat/lon
				a_lat, a_lon = sub_pnt[0].lat, sub_pnt[0].lon 
				i, j = lt.ll2ij( a_lat, a_lon, lats, lons )

				#_FOR_PERCENTILES_______________________________
				fcst_sp = sub_mn.values[:,j,i] 

				#_Masked arrays and Stats module are not happy
				if hasattr( fcst_sp, 'mask' ) and	\
				type( fcst_sp.mask ) != np.bool_:
					fcst_sp = fcst_sp[fcst_sp.mask==False]
				#_fcst_sp == values of ensemble members at pnt
	
				#_Create list to write 
				values = [ p ]
				values.append( fcst_sp.mean() )
				values.append( fcst_sp.std() )
				values.append( percentile( fcst_sp, 100) )
				values.append( percentile( fcst_sp, 90 ) )
				values.append( percentile( fcst_sp, 75 ) )
				values.append( percentile( fcst_sp, 50 ) )
				values.append( percentile( fcst_sp, 25 ) )
				values.append( percentile( fcst_sp, 10 ) )
				values.append( percentile( fcst_sp,  0 ) )

				values.append( stats.nobs )

				values.append( stats.rmse )
				values.append( stats.briers )
				values.append( stats.bias )
				values.append( stats.mae )
				values.append( stats.nrmse )
				values.append( stats.fgs )
				values.append( nens )
				values.append( model )
				values.append( mem_str )
				[ values.append( a ) for a in stats.ranks ]

				fmt = '%-8s '
				fmt += '%10.6f ' * 9	#_Perc
				fmt += '%10i '		#_Nobs 
				fmt += '%10.6f ' * 6 	#_Stats
				fmt += '%10i %12s %-'+mfmt+'s'
				fmt += '%7i' * (nens+1)	
				line = fmt % tuple(values)
				f.write( line + '\n' )		

			#_ENSEMBLE_MEMBER_SCORING_______________________________
			#_create dictionary of member scores by point
			stats = {}
			for p in points:
				stats[p] = {}
				for m in members:
					#_subset by member and score
					sub_mp = subset( sub_comp, code=p,
						model=m )

					stats[p][m] = la.aeronet_score( sub_mp )
					dbg(( m, p, sub_mp.size, fhr, variable,
						err, label, mode ), l=5 )

 
			#_Score individual members
			#_Loop over statistics and write member values
			for stat in ['rmse','bias','mae','nrmse','fgs']:
				#_Setup header
				labels = [ 'aeronet-' + stat ]
				[ labels.append( mem ) for mem in members ]
				fmt = '%-12s '		#_aeronetsite
				fmt += '%16s ' * nens	#_stat (member name)
	
				#_Write header
				header = fmt % tuple(labels)
				f.write( '\n' + header + '\n' ) 
				for p, st in stats.iteritems():

					#_Setup output format
					fmt = '%-12s '#_aeronet site,
					fmt += '%16.5f ' * nens	#_stat

					#_Setup output
					values = [ p ]
					[ values.append( 
						st[m].__getattribute__(stat))
						for m in members ]
		
					#_Setup formatted print
					line = fmt % tuple(values)
					f.write( line + '\n' )	

			#_Close point stats file
			f.close()
			lt.make_readable( file_out )

			#_______________________________________________________
			#_REGIONAL_SCORING______________________________________
			#_Setup output file name
			file_out = dir_out + '/' + '_'.join(( label, dtg_start, 
				dtg_end, 'f'+str(fhr).zfill(3), variable, mode,
				'reg', err, model.lower()+'.txt' ))
			dbg( file_out )
			f = open( file_out, 'w' )

			head = '%13s ' 			#_pnt 
			head += '%10s ' * 16 		#_Stats
			head += '%10s %12s %-'+mfmt+'s'#_nens model members
			head += '%7s' * (nens+1)
			labels = [ 'REGION', 'MEAN', 'STDEV', 'MAX', '90TH',
			'75TH', '50TH', '25TH', '10TH', 'MIN', 'NOBS', 'RMSE', 
			'BRIERS', 'BIAS', 'MAE', 'NMRSE', 'FGE', 'NMEMBERS', 
			'MODEL', 'MEMBERS' ]

			#_Add Headers for Ranks
			[ labels.append(str(i).zfill(2)+'_RNK')  \
				for i in np.arange(nens+1) ]

			header = head % tuple(labels)
			f.write( header + '\n' )

			#_Loop over regions
			for r in reg_dict:
				#_Pass components to calculate scores
				sub_reg = subset( sub_mod, code=pnt_reg[r] )
				stats = la.aeronet_score( sub_reg ) 

				#_If not enough data were available to score
				if sub_reg.size == 0: 
					dbg( 'Cannot Score '+r+' '+err, l=2 )
					continue

				#_FOR_PERCENTILES_______________________________
				#_Subset records by region
				fcst_sr = sub_region( sub_mn, region=r )
				fcst_mn = np.mean( fcst_sr.values, axis=ens_idx)
				fcst_mn = fcst_mn.copy().flatten()

				#_Percentile gets weird with masked values
				if hasattr( fcst_mn, 'mask' ) and	\
				type( fcst_mn.mask ) != np.bool_:
					fcst_mn = fcst_mn[fcst_mn.mask==False]
				
				# fcst_sr == recarray.values[nens, ny, nx]
				# fcst_mn == (ny, nx )

				#_Create list to write 
				values = [ r ]
				values.append( fcst_mn.mean() )
				values.append( fcst_mn.std() )
				values.append( percentile( fcst_mn, 100) )
				values.append( percentile( fcst_mn, 90 ) )
				values.append( percentile( fcst_mn, 75 ) )
				values.append( percentile( fcst_mn, 50 ) )
				values.append( percentile( fcst_mn, 25 ) )
				values.append( percentile( fcst_mn, 10 ) )
				values.append( percentile( fcst_mn,  0 ) )

				values.append( stats.nobs )

				values.append( stats.rmse )
				values.append( stats.briers )
				values.append( stats.bias )
				values.append( stats.mae )
				values.append( stats.nrmse )
				values.append( stats.fgs )
				values.append( nens )
				values.append( model )
				values.append( mem_str )
				[ values.append( a ) for a in stats.ranks ]

				fmt = '%-13s '
				fmt += '%10.6f ' * 9	#_Perc
				fmt += '%10i '		#_Nobs 
				fmt += '%10.6f ' * 6 	#_Stats
				fmt += '%10i %12s %-'+mfmt+'s'
				fmt += '%7i' * (nens+1)	
				line = fmt % tuple(values)
				f.write( line + '\n' )		

			#_ENSEMBLE_MEMBER_SCORING_______________________________
			#_Score individual members
			stats = {}
			for r in reg_dict:
				stats[r] = {}
				for m in members:
					#_Subset by member and score
					sub_mp = subset( sub_comp, model=m,  
						code=pnt_reg[r] )
					stats[r][m] = la.aeronet_score( sub_mp )
		
			#_write rank header	
			labels = [ 'tmp_header' ]
			[ labels.append( mem ) for mem in members ]
			fmt = '%-12s '#_aeronet site
			fmt += ('%'+str(10*nens)+'s ') * nens	#_stat
	
			#_write header
			dbg(fmt)
			dbg( labels )
			header = fmt % tuple(labels)
			f.write( '\n' + header + '\n' ) 
			for r, st in stats.iteritems():
				#_setup output format
				fmt = '%-12s '#_aeronet site
				fmt += '%10i ' * nens**2	#_stat

				#_setup output
				values = [ r ]
				for m in members:
					for mm in xrange(nens):
					    try:
						values.append( st[m].ranks[mm] )
					    except:
					        values.append( 0 )
##					if hasattr( st[m].ranks, '__iter__'):
##						[ values.append( rank ) \
##						for rank in st[m].ranks ]
##					else:
##						[ values.append( zz ) \
##						for zz in np.zeros( nens ) ]
	
				#_Setup formatted print
##				dbg( fmt )
##				dbg( values )
				line = fmt % tuple(values)
				f.write( line + '\n' )	

			#_Score individual members
			#_Loop over statistics and write member values
			for stat in ['rmse','bias','mae','nrmse','fgs']:
				#_Setup header
				labels = [ 'aeronet-' + stat ]
				[ labels.append( mem ) for mem in members ]
				fmt = '%-12s '	#_aeronetsite
				fmt += '%16s ' * nens	#_stat (member name)
	
				#_Write header
				header = fmt % tuple(labels)
				f.write( '\n' + header + '\n' ) 

				for r, st in stats.iteritems():

					#_Setup output format
					fmt = '%-12s '#_aeronet site
					fmt += '%16.5f ' * nens	#_stat

					#_Setup output
					values = [ r ]
					[ values.append( 
						st[m].__getattribute__(stat))
						for m in members ]
		
					#_Setup formatted print
					line = fmt % tuple(values)
					f.write( line + '\n' )	
			f.close()
			lt.make_readable( file_out )

	      #_Kill child process
	      os._exit(0) 

	  #_Wait for PID exit signals
	  for kid in children: os.waitpid( kid, 0 )

def read_location_stats( path='.', bins=None, **kwargs ):
	'''
	Globs all files in path that have the dtg start time
	dtgs	: str*10,	I'm limiting this to one dtgs per object
				because I know I'm going to accidentally
				score incompatible data soon
	path	: str,		Path to directory with stat ascii files
	'''
	import libclass as cl
	from glob import glob
	from multiprocessing import Process, Pipe
	path = '/'.join(( path, 'stat' ))
	dbg( path )
	
	#_get list of files to read
	files = glob( path + '/*.txt' )
	if len(files) == 0: raise IOError, 'no input files ' + path	

	#_initialize ensemble and member output
	tmp = cl.location_stats( files[0] )
	ens, mem = tmp.ensemble, tmp.member
	del files[0] 

	ens_objs, mem_objs = [], []

	#_setup multproc
	groups = lt.setup_groups( files, **kwargs )
	for group in groups:
		#_initialize temporary storage
		l 	= len( group )
		dbg((l, group[0].split('/')[-1]))
		thread	= [None]*l
                p 	= [None]*l
                c	= [None]*l
                for i in range(l):  #-Pipe init, open connections 
                        file		= group[i]
                        p[i],c[i]       = Pipe()
                        args		= ( file, )
                        kwargs.update( { 'pipe' : c[i] } )
                        thread[i]       = Process( target=cl.location_stats, \
                                        args=args, kwargs=kwargs)
                        thread[i].start()

		#_recv data, close connections
                for i in range(l):  
			e, m = p[i].recv()
			ens_objs.append( e )
			mem_objs.append( m )
                        thread[i].join()

	dbg('merging')
	for e in ens_objs: ens = lt.merge(( ens, e ))
	for m in mem_objs: mem = lt.merge(( mem, m ))

	dbg(( lt.unique( ens.model ), '<=', lt.unique( mem.model ) ))
	dbg( lt.unique( ens.location ))
	dbg( lt.unique( ens.label ))
	dbg( lt.unique( ens.fhr ))
	dbg( lt.unique( ens.error_model ))
	dbg(( lt.unique( ens.dtgs ),  lt.unique( ens.dtge )) )

	#_reduce bins kept 
	ens = subset( ens, label=bins )
	mem = subset( mem, label=bins )

	#_Return ensemble and member stats
	return ens, mem

def generateFineCoarse( records ):
	''' 
	temporary kludge to make fine/coarse mode values from records 
		
	THIS DOES NO REAL RECORD CHECKING AND ASSUMES UNIQUE DTG_VALD

	'''

	specs = { 	'coarse_aod'	: ['seasalt_aod', 'dust_aod'],
			'fine_aod' 	: ['smoke_aod','sulfate_aod'] }
	out = cl.model_object()
	o = ('dtg_vald',)
	for mode, species in specs.iteritems():
		spec0, spec1 = species
		spec0 = sort_rec( subset( records, variable=spec0 ), order=o )
		spec1 = sort_rec( subset( records, variable=spec1 ), order=o )
		refrec = spec0[0].copy()

		dat0 = join_values( spec0 )
		dat1 = join_values( spec1 )

		data = lt.masked_sum( dat0, dat1 ); del dat0; del dat1
	
		dtg_vald = lt.unique( spec0.dtg_vald )
		dtg_init = lt.unique( spec0.dtg_init )
		nt = len( lt.unique( dtg_vald ))
		recs = cl.model_object()
		recs.resize( nt )

		atn, atv = get_attr( refrec )

		recs.dtg_vald 	= spec0.dtg_vald
		recs.dtg_init 	= spec0.dtg_init
		recs.fhr 	= spec0.fhr
		recs.model 	= [refrec.model] * nt
		recs.ensemble 	= [refrec.ensemble] * nt
		recs.region 	= [refrec.region] * nt
		recs.dimname 	= [refrec.dimname] * nt
		recs.dimsize 	= [refrec.dimsize] * nt
		recs.units 	= [refrec.units] * nt
		recs.long_name 	= [refrec.long_name] * nt
		recs.variable 	= [mode] * nt
		for t in xrange(nt): 
			recs.values[t] = cl.var( data[t], attrv=atv, attrn=atn )
		out = lt.merge((out, recs))
	return out

def running_summation( records, mean, min=None, geomu=None, **kwargs ):
	''' 
	Used with cum mean
	mean	: cl.cum_mean() object
	records	: Records include in sum, sorted by fhr, variable 
	min	: flt,	Minimum value to sum.  The rest will be masked.
	geolog	: recarray,	When not NoneType, doesn't do arithmetic summ,
				but log of value over geomean 
	Sums records with matching forecast hours, variable and model names
	'''
	import libaeronet as la
	min = 1e-8 if min == None else min

	exp = expected_vars()

	local = mean.copy()
	for rec in records:
		atn, atv = get_attr( rec )
		val = rec.values.copy()
		mbr = rec.values.member.tolist()
		mod = rec.model
		dti = rec.dtg_init
		dtg = rec.dtg_vald	
		fhr = rec.fhr
		reg = rec.region
		vrb = rec.variable
		ens = rec.ensemble
		dmn = rec.dimname
		dsz = rec.dimsize
		unt = rec.units
		lng = rec.long_name

		#_pull out matching mean currently stored 
		tmp = subset( local, fhr=fhr, model=mod, variable=vrb )

		#_If requiring a threshold value
		if min != None: val = np.ma.masked_less( val, min )

		if geomu != None:
###			for member in mbr:
###				idx = mbr.index(member)
###				if vrb not in exp[member]:
###					dbg(( 'masking_summ', member ))
###					dbg(( vrb, exp[member] ))
###					val[idx,:,:].mask = True	

			#_zero substitution...
			mask = val.mask
			tidx = np.where( val.data < 0.001 )
			val[tidx] = 0.001
			val.mask = mask

			#_this whole section is full of assumptions 
			# I don't have time to correct
			nens, ny, nx = dsz
			mean_append = val.mean(axis=0).reshape(1,ny,nx)
			val = np.append( val, mean_append, axis=0 )
			gmu = subset( geomu,variable=vrb,fhr=fhr,unique=True )
			if (vrb == 'fine_aod' or vrb == 'coarse_aod') \
			and fhr == 12:
				dbg((vrb, 'X:', fhr, val[:,89,179]))
			val = np.log( val / gmu.values )**2

			if (vrb == 'fine_aod' or vrb == 'coarse_aod') \
			and fhr == 12:
				dbg((vrb, 'W:', fhr, val[:,89,179]))
				dbg((vrb, 'G:', fhr, gmu.values[:,89,179]))

			atv[0] = np.append( atv[0], mod )

		#_Create a ny,nx array of ones where we will be addings data
		n = np.ones( val.shape )

		if hasattr( val, 'mask' ) \
		and val.mask.shape == dsz:
			dbg( 'masking', l=7 )
			n[val.mask] = 0

		#_replace fill values with zeroes to not mess up summs
		val = np.ma.MaskedArray.filled( np.ma.masked_array(val), 0 )

		#_initialize this fhr/mod/variable entry
		if tmp.size == 0:
			dbg(( tmp.size, 'init' ), l=7 )
		
			#_add attribute for N to calculate mean
			val = cl.var( val, attrn=atn, attrv=atv )
			val.__setattr__('n', n)	

			local.resize( local.size + 1 )
			local[-1] = ( val, mod, dti, dti, [dti], [dtg], fhr,
				reg, vrb, ens, dmn, dsz, unt, lng ) 
	
		#_add Data if Present
		elif tmp.size == 1:
			dbg(( tmp.size, 'summing' ), l=7 )
			tmp = tmp[0]
		
			#_sum values
			tmp.values[:] = lt.masked_sum( tmp.values, val )
			val = cl.var( val, attrn=atn, attrv=atv )

			#_Keep list of DTGS averaged
			tmp.dtg_inits[:] = np.append( tmp.dtg_inits, dti ) 
			tmp.dtg_valds[:] = np.append( tmp.dtg_valds, dtg ) 

			#_Increment N where appropriate
			tmp.values.n += n

		#_Shouldn't hit here for a single forecast
		else:
			raise ValueError( 'Mod/Fhr/Vrb should be unique' )
	return local 

def running_product( records, prod, min=None, **kwargs ):
	''' 
	Used with cum mean
	prod	: cl.cum_mean() object
	records	: Records include in sum, sorted by fhr, variable 
	min	: flt,	Minimum value to sum.  The rest will be masked.

	Sums records with matching forecast hours, variable and model names
	'''
	import libaeronet as la
	min = 1e-8 if min == None else min
	local = prod.copy()
	for rec in records:
		val = rec.values.copy()
		atn, atv = get_attr( rec )
		mod = rec.model
		dti = rec.dtg_init
		dtg = rec.dtg_vald	
		fhr = rec.fhr
		reg = rec.region
		vrb = rec.variable
		ens = rec.ensemble
		dmn = rec.dimname
		dsz = rec.dimsize
		unt = rec.units
		lng = rec.long_name

		#_pull out matching prod currently stored 
		tmp = subset( local, fhr=fhr, model=mod, variable=vrb )

		#_if requiring a threshold value
		if min != None: val = np.ma.masked_less( val, min )

		#_create a ny,nx array of ones where we will be addings data
		n = np.ones( dsz )
		if hasattr( rec.values, 'mask' ) \
		and rec.values.mask.shape == dsz:
			dbg( 'masking', l=7 )
			n[rec.values.mask] = 0

		#_replace fill values with zeroes to not mess up products 
		val = np.ma.MaskedArray.filled( np.ma.masked_array(val), 1 )

		#_initialize this fhr/mod/variable entry
		if tmp.size == 0:
			dbg(( tmp.size, 'init' ), l=7 )

			#_add attribute for N to calculate prod
			val = cl.var( val, attrn=atn, attrv=atv )
			val.__setattr__('n', n)	

			local.resize( local.size + 1 )
			local[-1] = ( val, mod, dti, dti, [dti], [dtg], fhr,
				reg, vrb, ens, dmn, dsz, unt, lng ) 
	
		#_add Data if Present
		elif tmp.size == 1:
			dbg(( tmp.size, 'multiplying' ), l=7 )
			tmp = tmp[0]
		
			#_Sum values	
			tmp.values[:] = lt.masked_prod( tmp.values, val )
			val = cl.var( val, attrn=atn, attrv=atv )

			#_Keep list of DTGS averaged
			tmp.dtg_inits[:] = np.append( tmp.dtg_inits, dti ) 
			tmp.dtg_valds[:] = np.append( tmp.dtg_valds, dtg ) 

			#_Increment N where appropriate
			tmp.values.n += n
		
		#_Shouldn't hit here for a single forecast
		else:
			raise ValueError( 'Mod/Fhr/Vrb should be unique' )
	
	return local
 
def read_ensperiod( dtg_range, fhr=None, fcst_finc=24, variables=None,**kwargs):
	'''
	read in a period of ensemble data

	dtg_start	: str*10, 	start of read period
	dtg_end		: str*10, 	end of read period
	fcst_inc	: int,		hours between model initializations
	fhrs		: list||str,	fhrs to return 
	'''
	from multiprocessing import Process,Pipe
	import libaeronet as la
	import libicap as li

	dtg_start, dtg_end = dtg_range
	dbg(( dtg_start, dtg_end ))
	aod = cl.model_object()

	#_Get plist of dtgs to read in
	dtgs = []
	dtg_loop = dtg_start
	while dtg_loop <= dtg_end:	
		dtgs.append( dtg_loop )
		dtg_loop = lt.newdtg( dtg_loop, fcst_finc )

	#_READ_IN_DATA__________________________________________________________
	#_Scoring is now done on the fly, which requires a whole lot more
	# disk I/O.  Reading is now threaded, as it happens so damn many times
	groups = lt.setup_groups( dtgs, **kwargs )
	for group in groups:
		l 	= len( group )	#_Can't use nproc if grp%np!=0	
		t 	= [None]*l	#_List of processing threads
		pi	= [None]*l	#_List of pipes in
		po	= [None]*l	#_List of pipes out

		fcst_recs = [] 

		#_reading of ensemble data
		for n in np.arange( l ):
			dtg = group[n]	#_Name of point

			pi[n], po[n] = Pipe()
			args = ( dtg, )
			kwargs.update( { 'pipe'	: po[n] } )
			t[n] = Process( target=read_ensfcst,
					args=args, 
					kwargs=kwargs )
			t[n].start()
		
		#_collect multi-processing threads of forecast data 
		for n in np.arange( l ):
			dtg = group[n]

			#_store in temporary dictionary. Missing data loops.
			tmp_fcst, attrv = pi[n].recv()
			if type( tmp_fcst ) != np.recarray: 
				dbg(('error reading model', dtg), l=3)
				continue

			#_i am not good at python... my dim attributes aren't
			# being transferred through the pipe (since it pickles
			# them), so I reset them here. Can be fixed with 
			# Manager?
			for rec in tmp_fcst:
				rec.values = cl.var( rec.values, attrv=attrv )

			#_reduce to only variables desired (2013.06)
			tmp_fcst = subset( tmp_fcst, variable=variables )

##			#_if ICAP and strict, filter out records with 
##			tmp_fcst = li.filter( tmp_fcst, **kwargs )

			#_filter by FHR
			tmp_fcst = subset( tmp_fcst, fhr=fhr )

			#_If all records weren't filtered out, continue 
			if tmp_fcst.size > 0: 
				aod = lt.merge(( aod, tmp_fcst ))
			else: 
				dbg(( 'no model records', dtg ), l=3 )
			t[n].join( 60 )

##	aod = subset( aod, variable=variables )
	return aod

def write_period_musig( recs, path='.', label=None, **kwargs ):
	'''
	Given forecast records, generates a file with mean and
	std for each variable

	recs		: np.recarray()
	path		: str,	Where to dump these files
        '''
	import libaeronet as la
        from netCDF4 import Dataset
	from scipy import stats
	percentile = stats.mstats.mquantiles
	
        #_initialize lists to use for averaging
	dtg_valds = []
	dtg_inits = []

	#_get list of every possible dtg
	dtg_valds = lt.unique( recs.dtg_vald )
	dtg_inits = lt.unique( recs.dtg_init )
	dtg_start = dtg_inits[0]
	dtg_end = dtg_inits[-1]

	#_use first record as a metadata reference
	ref_rec = recs[0]

	#_Get name of ensemble
	model = lt.unique( recs.model, unique=True ) 

	is_ens = True if hasattr( ref_rec.values, 'member' ) else False

	dir_out = '/'.join(( path, model.upper(), dtg_start[:6], 'musig' )) 
	lt.mkdir_p( dir_out )
	if label == None:
	        file_out = dir_out + '/' + '_'.join(( 'musig', model.lower(),
			dtg_start, dtg_end+'.nc' ))
	else:
       		file_out = dir_out + '/' + '_'.join(( 'musig', model.lower(),
			label+'.nc' ))
        dbg( file_out )

        #_Initialize output file
        ncdf = Dataset( file_out, 'w', format='NETCDF3_CLASSIC')
	
	#_Get list of variables
	vars = lt.unique( recs.variable )

	#_Add dtgs used
        ncdf.__setattr__( 'dtg_valds', ','.join( dtg_valds ) )

        #_Initialize model dimensions
        lats = ref_rec.values.lat
        lons = ref_rec.values.lon
	fhrs = lt.unique( recs.fhr )
	nt = len( fhrs )
	ny, nx = lats.size, lons.size
	ncdf.createDimension( 'lat', ny )
	ncdf.createDimension( 'lon', nx )
	ncdf.createDimension( 'fhr', nt ) 

	#_Write coordinates
	cdf = ncdf.createVariable( 'lat', 'f4', ('lat',) )
	cdf[:] = lats
	cdf = ncdf.createVariable( 'lon', 'f4', ('lon',) )
	cdf[:] = lons
	cdf = ncdf.createVariable( 'fhr', 'f8', ('fhr',) )
	cdf[:] = fhrs 

	members = ref_rec.values.member
	members = np.append( members, model.upper() )
	nens = members.size
	ncdf.createDimension( 'member', nens )

	dimname = list(ref_rec.dimname)
	dimname.insert(0,'fhr')

	#_add member names as global attribute
	ncdf.member = ','.join( members )

	#_find index of members
	midx = dimname.index('member')

	exp = expected_vars()

	#_loop over variables and write mean to file
	for v in vars:
	    dbg( v )
	    mus = np.zeros((nt,nens,ny,nx))
	    sigmas = np.zeros((nt,nens,ny,nx))
	    for fhr in fhrs:
		t = fhrs.index(fhr)

		#_subset record, the sort by dtg_vald
		quart_var = subset( recs, variable=v, fhr=fhr )

		#_join into single array
		quart_stack = join_values( quart_var )
		quart_stackmu = quart_stack.mean( axis=midx )

		#_calculate arithemetic mean for members and ensemble
		quart = np.mean( quart_stack, axis=0 )
		quartmu = np.mean( quart_stackmu, axis=0 )

		#_calculate standard deviation for members and ensemble
		quart_std = np.std( quart_stack, axis=0 )
		quartmu_std = np.std( quart_stackmu, axis=0 )

		#_merge into 
		mu 	= np.append( quart, quartmu.reshape(1,ny,nx), axis=0 )
		sigma 	= np.append( quart_std, quartmu_std.reshape(1,ny,nx),
				axis=0 )

		#_add to array
		mus[t] = mu[:]
		sigmas[t] = sigma[:]

	    #_sort records, then join along time and member axes
	    cdf = ncdf.createVariable( v+'_mean','f4',dimname,fill_value=-9999)
	    cdf[:] = mus 
	    cdf = ncdf.createVariable( v+'_sigma','f4',dimname,fill_value=-9999)
	    cdf[:] = sigmas 

        ncdf.close()
	lt.make_readable( file_out )
	return file_out

def write_period_quant( recs, path='.', label=None, **kwargs ):
	'''
	summ		: cl.cum_mean() class, 	if made pairwise, it is 
						done elsewhere
	path		: str,	Where to dump these files

	(nmember, ny, nx) x nspec model variables, averaged in time
	(npoint) x nmode aeronet observations, averaged for same period 
       
	mean should be cl.cum_mean() object produced by running_summation()

	Any attempts at pairwise means have to be done in RUNNING_SUMMATION() 
        '''
	import libaeronet as la
        from netCDF4 import Dataset
	from scipy import stats
	percentile = stats.mstats.mquantiles
	quantiles = [0, 16, 25, 50, 75, 84, 100]
	q = np.array( quantiles ) / 100.

        #_initialize lists to use for averaging
	dtg_valds = []
	dtg_inits = []

	#_get list of every possible dtg
	dtg_valds = lt.unique( recs.dtg_vald )
	dtg_inits = lt.unique( recs.dtg_init )

	#_use first record as a metadata reference
	ref_rec = recs[0]

	#_Get name of ensemble
	model = lt.unique( recs.model, unique=True ) 

	is_ens = True if hasattr( ref_rec.values, 'member' ) else False

	#_File out
	dtg_start = dtg_inits[0]
	dtg_end = dtg_inits[-1]

	dir_out = '/'.join(( path, model.upper(), dtg_start[:6], 'quart' )) 
	lt.mkdir_p( dir_out )
	if label == None:
	        file_out = dir_out + '/' + '_'.join(( 'quart', model.lower(),
			dtg_start, dtg_end+'.nc' ))
	else:
       		file_out = dir_out + '/' + '_'.join(( 'quart', model.lower(),
			label+'.nc' ))
        dbg( file_out )

        #_Initialize output file
        ncdf = Dataset( file_out, 'w', format='NETCDF3_CLASSIC')
	
	#_Get list of variables
	vars = lt.unique( recs.variable )

	#_Add dtgs used
        ncdf.__setattr__( 'dtg_valds', ','.join( dtg_valds ) )

        #_Initialize model dimensions
        lats = ref_rec.values.lat
        lons = ref_rec.values.lon
	fhrs = lt.unique( recs.fhr )
	nhrs = len( fhrs )
	ny, nx = lats.size, lons.size
	nq = len( q )
	ncdf.createDimension( 'lat', ny )
	ncdf.createDimension( 'lon', nx )
	ncdf.createDimension( 'fhr', nhrs ) 
	ncdf.createDimension( 'quantile', nq ) 

	#_Write coordinates
	cdf = ncdf.createVariable( 'lat', 'f4', ('lat',) )
	cdf[:] = lats
	cdf = ncdf.createVariable( 'lon', 'f4', ('lon',) )
	cdf[:] = lons
	cdf = ncdf.createVariable( 'fhr', 'f8', ('fhr',) )
	cdf[:] = fhrs 
	cdf = ncdf.createVariable( 'quantile', 'f8', ('quantile',) )
	cdf[:] = q

	members = ref_rec.values.member
	members = np.append( members, model.upper() )
	nens 	= members.size
	ncdf.createDimension( 'member', nens )

	dimname = list(ref_rec.dimname)
	dimname.insert(0,'quantile')
	dimname.insert(0,'fhr')

	#_add member names as global attribute
	ncdf.member = ','.join( members )

	#_find index of members
	midx = dimname.index('member')
	exp = expected_vars()

	#_loop over variables and write mean to file
	for v in vars:
	    	quants = np.zeros((nhrs,nq,nens,ny,nx))
	        dbg(v)

		for fhr in fhrs:	
			t = fhrs.index(fhr)
	
			#_subset record, the sort by dtg_vald
			quant_var = subset( recs, variable=v, fhr=fhr )
			ntx = quant_var.size
	
			#_join into single array
			quant_stack = join_values( quant_var )
			quant_stackmu = quant_stack.mean( axis=midx-1 )

			#_add ensemble mean to array
			dbg((quant_stack.shape,quant_stackmu.shape,ny,nx,midx))
			quant = np.append( quant_stack,
				quant_stackmu.reshape(ntx,1,ny,nx), axis=1 )

			#_figure out better way of indexing this
			for i in xrange(nx):
			  for e in xrange(nens):
			    for j in xrange(ny):
				col = quant[:,e,j,i]
				quants[t,:,e,j,i] = percentile( col, prob=q )

		cdf = ncdf.createVariable( v, 'f4', dimname,
			fill_value=-9999)
		cdf[:] = quants[:] #[:,idx,:,:,:]

		#_sort records, then join along time and member axes
	##	for perc in quantiles:
	##		idx = quantiles.index(perc)
	##		name = '_'.join(( v, string(perc).zfill(3) ))
	##
	##		cdf = ncdf.createVariable( name, 'f4', dimname,
	##			fill_value=-9999)
	##		cdf[:] = quants[:,idx,:,:,:]

        ncdf.close()
	lt.make_readable( file_out )
	return file_out

def write_period_median( recs, path='.', label=None, **kwargs ):
	'''
	summ		: cl.cum_mean() class, 	if made pairwise, it is 
						done elsewhere
	path		: str,	Where to dump these files

	(nmember, ny, nx) x nspec model variables, averaged in time
	(npoint) x nmode aeronet observations, averaged for same period 
       
	mean should be cl.cum_mean() object produced by running_summation()

	Any attempts at pairwise means have to be done in RUNNING_SUMMATION() 
        '''
	import libaeronet as la
        from netCDF4 import Dataset
	from scipy import stats
	percentile = stats.mstats.mquantiles
	
        #_initialize lists to use for averaging
	dtg_valds = []
	dtg_inits = []

	#_get list of every possible dtg
	dtg_valds = lt.unique( recs.dtg_vald )
	dtg_inits = lt.unique( recs.dtg_init )

	#_use first record as a metadata reference
	ref_rec = recs[0]

	#_Get name of ensemble
	model = lt.unique( recs.model, unique=True ) 

	is_ens = True if hasattr( ref_rec.values, 'member' ) else False

	#_File out
	dtg_start = dtg_inits[0]
	dtg_end = dtg_inits[-1]

	dir_out = '/'.join(( path, model.upper(), dtg_start[:6], 'median' )) 
	lt.mkdir_p( dir_out )
	if label == None:
	        file_out = dir_out + '/' + '_'.join(( 'median', model.lower(),
			dtg_start, dtg_end+'.nc' ))
	else:
       		file_out = dir_out + '/' + '_'.join(( 'median', model.lower(),
			label+'.nc' ))
        dbg( file_out )

        #_Initialize output file
        ncdf = Dataset( file_out, 'w', format='NETCDF3_CLASSIC')
	
	#_Get list of variables
	vars = lt.unique( recs.variable )

	#_Add dtgs used
        ncdf.__setattr__( 'dtg_valds', ','.join( dtg_valds ) )

        #_Initialize model dimensions
        lats = ref_rec.values.lat
        lons = ref_rec.values.lon
	fhrs = lt.unique( recs.fhr )
	nt = len( fhrs )
	ny, nx = lats.size, lons.size
	ncdf.createDimension( 'lat', ny )
	ncdf.createDimension( 'lon', nx )
	ncdf.createDimension( 'fhr', nt ) 

	#_Write coordinates
	cdf = ncdf.createVariable( 'lat', 'f4', ('lat',) )
	cdf[:] = lats
	cdf = ncdf.createVariable( 'lon', 'f4', ('lon',) )
	cdf[:] = lons
	cdf = ncdf.createVariable( 'fhr', 'f8', ('fhr',) )
	cdf[:] = fhrs 

	members = ref_rec.values.member
	members = np.append( members, model.upper() )
	nens = members.size
	ncdf.createDimension( 'member', nens )

	dimname = list(ref_rec.dimname)
	dimname.insert(0,'fhr')

	#_add member names as global attribute
	ncdf.member = ','.join( members )

	#_find index of members
	midx = dimname.index('member')

	exp = expected_vars()

	#_loop over variables and write mean to file
	for v in vars:
	    dbg( v )
	    medians = np.zeros((nt,nens,ny,nx))
	    for fhr in fhrs:
		t = fhrs.index(fhr)

		#_subset record, the sort by dtg_vald
		quart_var = subset( recs, variable=v, fhr=fhr )

		#_join into single array
		quart_stack = join_values( quart_var )
		quart_stackmu = quart_stack.mean( axis=midx )

		quart = np.median( quart_stack, axis=0 )
		quartmu = np.median( quart_stackmu, axis=0 )

		quart = np.append( quart, quartmu.reshape(1,ny,nx), axis=0 )
		medians[t] = quart[:]

	    #_sort records, then join along time and member axes
	    cdf = ncdf.createVariable( v, 'f4', dimname, fill_value=-9999 )
	    cdf[:] = medians 

        ncdf.close()
	lt.make_readable( file_out )
	return file_out

def write_period_mean( summ, path='.', aeronet=None, **kwargs ):
	'''
	summ		: cl.cum_mean() class, 	if made pairwise, it is 
						done elsewhere
	path		: str,	Where to dump these files

	(nmember, ny, nx) x nspec model variables, averaged in time
	(npoint) x nmode aeronet observations, averaged for same period 
       
	mean should be cl.cum_mean() object produced by running_summation()

	Any attempts at pairwise means have to be done in RUNNING_SUMMATION() 
        '''
	dbg( "THESE AREN'T TRUE PAIRWISE AERONET PLOTS", l=3 )
	dbg( "THAT WOULD REQUIRE ALL AERONET SITES REPORTING ALWAYS", l=3 )
        from netCDF4 import Dataset
	import libaeronet as la
	
        #_Initialize lists to use for averaging
	dtg_valds = []
	dtg_inits = []

	mean = summ.copy()

	#_Get list of every possible dtg
	[[dtg_valds.append(d) for d in dtg_list] for dtg_list in mean.dtg_valds]
	[[dtg_inits.append(d) for d in dtg_list] for dtg_list in mean.dtg_inits]
	dtg_valds = lt.unique( dtg_valds )
	dtg_inits = lt.unique( dtg_inits )

	#_Use first record as a metadata reference
	ref_rec = mean[0]

	#_Get name of ensemble
	model = lt.unique( mean.model, unique=True ) 

	#_File out
	dtg_start = dtg_inits[0]
	dtg_end = dtg_inits[-1]

	dir_out = '/'.join(( path, dtg_start[:6], 'mean' )) 
	lt.mkdir_p( dir_out )
        file_out = dir_out + '/' + '_'.join(( 'mean', model.lower(),
		dtg_start, dtg_end+'.nc' ))
        dbg( file_out )

        #_Initialize output file
        ncdf = Dataset( file_out, 'w', format='NETCDF3_CLASSIC')
	
	#_Get list of variables
	vars = lt.unique( mean.variable )
	fhrs = lt.unique( mean.fhr )

##	#_Subset data to appropriate times (model is 
##	# limiting factor,  so don't bother using epochs) 
##	aero_sub = subset( aeronet, dtg=aero_dtgs )
	nt = len( fhrs ) #_Is this ok?

	#_Add dtgs used
        ncdf.__setattr__( 'dtg_valds', ','.join( dtg_valds ) )

        #_Initialize model dimensions
        lats = ref_rec.values.lat
        lons = ref_rec.values.lon
	ncdf.createDimension( 'lat', lats.size )
	ncdf.createDimension( 'lon', lons.size )
	ncdf.createDimension( 'fhr', nt ) 
 

	#_Write coordinates
	cdf = ncdf.createVariable( 'lat', 'f4', ('lat',) )
	cdf[:] = lats
	cdf = ncdf.createVariable( 'lon', 'f4', ('lon',) )
	cdf[:] = lons
	cdf = ncdf.createVariable( 'fhr', 'f8', ('fhr',) )
	cdf[:] = fhrs 

	if hasattr( ref_rec.values, 'member' ):
		members = ref_rec.values.member
		ncdf.createDimension( 'member', members.size )

		#_add member names as global attribute
		ncdf.member = ','.join( members )

	#_Loop over variables and write mean to file
	for v in vars:
		#_subset record, the sort by dtg_vald
		mean_var = subset( mean, variable=v )
		mean_srt = sort_rec( mean_var, order=('fhr',) )

		#_Get mean at each fhr, N CAN VARY DEPENDING ON THRESHOLDS
		N = None
		for rec in mean_srt: 
			tmp_N = rec.values.n.copy()
			atn, atv = get_attr( rec )
			tmp_mean = rec.values / tmp_N 
			rec.values = cl.var( tmp_mean, attrn=atn, attrv=atv )
		
			#_Initialize map of number of events
			if N == None:			
				shp = list( tmp_N.shape )
				N = tmp_N.reshape( shp.insert( 0, 1 ) ) 
			else:
				N = np.append( N, tmp_N, axis=0 )	
			
		#_Join into single array
		mean_stk = join_values( mean_srt, newattr='fhr', newattrv=fhrs )

		atn.insert( 0, 'fhr' ) 
		dimname = tuple( atn )
		
		#_Sort records, then join along time and member axes
	        cdf = ncdf.createVariable( v, 'f4', dimname, fill_value=-9999. )
		cdf[:] 		= mean_stk 
		cdf.units	= lt.unique( mean_var.units, unique=True ) 
		cdf.long_name 	= lt.unique( mean_var.long_name, unique=True ) 

		#_Write map of N
		var_N = v + '_n'
	        cdf = ncdf.createVariable( var_N, 'i4', dimname )
		cdf[:] = N

        ncdf.close()
	lt.make_readable( file_out )
	return file_out

def read_period_geostat( file ):
	'''
	read ensemble geometric mean produced by write_period_geomu

	file	: str,	path to input file
	'''
	from netCDF4 import Dataset

	ncdf = Dataset( file, mode='r', format='NETCDF3_CLASSIC' )

	model = ncdf.model
	members = ncdf.member.split(',')
	lats = ncdf.variables['lat'][:]
	lons = ncdf.variables['lon'][:]
	try:
		fhrs = ncdf.variables['fhr'][:].tolist()
	except KeyError:
		import re
		res = re.search('geostd_([\w\d-]+)_([\w\d]+)-f(\d{3}).nc',file)
		fhrs = [int(res.group(3))]
	dtgs = ncdf.dtg_valds

	dtg_str = dtgs[0]
	dtg_end = dtgs[-1]

	dims = []
	for dim in ncdf.dimensions: dims.append( dim )

	records = cl.model_object()
	aod_dict = lm.aod_vars()

	nens = len( members )
	nx = lons.size
	ny = lats.size

	for variable in ncdf.variables:
	    if variable in dims: continue
	    if variable not in aod_dict: continue
	    units = aod_dict[variable]['units']
	    long_name = aod_dict[variable]['long_name']
	    for fhr in fhrs:
	  	t = fhrs.index(fhr)
 
		#_make list of attribute names for cl.var()
		attrn = list( ncdf.variables[ variable ].dimensions )

	    	#_make list of attribute values for cl.var()
	    	attrv = []
	    	for dim in attrn:
			if dim == 'fhr':
				continue
			elif dim == 'member': 
				attrv.append( np.array( members ))
				midx = attrn.index( 'member' )
			else: 
				attrv.append( ncdf.variables[ dim ][:] )

		#_read data from file
		data = ncdf.variables[variable][t]
		data = cl.var( data, attrv=attrv, attrn=attrn )
		dimsize = data.shape
	
		fidx = attrn.index('fhr')
		del attrn[fidx]

		#_put into recarray
		records.resize( records.size + 1 )
		records[-1] = ( data, model, dtg_str, dtg_end, 
			fhr, 'none', variable, nens, attrn, dimsize,
			units, long_name )
	
	return records

def write_period_own( recs, path='.', label=None, **kwargs ):
	'''
	recs	: np.recarray,	records from read_ensperiod()
				records contain geometricstd values over period
	path	: str,		where to dump these files

	(nmember, ny, nx) x nspec model variables, averaged in time
	(npoint) x nmode aeronet observations, averaged for same period 
       
	mean should be cl.cum_mean() object produced by running_summation()

	Any attempts at pairwise means have to be done in RUNNING_SUMMATION() 
        '''
	import libaeronet as la
        from netCDF4 import Dataset

	fhrs = lt.unique( recs.fhr )
	if len(fhrs) != 2: raise ValueError, 'only takes two fhr values'
	
        #_initialize lists to use for averaging
	dtg_valds = []
	dtg_inits = []

	#_get list of every possible dtg
	dtg_valds = lt.unique( recs.dtg_vald )
	dtg_inits = lt.unique( recs.dtg_init )

	#_use first record as a metadata reference
	ref_rec = recs[0]

	#_get name of ensemble
	model = lt.unique( recs.model, unique=True ) 

	is_ens = True if hasattr( ref_rec.values, 'member' ) else False

	#_file out
	dtg_start = dtg_inits[0]
	dtg_end	= dtg_inits[-1]

	dir_out = '/'.join(( path, model.upper(), dtg_start[:6], 'own' )) 
	if label == None:
	        file_out = dir_out + '/' + '_'.join(( 'own', model.lower(),
			dtg_start, dtg_end+'.nc' ))
	else:
       		file_out = dir_out + '/' + '_'.join(( 'own', model.lower(),
			label+'.nc' ))
	lt.mkdir_p( dir_out )
        dbg( file_out )

        #_initialize output file
        ncdf = Dataset( file_out, 'w', format='NETCDF3_CLASSIC')
	
	#_add dtgs used
        ncdf.__setattr__( 'dtg_valds', ','.join( dtg_valds ) )

        #_initialize model dimensions
        lats = ref_rec.values.lat
        lons = ref_rec.values.lon
	nx = lons.size
	ny = lats.size
	ncdf.createDimension( 'lat', ny )
	ncdf.createDimension( 'lon', nx )

	#_write coordinates
	cdf = ncdf.createVariable( 'lat', 'f4', ('lat',))
	cdf[:] = lats
	cdf = ncdf.createVariable( 'lon', 'f4', ('lon',))
	cdf[:] = lons

	members = ref_rec.values.member
	members = np.append( members, model.upper() )
	nens = members.size
	ncdf.createDimension( 'member', nens )

	#_find index of members
	midx = ref_rec.dimname.index('member')

	#_add member names as global attribute
	ncdf.member = ','.join( members )
	ncdf.__setattr__( 'model', model )

	#_get list of variables
	vars = lt.unique( recs.variable )
	expected = expected_vars()

	#_loop over variables and write mean to file
	dimname = ('member','lat','lon')
	for v in vars:
		dbg(v)

		#_pull out variables and see what valid times are available
		var0 = subset( recs, variable=v, fhr=fhrs[0] )
		var1 = subset( recs, variable=v, fhr=fhrs[1] )
		vald0 = lt.unique(var0.dtg_vald)
		vald1 = lt.unique(var1.dtg_vald)
		dtg_avail = lt.intersection([ vald0, vald1 ])

		#_remove records where only one is available
		var0 = subset( var0, dtg_vald=dtg_avail )
		var1 = subset( var1, dtg_vald=dtg_avail )
		var0 = sort_rec( var0, order=('dtg_vald',) )
		var1 = sort_rec( var1, order=('dtg_vald',) )

		stk0 = join_values( var0 )
		stk1 = join_values( var1 )

		#_get number of vald times used
		n = stk0.shape[0]
		
		#_take mean of difference
		diff = (stk0 - stk1).mean(axis=0)

		#_take mean of difference of ensemble mean
		diffmu = (stk0.mean(axis=midx+1) 
			- stk1.mean(axis=midx+1)).mean(axis=0)
		own = np.append( diff, diffmu.reshape(1,ny,nx), axis=0 )

		#_mask things that shouldn't be
		for member in members:
			if v not in expected[member]:
				midx0 = members.tolist().index(member)
				own[midx0].mask = True

	    	#_put data into file
		cdf = ncdf.createVariable( v, 'f4', dimname,fill_value=-9999 )
		cdf[:] = own 
		cdf.n = n

        ncdf.close()
	lt.make_readable( file_out )
	return file_out

def write_period_geostd( recs, geomu, path='.', label=None, **kwargs ):
	'''
	recs	: np.recarray,	records from read_ensperiod()
				records contain geometricstd values over period
	path	: str,		where to dump these files

	(nmember, ny, nx) x nspec model variables, averaged in time
	(npoint) x nmode aeronet observations, averaged for same period 
       
	mean should be cl.cum_mean() object produced by running_summation()

	Any attempts at pairwise means have to be done in RUNNING_SUMMATION() 
        '''
	import libaeronet as la
        from netCDF4 import Dataset
	from scipy import stats
	percentile = stats.scoreatpercentile
	
        #_initialize lists to use for averaging
	dtg_valds = []
	dtg_inits = []

	#_get list of every possible dtg
	dtg_valds = lt.unique( recs.dtg_vald )
	dtg_inits = lt.unique( recs.dtg_init )

	#_use first record as a metadata reference
	ref_rec = recs[0]

	#_get name of ensemble
	model = lt.unique( recs.model, unique=True ) 

	is_ens = True if hasattr( ref_rec.values, 'member' ) else False

	#_file out
	dtg_start = dtg_inits[0]
	dtg_end = dtg_inits[-1]

	dir_out = '/'.join(( path, model.upper(), dtg_start[:6], 'geostd' )) 
	if label == None:
	        file_out = dir_out + '/' + '_'.join(( 'geostd', model.lower(),
			dtg_start, dtg_end+'.nc' ))
	else:
       		file_out = dir_out + '/' + '_'.join(( 'geostd', model.lower(),
			label+'.nc' ))
	lt.mkdir_p( dir_out )
        dbg( file_out )

        #_initialize output file
        ncdf = Dataset( file_out, 'w', format='NETCDF3_CLASSIC')
	
	#_add dtgs used
        ncdf.__setattr__( 'dtg_valds', ','.join( dtg_valds ) )

        #_initialize model dimensions
	fhrs = lt.unique( recs.fhr )
        lats = ref_rec.values.lat
        lons = ref_rec.values.lon
	nx = lons.size
	ny = lats.size
	nt = len( fhrs )
	ncdf.createDimension( 'lat', ny )
	ncdf.createDimension( 'lon', nx )
	ncdf.createDimension( 'fhr', nt )

	#_write coordinates
	cdf = ncdf.createVariable( 'lat', 'f4', ('lat',))
	cdf[:] = lats
	cdf = ncdf.createVariable( 'lon', 'f4', ('lon',))
	cdf[:] = lons
	cdf = ncdf.createVariable( 'fhr', 'i4', ('fhr',))
	cdf[:] = fhrs 

	members = ref_rec.values.member
	members = np.append( members, model.upper() )
	nens = members.size
	ncdf.createDimension( 'member', nens )

	#_find index of members
	midx = ref_rec.dimname.index('member')

	#_add member names as global attribute
	ncdf.member = ','.join( members )
	ncdf.__setattr__( 'model', model )

	#_get summs of the np.ln(values / geomu)**2
	#_initialize summation object
	summ = cl.cum_mean()
	summ = running_summation( recs, summ, geomu=geomu )

	dimname = list(summ.dimname[0])
	dimname.insert(0,'fhr')

	#_get list of variables
	vars = lt.unique( recs.variable )
	del recs

	expected = expected_vars()

	#_loop over variables and write mean to file
	for v in vars:
	    dbg(v)
	    geos = np.ma.masked_array( np.zeros((nt,nens,ny,nx)))
	    sumz = np.zeros((nt,nens,ny,nx))
	    root = np.zeros((nt,nens,ny,nx))
	    n    = np.zeros((nt,nens,ny,nx))
	    for fhr in fhrs:
		t = fhrs.index(fhr)

		#_subset record, the sort by dtg_vald
		sum_by_var = subset( summ, variable=v, fhr=fhr, unique=True )

		#_calculate geometric stdev
		data 	= sum_by_var.values
		sqrt 	= np.sqrt( data / data.n  )
		exp 	= np.exp( sqrt )

		sumz[t] = data
		root[t] = sqrt
		geos[t] = exp
		n[t]	= data.n

	    #_mask things that shouldn't be
	    for member in members:
		if v not in expected[member]:
			midx = members.tolist().index(member)
			geos[:,midx,:,:].mask = True

	    dbg(( fhrs ))
	    if v == 'fine_aod' or v == 'coarse_aod':
		dbg((v, 'G:',geos[t,:,89,179]))
		dbg((v, 'S:',sumz[t,:,89,179]))
		dbg((v, 'R:',root[t,:,89,179]))
		dbg((v, 'N:',n[t,:,89,179]))
		dbg((v, 'M:',members))

	    #_sort records, then join along time and member axes
	    cdf = ncdf.createVariable( v, 'f4', dimname, fill_value=-9999 )
	    cdf[:] = geos
	    cdf = ncdf.createVariable( v+'_n', 'f4', dimname, fill_value=-9999 )
	    cdf[:] = n 
	    cdf = ncdf.createVariable( v+'_sum','f4',dimname, fill_value=-9999 )
	    cdf[:] = sumz 
	    cdf = ncdf.createVariable( v+'_sqr','f4',dimname, fill_value=-9999 )
	    cdf[:] = root 

        ncdf.close()
	lt.make_readable( file_out )
	return file_out

def write_period_geomu( recs, path='.', label=None, **kwargs ):
	'''
	summ	: cl.cum_mean(),if made pairwise, it is 
				done elsewhere
	path	: str,		Where to dump these files

	(nmember, ny, nx) x nspec model variables, averaged in time
	(npoint) x nmode aeronet observations, averaged for same period 
       
	mean should be cl.cum_mean() object produced by running_summation()

	Any attempts at pairwise means have to be done in RUNNING_SUMMATION() 
        '''
	import libaeronet as la
        from netCDF4 import Dataset
	from scipy import stats
	percentile = stats.scoreatpercentile
	
        #_initialize lists to use for averaging
	dtg_valds = []
	dtg_inits = []

	#_get list of every possible dtg
	dtg_valds = lt.unique( recs.dtg_vald )
	dtg_inits = lt.unique( recs.dtg_init )

	#_use first record as a metadata reference
	ref_rec = recs[0]

	#_get name of ensemble
	model = lt.unique( recs.model, unique=True ) 

	is_ens = True if hasattr( ref_rec.values, 'member' ) else False

	#_file out
	dtg_start = dtg_inits[0]
	dtg_end = dtg_inits[-1]

	dir_out = '/'.join(( path, model.upper(), dtg_start[:6], 'geomu' )) 
	if label == None:
	        file_out = dir_out + '/' + '_'.join(( 'geomu', model.lower(),
			dtg_start, dtg_end+'.nc' ))
	else:
       		file_out = dir_out + '/' + '_'.join(( 'geomu', model.lower(),
			label+'.nc' ))
	lt.mkdir_p( dir_out )
        dbg( file_out )

        #_initialize output file
        ncdf = Dataset( file_out, 'w', format='NETCDF3_CLASSIC')
	
	#_get list of variables
	vars = lt.unique( recs.variable )

	#_add dtgs used
        ncdf.__setattr__( 'dtg_valds', ','.join( dtg_valds ) )

        #_initialize model dimensions
        lats = ref_rec.values.lat
        lons = ref_rec.values.lon
	fhrs = lt.unique( recs.fhr )
	nt = len( fhrs )
	nx = lons.size
	ny = lats.size
	ncdf.createDimension( 'lat', lats.size )
	ncdf.createDimension( 'lon', lons.size )
	ncdf.createDimension( 'fhr', nt ) 

	#_Write coordinates
	cdf = ncdf.createVariable( 'lat', 'f4', ('lat',) )
	cdf[:] = lats
	cdf = ncdf.createVariable( 'lon', 'f4', ('lon',) )
	cdf[:] = lons
	cdf = ncdf.createVariable( 'fhr', 'f8', ('fhr',) )
	cdf[:] = fhrs 

	members = ref_rec.values.member
	members = np.append( members, model.upper() )
	nens = members.size
	ncdf.createDimension('member',nens)

	#_find index of members
	dimname = list( ref_rec.dimname )
	dimname.insert(0,'fhr') 
	midx = dimname.index('member')

	ncdf.__setattr__( 'model', model )

	#_loop over variables and write mean to file
	for v in vars:
	    geomean = np.ma.masked_array( np.zeros((nt,nens,ny,nx)) )
	    for fhr in fhrs:
		dbg(( v, fhr ), l=5)
		t = fhrs.index(fhr)

		#_subset record, the sort by dtg_vald
		geomu_var = subset( recs, variable=v, fhr=fhr )

		#_join into single array
		tmp = join_values( geomu_var )

		#_zero substitution...
		mask = tmp.mask.copy()
		tidx = np.where( tmp.data < 0.001 )
		tmp[tidx] = 0.001
		tmp.mask = mask

		#_get ensemble mean 
		tmpmu = tmp.mean( axis=midx )

		#_calculate gemoemtric mean
		geomu = stats.mstats.gmean( tmp, axis=0 )
		geomumu = np.ma.masked_array(stats.mstats.gmean(tmpmu,axis=0 ))
	    	geomu = np.ma.MaskedArray.filled( geomu, -9999. )
		geomu = np.append( geomu, geomumu.reshape(1,ny,nx), axis=0 )
	
		geomean[t] = geomu[:]

	    #_sort records, then join along time and member axes
	    cdf = ncdf.createVariable( v, 'f4', dimname, fill_value=-9999. )
	    cdf[:] = geomean

	#_add member names as global attribute
	ncdf.member = ','.join( members )

        ncdf.close()
	lt.make_readable( file_out )
	return file_out

def plot_enaaps_mean( file, dir_out='.', cb=False ):
        '''
	GENERALIZE THIS BACK INTO PLOT PERIOD MEAN
        Take netcdf file written by write_period_mean and plot
        overlay of AERONET on MODEL

        file_out = dir_out + '/' + model.lower() + '_mean_' + dtg_start \
                        + '_' + dtg_end + '_f' + str(fhr).zfill(3) + '.nc'

	At some point, allow all to be separate files, or woult plot rec suffice
        '''
        import re
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        from netCDF4 import Dataset
	import libcmap as lc
	
        #_Create output directory
	dir_out = dir_out + '/overlays'
        lt.mkdir_p( dir_out )

        ncdf = Dataset( file, mode='r', format='NETCDF3_CLASSIC' )
        grid = [-80,81,-180,181]
        grid = [0,0,0,0]
        corn = [-90,90,-180,180]
        proj = 'cyl'
        delt = 40
        levs = [ .1, .2, .3, .4, .5, .6, .7 ]
        levs_std = [ 0., 0.0015, 0.0030, 0.0045, 0.0060, 0.0075, 0.0090 ]
        levs_stdn= [ 0., 0.015, 0.030, 0.045, 0.060, 0.075, 0.090 ]
        cmap = lc.rgbcmap( 'aeronet', N=256 )
        lats = ncdf.variables['lat'][:]
        lons = ncdf.variables['lon'][:]
        fhrs = ncdf.variables['fhr'][:]
        members = ncdf.members.split(',')
        nens = len( members )

        #_Get fhr, dtg_start, dtg_end, ensemble name from filename str
        re_file = re.compile('([\d\w-]+)_mean_(\d{10})_(\d{10})_([\w\d-]+).nc$')
        results = re_file.search( file )
        model   = results.group(1)
        dtg_str = results.group(2)
        dtg_end = results.group(3)
        label   = results.group(4)

        #_Read in AERONET data, total
        aero = ncdf.variables['aeronet_total'][:]
        aero_lat = ncdf.variables['aeronet_lat'][:]
        aero_lon = ncdf.variables['aeronet_lon'][:]

        #_Set max and min aeronet tau values and constrain outliers
        taumax = levs[-1]
        taumin = 0.0
        mind = np.where( aero < taumin )
        xind = np.where( aero > taumax )
        aero[mind] = taumin
        aero[xind] = taumax

        #_Make up levels once we make some headway
        std_levs = {
                'sulfate'       : [],
                'dust'          : [],
                'smoke'         : [],
                'seasalt'       : [],
                'total'         : []    }

        #_Loop over specs and output overlay plots
        for variable in ncdf.variables:
           #_Skip non-aod vars
           if not re.search( r'_aod', variable ): continue

           #_Read in forecast data to local variable 
           aod = ncdf.variables[ variable ][:]

           #_Mask missing fields (should be automatic with FillValue attr)
           aod = np.ma.masked_where( aod == -9999., aod )

           #_Limit forecast data to min and max tau values
           min = np.where( aod.data < taumin )
           max = np.where( aod.data > taumax )
           if len( min[0] ) > 0: aod.data[min] = taumin
           if len( max[0] ) > 0: aod.data[max] = taumax

           #_get name of dimensions
           dims = ncdf.variables[ variable ].dimensions
           ens_idx = dims.index( 'member' )
           fhr_idx = dims.index( 'fhr' )        #_Don't overc copmlicate this

           for fhr in fhrs:
                t = np.where( fhrs == fhr )[0]
                aod_tmp = aod[t].copy().squeeze()

                fhr = str( int( fhr )).zfill(3)

                #_Filename for ICAP MEAN and STDV plots 
                file_out = dir_out + '_'.join(( '/overlay', model, dtg_str, 
			dtg_end, 'f'+str(fhr), label, variable, label+'.'+ifmt))

                #_Calculate number of rows for two column plots
###             extra = 1 if (nens+2) % ncol else 0 #_Nens+1 because of stdv plt
###             nrow = int( (nens+2) / ncol ) + extra
		nrow = 3
                ncol = 1 

                #_Initialize plotting map
                m = draw_map( grid, corn, proj, delt )
                x, y = latlon_to_xy( m, lats, lons )

                #_Initialize figure
                fig = plt.figure()
                nfig = 1
                a = fig.add_subplot( nrow, ncol, nfig )

                #_Plot model data
                d = aod_tmp.mean( axis=0 )
                c = m.contourf( x, y, d.data, levels=levs, cmap=cmap )
		if cb:
	 		bar = plt.colorbar( c, orientation='horizontal', 
				shrink=0.4, aspect=40, pad=0.05 )
			[ t.set_fontsize(4) for t in bar.ax.get_xticklabels() ] 

                #_Plot aeronet data over mean
                m = draw_map( grid, corn, proj, delt )
                xs, ys = m( aero_lon, aero_lat )
                sc = m.scatter( xs, ys, s=3., c=aero, marker='s',
                        linewidth=0.25, cmap=cmap, zorder=10 )

                #_Set mean panel title
                title = model.upper() + ' ' + variable.upper() + ' ' + str(fhr)
                a.set_title( title, size='xx-small')

                #_For stdv plot
                nfig +=1

                #_Plot standard deviation
                a = fig.add_subplot( nrow, ncol, nfig )
                m = draw_map( grid, corn, proj, delt )
                d = aod_tmp.std( axis=0 )
		max_std = d.max()
		idx_max = np.where( d.data >= levs_std[-1] )
		d[idx_max] = levs_std[-1]
                cs = m.contourf( x, y, d, levels=levs_std )#, cmap=cm.jet )
		if cb:
	 		bar = plt.colorbar( cs, orientation='horizontal', 
				shrink=0.4, aspect=40, pad=0.05 )
			[ t.set_fontsize(4) for t in bar.ax.get_xticklabels() ] 

                #_Set member title
                title_str = 'STDV ' + variable.upper() + ' ' + str(fhr) \
                        + ' MAX: '
                fmt = "%20s %7.4f"
                title = fmt % ( title_str, max_std )
                a.set_title( title, size='xx-small' )

		#_For normalized stdv plot
                nfig +=1

                #_Plot Normalized standard deviation
                a = fig.add_subplot( nrow, ncol, nfig )
                m = draw_map( grid, corn, proj, delt )
                d = aod_tmp.std( axis=0 ) / aod_tmp.mean( axis=0 )
		max_std = d.max()
		idx_max = np.where( d.data >= levs_stdn[-1] )
		d[idx_max] = levs_stdn[-1]
                cs = m.contourf( x, y, d, levels=levs_stdn )#, cmap=cm.jet )
		if cb:
	 		bar = plt.colorbar( cs, orientation='horizontal', 
				shrink=0.4, aspect=40, pad=0.05 )
			[ t.set_fontsize(4) for t in bar.ax.get_xticklabels() ] 

                #_Set member title
                title_str = 'NORMALIZED STDV ' + variable.upper() + ' ' \
			+ str(fhr) + ' MAX: '
                fmt = "%20s %7.4f"
                title = fmt % ( title_str, max_std )
                a.set_title( title, size='xx-small' )

                dbg( file_out )
                plt.savefig( file_out, dpi=(400) )
                plt.close()

        #_Close ncdf file

def read_period_own( file ):
	'''
	Read in file generated by write_period_mean()
	file	: str,	Path to file
	'''
        from netCDF4 import Dataset
	import re

	dbg( file )

	#_Open ncdf file object
        ncdf = Dataset( file, mode='r', format='NETCDF3_CLASSIC' )

	#_Metadata for variables (because it's not in the file for some reason)
	aod_dict = lm.aod_vars()

	dims = {}
	for dim in ncdf.dimensions: dims[dim] = ncdf.variables	

        #_get fhr, dtg_start, dtg_end, ensemble name from filename str
        re_file = re.compile('own_([\d\w-]+)_([\w\d-]*)-f(\d{3})-(\d{3}).nc$')
        results = re_file.search( file )
        model   = results.group(1).upper()
        dtg_str = ncdf.dtg_valds.split(',')[0]
        dtg_end = ncdf.dtg_valds.split(',')[-1]
	label 	= results.group(2)
	fhr0	= results.group(3)
	fhr1	= results.group(4)

	#_check netcdf file for member for ensembles
	#_DOES NOT LOOK FOR MISSING ENSEMBLES FOR NENS HERE
        members = ncdf.member.split(',')
        nens 	= len( members )

	dtype = [('values', cl.var),
		('model', 'a16'),
		('dtg_init', 'a10'),
		('dtg_vald', 'a10'),
		('fhr0', 'i4'),
		('fhr1', 'i4'),
		('region', 'a10'),
		('variable', 'a20'),
		('ensemble', 'i4'),
		('dimname', tuple),
		('dimsize', tuple),
		('units', 'a30'),
		('long_name', 'a40') ]

	own = np.recarray((0,),dtype=dtype)

	#_loop over variables to put in records
        for variable in ncdf.variables:
		if variable in ncdf.dimensions: continue
		#_names in this file have _mean and _sigma appended,
		# remove them and use it as keys

		#_do not store dimensions as variables
	    	if variable in ncdf.dimensions: continue 

	    	#_make list of attribute names for cl.var()
	    	attrn = list( ncdf.variables[ variable ].dimensions )
		
	    	#_make list of attribute values for cl.var()
	    	attrv = []
	    	for dim in attrn:
			if dim == 'member': 
				attrv.append( np.array( members ))
			else: 
				attrv.append( ncdf.variables[ dim ][:] )

		#_put data in cl.var() class
	        field = ncdf.variables[ variable ]
		field = cl.var( field, attrv=attrv, attrn=attrn )

		dimsize = field.shape
		units	= aod_dict[variable]['units']
		long_name = aod_dict[variable]['long_name']
		
		#_add ot records
		own.resize( own.size + 1 )
		own[-1] = ( field, model, dtg_str, dtg_end, 
			fhr0, fhr1, label, variable, nens, attrn, dimsize,
			units, long_name )
	ncdf.close()
	return own 

def read_period_musig( file ):
	'''
	Read in file generated by write_period_mean()
	file	: str,	Path to file
	'''
        from netCDF4 import Dataset
	import re

	dbg( file )

	#_Open ncdf file object
        ncdf = Dataset( file, mode='r', format='NETCDF3_CLASSIC' )

	#_Metadata for variables (because it's not in the file for some reason)
	aod_dict = lm.aod_vars()

	dims = {}
	for dim in ncdf.dimensions: dims[dim] = ncdf.variables	

        #_get fhr, dtg_start, dtg_end, ensemble name from filename str
        re_file =re.compile('musig_([\d\w-]+)_([\w\d-]*).nc$')
        results = re_file.search( file )
        model   = results.group(1).upper()
        dtg_str = ncdf.dtg_valds.split(',')[0]
        dtg_end = ncdf.dtg_valds.split(',')[-1]
	label 	= results.group(2)

	#_check netcdf file for member for ensembles
	#_DOES NOT LOOK FOR MISSING ENSEMBLES FOR NENS HERE
        members = ncdf.member.split(',')
        nens = len( members )

	fhrs = ncdf.variables['fhr'][:].tolist()
	mus = cl.model_object()	
	sigmas = cl.model_object()	
	#_loop over variables to put in records
        for variable in ncdf.variables:
	    for fhr in fhrs:
		if variable in ncdf.dimensions: continue
		#_names in this file have _mean and _sigma appended,
		# remove them and use it as keys
		tmp = variable.split('_')
		specie = '_'.join((tmp[:2]))

		t = fhrs.index(fhr)

		#_do not store dimensions as variables
	    	if variable in ncdf.dimensions: continue 

	    	#_make list of attribute names for cl.var()
	    	attrn = list( ncdf.variables[ variable ].dimensions )
		didx = attrn.index('fhr'); del attrn[didx]
		
	    	#_make list of attribute values for cl.var()
	    	attrv = []
	    	for dim in attrn:
			if dim == 'member': 
				attrv.append( np.array( members ))
			else: 
				attrv.append( ncdf.variables[ dim ][:] )

		#_put data in cl.var() class
	        field = ncdf.variables[ variable ][t]
		field = cl.var( field, attrv=attrv, attrn=attrn )

		dimsize = field.shape
		units = aod_dict[specie]['units']
		long_name = aod_dict[specie]['long_name']
		
		#_add ot records
		if tmp[2] == 'mean':
			mus.resize( mus.size + 1 )
			mus[-1] = ( field, model, dtg_str, dtg_end, 
				fhr, label, specie, nens, attrn, dimsize,
				units, long_name )
		elif tmp[2] == 'sigma':
			sigmas.resize( sigmas.size + 1 )
			sigmas[-1] = ( field, model, dtg_str, dtg_end, 
				fhr, label, specie, nens, attrn, dimsize,
				units, long_name )
	ncdf.close()
	return mus, sigmas 

def read_period_median( file ):
	'''
	Read in file generated by write_period_mean()
	file	: str,	Path to file
	'''
        from netCDF4 import Dataset
	import re

	dbg( file )

	#_Open ncdf file object
        ncdf = Dataset( file, mode='r', format='NETCDF3_CLASSIC' )

	#_Metadata for variables (because it's not in the file for some reason)
	aod_dict = lm.aod_vars()

	dims = {}
	for dim in ncdf.dimensions: dims[dim] = ncdf.variables	

        #_get fhr, dtg_start, dtg_end, ensemble name from filename str
        re_file =re.compile('median_([\d\w-]+)_([\w\d-]*).nc$')
        results = re_file.search( file )
        model   = results.group(1).upper()
        dtg_str = ncdf.dtg_valds.split(',')[0]
        dtg_end = ncdf.dtg_valds.split(',')[-1]
	label 	= results.group(2)

	#_check netcdf file for member for ensembles
	#_DOES NOT LOOK FOR MISSING ENSEMBLES FOR NENS HERE
        members = ncdf.member.split(',')
        nens = len( members )

	fhrs = ncdf.variables['fhr'][:].tolist()
	records = cl.model_object()

	specs = {}
	#_loop over variables to put in records
	re_ensmean = re.compile( '_ensmean' )
        for variable in ncdf.variables:
	    for fhr in fhrs:
		t = fhrs.index(fhr)

		#_do not store dimensions as variables
	    	if variable in ncdf.dimensions: continue 

	    	#_make list of attribute names for cl.var()
	    	attrn = list( ncdf.variables[ variable ].dimensions )
		didx = attrn.index('fhr'); del attrn[didx]
		
	    	#_make list of attribute values for cl.var()
	    	attrv = []
	    	for dim in attrn:
			if dim == 'member': 
				attrv.append( np.array( members ))
			else: 
				attrv.append( ncdf.variables[ dim ][:] )
	  
	        data = ncdf.variables[ variable ][t]

		#_put data in cl.var() class
		data = cl.var( data, attrv=attrv, attrn=attrn )
		dimsize = data.shape
		units = aod_dict[variable]['units']
		long_name = aod_dict[variable]['long_name']
		
		#_add ot records
		records.resize( records.size + 1 )
		records[-1] = ( data, model, dtg_str, dtg_end, 
			fhr, label, variable, nens, attrn, dimsize,
			units, long_name )
	ncdf.close()
	return records

def expected_vars():
	''' return list of expected vars by model '''
	exp = {}
	for member in mod_dict:
		exp[member] = [s.lower() for s in mod_dict[ member ]['specs']]

		#_add coarse_aod to all model lists	
		if member != 'NGAC' and member != 'NMMB':
			exp[member].append('coarse_aod')
		
		#_add fine_aod to all but NGAC and NMMB	
		if member != 'NGAC' and member != 'NMMB':
			exp[member].append('fine_aod')
	return exp

def mask_members( records, regional=False ):
	''' 
	look at record array for members, return number present

	regional : bool,	If True, read in predefined mask 
				and apply it to data

	mask should be a netcdf file, by species, 4 dimensional,
	tuned to model strengths, perhaps seasonal. (nfhr, nmember, ny, nx)

	'''
	import re
	re_mem = re.compile( '_m\d+$' )
	exp = expected_vars()
	for rec in records:
		members = rec.values.member.tolist()
	
		for member in members:
			mkey = re_mem.sub( '', member )
			if rec.variable not in exp[mkey]:
				dbg(( exp[mkey], mkey, rec.variable ), l=5)
	
				#_mask models that aren't expected
				midx = members.index(member)
				if type(rec.values[midx].mask) != np.bool_:
					rec.values[midx].mask[:] = True
				else:
					rec.values[midx].mask = True
						
				rec.values[midx][:] = -9999.

				#_update ensemble count
				rec.ensemble = check_members( rec.values )

def check_members( values ):
	''' look at record array for members, return number present '''
	nens = 0
	for idx in xrange( values.shape[0] ):
		#_make zero array same shape as data
		miss = np.zeros( values[idx].shape )

		#_set missing to 1 and add them together
		miss[values[idx].mask] = 1
		miss = np.sum( miss )

		#_if more than 1/4 missing, do not count in ensemble
		size = float( values[idx].size )
		if miss / size < 0.25: nens += 1
	
	return nens	

def check_member_order( records ):
        ''' 
        completely throwaway function to do no good but give me piece of mind
        '''
        order = None 
        for rec in records:
                if order == None:
                        order = rec.values.member.tolist()
                else:
                        if order != rec.values.member.tolist():
                                raise RuntimeError, 'models out of order'


def find_nearest_dtg(epoch, dt=6, **kwargs):
	'''
	epoch	float,	unix time to find the nearest standard NAAPS dtg 
	dt		int,	timestep between NAAPS dtgs, assuming 00
	'''
	from libtools import epoch2dtg, dtg2epoch
	from time import localtime
	from numpy import arange

	#_00z for today
	to = localtime(epoch)
	arg = (to.tm_year, to.tm_mon, to.tm_mday)
	day = dtg2epoch('{0}{1:02d}{2:02d}00'.format(*arg))
	sec = arange(day, day+86401, dt*3600)
	return epoch2dtg(sec[abs(epoch - sec).argmin()], **kwargs)


def get_naapsaod_track(epoch, lat, lon, species='dust_aod', 
	dir_naaps=os.path.join(os.environ['PRODUCTS'], 'NRL', 'NVA_CLIMO1'),
	**kwargs):
	'''
	2016.01.19	WRS
	Pull out NAAPS AOD values along epoch/lat/lon
	'''
	from numpy import array, tile

	#_initialize output array
	aod = tile(None, epoch.size)
	fmt = os.path.join(dir_naaps, 'NAAPSAOD', '{0}/{1}_aod')
	dtg_last = 0

	#_loop over track
	for i, (e, y, x) in enumerate(zip(epoch, lat, lon)):

		#-get nearest dtg
		dtg = find_nearest_dtg(e)
		fname = fmt.format(dtg[:6], dtg)

		#_check if available
		if not os.path.exists(fname):
			raise IOError, 'Invalid filename {0}'.format(fname)

		#_read in file
		if dtg_last == dtg:
			pass
		else:
			naaps = read_naapsaod(fname, **kwargs)	
	
		#_select record
		naod = subset(naaps, dtg_vald=dtg, variable=species, unique=True)

		idx_y = abs(naod.values.lat - y).argmin()	
		idx_x = abs(naod.values.lon - x).argmin()	

		aod[i] = naod.values[idx_y, idx_x]

		dtg_last = dtg

	return aod
			

def get_naapsconc_xsect(epoch, lat, lon, species='dust', heights=False, 
	dir_naaps=os.path.join(os.environ['PRODUCTS'], 'NRL', 'NVA_CLIMO1'),
	valid_time=False, **kwargs):
	'''
	2016.01.19	WRS
	Pull out NAAPS crosssection along time/lat/lon
	epoch, lat, lon		ndarray,	Arrays containing coords
	dir_naaps			str,		Location of NAAPS data
	heights				bool,		Return vertical height coordinate?
	valid_time			bool,		Return valid time distance
	'''
	from numpy import array, cos, sqrt
	from libgeo import sig2pres, p2z
	from libtools import dtg2epoch

	#_initialize output field
	prof = []	
	ffmt = os.path.join(dir_naaps, 'NAAPS', '{0}/{1}_conc')
	dtg_last = None
	dt = []	#_distance from valid time
	dx = [] #_distance from nearest neighbor
	kpd = 111.2

	for i, (e, y, x) in enumerate(zip(epoch, lat, lon)):
		#_get nearest dtg for filename
		dtg = find_nearest_dtg(e)
		fname = ffmt.format(dtg[:6], dtg)
		
		#_check if it exists
		if not os.path.exists(fname):
			raise IOError, 'Invalid file name {0}'.format(fname)

		#_read in concentration file
		if dtg_last != dtg:
			conc = read_conc(fname, **kwargs)

		#_calculate distance to valid time
		dt.append(e - dtg2epoch(dtg))

		#_pull out lat/lon/concentrations
		varname = conc.variable.tolist()
		idx_lat = varname.index('lat')
		idx_lon = varname.index('lon')
		idx_con = varname.index('conc')
		idx_spe = varname.index('species')

		#_get coords
		lat = conc.values[idx_lat] 
		lon = conc.values[idx_lon] 

		#_pull out concentration values for desired specie
		spe = conc.values[idx_spe].tolist()
		con = conc.values[idx_con][spe.index(species)] 

		#_find the index for this location
		idx_y = abs(lat - y).argmin()
		idx_x = abs(lon - x).argmin()

		#_calc distance to gridpoint
		a = ((lat[idx_y]-y)*kpd)**2 + ((lon[idx_x]-x)*kpd*cos(lat[idx_y]))**2
		dx.append(sqrt(a))

		#_pull out column
		column = con[:, idx_y, idx_x]
		prof.append(column)

		dtg_last = dtg

	#_turn profile into one giant array
	prof = array(prof).T	

	#_return profile
	if heights and valid_time:
		idx_sig = varname.index('sigma')
		idx_sga = varname.index('sig_a')
		sigma = conc.values[idx_sig]
		p_top = conc.values[idx_sga].min()
		return prof, p2z(sig2pres(sigma, ptop=p_top)), array(dt), array(dx)
	elif heights:
		raise RuntimeError, 'Put this in later if needed'
	elif valid_time:
		raise RuntimeError, 'Put this in later if needed'
	else:
		return prof


def read_period_mean( file ):
	'''
	Read in file generated by write_period_mean()
	file	: str,	Path to file
	'''
        from netCDF4 import Dataset
	import re

	dbg( file )

	#_Open ncdf file object
        ncdf = Dataset( file, mode='r', format='NETCDF3_CLASSIC' )
	fhrs = ncdf.variables['fhr'][:]

	#_Metadata for variables (because it's not in the file for some reason)
	aod_dict = lm.aod_vars()

	dims = {}
	for dim in ncdf.dimensions: dims[dim] = ncdf.variables	

	#_check netcdf file for member for ensembles
	if hasattr( ncdf, 'member' ):
		#_DOES NOT LOOK FOR MISSING ENSEMBLES FOR NENS HERE
		is_ens = True
        	members = ncdf.member.split(',')
        	nens = len( members )

	else:
		is_ens = False
		nens = 0

        #_get fhr, dtg_start, dtg_end, ensemble name from filename str
        re_file =re.compile('mean_([\d\w-]+)_(\d{10})_(\d{10})_?([\w\d-]*).nc$')
        results = re_file.search( file )
        model   = results.group(1).upper()
        dtg_str = results.group(2)
        dtg_end = results.group(3)
	try: label = results.group(4)
	except AttributeError: label = 'NoLab'

	records = cl.model_object()

	specs = {}
	#_Loop over variables to put in records
        for variable in ncdf.variables:
	    #_Do not store dimensions as variables
	    if variable in ncdf.dimensions: continue 
	    if variable not in aod_dict: continue

	    #_Make list of attribute names for cl.var()
	    attrn = list( ncdf.variables[ variable ].dimensions )
	    idx = attrn.index( 'fhr' )
	    del attrn[idx]

	    #_Make list of attribute values for cl.var()
	    attrv = []
	    for dim in attrn:
		if dim == 'member': 
			attrv.append( np.array( ncdf.member.split(',') ))
		else: attrv.append( ncdf.variables[ dim ][:] )
	  

  	    #_Loop over fhr row and put into records
	    for fhr in fhrs:
		t = fhrs.tolist().index( fhr )
		fhr = int( fhr )

	        data = ncdf.variables[ variable ][t]

		#_Put data in cl.var() class
		data = cl.var( data, attrv=attrv, attrn=attrn )
		dimsize = data.shape
		units = aod_dict[variable]['units']
		long_name = aod_dict[variable]['long_name']

		#_Add ot records
		records.resize( records.size + 1 )
		records[-1] = ( data, model, dtg_str, dtg_end, 
			fhr, label, variable, nens, attrn, dimsize,
			units, long_name )
	ncdf.close()
	return records
