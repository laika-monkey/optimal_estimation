#############################################################################_80
#_NAME: 	libicap.py						       #
#_PURPOSE:	ICAP specific tools for plotting and metadata.  Some of this   #
#		overlaps with libnva.py, but with added items.  Please be      #
#		cautious						       #
#									       #
#_AUTHOR:	Walter R. Sessions, September 2011			       #
#									       #
#############################################################################_80
##########_GENERAL_#############################################################
################################################################################
import libnva 	as ln
import libtools as lt
import libmeta 	as lm
import libclass as cl
import numpy as np
import os

debug 		= lm.debug #_High number == Noisier logs
dir_prod 	= lm.dir_prod 

#_Metadata used for models and plotting
mod_dict = lm.models()
reg_dict = lm.regions()
plt_dict = lm.fields()
pnt_dict = lm.aeronet_points()

#############################################################################_80
##########_PLOTTING_############################################################
################################################################################

def plot_icap( records, **kwargs ):
	'''
	Used by ensemble_post to produce expected plot for ICAP website
	
	records	: cl.model_object(), recarray of model AODs
	dir_web	: str, root output directory
	nproc	: int, number of processors to use for plotting
	Passed a recarray, attempts to generate plots in accordance with ICAP
	'''
	plots_det = lm.fields_det()
	plots_ens = lm.fields_ens()

	#_Loop over plottable fields
	for field in plots_det:
		dbg(( 'deterministic', field ))
		#_Don't subset, allow all to plot

		#_Loop over each plotting regions 
		for reg in reg_dict:
			#_Get only the data for this region
			r_reg = ln.sub_region( records, region=reg )

			#_Deterministic plots (aod....and?)
	                ln.plot_forker( r_reg, field=field, **kwargs )

	#_Loop over plottable ensemble fields
	for field in plots_ens:
		dbg(( 'ensembles', field ))
		ens_idx = np.where( records.ensemble != 0 )	
		records_ens = records[ens_idx]

		#_Loop over each plotting regions 
		for reg in reg_dict:
			#_Get only the data for this region
			r_reg = ln.sub_region( records_ens, region=reg )

			#_Deterministic plots (aod....and?)
	                ln.plot_forker( r_reg, field=field, **kwargs ) 

def plot_availability( dtgs, dtge, fcst_finc=24, models=lm.current_icap(), 
	variables=['dust_aod'], path='.', **kwargs ):
	'''
	For a given period, creates a bar chart to easily reference when
	model forecasts are present or missing.

	Only check dust_aod for now
	'''
	import matplotlib.pyplot as plt

	lw = 30 

	#_Create list of forecast DTGs to check 
	dtg_list = []
	dtg = dtgs
	while dtg <= dtge:
		dtg_list.append( dtg )
		dtg = lt.newdtg( dtg, fcst_finc )
	nt = len( dtg_list )
	x = np.arange( nt )

	#_Create dictionary of model keys that correspond to a y values
	false_mask = [False]*nt
	y = {}
	models.insert( 0, 'ICAP' ) 
	for model in models:
		n = models.index( model ) 
		y[ model ] = np.ma.masked_array( [n]*nt )
		y[ model ].mask = false_mask 

	h, w = plt.figaspect( 10 )
	fig = plt.figure( figsize=(w,h) )

	#_Loop over dtgs and try to read forecasts
	for dtg in dtg_list:
		n = dtg_list.index( dtg )
		try: 	#_read in forecast
			fcst = read_icapfcst( dtg )
		except:	#_if the ICAP read function fails, assume all missing
			continue
		
		#_can expand to others later, but for now focus on dust
		for variable in variables:
			#_grab a single record to represent period
			fhr = 72	#_Use 72 hours as the testing sample
			rec = ln.subset( fcst, variable=variable, fhr=fhr, 
						unique=True )
			icap_mask = False
			#_Check for each model
			for member in rec.values.member:
				#_Assumes member index is 0
				m = rec.values.member.tolist().index( member )
				col = mod_dict[member]['color']

				mask = rec.values.mask[m]
				#_Test full masked array
				idx = np.where( mask == True )[0]
				if len(idx) < (mask.size/4):
					plt.hlines( [m+1],[n-.55],[n+.55],
						lw=lw, colors=col )
					continue

				#_If any models fail, mask ICAP
				icap_mask = True
			else:
				#_After all members done, plot ICAP	
				col = mod_dict['ICAP']['color']
				if not icap_mask: 
					plt.hlines( [0],[n-.55],[n+.55],
						lw=lw, colors='red' )


	x_max = 0
	lab_y = []
	for member in models:
		#_Some Y-labeling specifics
		if member == 'NGAC':
			lab_y.append( 'NGAC\n(dust only)')
                elif member == 'NMMB':
			lab_y.append( 'NMMB/BSC-CTM\n(dust and seasalt only)')
                else:
			lab_y.append( member )

	#_Limit size of plot area
	plt.xlim( 0, nt-1 )
        ax = plt.subplot( 111 )
        mjl = plt.MultipleLocator( nt / 5 )
        ax.xaxis.set_major_locator( mjl )
        mjl = plt.MultipleLocator( 90 )
        ax.yaxis.set_major_locator( mjl )
        ax.grid( True, zorder=-1 )

	#_make list of dtgs for x axis
        dtg_labels = dtg_list[ 0:nt: nt/5 ]
        tick_lab = lt.human_date( dtg_labels )

	#_make list of x values to put the dtg labels
        tick_loc = []
        [ tick_loc.append( dtg_list.index(p) ) for p in dtg_labels ]

        plt.xticks( tick_loc, tick_lab, size='x-small' ) 

	#_Setup y axis
        tick_loc = range( 0, len( models ), 1)
        plt.yticks( tick_loc, lab_y, size='x-small' )
        plt.ylim( -0.5, len( models )-0.5 )

	#_Write plot title
        plt.title( 'ICAP Forecast Availablity at NRL-MRY for '+ dtgs+'-'+dtge )

	#_Save file
	lt.mkdir_p( path )
	file_out = path + '/' + dtgs + '.availability.png'
        dbg( file_out )
        plt.savefig( file_out )
        plt.close()

#############################################################################_80
######_FILE-IO_#################################################################
################################################################################

def read_icapfcst( dtg, **kwargs ):
	kwargs['model'] = 'ICAP' 
	return ln.read_ensfcst( dtg, **kwargs )

def read_icapfcst_raw( dtg, remove_members=False, fhr=None,
	require=['total_aod','sulfate_aod','dust_aod','smoke_aod',
	'seasalt_aod'], members=lm.current_icap(), **kwargs ):
	'''
	Attempts to read in all centers' five day forecasts
	Will not merge without members in members having require specs

	fhr is a list of desired fhr
	'''
	dbg( members )
	aod = cl.model_object()
	if 'MACC' in members: #_READ MACC
		try:
			macc = read_maccfcst( dtg, require=require )
			if type( macc ) == np.recarray:
				aod = lt.merge(( aod, macc ))
		except:
			dbg( 'Error reading MACC', l=2 )

	if 'GEOS5' in members: #_READ GEOS5
		try:
	  		dtg_g = lt.newdtg( dtg, -2 )
	  		geos5 = read_geos5fcst( dtg_g, require=require )
			if type( geos5 ) == np.recarray:
				aod = lt.merge(( aod, geos5 ))
		except:
			dbg( 'Error reading GEOS5', l=2 )

	if 'MASINGAR' in members: #_READ MASINGAR
		try:
			masingar = read_masingarfcst( dtg, require=require )
			if type( masingar ) == np.recarray:
				aod = lt.merge(( aod, masingar ))
		except:
			dbg( 'Error reading MASINGAR', l=2 )

	if 'NGAC' in members: #_READ GFS NGAC-
		try:
			require_ngac = ['dust_aod']
			ngac = read_ngacfcst( dtg, require=require_ngac )
			if type( ngac ) == np.recarray:
				aod = lt.merge(( aod, ngac ))
		except:
			dbg( 'Error reading NGAC', l=2 )

	if 'NMMB' in members: #_READ GFS NGAC-
		try:
			require_nmmb = ['dust_aod','seasalt_aod']
			nmmb = read_nmmbfcst( dtg, require=require_nmmb )
			if type( nmmb ) == np.recarray:
				aod = lt.merge(( aod, nmmb ))
		except:
			dbg( 'Error reading NMMB', l=2 )
 
	if 'NAAPS' in members:	#_READ NAAPS DETERMINISTIC
		try:
			nva = ln.read_naapsfcst( dtg, **kwargs )
			if type( nva ) == np.recarray:
				aod = lt.merge(( aod, nva ))
		except:
			dbg( 'Error reading NAAPS', l=2 )

	if 'ENAAPS-NAV' in members: #_READ NAAPS MET ENSEMBLE
		try:
			ensemble = ln.read_ensfcst( dtg, model='ENAAPS-NAV' )
			if type( ensemble ) == np.recarray:
				aod = lt.merge(( aod, ensemble ))
		except:
			dbg( 'Error reading eNAAPS-NAV MET', l=2 )
	if 'dart' in members:	#_READ NAAPS DART ENSEMBLE	
		try:
			ensemble = read_ensfcst( dtg, model='ENAAPS-DART' )
			if type( ensemble ) == np.recarray:
				aod = lt.merge(( aod, ensemble ))
		except:
			dbg( 'Error reading eNAAPS-DART', l=2 )

	#_If subsetting by fhr
	aod = ln.subset( aod, fhr=fhr ) 
	try:
		icap = join_icap( aod, dtg, members=members, **kwargs )
		aod = lt.merge(( aod, icap ))
	except:
		dbg(( 'error merging ICAP_FCST', dtg ))

	#_return only ICAP
	if remove_members: aod = ln.subset( aod, model='ICAP' )

	return aod

def read_icapanal( *args, **kwargs ):
	kwargs.update({'model':'ICAP'}) 
	return ln.read_ensanal( *args, **kwargs )

def filter( records, strict_icap=True, members=['NAAPS','GEOS5','MACC',
	'MASINGAR','NGAC'], modes=False, **kwargs ):
	'''
	Builds recarray model_object() containing only ICAP records that
		1. Contain all models in members
		2. Contain all species of that member as defined in libmeta

	records	: model_object(),	np.recarray() of aod model data
	members	: list,			list of names to require to return

	'''
	dbg( records.size )
	tmp_model = lt.unique( records.model )
	if tmp_model != ['ICAP']:
		dbg(( 'filter for icap only', tmp_model ), l=3 )
		return records
	if not strict_icap:
		dbg(( 'icap not set to strict, returning' ), l=3 )
		return records
	
	#_Make expected species list for each model
	specs = ln.expected_vars()

	#_Initialize return object
	out = cl.model_object()

	#_REDUCE________________________________________________________________
	#_CLEANUP_______________________________________________________________
	#_Remove records lacking any of the required members
	#_Loop over each ICAP record
	for rec in records:
		#_take slice to check for masked members
		mask = rec.values[:,0,0].mask
		v = rec.variable

		#_loop over each model for this record, see if variable
		# is both expected and present
		desired = [] 
		for model in members:
			#_make list of expected species for each model
			idx = rec.values.member.tolist().index( model )
###			dbg(( rec.values.member, model, idx ))
			#_make list of indices to keep
			desired.append( idx )

			#_see if model is masked, and if so, break loop
			# leaving record out
			test = mask if type(mask) == np.bool_ else mask[idx]
			if test and v in specs[model]: 
				dbg(( rec.dtg_vald, v, 'filtered' ), l=1)
				dbg(( model, 'was the cause' ), l=1 )
				break

		#_if it makes it passed all members, add to return array
		else:
			#_Need to reducse attrv to plug var back into recarry
			atn, atv = ln.get_attr( rec )
			mem_idx = atn.index( 'member' )
			atv[mem_idx] = atv[mem_idx][desired]

			#_Update dimsize
			vals = rec.values.copy()
			vals_out = vals[desired,:,:]
			rec.dimsize = vals_out.shape

			#_Put back into original record
			rec.values = cl.var( vals_out, attrn=atn, attrv=atv )
			out = lt.merge(( out, rec ))	
	
	dbg( out.size )
	return out	

def join_icap( aod, fhr=120, fstrt=0, nt=None, finc=6, 
	members=lm.current_icap(), **kwargs ):
	'''
	Put all icap forecasts on common NAAPS grid for usage with ensemble 
	statistics
	
	require_all limits the returned values to timesteps when every model 
	present (at all) 
		can provide data.  So if MACC is in the mix, no 00z.
	'''
	if 'ICAP' in members: members.remove('ICAP')
	dbg( aod.size )

	#_Calculate last dtg
	species	= [s.lower() for s in mod_dict['ICAP']['specs']]
	nx = mod_dict['ICAP']['nx']
	ny = mod_dict['ICAP']['ny']
	lons = np.linspace( -179.5, 179.5, nx )
	lats = np.linspace( -89.5, 89.5, ny )
	finc = 6

	icap = cl.model_object()
	vars = lm.aod_vars()

	#_Create list of models with ANY data and icap models
	#_Loop over species, join models we have for specs
	dtg_valds = set( aod.dtg_vald )	#_list of unique dtgs
	dtg_init = lt.unique( aod.dtg_init )[0]

	#_Create array of for missing data
##	nan_2darray = np.empty((ny, nx))
##	nan_2darray[:] = NaN #_There's gotta ba shorthand for this
	nan_2darray = np.zeros(( ny, nx )) - 9999.
	
	nens_max = len( members )
	for spec in species:
		dbg( spec, l=2 )
		long_name = vars[spec]['long_name']
		for dtg_vald in dtg_valds:
			#_make recarray for one dtg, spec, but multiple models
			aod_sub = ln.subset( aod, variable=spec, model=members, 
				dtg_vald=dtg_vald )

			#_regrid models
			aod_rgd = np.empty(( nens_max, ny, nx ))
			aod_rgd[:] = -9999. 
			for e in np.arange( nens_max ):
				#_get model name and append it to dimension
				name = members[e]

				#_pull gridded data for specific model	
				tmp = ln.subset( aod_sub, model=name )
				if tmp.size == 1:	#_Should have single rec
					d = tmp.values[0]
					x = d.lon
					y = d.lat

					#_Regrid model data to icap x/y
					aod_rgd[e,:,:] = lt.regrid_field( x, y,
						 d, lons, lats ).transpose()
				elif tmp.size == 0:	#_Model data missing
					aod_rgd[e,:,:] = nan_2darray.copy() 
				else:
					print 'How did this happen?'
					return -1	

			#_Get indices that are non-physical
			neg_idx = np.where( aod_rgd < -1e-5 )
			aod_rgd[neg_idx] = -9999. #_SLOW
			aod_rgd = np.ma.masked_where( aod_rgd==-9999., aod_rgd )

			#_Convert to masked array and count present models
			nens = ln.check_members( aod_rgd )
###			miss = ( aod_rgd[:,0,0] == -9999. ).tolist().count(True)
###			nens = nens_max - miss

			data = cl.var( aod_rgd, attrv=(members, lats, lons,) )
			dimsize = data.shape
			dimname = ('member', 'lat', 'lon',)
			vhr = lt.find_runlength( dtg_init, dtg_vald ) / 3600

			icap.resize( icap.size + 1 )
			icap[-1] = ( data, 'ICAP', dtg_init, dtg_vald, vhr, 
				'global', spec, nens, dimname, dimsize, '',
				long_name ) 
	
	#_Limit to forecasts every finc hours
	idx_fhr = np.where( icap.fhr % finc == 0 )[0]
	icap = icap[idx_fhr]

	return icap

def read_ngacfcst( dtg, path=dir_prod + '/NGAC', require=['dust_aod'] ):
	'''
	dtg:	First valid time retrieved, dtg_vald if fstrt!=0
	fstrt:	How far into a forecast to start ( dtg - fstrt = dtg_init )
	nt:	Number of timesteps retrieved
	require:Species to not produce error 
	'''
	dbg( dtg )
	from netCDF4 import Dataset,Variable
	model = 'NGAC'
	species = [ s.lower() for s in mod_dict[model]['specs'] ]
	dtg_fcst = dtg 
	prefix = '/'.join(( path, dtg_fcst[:6], dtg_fcst ))
	aod = cl.model_object()
	vars = lm.aod_vars()

	finc = 3
	file_spec = prefix + '_aod_550_ngac.nc';
	for spec in species:
		long_name = vars[spec]['long_name']
		if not os.path.isfile(file_spec):
			error = file_spec + ' is missing'
			dbg( error, l=2 )
			if spec in require: raise IOError, error
		else:
			dbg( file_spec, l=2 )
			handle	= Dataset( file_spec, mode='r',
				format='NETCDF3_CLASSIC')
			null, ny, nx 	= handle.variables[spec].shape
			lons		= handle.variables['lon'][:]
			lats		= handle.variables['lat'][:]
			times_all	= handle.variables['time'][:]
			times 		= handle.variables['time'][:]
			dt		= (times[1]-times[0])/3600
			dtg_end		= lt.epoch2dtg(times[-1])
			nt		= len(times) 

			dtg_init = lt.epoch2dtg( times_all[0] )

			#_Make loop dict for time index
			dtg2t_ind = {}
			for t in np.arange(len(times_all)):
				dtg2t_ind[lt.epoch2dtg(times_all[t])] = t 

			for time in times: #np.arange(nt):
				t = times.tolist().index(time) 
				vhr = t*finc
				dtg_loop = lt.epoch2dtg( time )

				tmp = handle.variables[spec][t,:,:]

				#_Certain time periods need to be scaled
				if dtg_loop >= '2011071300' \
				and dtg_loop <= '2012012000':
					tmp = tmp / 10. 

				attrv = (lats, lons,) 
				data = cl.var( tmp, attrv=attrv )
				dimname = vars[spec]['dims'] 
				dimsize = data.shape
				aod.resize( aod.size + 1 )
				aod[-1] = ( data, model, dtg_init, dtg_loop, 
					vhr, 'global', spec, 0, dimname,
					dimsize, '', long_name )
			handle.close()

	#_Check for required species
	if aod.size == 0: raise IOError, 'missing fields'

	return aod

def read_maccfcst( dtg, path=dir_prod + '/MACC', 
	require=['biomassburning_aod','dust_aod','seasalt_aod','sulfate_aod',
	'total_aod'] ):
	"""
	  MACC files are read together - pass the dtg and a path 
		to build everything before _aod in the filename
	  aod[dtg][spec]	= read_maccfcst(dtg)
		dtg		= string, 		initial date
		fstrt		= integer,		start # hours after dtg
		finc		= integer,		timestep
	"""
	from netCDF4 import Dataset,Variable
	dbg( dtg )
	model = 'MACC'
	aod = cl.model_object()
	vars = lm.aod_vars()
	if dtg < '2012102600':
		  species = {	'biomassburning_aod'	: '210.210',
				'blackcarbon_aod'	: '211.210',
				'dust_aod'		: '209.210',
				'organicmatter_aod'	: '210.210',
				'seasalt_aod'		: '208.210',
				'sulfate_aod'		: '212.210',
				'total_aod'		: '207.210' 	}
	else:
		  species = {	'biomassburning_aod'	: 'bbaod550',
				'blackcarbon_aod'	: 'bcaod550',
				'dust_aod'		: 'duaod550',
				'organicmatter_aod'	: 'omaod550',
				'seasalt_aod'		: 'ssaod550',
				'sulfate_aod'		: 'suaod550',
				'total_aod'		: 'aod550' 	}
	

	dtg_init = dtg
	prefix = path + '/' + dtg_init[:6] + '/' + dtg_init

	#_Create list of dtgs to return
	

	for spec in species:
	   file_spec = prefix + '_' + spec + '_550_macc.nc';
	   if os.path.isfile(file_spec) == False:
		error = file_spec + ' is missing' 
		dbg( error, l=2 )
		if spec in require: raise IOError, error
	   else:
		dbg(( 'Reading', file_spec ), l=2 )
		code = species[spec]
		key = 'smoke_aod' if spec == 'biomassburning_aod' else spec
		long_name = vars[key]['long_name'] 
		handle = Dataset(file_spec,mode='r',format='NETCDF3_CLASSIC')
		null,ny,nx = handle.variables[code].shape
		lons = np.append(handle.variables['longitude'][nx/2:] \
			- 360., handle.variables['longitude'][:nx/2] )
		lats = handle.variables['latitude'][::-1]
		times = handle.variables['time'][:]

		for time in times:
			t = times.tolist().index( time )
			dtg_loop = lt.ecmwf_day2dtg( time )	

			#_Fix the orientation of dateline split
			tmp = handle.variables[code][t,:,:]
			values = np.append( tmp[::-1, nx/2:], tmp[::-1,:nx/2], 
				axis=1 )
			
			aod.resize( aod.size + 1 )
			data = cl.var( values, attrv=[ lats, lons ] )
			dimname = vars[key]['dims']
			dimsize = data.shape
			vhr = lt.find_runlength( dtg_init, dtg_loop ) / 3600
			aod[-1] = ( data, model, dtg, dtg_loop, vhr, 'global', 
				key, 0, dimname, dimsize, '', long_name ) 
		handle.close()

	if aod.size == 0: raise IOError, 'All macc fields missing' 

	return aod

def read_geos5fcst( dtg, require=['dust_aod','smoke_aod','seasalt_aod',
	'sulfate_aod','total_aod'], path=dir_prod + '/GEOS5'):
	  '''
	  Read GEOS-5 netcdf file
	  aod[dtg][spec] = read_geos5fcst(dtg)
		dtg	= string, 	start dtg
		finc	= integer,	step up to fhr
		nt	= integer,	number of timesteps to fhr, 
					normally diagnosed
		fhr	= integer,	maximum forecast hour from dtg
		path	= string,	directory of GEOS-5 files
	  '''
	  #_These are the species FROM GMAO.  
	  # Not the list getting PLOTTED. That
	  # is in the model dictionary
	  from netCDF4 import Dataset,Variable
	  dbg( dtg )
	  aod = cl.model_object()
	  vars = lm.aod_vars()
	  model = 'GEOS5'
	  finc = 6
	  fhr = 120
	  species = ['dust_aod','blackcarbon_aod',
			'organiccarbon_aod','seasalt_aod','sulfate_aod']
	  file		= path + '/' + dtg[:6] + '/' + dtg + '_aod_geos5.nc'
	  dbg( 'Reading ' + file, l=2 )

	  if os.path.isfile(file) == False:
	    raise IOError, file + ' is missing.' 
	  else:
	    handle	= Dataset( file, mode='r', format='NETCDF3_CLASSIC' )
	    lons	= handle.variables['lons'][:]
	    lats	= handle.variables['lats'][:]
	    times	= handle.variables['times'][:]
	    dtg_init	= lt.gsfc_day2dtg( times[0] ) #_Brings dtg to 00z, not 22
	    dtg_fcst	= lt.newdtg( dtg_init, fhr )
	    times	= times[::finc]
	    nt		= len(times)

	    for spec in species: # handle.variables: 
	      spec_geos = spec.replace( '_aod','' )
	      tmp = handle.variables[spec_geos]
	      long_name = vars[spec]['long_name']
	      for t in np.arange( nt ):
		days = times[t]
		dtg_loop = lt.gsfc_day2dtg( days )
		if dtg_loop > dtg_fcst: 
			break

		vhr = lt.find_runlength( dtg_init, dtg_loop ) / 3600
		data = cl.var( tmp[t*finc,:,:], attrv=[lats,lons])

		dimname = vars[spec]['dims'] 
		dimsize = data.shape				

		aod.resize( aod.size + 1 )
		aod[-1] = ( data, model, dtg_init, dtg_loop, vhr, 'global', 
			spec, 0, dimname, dimsize, '', long_name )

		#_Replace this with masked_outside
		if len(np.where( tmp > 1e3 ))  > 0:
			raise ValueError, 'GEOS-5 data out of range for ' + \
				dtg_loop + ' ' + spec
		
	    handle.close()
	
	    #_If we have bc and oc, sum them to generate smoke specie
	    spec_smoke	= ['blackcarbon_aod','organiccarbon_aod']
	    dtg_pres = aod.dtg_vald[0]	
	    for spec in spec_smoke: 
		test = ln.subset( aod, variable=spec )
		if test.size > 0:
			pass
		else: 	
			break
	    else: #-If all species present, sum and generate 'total'i WALTER
		long_name = vars['smoke_aod']['long_name']
		bc = ln.subset( aod, variable='blackcarbon_aod' )
		oc = ln.subset( aod, variable='organiccarbon_aod' )

		#_Get arrays of dtgs available for species
		dtg_bc = bc.dtg_vald
		dtg_oc = oc.dtg_vald

		#_Find where we have both values
		dtg_smoke = lt.intersection( [dtg_bc, dtg_oc] )
	
		for dtg_loop in dtg_smoke:
			if dtg_loop > dtg_fcst: break	#_KLUUUUUUUUDGE
		
			vhr = lt.find_runlength( dtg_init, dtg_loop ) / 3600
			bc_loop = ln.subset( bc, dtg_vald=dtg_loop ).values[0]
			oc_loop = ln.subset( oc, dtg_vald=dtg_loop ).values[0]

			attrv = [ lats, lons ]
			tmp = np.sum([ bc_loop, oc_loop,], axis=0 )
			data = cl.var( tmp, attrv=attrv )
			dimname = vars['smoke_aod']['dims'] 
			dimsize = data.shape				
			aod.resize( aod.size + 1 )
			aod[-1] = ( data, model, dtg_init, dtg_loop, vhr, 
				'global', 'smoke_aod', 0, dimname, dimsize,
				'', long_name )

		
			#_If we have all necessary species, sum total
			spec_total = ['dust_aod','seasalt_aod','sulfate_aod',
				'smoke_aod']
			for spec in spec_total:
				test = ln.subset( aod, variable=spec, 
					dtg_vald=dtg_loop )
				if test.size > 0:
					pass
				else: 
					dbg( 'Cannot calc total AOD ' 
						+ dtg_loop+' '+spec, l=2 )
					break
			else:	#-If bc and oc present, generate 'smoke' specie
				long_name = vars['total_aod']['long_name']
				aod_loop = ln.subset( aod, dtg_vald=dtg_loop )
				dust_loop = ln.subset( aod_loop,  
					variable=['dust_aod'] ).values[0]
				salt_loop = ln.subset( aod_loop,  
					variable=['seasalt_aod']).values[0]
				sulf_loop = ln.subset( aod_loop,  
					variable=['sulfate_aod']).values[0]
				smoke_loop = ln.subset( aod_loop, 
					variable=['smoke_aod']).values[0]
				tmp = np.sum([dust_loop, salt_loop, sulf_loop, 
					smoke_loop ], axis=0)
				attrv = [ lats,lons ]
				data = cl.var( tmp, attrv=attrv )
			
				dimname = vars['total_aod']['dims'] 
				dimsize = data.shape	
			
				aod.resize( aod.size+1 )	          
				aod[-1] = ( data, model, dtg_init, dtg_loop, 
					vhr, 'global', 'total_aod', 0,
					dimname, dimsize, '', long_name )

	  #_Check for required species
	  for s in require:
		size = ln.subset( aod, variable=s ).size
		if size == 0:
			raise RuntimeError, 'Not all GEOS-5 species available'
	  return aod

def read_masingarfcst( dtg,  
	require=['sulfate_aod','dust_aod','smoke_aod','seasalt_aod',
	'total_aod'], path=dir_prod+'/MASINGAR' ):
	'''
	  Read MASINGAR netcdf file
	  aod[dtg][spec] = read_geos5fcst(dtg)
		dtg	= string, 	start dtg
		finc	= integer,	step up to fhr
		nt	= integer,	number of timesteps to fhr, 
					normally diagnosed
		fhr	= integer,	maximum forecast hour from dtg
		path	= string,	directory of MASINGAR files
	'''
	from netCDF4 import Dataset,Variable
	dtg_fcst	= dtg
	file		= path + '/' + dtg_fcst[:6] + '/' + dtg_fcst \
			+ '_aod_masingar.nc'
	dbg( dtg )

	model = 'MASINGAR'
	aod = cl.model_object()
	vars = lm.aod_vars()

	fhr = 120
	specs = [ s.lower() for s in mod_dict[model]['specs'] ]
        if os.path.isfile(file) == False:
		raise IOError, 'WARNING: ' + file + ' is missing.' 
	else:
	    dbg( file, l=2 )
	    handle	= Dataset( file, mode='r', format='NETCDF3_CLASSIC' )
  	    null,ny,nx	= handle.variables['total'].shape
	    lons	= handle.variables['lon'][:]
	    lats	= handle.variables['lat'][:]
	    times	= handle.variables['time'][:]
	    dtg_init	= lt.epoch2dtg(times[0]) #_Brings dtg to 00z instead of 22
	    dtg_fhr	= lt.newdtg( dtg_fcst, fhr )
	    nt		= len(times) 

	    #_Loop over each variable, store in dictionary
	    for spec in specs: #handle.variables:
	      if spec == 'lat' or spec == 'lon' or spec == 'time': 
			continue
	      long_name = vars[spec]['long_name']
	      spec_mas = spec.replace( '_aod', '' )
	      tmp = handle.variables[spec_mas]
	      for t in np.arange(nt):
	        s = times[t]
	        dtg_loop = lt.epoch2dtg(s)
		if dtg_loop > dtg_fhr: break

		vhr = lt.find_runlength( dtg_init, dtg_loop ) / 3600

		attrv = (lats, lons,) 
  	        data = cl.var( handle.variables[spec_mas][t,:,:], attrv=attrv )

		dimname = vars[spec]['dims']
		dimsize = data.shape

		aod.resize( aod.size + 1 )
		aod[-1] = ( data, model, dtg_init, dtg_loop, vhr, 'global', 
			spec, 0, dimname, dimsize, '', long_name )
		
	    handle.close()

	#_Check for required species
	for s in require:
		size = ln.subset( aod, variable=s ).size
		if size == 0:
			err = 'Not all '+model+' species available'
			raise RuntimeError, err 

	return aod

def read_nmmbfcst( dtg, path=dir_prod+'/NMMB', 
	require=['dust_aod','seasalt_aod'] ):
	'''
	  Read NMMB netcdf file
	  aod[dtg][spec] = read_geos5fcst(dtg)
		dtg	= string, 	start dtg
		finc	= integer,	step up to fhr
		nt	= integer,	number of timesteps to fhr, 
					normally diagnosed
		fhr	= integer,	maximum forecast hour from dtg
		path	= string,	directory of MASINGAR files
	'''
	from netCDF4 import Dataset,Variable
	dtg_fcst = dtg
	file = '/'.join(( path, dtg_fcst[:6], dtg_fcst+'-NMMB_BSC_CTM-ICAP.nc'))
	dbg( dtg )
	model = 'NMMB'
	aod = cl.model_object()
	vars = lm.aod_vars()
 
	species = { 	'dust_aod550' : 'dust_aod', 
			'salt_aod550' : 'seasalt_aod' }

	#_BSC added seasalt on this date
	if dtg < '2012112000': del species['salt_aod550']

        if not os.path.isfile(file):
	  	error = file + ' is missing'
	  	dbg( error, l=2 ) 
	else:
		dbg( file, l=2 )
	    	handle	= Dataset(file,mode='r',format='NETCDF3_CLASSIC')
  	    	x,ny,nx	= handle.variables['dust_aod550'].shape
	    	lons	= handle.variables['lon'][nx/2,1:]
	    	lats	= handle.variables['lat'][:,ny/2]
	    	times	= handle.variables['time'][:]
	    	dtg_init= handle.dtg_init #_Brings dtg to 00z instead of 22
	    	nt	= len(times) 

		#_Loop over species dictionary to read in ncdf variables
		for spec in species:
			spec_loc = species[spec]
	      		long_name = vars[spec_loc]['long_name']
	      		tmp = handle.variables[spec]
	      		for t in np.arange(nt):
	        		s = times[t]
	        		dtg_loop = lt.newdtg( dtg_init, s )
				vhr = lt.find_runlength( dtg_init, dtg_loop ) \
							/ 3600 
			
				#_Setup attributes and add to recarray
				attrv = (lats, lons,)
  	        		data = cl.var( tmp[t,:,1:], attrv=attrv )
				dimname = vars[spec_loc]['dims']
				dimsize = data.shape
				aod.resize( aod.size + 1 )
				aod[-1] = ( data, model, dtg_init, dtg_loop, 
					vhr, 'global', spec_loc, 0, 
					dimname, dimsize, '', long_name )
    		handle.close()

	#_Check for required species
	for s in require:
		size = ln.subset( aod, variable=s ).size
		if size == 0:
			err = 'Not all '+model+' species available'
			raise RuntimeError, err 
	return aod

def cronlog( dtg, path='/'.join(( dir_prod, 'status' )), job={} ):
	'''
	Looks into $PRODUCT/status directory and uses True/False files
	to figure out the current progress of runs, downloads, and plotting

	dtg	: str*10, 	DTG for log file
	path	: str,		Path to output directory
	job	: dict,		keys == jobnames, values are status
	'''
	from time import strftime
	from glob import glob 
	import re
	dbg(( 'passed', dtg, job ))
	#_get current time
	current_time = strftime( '%Y-%m-%d_%H:%M:%S' )

	#_set log filename
	dir_out = '/'.join(( path, dtg[:6] ))
	lt.mkdir_p( dir_out )
	file = '/'.join(( dir_out, dtg+'_status' ))

	#_setup format for log writes
	fmt = '%-20s %-10s %-22s %-22s\n'

	#_Initialize log file if it does not exist
	if not os.path.exists( file ):
		with filelock( file, mode='w' ) as f:
			#_write header
			header = ( 'job', 'status', 'start-time', 'end-time' )
			f.write( fmt % header )

			#_write model statuses
			for member in lm.current_icap():
				line = ( member, 'UNKNOWN','UNKNOWN','UNKNOWN' )
				dbg( line, l=9 )
				f.write( fmt % line )
	
			#_add e-naaps to file
			line = ( 'ENAAPS-NAV', 'UNKNOWN','UNKNOWN','UNKNOWN' )
			dbg( line, l=9 )
			f.write( fmt % line )

			#_add plotting
			line = ( 'PLOTTING', 'UNKNOWN','UNKNOWN','UNKNOWN' )
			dbg( line, l=9 )
			f.write( fmt % line )

####	#_look for status files	
####	status_regex = re.compile( '(\w+)_\d{10}_status_(\w+)' )
####	glob_search = '/'.join(( dir_prod, 'status', dtg[:6], '*'+dtg+'*' ))
####	for statusfile in glob( glob_search ):
####		s = status_regex.search( statusfile )
####
####		#_if the file matches, get status
####		if s:
####			result = s.group(1) 
####			jobname = s.group(2)
####			if result == 'True':
####				pass
####	##			#_LEAVE THIS COMMENTED OUT
####	##			#_leave False files to trigger backfill() 
####	##			os.unlink( statusfile )
####				status = 'COMPLETE' 
####			else: 
####				status = 'FAILED' 
####
####			#_update job dictionary
####			if jobname not in job:
####				job.update({ jobname : status })
####				
	dbg(( 'writing', dtg, job ))
	#_loop over job dictionary and update main log file
	for jobname, status in job.iteritems():
		#_read in old log
		with filelock( file, 'r' ) as f: lines = f.readlines()

		#_write updated log
		with filelock( file, 'w' ) as f:
			for line in lines:
				values = line.split()
				
				#_if reporting on this job, update entry
				if values[0] == jobname:
					#_check to see what is being	
					# reported, then set the appropriate
					# value.  Times could be gathered
					# from the filestats of the status 
					# files, but that's for another day
					if status == 'STARTED':
						start = current_time 
					else:
						start = values[2]
			
					if status == 'COMPLETE': 
						end = current_time 
					elif status == 'FAILED':
						end = current_time
					else:
						end = 'UNKNOWN'

					line = ( jobname, status, start, end )
					f.write( fmt % line )

				#_otherwise, don't change
				else:
					f.write( line ) 
					
class filelock():
	'''
	Class to be used in 'with' structure. 
	Open a formatted file using mode, and creates a lockfile
	delaying any other operations attempting to read the same file.

	file	: str,	full path to file
	mode	: char,	w, r, a

	with filelock( '/tmp/test.txt', mode='w' ) as f:
		f.write( 'testing\n' )

	When 'with' block exits, file closed
	'''
	def __init__( self, file, mode='r' ):
		dbg( file, l=9 )

		#_gerate lockfile to prevent file from being corrupted
		lockfile( file, option='create' )

		#_set values as attributes to self, open file
		self.f = open( file, mode )

		self.filename = file
		self.mode = mode

	def __enter__( self ):
		#_set object that 'with' passes 'as'
		return self.f

	def __exit__( self, type, value, traceback ):
		#_executed regardless of outcome of 'with' block
		dbg(( type, value, traceback ), l=5 )

		#_clear off lockfile and close
		lockfile( self.filename, option='remove' )
		self.f.close()

		#_make readable by everyone
		if self.mode == 'w': lt.make_readable( self.filename )

def lockfile( file, option='check', **kwargs ):
        '''
        The lockfiles purpose is to prevent log from being opened by two jobs 

        option  : string, 'check','create','remove'
                        check   : returns a boolean indicating current status
                        create  : generate lockfile for model/dtg
                        remove  : delete lockfile for model/dtg
        '''
	import time
        lockfile = '.'.join(( file, 'lock' ))
	dbg( file, l=7 )
	dbg( lockfile, l=7 )

        if option == 'create':
                while os.path.exists( lockfile ): 
			dbg(( 'waiting; lockfile already exists', lockfile ))
			time.sleep( 3 )
                res = open( lockfile, 'w' ).close()
                return res
        elif option == 'remove':
                if not os.path.exists( lockfile ): return 0 
###			raise OSError, 'lockfile doesn\'t exists ' + lockfile 
                return os.unlink( lockfile )
        elif option == 'check':
                return os.path.exists( lockfile )
 
def dbg( msg, l=1 ):
	''' if global debug is set to true, be more verbose '''	
	import inspect
	msg = lt.to_string( msg )
        if hasattr( msg, '__iter__' ): msg = ' '.join( msg )

	if debug >= l:
		curf = inspect.currentframe()
		calf = inspect.getouterframes( curf, 2 )
		file, line, method = calf[1][1:4]
		file = file.split('/')[-1]
		print '[%s.%s.%i] %s' % ( file, method, line, msg )
