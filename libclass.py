#############################################################################_80
#_NAME: 	libclass.py						       #
#_PURPOSE:	Classes for NVA						       #
#_USAGE:        Not a stand alone script.  Should only be imported, not run    #
#_AUTHOR:	Walter R. Sessions, September 2011			       #
#############################################################################_80
##########_GENERAL_#############################################################
################################################################################

import libtools as lt
import libmeta as lm
import numpy as np
from numpy import ma

NaN = np.nan

class var( ma.masked_array ):
	''' Setup to have basic expected information about 3d variables '''	
	def __new__( cls, 
		input, 
		attrv=None, 
		attrn=[ 'time','member','lat','lon' ] ):
		''' 
		attrv: attribute values, generally coordinates
		attrn: attribute names, lat/lon/etc

		attrn is is often not specified for AOD files, however
			it must be for concentration to replace 'member'
			with 'sigma'

		THIS IS MASKING VALUES BASED UPON PHYSICAL PROPERTIES OF AOD
		YET IS OCCASIONALLY USED FOR NON-AOD VARIABLES!!!
		'''
		aod_range = [ -9.e6, 9.e6 ] 

		#_Check number of dimensions, then limit attributes to be made
		ndim = input.ndim
		attrn = attrn[ -ndim: ]

		#_Initialize object
		input = ma.masked_where( np.isnan( input ), input ) 
					#_leave in for legcy files
		obj = ma.masked_array( input, fill_value=-9999. )
		idx = np.where( obj.data < aod_range[0] )
		obj.data[idx] = -9999. 
		idx = np.where( obj.data > aod_range[1] )
		obj.data[idx] = -9999.
		obj = ma.masked_where( obj == -9999., obj )

		if type( obj.mask ) == np.bool_:
			mask = np.empty( obj.shape )
			mask[:] = obj.mask
			obj.mask = mask

		#_Set dimensional attributes
		for i in np.arange( ndim ):
			attr = str( attrn[i] )
			value = np.array( attrv[i] )
			obj.__setattr__( attr, value )
			if input.shape[i] != value.size:
				err = 	'\nERROR: Dim size mismatch!\n' +\
				'ERROR: ' + attr + ' ' + str(i) +'\n'+ \
				'ERROR_VSIZE: ' + str(input.shape[i]) + '\n' \
				'ERROR_DSIZE: ' + str(value.size) + '\n' 
				raise AttributeError( err ) 

		return obj
	def __array_finalize__(self,obj):
		if obj is None: return

class model_object( np.recarray ):
	''' 
	Creates an empty record array for model data, 
	similar to an IDL structure

	INTIALIZING ARRAY OBJECT:
		model_records = model_object()

	Every record will have the follow fields:
		values		: class.var(), np.ndarray.  
					Will have attributes that list 
					dimensional coordinates.
					Prefer NCDF ordering
		model		: string, name.  Used as key for
					dictionaries and plotting
		dtg_vald	: string*10, valid time of data in values
		dtg_init	: string*10, initialization time of model
		fhr		: int, hours between dtg_vald, dtg_init
		region		: string*10, some regions are prefined,
					This should never be used as a key,
					but will be used to subset data. 
					Most start as 'global,' then are reset
					by sub_region()
		variable	: string*20, name of variable in values.
					Occasionally used as dictionary key,
					with headachy results
		ensemble	: bool, are the data in values part of an 
					ensemble.  Use to check to see if 
					'members' was in dimname, but it 
					comes up enough to warrant this
		dimname		: tuple, len*ndim, tuple listing names of
					dimensions of values.  Will correspond
					to attributes containing those coords
		dimsize		: tuple, values.shape
		units		: string*30, physical units of values
		long_name	: string*40, more descriptive name of variable 

	ADDING A RECORD:
		Add an empty record to object:
		model_records.resize( model_records.size + 1 )

		Fill it with data:
		model_records[-1] = ( values, model, dtg_vald, dtg_init, fhr,
			region, variable, ensemble, dimname, dimsize, units,
			long_name, [file] )

	This style recarray is not recommended for extremely large datasets where
	it becomes more prudent to use look up tables for repeated values such as
	dimensions in .values
	'''
	def __new__( self ):
		dtype = [ 	('values', var), 
				('model', 'a16'), 
				('dtg_init', 'a10'),
				('dtg_vald', 'a10'),
				('fhr', 'i4'),
				('region', 'a10'),
				('variable', 'a20'),
				('ensemble', 'i4'), 
				('dimname', tuple), 
				('dimsize', tuple),
				('units', 'a30'),
				('long_name', 'a40') ]
		obj = np.recarray( (0,), dtype=dtype )
		return obj

class sat_object(np.recarray):
	def __new__( self ):
		dtype = [ 	('tau', 'f4'), 
				('tausd', 'f4'), 
##				('satellite', 'a16'), 
##				('sensor', 'a16'),
				('lat', 'f4'),
				('lon', 'f4'),
				('sensor', 'a20'),
				('dtg', 'a10'),
				('wavelength', 'i4'),
				('level', 'i4'),
				('units', 'a30'),
				('long_name', 'a40') ]
		obj = np.recarray( (0,), dtype=dtype )
		return obj

class location_stats( object ):
	''' reads in ensemble stat file produced by write_ensemble_stats() '''
	def __init__( self, file, pipe=None, **kwargs ):
		import re
		#_Get metadata from filename
		fn = file.split('/')[-1]
		print fn
		match_string = '([\w\d-]+)_(\d{10})_(\d{10})_f[day]*(\d+)' \
			+ '_(\w+_aod)_(\w+)_\w{3}_(\w+)_[\w\d-]+.txt' 
		fn_reg = re.search( match_string, fn )

		label 	= fn_reg.group(1)
		dtgs	= fn_reg.group(2)
		dtge	= fn_reg.group(3)
		fhr	= fn_reg.group(4)
		variable= fn_reg.group(5)
		mode	= fn_reg.group(6)
		error	= fn_reg.group(7)

		#_Open ASCII file, save lines, close
		f = open( file, 'r' )		
		lines = f.readlines()
		f.close()

		dtype_ens = [ 	('location', 'a16'),
				('variable', 'a16'),
				('dtgs', 'a10'),
				('dtge', 'a10'),
				('fhr', 'i4'),
				('mode', 'a12'),
				('error_model', 'a12'),
				('label', 'a20'),
				('mean', 'f8'),
				('stdev', 'f8'),
				('max', 'f8'),
				('p90', 'f8'), 
				('p75', 'f8'), 
				('p50', 'f8'), 
				('p25', 'f8'), 
				('p10', 'f8'), 
				('min', 'f8'), 
				('nobs', 'i4'),
				('rmse', 'f8'), 
				('briers', 'f8'), 
				('bias', 'f8'), 
				('mae', 'f8'), 
				('nrmse', 'f8'), 
				('fge', 'f8'), 
				('nens', 'i4'), 
				('model', 'a10'),
				('members', np.ndarray ),
				('ranks', np.ndarray ) ]
		dtype_mem = [	('location', 'a16'),
				('variable', 'a16'),
				('dtgs', 'a10'),
				('dtge', 'a10'),
				('fhr', 'i4'),
				('mode', 'a12'),
				('error_model', 'a12'),
				('label', 'a20'),
				('rmse', 'f8'),
				('bias', 'f8'), 
				('mae', 'f8'), 
				('nrmse', 'f8'), 
				('fge', 'f8'), 
				('nobs', 'i4'),
				('model', 'a30')	]
		ens = np.recarray( (0,), dtype=dtype_ens )
		mem = np.recarray( (0,), dtype=dtype_mem )

		regex_stat = re.compile( 'aeronet-(\w+)' )
		loc 	= {} 	#_temp storage for reg[ mem_name ][ stat_name ]
		mem_idx = {} 	#_tmp lookup for member indices
		nobs	= {}	#_keep track of n for regions in ens section

		#_Loop over lines, put in record
		for l in lines:
			#_Split columns by whitespace into list
			c = l.split()

			#_Check if a header line
			if len( c ) == 0:		#_Spacer
				continue
			elif c[0] == 'REGION' or c[0] == 'AERONET':#_Ens Header
				loc_ind = 0 
				mu_ind 	= c.index( 'MEAN' )
				std_ind = c.index( 'STDEV' )
				max_ind = c.index( 'MAX' )
				p90_ind = c.index( '90TH' )
				p75_ind = c.index( '75TH' )
				p50_ind = c.index( '50TH' )
				p25_ind = c.index( '25TH' )
				p10_ind = c.index( '10TH' )
				min_ind = c.index( 'MIN' )
				nbs_ind = c.index( 'NOBS' )
				rms_ind = c.index( 'RMSE' )
				bri_ind = c.index( 'BRIERS' )
				bia_ind = c.index( 'BIAS' )
				mae_ind = c.index( 'MAE' )
				nrl_ind = c.index( 'NMRSE' )
				fge_ind = c.index( 'FGE' )
				ens_ind = c.index( 'NMEMBERS' )
				mod_ind	= c.index( 'MODEL' )
				mem_ind = c.index( 'MEMBERS' )

				#_Find rank histogram headings
				rex = re.compile( '[\w\d]+_RNK' )
				rnk_head = [ rex.search(b).group(0) for b in c\
					if rex.search(b) != None ]
				rnk_ind = [ c.index(b) for b in rnk_head ]

				section = 'ensemble'
				continue

			#_When it hits the headers for Member Sections
			elif regex_stat.search( c[0] ): #== '':	#_Member Header
				stat = regex_stat.search( c[0] ).group(1)
				loc_ind = 0 

				#_Find Model stat headings (gets reset a lot) 
				[mem_idx.update({m:c.index(m)}) \
					for m in members ]
				section = 'member'
				continue

			#_Else, add data to record
			#_ENSEMBLE
			if section == 'ensemble':
				members = tuple( c[mem_ind].split(',') )
				#ranks = tuple( np.array(c)[rnk_ind] )
				ranks = [ int(a) for a in np.array(c)[rnk_ind]] 
				vals = ( c[loc_ind], 
					variable, 
					dtgs, 
					dtge, 
					fhr,
					mode, 
					error, 
					label, 
					c[mu_ind], c[std_ind],
					c[max_ind], c[p90_ind], c[p75_ind],
					c[p50_ind], c[p25_ind], c[p10_ind],
					c[min_ind], c[nbs_ind], c[rms_ind], 
					c[bri_ind], c[bia_ind], c[mae_ind], 
					c[nrl_ind], c[fge_ind], c[ens_ind], 
					c[mod_ind], members, ranks )
				ens.resize( ens.size + 1 )
				ens[-1] = vals 
				nobs[c[loc_ind]] = c[nbs_ind]
		
			#_MEMBERS
			elif section == 'member':
				#_And this just gets rewritten ad nauseum
				r = c[loc_ind]
				if r not in loc: loc[ r ] = {}
				if r not in nobs: continue
				loc[r]['nobs'] = nobs[r]

				#_Get member statistics in tmp dict
				for (m, idx) in mem_idx.items():
					#_Initialize member
					if m not in loc[r]:
						loc[r][m] = {}

					#_Set value
					loc[r][m][stat] = c[idx]

		#_Add member statistic records
		for r in loc:
			for m in loc[ r ]:
				if m == 'nobs': continue
				vals = ( r, variable, dtgs, dtge, fhr, 
					mode, error, label, 
					loc[r][m]['rmse'], loc[r][m]['bias'], 
					loc[r][m]['mae'], loc[r][m]['nrmse'], 
					loc[r][m]['fgs'], loc[r]['nobs'], m )
				mem.resize( mem.size + 1 )
				mem[-1] = vals

		self.ensemble = ens
		self.member = mem
	
		if pipe != None:
			pipe.send(( ens, mem ))
			pipe.close()

class cum_mean( object ):
	'''
	Used with running_summation() to hold values to
	eventually calculate the cummulative mean
	'''
	def __new__( self ):
		dtype = [ 	('values', var),
				('model', 'a16'), 
				('dtg_vald', np.ndarray),  
				('dtg_init', np.ndarray),  
				('dtg_valds', np.ndarray),  
				('dtg_inits', np.ndarray),  
				('fhr', 'i4'),
				('region', 'a10'),
				('variable', 'a20'),
				('ensemble', 'i4'), 
				('dimname', tuple), 
				('dimsize', tuple),
				('units', 'a30'),
				('long_name', 'a40') ]
		obj = np.recarray( (0,), dtype=dtype )
		#_Create layer of n for each pointin running_summation()	
		return obj

class lidar( np.recarray ):
	'''
	Used to contain multiple lidar observations
	Initialized by passing an array of lidar() objects
	'''
	def __new__( self ):
		dtype = [	('tau','f8'), 
				('epoch','i8'),
				('dtg','a10'),
				('code', 'a6'),
				('nrl_name', 'a20'), 
				('long_name','a20'), 
				('wavelength','i4'), 
				('mode','a10'), 
				('lat','f8'), 
				('lon','f8')	]
		obj = np.recarray( (0,), dtype=dtype )	
		return obj	

class aeronet_stats( object ):
	def __init__( self, stats ):
		self.nobs 	= stats[0] 	#_Number of Observations
		self.bias 	= stats[1] 	#_Bias (NOT IMPLEMENTED)
		self.rmse	= stats[2] 	#_Root Mean Square Error
		self.rmsd	= stats[3] 	#_Root Mean Square Dev 
		self.nrmse 	= stats[4] 	#_Normalized RMSE
		self.mae 	= stats[5] 	#_Mean Absolute Error
		self.fgs 	= stats[6] 	#_Fractional Gross Error
		self.briers 	= stats[7] 	#_Half-Briers Score
		self.ranks 	= stats[8]	#_Rank array of nens+1 

class aeronet_compz( np.recarray ):
	def __new__( self ):
		dtype = [	('code', 'a6'),
				('lat', 'f4'),
				('lon', 'f4'),
				('epoch', 'f8'),
				('mode', 'a10'),
				('variable', 'a16'),
				('model', 'a12'),
				('fhr', 'i4'),
				('residual','f4'),
				('gross', 'f4'),
				('nrml_res', 'f4'),
				('briers_res', 'f4'),
				('rank', np.ndarray),
				('ensemble', bool),
				('members', list),
				('label', 'a30'),]
		obj = np.recarray( (0,), dtype=dtype )
		return obj

class aeronet_comps( np.recarray ):
	def __new__( self ):
		dtype = [	('code', 'a6'),
				('lat', 'f4'),
				('lon', 'f4'),
				('epoch', 'f8'),
				('mode', 'a10'),
				('variable', 'a16'),
				('model', 'a12'),
				('fhr', 'i4'),
				('residual','f4'),
				('gross', 'f4'),
				('nrml_res', 'f4'),
				('briers_res', 'f4'),
				('rank', np.ndarray),
				('error_model', 'a16'),
				('ensemble', bool),
				('members', list),
				('label', 'a30')	]
		obj = np.recarray( (0,), dtype=dtype )
		return obj

