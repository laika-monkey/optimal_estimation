
def read_cod(flight, dtg, fidx=None, out_label=None, thresh=0.3, 
	dir_lbl=os.path.join(DIR_PROD, 'LBL-RTM_hs3'), experiment='HS3', **kwargs):
	''' 
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
	dir_dtg = '{0}_{1}.{2}_fov{{0:06d}}'.format(experiment, d0, d1)

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

		#_read in layers
		lbl = Lbldis_input(lbl_input)

		cod = 0.
		for j, layer in enumerate(lbl.layers):
			#_pull out current layers ssp
			ssp = lbl.ssp_db[layer['dbnum']]
			
			#_if it's the water cloud, include its OD in total
			if reg.search(ssp):
				cod += sum(layer['tau'])

		#_add values for cod to record array for SHIS and CPL
		data[i] = (epoch, f, cod, _get_cpl_cod(flight, f))	

	return data
