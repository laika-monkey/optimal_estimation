#!/usr/bin/env python

'''
KEEP TRACK OF DIFFERENCES BASED ON CLOUD TYPE

'''


import os, sys


def run_main(**kwargs):
	from libtools import subset, unique
	from lblrtm_utils import airs_channels
	from numpy import array
	wavenumbers = airs_channels()
	#_run this to initially create a matched dataset to test
	write_matched_data(read_data_raw())

	#_read that in with this
	data = read_data_matched()
	
	#_remove -8888 clouds	
	data = subset(data, cloud=-101.0)

	#_latitude ranges
	ranges = {	'lat_hgh' : [60, 90],
				'lat_mid' : [20, 60],
				'lat_low' : [0, 20]	}
	for label, lat_range in ranges.iteritems():
		idx0 = abs(data.latitude) >  lat_range[0]
		idx1 = abs(data.latitude) <= lat_range[1]

		#_subset the data by indices
		idx = array([idx0, idx1]).all(axis=0)
		dat = data[idx].copy()

		bias_histogram(dat, label=label, **kwargs)

	#_longitude ranges
	ranges = {	'lon_west' : [-180, 0],
				'lon_east' : [0,  180]	}
	for label, lon_range in ranges.iteritems():
		idx0 = abs(data.longitude) >  lon_range[0]
		idx1 = abs(data.longitude) <= lon_range[1]

		#_subset the data by indices
		idx = array([idx0, idx1]).all(axis=0)
		dat = data[idx].copy()

		bias_histogram(dat, label=label, **kwargs)

###	#_cloud mask
###	for label in unique(data.cloud):
###		#_subset the data by indices
###		dat = subset(data, cloud=label)
###		bias_histogram(dat, label=str(label), **kwargs)

	#_wavenumber ranges
	ranges = {	'co2'	: [600, 740],
				'win0'	: [740, 990],
				'o3'	: [990, 1080],
				'wv'	: [1080, 1600],
				'ch4'	: [1600, 9999] }	
	for label, wvn_range in ranges.iteritems():
		idx0 = wavenumbers >=  wvn_range[0]	#_inside == True
		idx1 = wavenumbers < wvn_range[1]
		idx = array([idx0, idx1]).all(axis=0)
##		idx0 = wavenumbers <=  wvn_range[0] #_outside == True
##		idx1 = wavenumbers > wvn_range[1]
##		idx = array([idx0, idx1]).any(axis=0)

		#_subset the data by indices
		bias_histogram(data, label=label, mask=idx, **kwargs)
	
	return 0


################################################################################
################################################################################
################################################################################
def plot_map(data, **kwargs):
	from libnva import draw_map
	pass


def bias_histogram(data, nbin=50, label='label', mask=None, **kwargs):
	'''
	produce histograms of bias size between LBL-DIS and AIRS

	Subset outside of this based upon latitude, cloud mask value,
	wavenumber, et cetera
	'''
	import matplotlib.pyplot as plt
	from numpy import tile, vstack, linspace
	from numpy.ma import masked_where

	if not data.size:
		print 'not enough data in ', label
		return

	#_histogram spacing
	nbin = linspace(-0.05, .05, 101)

	#_get differences between passed data
	stk_l = vstack(data.copy().lbldis)
	stk_a = vstack(data.copy().airs)

	#_mask arbitrary
	if mask is not None:
		mask = tile(mask, (stk_l.shape[0], 1))
		stk_l = stk_l[mask]
		stk_a = stk_a[mask]
##		stk_l = masked_where(mask, stk_l)  
##		stk_a = masked_where(mask, stk_a)  
	stk_l = stk_l.flatten()
	stk_a = stk_a.flatten()

	#_mask missing data
	idx = stk_a != -9999
	stk_a = stk_a[idx]
	stk_l = stk_l[idx]
	idx = stk_a > 1e-6 
	stk_a = stk_a[idx]
	stk_l = stk_l[idx]

	reslt = stk_l - stk_a
	plt.title(	label.upper().replace('_',' ') + 
				' LBL-DIS minus AIRS nchan=' + str(reslt.size))
	n, bins, patches = plt.hist(reslt, bins=nbin) #_log=True)
	plt.xlim(-0.1,0.1)
	plt.xlabel('W/m$^2$/cm$^{-1}$/sr')
	plt.ylim(0, 10000)

	pname = '.'.join(('bias', label, 'png'))
	print pname
	plt.savefig(pname)
	plt.close()


class Matched_data(object):
	def __new__(self, size):
		from numpy import recarray, ndarray
		dtype = [	('cloud', 'f4'), 
					('lbldis', ndarray),	('airs', ndarray),
					('lbldis_bt', ndarray),	('airs_bt', ndarray),
					('latitude','f8'),		('longitude', 'f8'),
					('time', 'f8'),	]
		data = recarray((size,), dtype=dtype)
		return data


def read_data_matched(fname='matched_data_live.nc', mask=True, **kwargs):
	''' read in file produced by write_matched_data() '''
	from netCDF4 import Dataset
	from numpy import zeros
	from libclass import var
	from libgeo import planck_inv
	from lblrtm_utils import AIRS_ancillary
	from numpy.ma import masked_where

	cdf = Dataset(fname, 'r')

	nw = len(cdf.dimensions['wavenumber'])
	nt = len(cdf.dimensions['time'])

	data = Matched_data(size=nt)

	wvn = cdf.variables['wavenumber'][:]
	cld = cdf.variables['cloud']
	dis = cdf.variables['lbldis']
	air = cdf.variables['airs']
	lat = cdf.variables['latitude']
	lon = cdf.variables['longitude']
	tim = cdf.variables['time']

	arg = {'attrn' : ['wavenumber'], 'attrv' : [wvn]}
	mask_l2 = AIRS_ancillary().l2_ignore

	for n in xrange(nt):
		lbldis	= dis[n,:]
		airs	= air[n,:]

		#_blank out bad channels
		if mask:
			lbldis	= masked_where(mask_l2, lbldis)	
			airs	= masked_where(mask_l2, airs)	

			lbldis	= masked_where(lbldis == -9999, lbldis)	
			airs	= masked_where(airs == -9999, airs)	
			
		lbldis	= var(lbldis, **arg)
		airs	= var(airs, **arg)
		lbldis_tb	= var(planck_inv(lbldis, wvn, domain='wavenumber'), **arg)
		airs_tb		= var(planck_inv(airs, wvn, domain='wavenumber'), **arg)

		data[n] = (	cld[n], lbldis, airs, lbldis_tb, airs_tb, 
					lat[n], lon[n], tim[n])

	return data


def read_data_raw(dir_gdas=os.path.join(os.environ['PRODUCTS'], 'gdas'),
	ngdas=99999, **kwargs):
	from numpy import ndarray, recarray
	from lblrtm_utils import read_lbldis_out, read_airs, convolve2airs
	from glob import glob
	'''
	reads in all gdas files in dir_gdas, then tries to find 
	LBL-DIS output directories that match the timestamp/naming

	returns: all dat noise in a recarray
	'''

	#_get list of directories to read in
	lbldis_dirs, files_gdas = get_lbldis_directories()
	nrun = len(lbldis_dirs)				#_number of scenes
	nair = 2378							#_number of airs channels

	#_testing purposes
	files_gdas = files_gdas[:ngdas]

	#_declare recarray
	dtype = [	('wavenumber', 'f8'),	('cloud', 'f4'), 
				('lbldis', 'f8'),		('airs', 'f8'),
				('latitude','f8'),		('longitude', 'f8'),
				('time', 'f8'),	]
	data = recarray((nrun*nair,), dtype=dtype)

	n = 0	
	for file_gdas in files_gdas:
		airs = read_airs(file_gdas)
##		f0 = airs_channels.max() 
##		f1 = airs_channels.min() 

		dirs_lbldis = glob(file_gdas[:-4].replace('gdas', 'LBL-RTM') + '/0*')
		for dir_lbldis in dirs_lbldis: 	
			#_read lblrtm output
			fname	= os.path.join(dir_lbldis, 'lbldis_out.cdf')
			lbldis	= read_lbldis_out(fname)
			f0 = lbldis.wavenumber.min()
			f1 = lbldis.wavenumber.max()
			lblcon	= convolve2airs(lbldis, f0=f0, f1=f1)

			#_get time index of collocated scene and where to start/end	
			idx		= abs(airs.time - lblcon.time).argmin()
			idxs	= n*nair
			idxe	= idxs + nair
	
			#_get differences (with wnum, T, lat, lon, cloud, time attached)
			data.wavenumber[idxs:idxe]	= lblcon.wavenumber
			data.cloud[idxs:idxe]		= lblcon.cloud
			data.latitude[idxs:idxe]	= lblcon.latitude
			data.longitude[idxs:idxe]	= lblcon.longitude
			data.time[idxs:idxe]		= lblcon.time
	
			data.lbldis[idxs:idxe]		= lblcon.radiances
			data.airs[idxs:idxe]		= airs[idx].spectrum.radiances

			print nrun*nair, idxs, idxe
			n += 1

	return data


def write_matched_data(data):
	from netCDF4 import Dataset
	from libtools import unique, subset
	from libnva import sort_rec
	from numpy import argsort
	raise RuntimeError, 'USE LIVE_WRITE'
	cdf = Dataset('matched_data.nc', 'w', format='NETCDF3_CLASSIC')

	#_define dimensions
	time = unique(data.time)
	wnum = unique(data.wavenumber)
	nt	= len(time)
	nw	= len(wnum)
	cdf.createDimension('wavenumber', nw)
	cdf.createDimension('time', nt)

	cdf.createVariable('wavenumber', 'f8', ('wavenumber',))
	cdf.createVariable('time', 'f8', ('time',))
	cdf.variables['wavenumber'][:] = wnum
	cdf.variables['time'][:]	= time 

	cdf.createVariable('cloud', 'f4', ('time',))
	cdf.createVariable('latitude', 'f8', ('time',))
	cdf.createVariable('longitude', 'f8', ('time',))

	cdf.createVariable(	'radiances_lbldis', 'f8',
						('time', 'wavenumber'), fill_value=-9999)
	cdf.createVariable(	'radiances_airs', 'f8',
						('time', 'wavenumber'), fill_value=-9999)
	cdf.createVariable(	'tb_lbldis', 'f8',
						('time', 'wavenumber'), fill_value=-9999)
	cdf.createVariable(	'tb_airs', 'f8',
						('time', 'wavenumber'), fill_value=-9999)

	for tidx in xrange(nt):
		t = time[tidx]

		#_pull out this time
		data_t	= subset(data, time=t) 

		#_pull out thiiiisss
		cdf.variables['cloud'][tidx]	= unique(data_t.cloud, unique=True)
		cdf.variables['latitude'][tidx]	= unique(data_t.latitude, unique=True)
		cdf.variables['longitude'][tidx]= unique(data_t.longitude, unique=True)

		#_sort these
		idx = argsort(data_t.wavenumber)

		cdf.variables['lbldis'][tidx]	= data_t.lbldis[idx]
		cdf.variables['airs'][tidx]		= data_t.airs[idx]

	cdf.close()


def get_lbldis_directories(
	dir_lbldis_out=os.path.join(os.environ['PRODUCTS'], 'LBL-RTM'), **kwargs):
	'''get list of all lbldis directories within dir_lbldis_out'''
	from glob import glob

	lblrtm_dirs = []
	gdas_dirs	= glob(dir_lbldis_out + '/COL*')
	[lblrtm_dirs.extend(glob(f + '/*')) for f in gdas_dirs]
	files_gdas = ['.'.join((g.replace('LBL-RTM', 'gdas'), 'hdf')) \
		for g in gdas_dirs]

	#_get number of input collocated files
	files_gdas.sort()
	lblrtm_dirs.sort()	
	return lblrtm_dirs, files_gdas	


if __name__ == '__main__': run_main()
