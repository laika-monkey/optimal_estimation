#!/usr/bin/env python

'''
KEEP TRACK OF DIFFERENCES BASED ON CLOUD TYPE

'''
import os, sys


def run_main():
	#_run this to initially create a matched dataset to test
	write_data_live()

	#_read that in with this
##	data = read_data_matched()

	#_slice, dice
	
	return 0


################################################################################
################################################################################
################################################################################


class matched_data(object):
	def __new__(self, size):
		dtype = [	('wavenumber', 'f8'),	('cloud', 'f4'), 
					('rads_lbldis', 'f8'),	('rads_airs', 'f8'),
					('tb_lbldis', 'f8'),	('tb_airs', 'f8'),
					('latitude','f8'),		('longitude', 'f8'),
					('time', 'f8'),	]
		data = recarray((size,), dtype=dtype)
		return data


def read_data_matched(fname, **kwargs):
	''' read in file produced by write_matched_data() '''
	from netCDF4 import Dataset

	cdf = Dataset(fname, 'r')

	nw = len(cdf.dimension['wavenumber'])
	nt = len(cdf.dimension['time'])

	data = Matched_data(size=nw*nt)

	data.wavenumber[:] = cdf.variables['wavenumber'][:]
	data.cloud[:] = cdf.variables['cloud'][:]

	data.lbldis[:] = cdf.variables['lbldis'][:]
	data.airs[:] = cdf.variables['airs'][:]

	data.latitude[:] = cdf.variables['latitude'][:]
	data.longitude[:] = cdf.variables['longitude'][:]

	return data


def write_data_live(dir_gdas=os.path.join(os.environ['PRODUCTS'], 'gdas'),
	ngdas=99999, **kwargs):
	from numpy import ndarray, recarray
	from lblrtm_utils import read_lbldis_out, read_airs, convolve2airs
	from lblrtm_utils import airs_channels 
	from hs3_utils import Flight_data, FLIGHTS #_read_shis
	from glob import glob
	from libgeo import planck_inv

	'''
	reads in all gdas files in dir_gdas, then tries to find 
	LBL-DIS output directories that match the timestamp/naming

	returns: all dat noise in a recarray
	'''

	#_get list of directories to read in
	dirs_lbldis, files_gdas = get_lbldis_directories()
	nrun = len(dirs_lbldis)		

	
	wave = airs_channels()[:].copy()
	nair = len(wave)

	#_testing purposes
	files_gdas = files_gdas[:ngdas]

	#_initialize output file
	cdf = initialize_cdf(nrun, nair)

	vvv = cdf.variables
	vvv['wavenumber'][:]= wave 
 
	#_conversion args
	args = (wave,)
	kwargs.update({'domain' : 'wavenumber'})
	
	n = 0	
	for file_gdas in files_gdas:
		airs = read_airs(file_gdas)
			
		#_redundancy!  This limits the case to this gdas file
		dirs_lbldis = glob(file_gdas[:-4].replace('gdas', 'LBL-RTM') + '/0*')
		dirs_lbldis.sort()
		for dir_lbldis in dirs_lbldis: 	
			#_read lblrtm output
			fname	= os.path.join(dir_lbldis, 'lbldis_out.cdf')
			lbldis	= read_lbldis_out(fname)
			f0		= lbldis.wavenumber.min()
			f1		= lbldis.wavenumber.max()
			lblcon	= convolve2airs(lbldis, f0=f0, f1=f1)

			#_get time index of collocated scene and where to start/end	
			idx	= abs(airs.time - lblcon.time).argmin()
	
			#_get differences (with wnum, T, lat, lon, cloud, time attached)
			vvv['cloud'][n]		= lblcon.cloud
			vvv['latitude'][n]	= lblcon.latitude
			vvv['longitude'][n]	= lblcon.longitude
			vvv['time'][n]		= lblcon.time
	
			vvv['rads_lbldis'][n,:]	= lblcon.radiances
			vvv['rads_airs'][n,:]	= airs[idx].spectrum.radiances
			vvv['tb_lbldis'][n,:]	= planck_inv(lblcon.radiances, 
												*args, **kwargs)
			vvv['tb_airs'][n,:]		= planck_inv(airs[idx].spectrum.radiances,
												*args, **kwargs)

			print nrun, n 
			n += 1

	cdf.close()


def initialize_cdf(nt, nw):
	from netCDF4 import Dataset
	from libtools import unique, subset
	from libnva import sort_rec
	from numpy import argsort

	cdf = Dataset('matched_data_live_tb.nc', 'w', format='NETCDF3_CLASSIC')
	cdf.createDimension('wavenumber', nw)
	cdf.createDimension('time', nt)
		
	cdf.createVariable('wavenumber', 'f8', ('wavenumber',))
	cdf.createVariable('time', 'f8', ('time',))
	cdf.createVariable('cloud', 'f4', ('time',))
	cdf.createVariable('latitude', 'f8', ('time',))
	cdf.createVariable('longitude', 'f8', ('time',))
	cdf.createVariable(	'rads_lbldis', 'f8', ('time', 'wavenumber'),
						fill_value=-9999)
	cdf.createVariable(	'rads_airs', 'f8', ('time', 'wavenumber'),
						fill_value=-9999)
	cdf.createVariable(	'tb_lbldis', 'f8', ('time', 'wavenumber'), 
						fill_value=-9999)
	cdf.createVariable(	'tb_airs', 'f8', ('time', 'wavenumber'),
						fill_value=-9999)

	return cdf


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
