#!/usr/bin/env python
############################################################################_80
# First whack at trying to play with S-HIS data	     			      #
############################################################################_80
#
import libtools as lt

def run_main():
    import matplotlib.pyplot as plt
    
    file='/Users/wsessions/Work/test/shis_130820_filtered_sci.nc'
    data=read_radiances(file)
    
    print lt.epoch2dtg(data.epoch)
        
class structure(object):
    import numpy as np
    
    #_initialize variable structure
    def __init__(self):
        self.values=None
        self.epoch=None
        self.wavenumber=None
    
def read_radiances(file, var_name='radiance'):
    '''Read in radiances from netcdf files and store in recarray'''
    from netCDF4 import Dataset
    import numpy as np
        
    #_initialize netcdf input
    ncdf=Dataset(file, mode='r')
    
    #_get number of observations
    nobs=len(ncdf.dimensions['time'])
    
    #_initialize structure
    data=structure()
    
    #_get dimensions
    epoch_init=ncdf.variables['base_time'][0]
    epoch_offset=ncdf.variables['time_offset'][:]
    #epoch=ncdf.variables['time_offset'][:]+epoch_init
    epoch=np.array([epoch_init]*epoch_offset.size)+epoch_offset
    #epoch=epoch_offset+99999999
    #print epoch[0], epoch_offset[0]
    #print epoch_init+epoch_offset
    
    #_add to recarray
    data.epoch=epoch
    data.wavenumber=ncdf.variables['wavenumber'][:]
    data.values=ncdf.variables[var_name][:]
    
    #_close dataset
    ncdf.close()
    return data

class interferogram(object):
	'''
	igmH0	sweep direction 0 HBB reference interferogram (hot?)	
	igmA0	sweep direction 0 ABB reference interferogram (ambient?)
	igmO0	sweep direction 0 sky interferogram
	TH		H black body temperature
	TA		A blackbody temperature
	TH_bg	HBB reflected temperature (background)
	TA_bg	ABB reflected temperature (background)
	wnlaser	metrology laser wavenumber
	'''
	def __init__(self, mfile):
	    if mfile[-4:] == '.mat':
		from scipy.io import loadmat

		matlab = loadmat(mfile)
		self.__header__ 	= matlab['__header__']
		self.__version__ 	= matlab['__version__']
		
		self.TH 	= matlab['TH'][0][0]
		self.TA 	= matlab['TA'][0][0]
		self.TH_bg 	= matlab['TH_bg'][0][0]
		self.TA_bg 	= matlab['TA_bg'][0][0]
		self.igmH0 	= matlab['igmH0'][:,0]
		self.igmA0 	= matlab['igmA0'][:,0]
		self.igmO0 	= matlab['igmO0'][:,0]
		self.wnlaser = matlab['wnlaser'][0][0]

if __name__ == '__main__': run_main()
