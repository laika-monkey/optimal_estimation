#!/usr/bin/env python

def run_main():
	from netCDF4 import Dataset
	from numpy import append, array
	from glob import glob
	from lblrtm_utils import read_caliop, read_gdas
	import os
	from libnva import draw_map
	import matplotlib.pyplot as plt
	
	path_in = os.path.join(os.environ['PRODUCTS'], 'gdas')
	files = glob(path_in + '/*.hdf')
	lats = array([])
	lons = array([])

	for fname in files:
		with Dataset(fname, 'r', format='HDF4') as hdf:
			caliop = read_caliop(hdf)
			gdas = read_gdas(hdf)

##		print '{0:8.2f} {1:8.2f}'.format(caliop.longitude.max(), caliop.longitude.min())
		#_locations
		idx = caliop.latitude >= -90 
		gdas = gdas[idx]
		caliop = caliop[idx]
		idx = caliop.latitude <= 90 
		gdas = gdas[idx]
		caliop = caliop[idx]
		if gdas.size == 0:
			continue

		lats = append(lats, caliop.latitude.copy())
		lons = append(lons, caliop.longitude.copy())

		#_cloud free
		idx = caliop.cloud == -101
		gdas = gdas[idx]
		caliop = caliop[idx]
		if gdas.size == 0:
			continue
	
		#_ocean only
		idx = gdas.land_fraction == 0
		gdas = gdas[idx]
		caliop = caliop[idx]
		if gdas.size == 0:
			continue

		#_lat/lon requirements
		lat0_idx = caliop.latitude > 10 
		lat1_idx = caliop.latitude < 20 
##		lon0_idx = caliop.longitude > -46 
##		lon1_idx = caliop.longitude < -18

		def sing(a): return a.reshape(1,a.size)		
		idx = append(sing(lat0_idx), sing(lat1_idx), axis=0)
##		idx = append(idx, sing(lon0_idx), axis=0)
##		idx = append(idx, sing(lon1_idx), axis=0)
		idx = idx.all(axis=0)
		gdas = gdas[idx]
		caliop = caliop[idx]

		if gdas.size > 0:
			print fname
			print gdas.time
			print gdas.latitude, gdas.longitude

	with open('coords', 'w') as f:
		for n in xrange(lats.size):
			f.write('{0:8.5f} {1:8.5f}\n'.format(lons[n], lats[n]))

	m = draw_map()
	m.scatter(lons, lats, latlon=True)
	plt.savefig('map.png')

if __name__ == '__main__': run_main()
