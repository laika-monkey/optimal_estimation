def read_modis_aod_xy(fname, x, y,
	variable='Effective_Optical_Depth_Average_Ocean',
#	variable='Optical_Depth_Land_And_Ocean',
#	variable='Effective_Optical_Depth_Average_Ocean',
#	variable='Effective_Optical_Depth_Average_Ocean',
	 **kwargs):
	'''
	Read nearest MODIS AOD within modis_thresh distance

	fname			str,	path to modis file
	lat				flt,	latitude
	lon				flt,	longitude
	modis_thresh	flt,	distance radius from lat/lon in km

	2016.02.17
	'''
	from netCDF4 import Dataset
	from libgeo import great_circle
	from numpy import array
	from libtools import modis2epoch

	cdf = Dataset(fname, 'r')
	try:
		vars = cdf.variables[variable]
	except KeyError:
		raise RuntimeError, 'File does not contain variable: {0} {1}'.format(
			fname, variable) 

	return cdf.variables[variable][1, x, y]


def read_modis_aod_latlon(fname, lat, lon, epoch, thresh_space=2e4,
	thresh_time=3600*6, 
	variable='Effective_Optical_Depth_Average_Ocean',
#	variable='Optical_Depth_Land_And_Ocean',
#	variable='Effective_Optical_Depth_Average_Ocean',
#	variable='Effective_Optical_Depth_Average_Ocean',
	 **kwargs):
	'''
	Read nearest MODIS AOD within modis_thresh distance

	fname			str,	path to modis file
	lat				flt,	latitude
	lon				flt,	longitude
	modis_thresh	flt,	distance radius from lat/lon in km

	2016.02.17
	'''
	from netCDF4 import Dataset
	from libgeo import great_circle
	from numpy import array
	from libtools import modis2epoch

	cdf = Dataset(fname, 'r')
	try:
		vars = cdf.variables[variable]
	except KeyError:
		raise RuntimeError, 'File does not contain variable: {0} {1}'.format(fname, variable) 

	#_match by lat/lon
	lats = cdf.variables['Latitude'][:].flatten()
	lons = cdf.variables['Longitude'][:].flatten()
	time = cdf.variables['Scan_Start_Time'][:].flatten()
	xx, yy = cdf.variables['Latitude'][:].shape

	#_do a quick check to see if any obs are within time window
	if min(abs(epoch - modis2epoch(time))) > thresh_time:
		print 'Skipping {0}: time'.format(fname)
		return False, 9e9, 9e9

	#_do a quick check to see if it's worth looking for the closest
	if min(abs(lat - lats)) > 3 or min(abs(lon - lons)) > 3:
		print 'Skipping {0}: distance'.format(fname)
		return False, 9e9, 9e9

	#_if multiple lat/lons passed, loop and find nearest spot for each
	req = (lat, lon); d = []
	[d.append(great_circle(req, mod)) for mod in zip(lats, lons)]
	ind = array(d).argmin()
	idx, idy = ind / yy, ind % yy

	if d[ind] > thresh_space:
		print 'None within acceptable limits, min distance {0:7.4f}'.format(d[ind])
		return False, 9e9, 9e9

	print cdf.variables[variable][:].shape, idx, idy 
	return cdf.variables[variable][1, idx, idy], d[ind], modis2epoch(time[ind]) 



