#!/usr/bin/env python
# PURPOSE	Create netcdf files from xdr files and renames them based upon
#			time range. Relies upon Bob Holz's CPL_OPT to do actual conversion,
#			then renames.
# AUTHOR	Walter R Sessions
# DATE		2015.06.16


def main(indir='.', script='CPL_OPT'):
	from glob import glob
	from subprocess import call
	from os import rename

	#_get list of CPL files
	xdr_files = glob(indir + '/*.xdr')	

	for fxdr in xdr_files:
		#_create initial output filename
		fcdf = fxdr.replace('.xdr', '.nc')		
	
		#_convert file
		if call([script, fxdr, fcdf]):
			raise RuntimeError, 'Failed to convert file {0}'.format(fxdr)

		#_rename file
		fname = _gen_final_name(fxdr, fcdf)
		print fxdr, '=>', fname	
		rename(fcdf, fname)


def _gen_final_name(fxdr, fcdf):
	'''
	Given the original xdr filename and the dates from the 
	converted netcdf, return the final file name of the netcdf.
	'''
	import re

	#_open netcdf file, get start and end times
	yymmdd, shsmss, ehemes = _get_time_bounds(fcdf)

	#_get julian date from xdr file
	res = re.search('OP_(\d{5}\w?)_', fxdr)
	yyjjj = res.group(1)

	return 'OP_{0}_{1}_{2}_{3}.nc'.format(yyjjj, yymmdd, shsmss, ehemes)


def _get_time_bounds(fcdf):
	''' open netcdf file and get start and end times '''
	from netCDF4 import Dataset
	from libtools import julian2dtg, dtg2epoch, epoch2dtg

	with Dataset(fcdf, 'r') as cdf:
		jday = cdf.variables['Jday'][0]
		year = cdf.variables['Year'][0]
		time = cdf.variables['Time'][:]

		#_convert from julian to dtg
		dtg = julian2dtg(year, jday)
		ymd = dtg[2:-2]		

		#_get start and end times
		epoch = dtg2epoch(dtg) + (time - jday) * 86400
		start = epoch2dtg(epoch.min(), full=True)[-6:]
		final = epoch2dtg(epoch.max(), full=True)[-6:]

	return ymd, start, final 


if __name__ == '__main__':
	main()

