#!/usr/bin/env python
import sys, os


DIR_PROD = os.environ['PRODUCTS']

################################################################################
################################################################################
################################################################################
'''
When updating for other sensors, it will need to be passed
	sensor_source	-> tells it which oe_utils.read_sensor method to use
					-> format of collocation file
					-> height of sensor
					-> lat/lon/epoch attributes to pull out
	updated pickle with collocation <-> observation file lookup


	Search file for WRS
'''

def run(idx_cross, idx_along, experiment, rerun=False, 
	pk_file='/home/wsessions/data/updated_airscal.pk',#_WRS
	dir_airs=os.path.join(os.environ['PRODUCTS'], 'sips/airs'),#_WRS
	dir_lbl=os.path.join(DIR_PROD, 'LBL-RTM'),
	file_col=None, **kwargs):
	''' idx defines which profile index to run simulation on '''
	from lblrtm_utils import run_lblrtm_profile
	import re
	from pickle import load
	from oe_utils import read_profile_gdas, read_airs
	from libtools import epoch2dtg, remote_copy, dbg

	#_read in run options
	kw = {	'iatmm': 1, 
			'v1': 600.0,
			'v2': 2500.0,
			'ioddd': 0,
			'dummy_levs': 0.5,
			'dir_lblrtm_fmt': '/data/wsessions/LBL-RTM/{0}.{1}', #_experiment, dtg
			'cntnm': 5,
			'merge': 1, 
			'ipunch': 1,
			'iemit': 0,
			'nproc': 1,
			'dv': 0.5	} 
	#'clddef': [{'z': 2.0, 'z_top': 4.0}, {'z': 3.0, 'z_top': 5.0}, {'z': 4.0, 'z_top': 6.0}]}

	#_get press/prof data
	file_key = load(open(pk_file, 'rb'))[file_col.split('/')[-1]]['airs']#_WRS
	file_key = os.path.join(dir_airs, file_key.split('/')[-1])#_WRS

	#_copy over airs/gdas files

	gdas = read_profile_gdas(file_key, idx_cross, idx_along)
	airs, wave = read_airs(file_key, idx_cross, idx_along)#_WRS
	profile = gdas[0]
	radfile = airs[0] #_WRS

	#_open pickle, get out collocation files, use file_col to find associated
	# airs file, which will dictate the GDAS files yes this is dumb.	
	pressure	= profile.GDAS_pressure
	latitude	= radfile.AIRS_latitude #_WRS
	longitude	= radfile.AIRS_longitude #_WRS
	epoch		= radfile.AIRS_epoch #_WRS
##	latitude	= profile.GDAS_latitude
##	longitude	= profile.GDAS_longitude
	sensor_altitude = 7077.75 #e3	#_shoudl this be in km? 

	#_output and run directory
	dtg = epoch2dtg(epoch, full=True) 
	label = '{0}_{1}'.format(experiment, dtg)
	path_out = os.path.join(dir_lbl, label)
	src	= kwargs.get('profile_source', 'GDAS')
	kwargs.update({	'label'				: label, 
					'dir_lblrtm_out'	: path_out,
					'pressure'			: pressure,
					'h1'				: sensor_altitude,
					'latitude'			: latitude,
					'longitude'			: longitude,
					'epoch'				: epoch		})
	kwargs.update(kw)

	#_do not rerun if unnecessary
	if not rerun and os.path.exists(path_out):
		dbg('LBL-RTM already processed for {0}'.format(path_out))
		return 0

	#_run lblrtm from profile
	run_lblrtm_profile(profile, **kwargs)
	return 0	


################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
	#_first argument is the file number or profile index
	try:
		script, out_label, idx_cross, idx_along, file_col = sys.argv
		run(int(idx_cross), int(idx_along), out_label, file_col=file_col,
			lbldis=True)

	except:
		print '''		./run_lblrtm.py label idx_Xtrack, idx_Atrack file_col '''
		os._exit(0)

