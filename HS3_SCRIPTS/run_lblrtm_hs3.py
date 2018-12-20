#!/usr/bin/env python
import sys, os


DIR_PROD = os.environ['PRODUCTS']

#_dummy_levs controls minimum space between vertical levels
kwargs = {'dummy_levs' : 0.5}

################################################################################
################################################################################
################################################################################


def run(fov_idx, dir_lbl=os.path.join(DIR_PROD, 'LBL-RTM'), file_seg=None,
	label='hs3', **kwargs):
	''' idx defines which profile index to run simulation on '''
	from lblrtm_utils import run_lblrtm_profile
	from hs3_utils import Flight_segment
	import pickle

	#_read in run options
	file_kwarg = os.path.join(dir_lbl, 'kwarg_lblrtm.pk')
	kw = pickle.load(open(file_kwarg, 'rb'))

	#_read in flight data
	flight = Flight_segment(file_seg=file_seg, **kwargs)
	
	#_output and run directory 
	label = '{3:s}_{1:s}.{2:s}_fov{0:06d}'.format(fov_idx, 
		flight.dtg0, flight.dtg1, label)
	path_out = os.path.join(dir_lbl, label)
	src	= kwargs.get('profile_source', 'GDAS')
	pressure = flight.__getattribute__('{0}_pressure'.format(src)) 
	kwargs.update({	'label'				: label, 
					'dir_lblrtm_out'	: path_out,
					'pressure'			: pressure })
	kwargs.update(kw)

	#_run lblrtm from profile
	profile = flight[fov_idx]
	run_lblrtm_profile(profile, h1=profile.SHIS_altitude, **kwargs)
	return 0	


################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
	#_first argument is the file number or profile index
	try:
		label = sys.argv[1]
		idx_fov = int(sys.argv[2])
		run(idx_fov, file_seg=sys.argv[3], lbldis=True, label=label, **kwargs)

	except:
		print '''		./run_lblrtm.py hs3 fov_idx file_flight '''
		os._exit(0)

