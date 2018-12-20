#!/usr/bin/env python
import sys, os


DIR_PROD = os.environ['PRODUCTS']
DIR_HS3 = os.path.join(DIR_PROD, 'LBL-RTM')


def lbldis_hs3(fov_idx, kw_idx, fname, kname, path_hs3=DIR_HS3,
	experiment='hs3', **kwargs):
	''' 
	submission script for running with defined r_eff, optical 
	depth and layer height
	'''
	from lblrtm_utils import run_lbldis
	from hs3_utils import Flight_segment
	import pickle

	#_read in flight options
	flight = Flight_segment(file_seg=fname) 

	#_read in run options
	file_kwarg = os.path.join(path_hs3, kname)
	kw = pickle.load(open(file_kwarg, 'rb'))[kw_idx]

    #_pull out cloud info if present for labeling
	if 'clddef' in kw:
		clddef = kw['clddef']
		ref = clddef[0]['r_eff']
		aod = clddef[0]['tau'][0]
		zkm = [clddef[0]['z'], clddef[0]['z_top']]
	else:
		ref, aod, zkm = 0., 0., [0., 0.]

	#_output and run directory 
	lblrtm_out = kw['dir_lblrtm_fmt'].format(fov_idx, flight.dtg0, flight.dtg1,
												experiment)
	run_lbldis(lblrtm_out, **kw) 
###_run_lbldis(lblrtm_out, out_lbldis=out_lbldis, **kw)

	return 0	


################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
	#_HS3 2013 run
	label = sys.argv[1]
	#_./run_lbldis.py experiment_label fov_idx kw_idx fname kname
	f, k, fname, kname  = int(	sys.argv[2]),	int(sys.argv[3]), \
								sys.argv[4],	sys.argv[5]
	sys.exit(lbldis_hs3(f, k, fname, kname))	

