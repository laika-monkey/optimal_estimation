#!/usr/bin/env python


def run_main():
	#_count number of fields of view currently being used, 
	# how many there are total possible
	# And report on current usage and potential
	import subprocess as s
	from glob import glob
	import numpy as np
	import os

	dir_hs3 = '/data/wsessions/hs3'
	dir_lbl = '/data/wsessions/LBL-RTM_hs3'
	sizes = []
	files = glob(dir_lbl + '/HS3*fov*')
	print dir_lbl
	mods = { 'M' : 1, 'G' : 10e3, 'T' : 10e6 }
	for f in files:
		#_start process, capture output
		out = s.Popen(['du', '-hs', f], stdout=s.PIPE).communicate()[0]

		#_read in size 
		mod = out.split('\t')[0][-1:]
		size = float(out.split('\t')[0][:-1])
		size *= mods[mod]
	
		sizes.append(size)	

	#_convert to math friendly
	sizes = np.array(sizes)

	#_print stats about current scenario
	arg = (sizes.sum(), sizes.mean(), sizes.max(), sizes.min())
	print 'total: {0}, mean: {1}, max: {2}, min: {3}\n'.format(*arg) 

	#_get all potential fields of view
	nfovs = []
	hs3_2013_files = glob(dir_hs3 + '/SHIS.CPL.GDAS.COLLOC.13*nc')
	for f in hs3_2013_files:
		from netCDF4 import Dataset as D
		#_open file, get size of fov
		nfov = len(D(f, 'r').dimensions['fov'])

		nfovs.append(nfov)

	nfovs = np.array(nfovs)
	tot_poss = nfovs.sum()
	tot_size = tot_poss * sizes.mean()
	tot_diff = tot_size - sizes.sum() 
	arg = (tot_poss, tot_size, tot_diff)
	print 'total_possible: {0}, total_size: {1}, total_diff: {2}\n'.format(*arg)

if __name__ == '__main__':
	run_main()
