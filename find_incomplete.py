#!/usr/bin/env python
# some lblrtm runs weren't completed, find them


def main():
	import os
	import re
	from lblrtm_utils import read_tape5
	from glob import glob
	from libtools import epoch2dtg, dtg2epoch
	from pickle import dump

	#_directory with all LBLRTM runs
	dir_lbl = os.path.join(os.environ['PRODUCTS'], 'LBL-RTM') 
	dirs = glob(dir_lbl + '/*')

	#_build regular expression to pull out dates
	re_dtg = re.compile('\w+_(\d{14}).\d{14}_fov(\d{6})')

	#_initialize dictionary of missing fovs
	missing = {}

	#_loop over each dir, see if they have TAPE5
	for directory in dirs:
		print 'Checking {0}'.format(directory)

		ftape5 = os.path.join(directory, 'TAPE5')
		if not os.path.exists(ftape5):
			continue

		#_get number of levels
		nz = int(read_tape5(ftape5)[2]['ibmax'])

		#_check that all output levels exist
		for i in range(nz):
			fOD = 'ODdeflt_{0:03d}'.format(i+1)
			pOD = os.path.join(directory, fOD)

			if not os.path.exists(pOD):
				#_pull out a date and fov
				regex = re_dtg.search(directory)
				dtg, fov = regex.group(1), int(regex.group(2))			

				#_write missing fov to file
				if dtg not in missing:
					missing.update({dtg : [fov]})
				else:
					missing[dtg].append(fov)

				#_go ahead and rerun
		##		os.chdir(directory)
		##		os.system('lblrtm')
			
				#_abandon this dir
				break

		else:
			#_all files accounted for
			continue

	#_dump missing into pickle
	f = dump(missing, open('lbl_incomplete.pk', 'wb'))
	
	#_then print them
	for dtg, fovs in missing.iteritems():
		print dtg, len(fovs), fovs
	

if __name__ == '__main__':
	main()
