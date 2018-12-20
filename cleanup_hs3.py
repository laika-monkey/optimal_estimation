#!/usr/bin/env python


def main():
	import re
	import os
	from glob import glob
	from shutil import rmtree
 
	dir_hs3 = os.path.join(os.environ['PRODUCTS'], 'LBL-RTM')

	re_hs3 = re.compile('HS3_(\d{14}).(\d{14})_fov\d{6}')

	#_get list of dtgs to keep
	dtgs = get_list_dtgs()

	keep, kill = [], []
	for directory in glob('{0}/HS3_*'.format(dir_hs3)):
		res = re_hs3.search(directory)

		if res:
			if res.group(1) in dtgs:
				keep.append(directory)
			else:
				kill.append(directory)

	print len(keep), len(kill)

	for p in kill:
		print p
		rmtree(p)


def get_list_dtgs():
	from flight_namelists import experiments
	from hs3_utils import Flight_segment as F
	segments = experiments['HS3']
	dtgs = []
	for dtg, opts in segments.iteritems():
		dtgs.append(F(dtg=dtg).dtg0)

	return dtgs


if __name__ == '__main__':
	main()
