#!/usr/bin/env python
#_wrapper to plot all flight segments with COLLOC files

import os
from glob import glob
from hs3_utils import Flight_segment
import re

re_dtg = re.compile('COLLOC.(\d{12}).')

#dir_col = os.path.join(os.environ['PRODUCTS'], 'tc4')
dir_col = os.path.join(os.environ['PRODUCTS'], 'hs3')
files = glob(dir_col + '/*COLLOC*nc')

kwargs = {'nproc' : 20, 'out_label' : 'FLIGHT-ONLY-JULY01'}
files.sort()
for fname in files:
	dtg = '20{0}'.format(re_dtg.search(fname).group(1))
##	if dtg < '20130825000000':
##		print 'ALREADY DONE', dtg
##		continue

	flight = Flight_segment(file_seg=fname)

	try:
		flight.plot_flight(dir_out=dir_col, **kwargs)
	except:
		os.system('echo {0} >> PLOT_ALL.err'.format(fname))
