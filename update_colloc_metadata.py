#!/usr/bin/env python

#_for some reason, filenames not in collocation thing.  This fixes that

from pyhdf.SD import SD, SDC
from pickle import load
import os

pname = '/home/wsessions/data/updated_airscal.pk'
pickle = load(open(pname, 'rb'))

#_local directories
dir_colloc = os.path.join(os.environ['WORK'], 'colloc')
dir_airs = os.path.join(os.environ['PRODUCTS'], 'sips/airs')
dir_caliop = os.path.join(os.environ['PRODUCTS'], 'sips/caliop')

for file_col, file_dict in pickle.iteritems():

	#_fix local location
	file_col = os.path.join(dir_colloc, file_col.split('/')[-1])

	#_pull out airs and caliop files and fix them, too
	file_airs = file_dict['airs']
	file_caliop = file_dict['caliop']
	file_airs = os.path.join(dir_airs, file_airs.split('/')[-1]) 
	file_caliop = os.path.join(dir_caliop, file_caliop[0].split('/')[-1]) 

	#_write these to the file_col
	hdf = SD(file_col, SDC.WRITE)
	att = hdf.attr('fname_AIRS')
	att.set(SDC.CHAR8, file_airs)
	att = hdf.attr('fname_CALIOP')
	att.set(SDC.CHAR8, file_caliop)

	hdf.end()

	#_test output
	hdf = SD(file_col, SDC.READ)
	print hdf.attributes()
	hdf.end()
