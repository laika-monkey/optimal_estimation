#!/usr/bin/env python


import os
import sys


def main(fname, fkwargs):
	from hs3_utils import read_collocated, write_collocated
	import pickle
	from os import unlink

	#_get kwargs
	kwargs = pickle.load(open(fkwargs,'rb'))
	print kwargs

	#_read in shis, cpl, gdas
	s, c, p = read_collocated(fname, **kwargs)

	#_write segment file
	write_collocated(s, c, p, fname, **kwargs)

	if fkwargs[-2:] == 'pk':
		unlink(fkwargs)


if __name__ == '__main__':
	if os.environ['HOSTNAME'] != 'iris.ssec.wisc.edu':
		fname, fkwargs = sys.argv[1:]
		main(fname, fkwargs)
	else:
		#_run sbatch'd version
		pass
