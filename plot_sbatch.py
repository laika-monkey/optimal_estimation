#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

def run_main(dtg='', **kwargs):
	from libtools import dummy_sbatch
	from hs3_utils import Flight_segment
	
	flight = Flight_segment(dtg)
	flight.plot(**kwargs)


if __name__ == '__main__':
	from pickle import load
	from os import unlink
	import sys
	script, kwarg_file = sys.argv
	kwargs = load(open(kwarg_file, 'rb'))
	run_main(**kwargs)
	unlink(kwargs_file)
