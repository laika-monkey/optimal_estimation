#!/usr/bin/env python
'''
meant to be used with bulk submission of jobs for CALZONE.ssec.wisc.edu

swap out which program is run run_main() 
'''


import os, sys, imp
n = int(sys.argv[1]) if len(sys.argv) > 1 else 99999


#############################################################################_80
#############################################################################_80
#############################################################################_80



#############################################################################_80
#############################################################################_80
#############################################################################_80


def lblrtm(nfile=999999):
	'''	loop over lblrtm files '''
	from glob import glob

	path_in = os.path.join(os.environ['PRODUCTS'], 'gdas')
	files = glob(path_in + '/*.hdf')[:nfile]

	#_setup submission command
	env = ','.join(['='.join((var, os.environ[var])) \
		for var in ['PYTHONPATH', 'PRODUCTS', 'PATH', 'LOG']])
	qsub = ' '.join((	'qsub -v', env, '-o', DIR_LOG, '-e', DIR_LOG,
						'-cwd -S /opt/ShellB3/bin/python'))
	srpt = os.path.expanduser('~/lib/run_lblrtm.py')

	#_loop over collocated files, run lblrtm	
	for file in files: 
		cmd = ' '.join((qsub, srpt, file))
		os.system(cmd)  


def lbldis(verbose='1', nscene=99999):
	''' loop over lblrtm output directories and run lbldis '''
	from glob import glob

	#_generate list of all output directories following COL* naming
	lblrtm_dirs = []
	path_in = os.path.join(os.environ['PRODUCTS'], 'LBL-RTM')
	for file in glob(path_in + '/COL*'):
		lblrtm_dirs.extend(glob(file + '/*'))

	#_setup submission command
	env = ','.join(['='.join((var, os.environ[var])) \
		for var in ['PYTHONPATH', 'PRODUCTS', 'PATH', 'LOG']])
	env += ',LD_LIBRARY_PATH=/opt/netcdf-4.1.2/lib'
	qsub = ' '.join((	'qsub -v', env, '-o', DIR_LOG, '-e', DIR_LOG, 
						'-cwd -S /opt/ShellB3/bin/python')) #/bin/bash'))
	srpt = os.path.expanduser('~/lib/run_lbldis.py')

	#_loop over each lblrtm output and make input file, submit
	for lblrtm_out in lblrtm_dirs[:nscene]:
		cmd = ' '.join((qsub, srpt, lblrtm_out))	
		os.system(cmd)


def lblrtm_hs3(nfov=99999, file_flight='hs3_flight.pk', **kwargs):
	''' run model on hs3 profiles '''
	from hs3_utils import Flight_data
	from libtools import combination, mkdir_p
	import pickle
 
	#_read profiles from file 
	fname = os.path.join(DIR_PROD, 'LBL-RTM', 'hs3', file_flight) 
	flt = pickle.load(open(fname, 'rb'))

	#_options to loop over (for lbldis)	
	idx = range(flt.shis.size)	

	#_had to do all this to get the size of the collocated/avg
	#_setup submission command
	env = ','.join(['='.join((var, os.environ[var])) \
		for var in ['PYTHONPATH', 'PRODUCTS', 'PATH', 'LOG']])
	qsub = ' '.join((	'qsub -v', env, '-o', DIR_LOG, '-e', DIR_LOG,
						'-cwd -S /opt/ShellB3/bin/python'))
	srpt = os.path.expanduser('~/lib/run_lblrtm.py')

	#_loop over collocated files, run lblrtm
	for fov_idx in idx:
		#_write run options
		cmd = ' '.join((qsub, srpt, 'hs3', str(fov_idx), file_flight))
		os.system(cmd)  


def lbldis_hs3(nfov=99999, file_kw='kwarg_lbldis.pk', 
	file_flight='hs3_flight.pk', **kwargs):
	''' run model on hs3 profiles '''
	from libtools import combination
	import pickle
 
	#_read lbldis run options from file
	file_kwarg = os.path.join(DIR_PROD, 'LBL-RTM', 'hs3', file_kw)
	kw = pickle.load(open(file_kwarg, 'rb'))

	#_read profiles from file 
	file_flight= os.path.join(DIR_PROD, 'LBL-RTM', 'hs3', file_flight)
	flt = pickle.load(open(file_flight, 'rb'))

	#_had to do all this to get the size of the collocated/avg
	#_setup submission command
	env = ','.join(['='.join((var, os.environ[var])) \
		for var in ['PYTHONPATH', 'PRODUCTS', 'PATH', 'LOG']])
	qsub = ' '.join((	'qsub -v', env, '-o', DIR_LOG, '-e', DIR_LOG,
						'-cwd -S /opt/ShellB3/bin/python'))
	srpt = os.path.expanduser('~/lib/run_lbldis.py')

	#_loop over collocated files, run lblrtm
	runs = combination([xrange(flt.shis.size), xrange(len(kw))])
	for fov_idx, kw_idx in runs:
		#_write run options
		cmd = ' '.join((qsub, srpt, 'hs3', str(fov_idx), str(kw_idx)))
		os.system(cmd)  


def lblrtm_sbatch(fidx=[], rerun=True, experiment='hs3', **kwargs):
	'''
	2016.04.20
	Run LBL-RTM on profiles defined
	'''
	from libtools import combination, mkdir_p
	from hs3_utils import Flight_segment as read_flight_file
 
	#_had to do all this to get the size of the collocated/avg
	#_setup submission command
	env = ','.join(['='.join((var, os.environ[var])) \
		for var in ['PYTHONPATH', 'PRODUCTS', 'PATH', 'LOG']])
	qsub = ' '.join((	'qsub -v', env, '-o', DIR_LOG, '-e', DIR_LOG,
						'-cwd -S /opt/ShellB3/bin/python'))
	srpt = os.path.expanduser('~/lib/run_lblrtm.py')

	#_if none passed, do them all
	if len(fidx) == 0:
		flight = read_flight_file(file_seg=file_flight)
		fidx = range(flight.size)

	#_loop over collocated files, run lblrtm
	flight = read_flight_file(file_seg=file_flight)	#_DELETE
	for fov_idx in fidx:
		#_KLUDDDGE
		dtgfov = (flight.dtg0, flight.dtg1, fov_idx, experiment)
		dir_lblrtm = kwargs.get('dir_lblrtm_fmt').format(*dtgfov)
	
		#_skip already run fovs
		if os.path.exists(dir_lblrtm) and not rerun:
			print 'BAILING OUT, NO RERUN', dir_lblrtm
			continue

		print 'running {0}'.format(dir_lblrtm)

		#_write run options
		cmd = ' '.join((qsub, srpt, experiment, str(fov_idx), file_flight))
		print cmd
		os.system(cmd)  


def lblrtm_hs3_lim(fidx=[], file_flight='hs3_flight.nc', rerun=True,
	experiment='hs3', **kwargs):
	''' run model on hs3 profiles '''
	from libtools import combination, mkdir_p
	from hs3_utils import Flight_segment as read_flight_file
 
	#_had to do all this to get the size of the collocated/avg
	#_setup submission command
	env = ','.join(['='.join((var, os.environ[var])) \
		for var in ['PYTHONPATH', 'PRODUCTS', 'PATH', 'LOG']])
	qsub = ' '.join((	'qsub -v', env, '-o', DIR_LOG, '-e', DIR_LOG,
						'-cwd -S /opt/ShellB3/bin/python'))
	srpt = os.path.expanduser('~/lib/run_lblrtm.py')

	#_if none passed, do them all
	if len(fidx) == 0:
		flight = read_flight_file(file_seg=file_flight)
		fidx = range(flight.size)

	#_loop over collocated files, run lblrtm
	flight = read_flight_file(file_seg=file_flight)	#_DELETE
	for fov_idx in fidx:
		#_KLUDDDGE
		dtgfov = (flight.dtg0, flight.dtg1, fov_idx, experiment)
		dir_lblrtm = kwargs.get('dir_lblrtm_fmt').format(*dtgfov)
	
		#_skip already run fovs
		if os.path.exists(dir_lblrtm) and not rerun:
			print 'BAILING OUT, NO RERUN', dir_lblrtm
			continue

		print 'running {0}'.format(dir_lblrtm)

		#_write run options
		cmd = ' '.join((qsub, srpt, experiment, str(fov_idx), file_flight))
		print cmd
		os.system(cmd)  

def lbldis_hs3_lim(kidx=[0], fidx=[], file_kw='kwarg_lbldis.pk', 
	file_flight='hs3_flight.nc', experiment='hs3', **kwargs):
	''' run model on hs3 profiles '''
	from libtools import combination
	from hs3_utils import Flight_segment as read_flight_file
 
	#_had to do all this to get the size of the collocated/avg
	#_setup submission command
	env = ','.join(['='.join((var, os.environ[var])) \
		for var in ['PYTHONPATH', 'PRODUCTS', 'PATH', 'LOG']])
	qsub = ' '.join((	'qsub -v', env, '-o', DIR_LOG, '-e', DIR_LOG,
						'-cwd -S /opt/ShellB3/bin/python'))
	srpt = os.path.expanduser('~/lib/run_lbldis.py')

	#_if none passed, do them all
	if len(fidx) == 0:
		flight = read_flight_file(file_seg=file_flight)
		fidx = range(flight.size)

	runs = combination([fidx, kidx])
	print fidx, kidx
	for f, k in runs:
		#_write run options
		cmd = ' '.join((qsub,srpt,experiment,str(f),str(k),file_flight,file_kw))
		os.system(cmd)  
		print cmd


################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
	sys.exit(run_main(n=n))
