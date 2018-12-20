#!/usr/bin/env python
'''
meant to be used with bulk submission of jobs for CALZONE.ssec.wisc.edu

swap out which program is run run_main() 
'''


import os, sys, imp
n = int(sys.argv[1]) if len(sys.argv) > 1 else 99999


#_SGE log directory	
# DIR_LOG	= os.path.expanduser('~/qsub_logs/')
DIR_LOG	= os.environ['LOG'] 
DIR_PROD = os.environ['PRODUCTS']
DIR_WORK = os.environ['WORK']
DIR_HS3 = os.path.join(DIR_PROD, 'LBL-RTM', 'hs3')
DIR_TMP = os.path.join(DIR_WORK, 'kwargs')

#############################################################################_80
#############################################################################_80
#############################################################################_80


def run_main(n=1):
	return 0


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


def lblrtm_sbatch(file_col, fidx=None, experiment='test',
	dir_tmp=os.environ['WORK'], dir_batch=os.path.join(os.environ['WORK'],
	'batch_scripts'), **kwargs):
	'''
	2016.04.20
	Run LBL-RTM on profiles defined

	fidx limits the indices run
	'''
	from libtools import combination, dbg, mkdir_p, capture_stdout
	from pyhdf.SD import SD, SDC
	import re

	#_had to do all this to get the size of the collocated/avg
	#_setup submission command
	scpt = os.path.expanduser('~/lib/run_lblrtm.py')

	#_open collocation file
	hdf = SD(file_col, SDC.READ)
	idx = range(hdf.select('Master_Index_1')[:].size) if fidx is None else fidx

	#_loop over collocated files, run lblrtm
	for xxx, fov_idx in enumerate(idx):
	
		#_pull out cross and along track indices
		idx_x = hdf.select('Master_Index_1')[fov_idx] - 1	
		idx_a = hdf.select('Master_Index_2')[fov_idx] - 1

		#_create job label for sbatch
		job_label='{0}_{1:04d}_{2:04d}'.format(file_col.split('/')[-1][:-4], idx_x, idx_a) 

		#_initialize output for batch script
		out = '#!/bin/bash\n'
		out += '#SBATCH --job-name={0}\n'.format(job_label)
		out += '#SBATCH --partition=all\n'
		out += '#SBATCH --share\n'
		out += '#SBATCH --time=2:00:00\n'
		out += '#SBATCH --ntasks=1\n'
		out += '#SBATCH --cpus-per-task=1\n'
		out += '#SBATCH --nice=1000\n'
		out += '#SBATCH --output=/odyssey/scratch/%u/logs/{0}-%A.txt\n'.format(job_label)
		out += 'module purge\n'
		out += 'module load license_intel\n'
		out += 'module load impi\n'
		out += 'module load intel/15.0-2\n'
		out += 'module load hdf5/1.8.14\n'
		out += 'module load netcdf4/4.3.3\n'
		out += 'module load anaconda27/base_2015a_sharppy13\n'
		out += 'export TMPDIR={0}/${{SLURM_JOB_NAME}}.${{SLURM_JOB_ID}}\n'.format(dir_tmp)
		out += 'mkdir -p $TMPDIR\n'
	
		#_add environment variables to script
		out += '\n'.join(['export {0}={1}'.format(var, os.environ[var]) \
            for var in ['PYTHONPATH', 'PRODUCTS', 'WORK', 'PATH', 'LOG']])

		#_DO THIS ONCE WE USE OTHER SENSORS WRS
		'''
		#_generate kwarg file that scrpt will load to get settings
        fmt = 'kwgs.{4}.{2}.{3}.{0:04d}.{1:07d}'
        fmt_args = (i, pid, flight.dtg0, flight.dtg1, out_label)
        file_kwarg = os.path.join(DIR_TMP, fmt.format(*fmt_args))
        pickle.dump(kwargs, open(file_kwarg, 'wb'))	
		'''
   
		#_build actual script call string
		out += 'source activate my_root\n'
		out += 'echo `which python`\n'
		out += ' '.join((scpt, experiment, str(idx_x), str(idx_a), file_col)) + '\n'
		out += 'rm -rf $TMPDIR\n'

		#_write batch script
		sname = '{0}/sbatch.lblrtm.{1}'.format(dir_batch, job_label)
		with open(sname, 'w') as f:
			print out
			f.write(out)
		        #_put it all together and submit
		print sname

        #_GET BATCH ID AND PASS TO NEXT AS OK IF PERSISTENCE SET
		cmd = ' '.join(('sbatch', sname))
		dbg(cmd)
		stdout = capture_stdout(cmd)

		#_get ID
		last_jid = int(re.search('Submitted batch job (\d+)', stdout[0]).group(1))



def lbldis_sbatch(file_col, fidx=None,  
	dir_tmp=os.environ['WORK'], dir_batch=os.path.join(os.environ['WORK'],
	'batch_scripts'), hostname=False, **kwargs):
	'''
	2016.04.20
	Run LBL-RTM on profiles defined

	fidx limits the indices run

	purpose of thise level is to generate appropriate sbatch and kwargs option files
	'''
	from libtools import combination, dbg, mkdir_p, capture_stdout
	from pyhdf.SD import SD, SDC
	from pickle import dump
	import re

	#_had to do all this to get the size of the collocated/avg
	#_setup submission command
	scpt = os.path.expanduser('~/lib/run_oe.py')

	#_open collocation file
	hdf = SD(file_col, SDC.READ)
##	idx = range(hdf.select('Master_Index_1')[:].size) if fidx is None else fidx
	idx = fidx

	#_get out_label
	out_label = kwargs.get('out_label')

	#_loop over collocated files, run lblrtm
	for xx, fov_idx in enumerate(idx):
	
		#_create job label for sbatch
		job_label='{0}_{1:05d}'.format(file_col.split('/')[-1][:-4], fov_idx) 

		#_initialize output for batch script
		out = '#!/bin/bash\n'
		out += '#SBATCH --job-name={0}\n'.format(job_label)
		out += '#SBATCH --partition=all\n'
		out += '#SBATCH --share\n'
		out += '#SBATCH --time=2:00:00\n'
		out += '#SBATCH --ntasks=1\n'
		out += '#SBATCH --cpus-per-task=1\n'
		out += '#SBATCH --nice=1000\n'
		out += '#SBATCH --output=/odyssey/scratch/%u/logs/{0}-%A.txt\n'.format(job_label)
		out += 'module purge\n'
		out += 'module load license_intel\n'
		out += 'module load impi\n'
		out += 'module load intel/15.0-2\n'
		out += 'module load hdf5/1.8.14\n'
		out += 'module load netcdf4/4.3.3\n'
		out += 'module load anaconda27/base_2015a_sharppy13\n'
		out += 'export TMPDIR={0}/${{SLURM_JOB_NAME}}.${{SLURM_JOB_ID}}\n'.format(dir_tmp)
		out += 'mkdir -p $TMPDIR\n'
	
		#_add environment variables to script
		out += '\n'.join(['export {0}={1}'.format(var, os.environ[var]) \
            for var in ['PYTHONPATH', 'PRODUCTS', 'WORK', 'PATH', 'LOG']])

		#_if needing to remote copy, add hostname to front
		args = ('.'.join(file_col.split('/')[-1].split('.')[1:3]),
			fov_idx, out_label)
		kname = 'kwgs.{0}.{1:04d}.{2}.pk'.format(*args)
		file_kwarg = os.path.join(DIR_TMP, kname)
		if hostname:
			file_tmp = '{0}:{1}'.format(hostname, file_col)
		else:
			file_tmp = file_col

		sname = '{0}/sbatch.lbldis.{1}'.format(dir_batch, job_label)
		kwargs.update({	'file_col'		: file_tmp,
						'sbatch_file'	: sname	})	

		#_DO THIS ONCE WE USE OTHER SENSORS WRS
		#_generate kwarg file that scrpt will load to get settings
		dump(kwargs, open(file_kwarg, 'wb'))	
		if hostname:
			file_tmp = '{0}:{1}'.format(hostname, file_kwarg)  
		else:
			file_tmp = file_kwarg

		#_build actual script call string
		out += 'source activate my_root\n'
		out += 'echo `which python`\n'
		out += ' '.join((scpt, str(fov_idx), file_tmp)) + '\n'
		out += 'rm -rf $TMPDIR\n'

		#_write batch script
		with open(sname, 'w') as f:
		##	print out
			f.write(out)
		        #_put it all together and submit
		print sname

        #_GET BATCH ID AND PASS TO NEXT AS OK IF PERSISTENCE SET
	##	cmd = ' '.join(('sbatch', dependency, sname))
		cmd = ' '.join(('sbatch', sname))
		dbg(cmd)
		stdout = capture_stdout(cmd)

		#_get ID
		last_jid = int(re.search('Submitted batch job (\d+)', stdout[0]).group(1))


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
