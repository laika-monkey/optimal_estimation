#!/usr/bin/env python
# script	find_clear_sky_cases.py
# author	walter r sessions
# purpose	find field of views without cpl aerosols or clouds reported

#_2015.05.14	THIS NEEDS TO BE RUNRUN FROM MAIN ONWARD
#				The old clear sky files are no longer named the same,
#				which I find confusing, but whatever.  COLLOCs are now different

#_

import os
DIR_PLOT = os.path.join(os.environ['PRODUCTS'], 'plots')
namelist = {
	'experiment'	: 'HS3',
	'out_label'		: 'clear-sky', 
	'surf_temp_src'	: 'GDAS',

	'nproc'			: 20,

	'lblrtm'		: False, 
	'lbldis'		: True,

	#_use mw
	'dv'			: -27,

	#_dir plot
	'dir_plot'	: os.path.join(DIR_PLOT, 'clear_sky_fovs'),
	'dir_out'	: os.path.join(DIR_PLOT, 'clear_sky_fovs'),
	}


def main(dir_hs3=os.path.join(os.environ['PRODUCTS'], 'hs3')):
	'''
	Auto generates dictionary of cases with no AOD reported in cpl.

	Hand filter afterward because it is full of crap
	'''

	from hs3_utils import Flight_segment as F
	from glob import glob
	import re
	from numpy import arange, vstack

	re_colloc = re.compile('SHIS.CPL.GDAS.COLLOC.(\d{12}).(\d{12}).nc')

	#_initialize output dictionary
	output = {}

	#_get all flight segment files
	for fname in glob(dir_hs3 + '/*'):
		res = re_colloc.search(fname)

		#_skip if sile doesn't match re_colloc
		if not res:
			continue

		#_read in flight segment
		flight = F(file_seg=fname)
		nfov = flight.size

		#_check cpl data for clear skies
		cpl_tau = vstack(flight.CPL_tau_532)
		cpl_typ = vstack(flight.CPL_tau_type)

		#_look for clear
		clr_typ = (cpl_typ == 0).all(axis=1)
		clr_tau = (cpl_tau <= 0).all(axis=1)

		#_where both are true
		fidx = arange(nfov)[clr_tau]

		#_exit after first just because
		output.update({ fname : { 'fidx' : fidx.tolist() } })
		

	#_write dictionary in text file
	with open('clear_sky_fovs.py', 'w') as f:
		f.write("clear_sky = {\\\n")

		#_loop over
		for fname, opts in output.iteritems():
			f.write("'{0}' : {1},\n".format(fname, str(opts)))
	
		f.write("}\n")


def run_clear_skies(lblrtm=False, lbldis=False,
	fmt = '/data/wsessions/LBL-RTM_hs3/{3}_{0}.{1}_fov{2:06d}',
	ecmwf_file='/data/wsessions/era-interim/hs3_2013.nc', **kwargs):
	''' read in field of views and run cases through lblrtm and lbldis '''
	from clear_sky_fovs import clear_sky as cs
	from hs3_utils import Flight_segment as F
	from qsubmissions import lblrtm_hs3_lim, lbldis_hs3_lim
	from optimal_estimation import clear_sky_radiance as csr
	from libtools import setup_groups
	import os, sys

	if lblrtm and lbldis:
		raise RuntimeError, 'Submit LBLRTM and LBLDIS separate, not synced'

	#_run lblrtm for all
	if lblrtm:
	  for file_seg, opts in cs.iteritems():

		#_check if files exist
		if not os.path.exists(file_seg):
			raise IOError, 'Missing file {0}'.format(file_seg)

		#_read in flight
		flight = F(file_seg=file_seg)
		kwargs.update({'file_flight' : file_seg})
		kwargs.update(opts)

		#_run lblrtm
		lblrtm_hs3_lim(dir_lblrtm_fmt=fmt, **kwargs)

	#_DO EITHER LBLRTM OR LBLDIS, NOT BOTH
	if lbldis:
	  experiment = kwargs.get('experiment')
	  for file_seg, opts in cs.iteritems():

		#_check if files exist
		if not os.path.exists(file_seg):
			raise IOError, 'Missing file {0}'.format(file_seg)

		flight = F(file_seg=file_seg)
	
		#_run lbldis for all
		groups = setup_groups(opts['fidx'], **kwargs)
		for group in groups:
			children = []
			for i in group:
				pid = os.fork()
				if pid != 0:
					children.append(pid)
				else:
					surf_temp = get_surf_temp(flight[i], **kwargs)
					kwargs.update({'surf_temp':surf_temp})

					#_run lbldis
					args =  (flight.dtg0, flight.dtg1, i, experiment) 
					dir_lblrtm = fmt.format(*args)

					#_run clear sky case
					csr(dir_lblrtm, **kwargs)

					os._exit(0)

			for kid in children:
				os.waitpid(kid, 0)
			

def calc_bias(experiment='hs3', dv=None, dir_out=DIR_PLOT, 
	fmt = '/data/wsessions/LBL-RTM/{3}_{0:s}.{1:s}_fov{2:06d}', **kwargs):
	''' read in all the output from run_clear_skies and get average bias '''
	from clear_sky_fovs_hand_filtered import clear_sky as cs
	from hs3_utils import Flight_segment as F
	import os
	from lblrtm_utils import read_lbldis_out, microwindow_average
	from numpy import append
	from libgeo import planck_inv
	
	#_open a csv to put all the values into
	fname = 'shis-lbl.{0}.csv'.format(experiment)
	fname = os.path.join(dir_out, fname)
	f = open(fname, 'w')
	f.write('SHIS - LBL\n')

	first = True
	for file_seg, opts in cs.iteritems():
		#_read in flight
		flight = F(file_seg=file_seg)

		for i in opts['fidx']:
			#_construct directory name
			dir_lblrtm = fmt.format(flight.dtg0, flight.dtg1, i, experiment)
			fname = os.path.join(dir_lblrtm, 'lbldis_output.clear.cdf')

			print fname	
			#_read in data
			clear = read_lbldis_out(fname)			

			rad_shs = flight[i].SHIS_radiances
			wvn_shs = flight.SHIS_wavenumber
			rad_lbl = clear.radiances
			wvn_lbl = clear.wavenumber
		
			#_if run on microwindows, only average the obs
			if dv < 0:
				rad_shs, wvn = microwindow_average(rad_shs, wvn_shs, dv)
 
			else:	
				rad_shs, wvn = microwindow_average(rad_shs, wvn_shs, -7) 
				rad_lbl, wvn = microwindow_average(rad_lbl, wvn_lbl, -7) 

			tb_shs = planck_inv(rad_shs, wvn*100, domain='wavenumber')
			tb_lbl = planck_inv(rad_lbl, wvn*100, domain='wavenumber')

			#_write wavenumbers first time
			if first:
				fmt0 = ','.join(['{{{0:d}:10.5f}}'.format(j) 
					for j in range(wvn.size)])
				f.write(fmt0.format(*wvn) + '\n')

			#_put into one thinger	
			rads = append(tb_shs, tb_lbl).tolist()
			fmt0 = ','.join(['{{{0:d}:10.5f}}'.format(j) 
				for j in range(wvn.size*2)])
			fmt0 += '\n'
		
			f.write(fmt0.format(*rads))
			first = False
	
	#_write
	f.close()
	print fname	
	return fname


def plot_bias(fname=None, out_label='', experiment='hs3', dir_out='.', 
	dir_plot='.', **kwargs): 
	''' read in output from calc bias, plot '''
	import numpy as np
	import matplotlib.pyplot as plt
	from libtools import shrink_ticks

	#_
	if fname is None:
		fname = 'shis-lbl.{0}.csv'.format(experiment)
		fname = os.path.join(dir_out, fname)

	#_open file
	f = open(fname, 'r')

	#_get wavenumbers
	f.readline()
	wvn = f.readline().split(',')
	wvn = [float(a) for a in wvn]

	data = []
	for line in f.readlines():
		data.append([float(a) for a in line.split(',')])

	data = np.array(data)
	nfov, nchan2 = data.shape

	#_shis is [:,0,:], lbldis is [:,1,:]
	data = data.reshape(nfov, 2, nchan2/2)

	diff = data[:,1,:] - data[:,0,:]
	bias_mu = diff.mean(0)
	bias_sd = diff.std(0)

	h, w = plt.figaspect(6)
	fig = plt.figure(figsize=(w,h))
	ax = fig.add_subplot(111)
	#fig, ax = plt.subplots(1)

	with open('/home/wsessions/bias.py', 'w') as f:
		for i, wave in enumerate(wvn):
			f.write("{0:7.3f}\t:\t{1:8.5f}\n".format(wave, bias_mu[i]))

	arg = (experiment, diff.shape[0])
	ax.set_title('Clear Sky Bias (LBLDIS minus S-HIS)')
#	ax.set_title('LBLDIS minus S-HIS, experiment={0}, N={1:d}'.format(*arg))
	ax.scatter(wvn, bias_mu)
	ax.errorbar(wvn, bias_mu, yerr=bias_sd, linestyle='None')
	ax.set_ylim(-3, 3)
	ax.set_xlim([800, 1200])
	ax.grid(True)
	ax.set_xlabel('wavenumber')
	ax.set_ylabel('$\Delta B_{\\nu}T$ (K)')
	shrink_ticks(ax)
	pname = 'bias.{0}.png'.format(out_label)
	pname = os.path.join(dir_plot, pname)
	print pname
	plt.savefig(pname)

	#_close up file
	f.close()


################################################################################
################################################################################
################################################################################


def get_surf_temp(fov, surf_temp_src='GDAS', varname='skt', #_sst || skt
	ecmwf_file='/data/wsessions/era-interim/hs3_2013.nc', **kwargs):
	'''
	return surf temperature depending on source

	The default behavior of the code is to use GDAS
	or LBL-ATM based surface temperatures.  If the
	surf_temp_src str is passed as ECMWF, this method
	returns a reaaalllly poorly done nearest neighbor
	interp from the ERA Interim data.
	'''
	from netCDF4 import Dataset
	from libtools import ecmwf_day2dtg as d2d
	from libtools import dtg2epoch as d2e 
	from libtools import epoch2dtg as e2d 
	from numpy import array

	#_only do this if ECMWF	
	if surf_temp_src == 'GDAS':
		return fov.GDAS_sfc_temperature
	elif surf_temp_src == 'LBL-ATM':
		return -1
	else:
		pass

	#_open netcdf file
	ncdf = Dataset(ecmwf_file, 'r')

	#_get indices
	yidx = abs(fov.SHIS_latitude - ncdf.variables['latitude']).argmin() 
	xidx = abs(fov.SHIS_longitude - ncdf.variables['longitude']).argmin() 
	epoch = array([d2e(d2d(t)) for t in ncdf.variables['time']])
	tidx = abs(fov.SHIS_epoch - epoch).argmin()

	print d2d(ncdf.variables['time'][tidx]), ncdf.variables['latitude'][yidx], \
		ncdf.variables['longitude'][xidx]
	print e2d(fov.SHIS_epoch), fov.SHIS_latitude, fov.SHIS_longitude

	surf_temp = ncdf.variables[varname][tidx,yidx,xidx]
	ncdf.close()

	return surf_temp


def plot(**kwargs):
	''' plotting clear sky case dashboards '''
	from clear_sky_fovs import clear_sky as cs
	from hs3_utils import Flight_segment as F

	for file_seg, opts in cs.iteritems():
		#_read in flight
		flight = F(file_seg=file_seg)

		#_get indices
		fidx = opts['fidx']
	
		#_plot just clear sky instances
		flight.plot_flight(fidx=fidx, **kwargs)

	
if __name__ == '__main__':
	#_generate list, filter for use in later tests
##	main()				#_find indices and produce text file of them
##	plot(**namelist)	#_take those and plot out flights to filter by sight

	#_to rerun
##	run_clear_skies(**namelist)		#_run lblrtm and lbldis over clear sky cases
##	calc_bias(**namelist)			#_read in output of rcs() and get avg bias
	plot_bias(**namelist)
