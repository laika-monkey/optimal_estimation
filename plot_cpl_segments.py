#!/usr/bin/env python
#
#    NAME:   plot_all_cpl.py
#  AUTHOR:   Walter R. Sessions
#    DATE:   2015.Early
# PURPOSE:   Create multi-panel image showing backscatter with AOD overlaid 
#            from CPL data.
#			THIS VERSION USES FLIGHT FILES GENERATED FROM COLLOCATED HS3
#			FIELDS.  TO PLOT FROM XDR/NC, USE plot_cpl_all.py

import os, sys

namelist = {
	#_if a subset is desired
	'start_dtg'	: '20130101000000',
	'end_dtg'   : '20140101000000',

	#_cloud or aerosol?
	'type_coa' : 'a',

	#_plotting directory
	'dir_plot'	: os.path.join(os.environ['PRODUCTS'], 'plots','cpl'),
}

def main(dir_cpl='{0}/cpl'.format(os.environ['PRODUCTS']), start_dtg=None,
	dir_plot='.', end_dtg=None, type_coa='a', **kwargs):
	'''
	start_dtg   str{14},    Starting date-time-group for subsetting
	end_dtg     str{14},    Can you solve the mystery?	

	'''
	from glob import glob
	from hs3_utils import Flight_segment 
	import matplotlib.pyplot as plt
	from libtools import shrink_ticks, epoch2iso, dtg2epoch, mkdir_p
	from numpy import arange, array, ceil, linspace
	from numpy.ma import mean
	from libcmap import rgbcmap
	from flight_namelists import experiments

	mkdir_p(dir_plot)

	#_get segments
	segments = experiments['HS3']

	#_convert bounding limits
	if start_dtg is not None:
		ep0 = dtg2epoch(start_dtg, full=True)
		ep1 = dtg2epoch(end_dtg, full=True)
	else:
		ep0 = '19000101000000'
		ep1 = '99999999999999'

	#_image setup
	smooth = False 
	label = 'nosmooth'
	files = glob('{0}/nc/*nc'.format(dir_cpl))
	files.sort()
	nrow = 6
	ncol = 2
	npag = ceil(len(files) / (nrow/ncol))
	calipso = rgbcmap('calipso')

	dtgs = segments.keys()
	dtgs.sort()
	#_loop over all cpl files (netcdf)
#	for i, fname in enumerate(files):
	for i, dtg in enumerate(dtgs):

		#_start new page
		if not i % (nrow*ncol):
			fig = plt.figure()
			fig.suptitle('CPL 532 Backscatter', size='small')
		#_calc column index

		#_calc row index

		ax_backscatter = fig.add_subplot(nrow, ncol, i%(nrow*ncol)+1)
		ax = ax_backscatter.twinx() 

		#_read in cpl data		
	#	cpl = read_cpl(f, None)
		cpl = Flight_segment(dtg=dtg)

		#_check if file is within limits
		if cpl.CPL_epoch.min() > ep1 and cpl.CPL_epoch.max() < ep0:
			continue

		#_get values of just aerosl
		aod = calc_aod(cpl, type_coa)
		bck	= calc_backscatter(cpl)
		if smooth:
			aod = smooth_aod(aod)

		#_generate list of times
		time = [epoch2iso(e) for e in cpl.CPL_epoch]

		#_get number of fovs
		nfov = aod.size
		nt = nfov / 2 

		ax.set_xticks(arange(aod.size)[nt:-nt:nt])
		ax.set_xticklabels(time[nt:-nt:nt])
		ax.xaxis.set_visible(False)

		props = dict(boxstyle='round', facecolor='white', alpha=.5)
	#	ax.text(nt/4., 3, f, color='k', size='x-small', bbox=props)
		f = '{0} - {1}'.format(time[0], time[-1]) 
		ax.text(nt/4., .75, f, color='k', size='xx-small', bbox=props)
		ax.set_ylim(0, 1)

		#_plotbackscatter
		cb=ax_backscatter.pcolormesh(bck,vmin=-4e-7,vmax=1e-4,cmap=calipso,
			zorder=0)
		ax.plot(aod, linewidth=0.2, zorder=1, color=(1,1,1,.5))

		if i % ncol:
		#	ax.set_ylabel('aod', size='xx-small')
			ax_backscatter.yaxis.set_visible(False)
		elif not i % ncol:
			ax.yaxis.set_visible(False)
			
		xlim=ax_backscatter.xaxis.get_data_interval()
		ax_backscatter.set_xlim(xlim)
		ax_backscatter.set_yticks(linspace(170,900,9))
		ax_backscatter.set_yticklabels(linspace(0,20,9)) 
		ax_backscatter.set_ylim(170,900)

		shrink_ticks(ax)
		shrink_ticks(ax_backscatter)

		#_I think this is bailing
		if (i % (nrow*ncol)) == ((nrow*ncol) - 1):
			page_num = i / (nrow*ncol)
			pname = 'page_{1}_{0:02d}_{2}OD_segments.png'.format(page_num,
				label, type_coa)
			pname = os.path.join(dir_plot, pname)
			print pname

			#_make a common y label
			fig.text(0.04, 0.5, 'altitude (km)', va='center',
				rotation='vertical', size='x-small')
			fig.text(0.96, 0.5, '{0}OD'.format(type_coa.upper()), va='center',
				rotation=270, size='x-small')

			fig.savefig(pname)
			plt.close()

	#_after everything is done	
	else:
		page_num = i / (nrow*ncol)
		pname = 'page_{1}_{0:02d}_{2}OD_segments.png'.format(page_num,
			label, type_coa)
		pname = os.path.join(dir_plot, pname)
		print pname

		#_make a common y label
		fig.text(0.04, 0.5, 'altitude (km)', va='center',
			rotation='vertical', size='x-small')
		fig.text(0.96, 0.5, 'AOD', va='center',
			rotation=270, size='x-small')

		fig.savefig(pname)
		plt.close()


def smooth_aod(a):
	from numpy import array, mean
	from libfilt import onethreeone as oto
##	a = array([ mean(a[max([n-1,0]):n+2]) for n in range(a.size) ])
##	return a
	return oto(a)


def calc_aod(c, type_coa):
	from numpy import vstack

	#_find locations with PBL/AEROSOL
	tau = vstack(c.CPL_tau_532)
	typ = vstack(c.CPL_tau_type)

	#_find where not missing
	pos = tau >= 0
	if type_coa == 'a':

		#_for aerosols:
		pbl = tau * (typ == 1)	* pos	#_PBL
		ele = tau * (typ == 2)	* pos	#_AEROSOL

		#_get column optical depths
		pbl = pbl.sum(axis=1)
		ele = ele.sum(axis=1)
		tot = pbl + ele

	elif type_coa == 'c':
		tot = tau * (typ == 0)	* pos	#_CLD
		
	else:
		raise RuntimeError, 'Invalid optical depth type. c | a.'

	return tot


def calc_backscatter(c):
	from numpy import vstack
	return vstack(c.CPL_ext_532).T


if __name__ == '__main__':
	main(**namelist)
