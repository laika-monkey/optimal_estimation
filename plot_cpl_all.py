#!/usr/bin/env python
#
#    NAME:   plot_all_cpl.py
#  AUTHOR:   Walter R. Sessions
#    DATE:   2015.Early
# PURPOSE:   Create multi-panel image showing backscatter with AOD overlaid 
#            from CPL data.

import os, sys

namelist = {
	#_if a subset is desired
	'start_dtg'	: '20130820000000',
	'end_dtg'   : '20130902000000',

	'dir_plot'	: os.path.join(os.environ['PRODUCTS'], 'plots', 'cpl'),
}

def main(dir_cpl='{0}/cpl'.format(os.environ['PRODUCTS']), start_dtg=None,
	dir_plot='.', end_dtg=None, **kwargs):
	'''
	start_dtg   str{14},    Starting date-time-group for subsetting
	end_dtg     str{14},    Can you solve the mystery?	

	'''
	from glob import glob
	from hs3_utils import read_cpl
	import matplotlib.pyplot as plt
	from libtools import shrink_ticks, epoch2iso, dtg2epoch, mkdir_p
	from numpy import arange, array, ceil, linspace
	from numpy.ma import mean
	from libcmap import rgbcmap

	mkdir_p(dir_plot)

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
	npag = ceil(len(files) / nrow)
	calipso = rgbcmap('calipso')

	#_loop over all cpl files (netcdf)
	i = 0
	for q, fname in enumerate(files):

		f = fname.split('/')[-1]

		#_read in cpl data		
		cpl = read_cpl(f, None)

		#_check if file is within limits
		if cpl.epoch.min() > ep1 or cpl.epoch.max() < ep0:
			continue

		#_start new page
		if not i % nrow:
			fig = plt.figure()

		ax_backscatter = fig.add_subplot(nrow, 1, i%nrow+1)
		ax = ax_backscatter.twinx() 

		print 'MADE IT'

		#_get values of just aerosl
		aod = calc_aod(cpl)
		bck	= calc_backscatter(cpl)
		if smooth:
			aod = smooth_aod(aod)

		#_generate list of times
		time = [epoch2iso(e) for e in cpl.epoch]

		#_get number of fovs
		nfov = aod.size
		nt = nfov / 5

		ax.set_xticks(arange(aod.size)[nt:-nt:nt])
		ax.set_xticklabels(time[nt:-nt:nt])

		props = dict(boxstyle='round', facecolor='white', alpha=.5)
		ax.text(nt/4., 3, f, color='k', size='x-small', bbox=props)
		ax.set_ylim(0, 4)

		#_plotbackscatter
		cb=ax_backscatter.pcolormesh(bck,vmin=-4e-7,vmax=1e-4,cmap=calipso,
			zorder=0)
		ax.plot(aod, linewidth=0.2, zorder=1, color=(1,1,1,.5))
		ax.set_ylabel('aod', size='xx-small')
		xlim=ax_backscatter.xaxis.get_data_interval()
		ax_backscatter.set_xlim(xlim)
	##	ax_backscatter.set_ylim(0,900)
		ax_backscatter.set_yticks(linspace(170,900,11))
		ax_backscatter.set_yticklabels(linspace(0,20,11)) 
		ax_backscatter.set_ylim(170,900)

		shrink_ticks(ax)
		shrink_ticks(ax_backscatter)

		i += 1
		if (i % nrow) == (nrow - 1):
			page_num = i / nrow
			pname = 'page_{1}_{0:02d}_0.png'.format(page_num, label)
			pname = os.path.join(dir_plot, pname)
			print pname

			#_make a common y label
			fig.text(0.04, 0.5, 'altitude (km)', va='center',
				rotation='vertical', size='x-small')

			fig.savefig(pname)
			plt.close()
	
	else:
		page_num = i / nrow
		pname = 'page_{1}_{0:02d}_0.png'.format(page_num, label)
		pname = os.path.join(dir_plot, pname)
		print pname

		#_make a common y label
		fig.text(0.04, 0.5, 'altitude (km)', va='center',
			rotation='vertical', size='x-small')

		fig.savefig(pname)
		plt.close()


def smooth_aod(a):
	from numpy import array, mean
	from libfilt import onethreeone as oto
##	a = array([ mean(a[max([n-1,0]):n+2]) for n in range(a.size) ])
##	return a
	return oto(a)

def calc_aod(c):
	from numpy import vstack

	#_find locations with PBL/AEROSOL
	tau = vstack(c.tau_532)
	typ = vstack(c.tau_type)

	#_find where not missing
	pos = tau >= 0
	pbl = tau * (typ == 1)	* pos	#_PBL
	ele = tau * (typ == 2)	* pos	#_AEROSOL

	#_get column optical depths
	pbl = pbl.sum(axis=1)
	ele = ele.sum(axis=1)
	tot = pbl + ele

	return tot


def calc_backscatter(c):
	from numpy import vstack
	return vstack(c.ext_532).T


if __name__ == '__main__':
	main(**namelist)
