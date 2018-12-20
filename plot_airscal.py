#!/usr/bin/env python
# AUTHOR	W R Sessions
# DATE		2016.04.28
# PURPOSE	Look at some of this shit.

import os
import sys
from libtools import dtg2epoch, dbg 

if 'DISPLAY' not in os.environ:
	import matplotlib
	matplotlib.use('Agg')
	 
namelist = {
	'out_label' : 'NVACLIMO1_tau-ref_M2_DV28_CLD_NOCTP_ATP-NAAPS_3KM',

	#_start and end times for period
	'e0'	: dtg2epoch('2013082400'),
	'e1'	: dtg2epoch('2013082500'),

	#_location of collocation files
	'dir_colloc'	: os.path.join(os.environ['PRODUCTS'], 'colloc'),

	#_name of experiment used in determining lblrtm location
	'dir_lblrtm_fmt'	: os.path.join(os.environ['PRODUCTS'],
								'LBL-RTM','{0}_{1}'),
	'experiment'		: 'AIRSCAL',

	#_bounding box of interest
	'latlon'	: [-30, 30, -90, 20],


	#_What to put into comparison
	'plot_oe'		: True,
	'plot_caliop'	: True,
	'plot_map'		: True, 
	}


################################################################################
################################################################################
################################################################################


def run_main(**kwargs):

	#_read fields we want to plot into dictionary
	data = read_period(**kwargs)
	dbg(('DONE READING'))

	#_plot time frame
	#_loop over each thing that is more than ~an hour apart
	windows = get_windows(data)
	#for e0, e1 in windows:
##	for e0, e1 in [[0, 9e9]]:
		kwargs.update({'e0' : e0, 'e1' : e1})	
		plot_period(data, **kwargs)

		#_remove what we just used up
		idx = (data.epoch >= e0) * (data.epoch <= e1)
		data = data[idx == False]
		dbg('done for one')


def get_windows(d):
	from numpy import diff, arange
	from libtools import epoch2iso as e2i

	#_pull out just optimal estimation times and put them in order
	oe = d[d.name == 'oe']
	oe = oe[oe.epoch.argsort()]
	e = oe.epoch

	#_make a list of where more than an hour passes between retrievals
	idx = arange(e.size - 1)[diff(e) > 3600.]
	es = []
	time_str = e[0] 
	for ii, (x, y) in enumerate(zip(e[idx], e[idx+1])):
		time_end = y+1
		es.append([time_str, time_end])
		time_str = y+2

	es.append([time_str, e[-1]])

	for x, y in es:
		dbg(e2i((x,y)))

	return es 


def read_period(plot_oe=False, plot_caliop=False, **kwargs):
	'''
	plot CALIOP and retrieval for time frame
	return record array with
	dtype = data, epoch, lat, lon, source
	'''
	from glob import glob
	from oe_utils import read_airs, read_caliop, get_airs_times, read_oe
	from oe_utils import get_collocation_indices, precheck
	from pyhdf.SD import SD, SDC
	from libtools import epoch2dtg as e2d
	from numpy import array, append, vstack, recarray, ndarray

	#_pull out some vars
	e0 = kwargs.get('e0', 0)
	e1 = kwargs.get('e1', 9e9)
	latlon = [-90,90,-180,180]

	#_read in previously generated file list
	files_colloc, pre_create = precheck(**kwargs)

	#_create recarray for these things.
	dtype = [('values', ndarray), ('lat','f4'), ('lon','f4'), ('epoch', 'f8'),
			('name','a20')]
	data = recarray((0,), dtype)

	#_loop over collocation times and read in crap
	all_oe = []
	for ii, file_colloc in enumerate(files_colloc):
		#_kklluuddggee
		file_colloc = file_colloc.split('/')[-1]
		file_colloc = os.path.join(os.environ['PRODUCTS'], 'colloc', file_colloc)
		dbg(file_colloc)

		#_get name of AIRS file
		hdf = SD(file_colloc, SDC.READ)
		file_airs = hdf.attributes()['fname_AIRS'] 
		file_caliop	= hdf.attributes()['fname_CALIOP']	

		#_check if ANYTHING is within the period
		e0_airs, e1_airs = get_airs_times(file_airs)
		
		#_if completely outside window, skip
		if (e1_airs < e0) or (e0_airs > e1):
			continue

		#_get indices
		idx_x, idx_l, idx_c = get_collocation_indices(file_colloc)
	
		#_read in airs and caliop fields
		#_there's no time data in colloc file, so use AIRS to filter
		airs, wave = read_airs(file_airs, idx_x, idx_l, **kwargs)

		#_check airs times
		idx = (airs.AIRS_epoch >= e0) * (airs.AIRS_epoch < e1)
		lat0, lat1, lon0, lon1 = latlon
		lidx = (airs.AIRS_latitude >= lat0) * (airs.AIRS_latitude < lat1)
		nidx = (airs.AIRS_longitude >= lon0) * (airs.AIRS_longitude < lon1)
		idx = idx * lidx * nidx	

		#_skip if none match
		if idx.sum() == 0:
			continue

		#_keep track of which files for future runs
		if pre_create:
			pre_create.write(file_colloc + '\n')	

		#_limit to 4d box
		airs = airs[idx]
		idx_c = idx_c[idx]
		idx_x = idx_x[idx]
		idx_l = idx_l[idx]

		#_read oe output
		if plot_oe:
			dbg('READING AIRS')
			oe = read_oe(airs.AIRS_epoch, **kwargs)
			new = len(oe)

			#_add to recarray
			data.resize(data.size + new)
			data.epoch[-new:] = airs.AIRS_epoch[:]
			data.lat[-new:] = airs.AIRS_latitude[:]
			data.lon[-new:] = airs.AIRS_longitude[:]
			data.epoch[-new:] = airs.AIRS_epoch[:]
			data.name[-new:] = 'oe'	
			for nnn, v in enumerate(oe):
				data.values[-new+nnn] = v #oe[nnn]

		#_do caliop
		if plot_caliop:
			idx_c = idx_c[:,0][:,None]
			caliop = read_caliop(file_caliop, idx_c, **kwargs)

			#_flatten caliop and drop missing fields
			mask = idx_c.flatten() != -1
			backscatter = vstack(caliop.CALIOP_total_attenuated_backscatter_532)[mask]
			longitude = vstack(caliop.CALIOP_longitude).flatten()[mask]
			epoch = vstack(caliop.CALIOP_epoch).flatten()[mask]
			latitude = vstack(caliop.CALIOP_latitude).flatten()[mask]
			
			new = mask.sum() #caliop.size
			data.resize(data.size + new)
			data.epoch[-new:] = epoch[:] 
			data.lat[-new:] = latitude[:] 
			data.lon[-new:] = longitude[:] 
			data.name[-new:] = 'caliop'	
			for nnn in range(new):
				data.values[-new+nnn] = backscatter[nnn]

		##	#_flatten caliop and drop missing fields
		##	mask = idx_c != -1
		##	longitude = vstack(caliop.CALIOP_longitude).flatten()[mask]
		##	latitude = vstack(caliop.CALIOP_latitude).flatten()[mask]
		##	epoch = vstack(caliop.CALIOP_epoch).flatten()[mask]
		##	backscatter = vstack(caliop.CALIOP_total_attenuated_backscatter_532)[mask]
		##	
		##	new = idx_c.size
		##	data.resize(data.size + new)
		##	data.epoch[-new:] = epoch[:] 
		##	data.lat[-new:] = latitude[:] 
		##	data.lon[-new:] = longitude[:] 
		##	data.name[-new:] = 'caliop'	
		##	for nnn in range(new):
		##		data.values[-new+nnn] = backscatter[nnn]

	#_close precreated file
	if pre_create:
		pre_create.close()

	#_send recarray back to sender
	return data


def plot_period(data, e0, e1, plot_oe=True, plot_caliop=False,
	cmfile='/home/wsessions/lib/cmaps/calipso-backscatter.cmap', **kwargs):
	''' plot all the crep in data '''
	#_read in retrieval data for period
	##	obs = read_oe(airs.AIRS_epoch, **kwargs)
	import matplotlib.pyplot as plt
	from libcmap import rgbcmap
	from libtools import epoch2iso as e2i
	from numpy import array, vstack, arange, recarray, append
	from libtools import shrink_ticks
	from numpy.ma import masked_where as mw
	import ccplot_cmap as ccc

	nplot = sum([plot_oe, plot_caliop])
	fig, ax = plt.subplots(nplot)

	try:
		ax_oe, ax_cal = ax
	except:
		ax_oe = ax

	dbg(e2i([e0, e1]))

	#_plot retrievals from optimal estimation
	#_AIRS_____#
	def fix_oe_crap(d, **kwargs):
		''' get the oe values out of stupid dictionaries '''
		example = array([f is not None for f in d.values]) 
		example = d[example.argmin()]
		ssp_dbs = example.values.ssp_db
		names = [n.split('.')[-2] for n in ssp_dbs]

		#_build recarray to store these
		nssp = len(ssp_dbs)
		nt = d.size
		dtype = [('tau','f4'), ('lat','f4'), ('lon','f4'), ('epoch', 'f8'),
				('name','a20'), ('habit', 'a20'), ('ref_wvnum', 'f4'),
				('z_top', 'f4'), ('z', 'f4'), ('ref', 'f4')]  
		data = recarray((nssp*nt), dtype)

		#_loop over layers and pull all this crap out
		for i, vals in enumerate(d.values):
			for ii, l in enumerate(vals.layers):
				idx = ii + nssp*i

				habit = vals.ssp_db[l['dbnum']].split('.')[-2]
				arg = (l['tau'][0], d.lat[i], d.lon[i], d.epoch[i], d.name[i],
					habit, l['ref_wvnum'], l['z_top'], l['z'], l['ref'])
				data[idx] = arg

		data = data[data.tau != -9999]		
		return data

	#_pull out oe
	oe = data[data.name == 'oe']
	idx = array([d is not None for d in oe.values])
	oe = oe[idx]
	oe = oe[oe.epoch.argsort()]
	oe = oe[(oe.epoch >= e0) * (oe.epoch <= e1)]

	#_pull out optical depths and put into dictionary by specie
	oe = fix_oe_crap(oe, **kwargs)
	ssps = set(oe.habit)
	for habit in set(oe.habit):
		if habit == 'mie_wat':
			continue
	
		#_pull out this ssp habit
		oe_habit = oe[oe.habit == habit]
		oe_habit = oe_habit[oe_habit.epoch.argsort()]

		#_
		x = oe_habit.epoch
		y = oe_habit.tau
		nt = x.size
		max_oe_epoch, min_oe_epoch = x.max(), x.min()
		
		ax_oe.plot(x, y, label=habit, linewidth=0.5)
		xticks = append(x[::nt/5], x[-1])
		if xticks[-1] - xticks[-2] < nt/20:
			xticks[-2] = xticks[-1]
		xticklabels = [tmp[11:19] for tmp in e2i(xticks)]
		ax_oe.set_xticks(xticks)
		ax_oe.set_xticklabels(xticklabels)
		shrink_ticks(ax_oe)

		#_crop plotting area.
		ax_oe.set_xlim(x.min(), x.max())
		ax_oe.set_ylim(0, 3.)

	#_drop box, make smaller
	ax_oe.legend()

	##########
	#_CALIOP_#
	if plot_caliop:
		import numpy as np
	
		#_ pull out caliop data
		caliop = data[data.name == 'caliop']
		idx_c = (caliop.epoch >= e0) * (caliop.epoch <= e1)
		caliop = caliop[(caliop.epoch >= min_oe_epoch) * (caliop.epoch <= max_oe_epoch)]
	##	caliop = caliop[(caliop.epoch >= e0) * (caliop.epoch <= e1)]
		caliop = caliop[caliop.epoch.argsort()]
		x = caliop.epoch

		#_put into ONE BIG OL 2d array 
		caliop = vstack(caliop.values)
		caliop = mw(caliop == -9999, caliop)
		calipo = np.ma.masked_invalid(caliop)
	##	caliop[caliop < 0] = 0.
	##	caliop[caliop > 1.5] = 1.5

		#_coords
		y = arange(caliop.shape[1])

		dbg(('REMOVE TRUNCATION'))
		#_load up backscatter colormaps
		cmap, norm, ticks = ccc.loadcolormap(cmfile, 'caliop')
	##	cb = ax_cal.imshow(caliop.T, cmap=cmap, norm=norm, interpolation='nearest')
		im = ax_cal.pcolormesh(x, y, caliop.T, cmap=cmap, norm=norm)#, vmin=0, vmax=1.5)
		ax_cal.set_xlabel('{0:e}, {1:e}'.format(caliop.max(), caliop.min()))

		#_label bottom
		xticks = append(x[::x.size/5], x[-1])
		if xticks[-1] - xticks[-2] < 100:
			xticks[-2] = xticks[-1]
		xticklabels = [tmp[11:19] for tmp in e2i(xticks)]

		ax_cal.set_xticks(xticks)
		ax_cal.set_xticklabels(xticklabels)
		ax_cal.invert_yaxis()
		ax_cal.get_yaxis().set_visible(False)
		ax_cal.set_xlim(x.min(), x.max())
		shrink_ticks(ax_cal)

	##	cb = fig.colorbar(im, ax=axes, cax=cbaxes, orientation="vertical",
    ##                  extend="both", ticks=ticks, norm=norm,
    ##                  format=SciFormatter())

	##	cb.set_label(name)

	pname = 'oe_{0}_{1}.png'.format(e2i(e0), e2i(e1))
	pname = pname.replace(':', '-')
	dbg(pname)
	plt.savefig(pname)


################################################################################
################################################################################
################################################################################


if __name__ == '__main__':
	run_main(**namelist)
