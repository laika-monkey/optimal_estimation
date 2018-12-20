#!/usr/bin/env python
#############################################################################_80
#	HS-OE_fig_01.py
#	Create figures for hypserspectral optimal estimation paper
################################################################################
################################################################################


import re
import os
import sys
from optimal_estimation import *
from numpy import arange


# DIR_LOG  = os.path.expanduser('~/qsub_logs/')
DIR_LOG  = os.environ['LOG']
DIR_PROD = os.environ['PRODUCTS']
DIR_TMP  = '/data/wsessions/TMP'
DIR_NRL  = os.path.join(DIR_PROD, 'NRL', 'NVA_CLIMO1misr')
DIR_PLOT = os.path.join(DIR_PROD, 'plots')
DIR_SHIS = os.path.join(DIR_PROD, 'hs3')
DEBUG    = 1


#_used instead of config_file
namelist = {

	'fnames' : ['SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.nc',
				'SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.NOISE.0004485.nc', 
				'SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.NOISE.0021783.nc',
				'SHIS.CPL.GDAS.COLLOC.SIMULATED.tau_ref_z.NOISE.0030252.nc'],
	'vnames'	: [	'ref', 'z', 'thickness' , 'surf_emis', 
					'instrument', 'surf_temp'],

	#_ok, see if this breaks everything
	'state_vars'	: ['tau'],
	'uncert_vars'	: ['ref','z','thickness','surf_emis','surf_temp'],#'rel_hum'

	#_location of input radiances and output directory for fov plots
#	'dir_out'		: os.path.join(DIR_PROD, experiment.lower()), 

	#_directory of simulated run static case
	'dir_lblrtm'	: '/data/wsessions/LBL-RTM/test_fov',

	#_top directory to dump plots
	'dir_plot'		: os.path.join(DIR_PROD, 'plots'),
	
	#_directory with campaign data
	'dir_hs3'		: os.path.join(os.environ['PRODUCTS'],'hs3'),

#	'delta_ts'		: arange(-5, 5.1, .25),
	'delta_ts'		: arange(11) - 5, 
	'by_var'		: False, #_False, z, ref, tau

	#_threshold for optical depth variance
	'threshold_aod'	: 0.1, 

	#_which naaps species to plot
	'species'	: ['dust_aod', 'smoke_aod', 'total_aod'],
	'region'	: 'hs3_2013c',	
	}


################################################################################
################################################################################
################################################################################


def run_main(vnames=None, **kwargs):
	''' plot all figures '''

	from flight_namelists import experiments
	from hs3_utils import Flight_segment as F
	from libtools import epoch2dtg, newdtg

	#_plot only ones during flight times
	fig_06(**kwargs)


################################################################################
################################################################################
################################################################################



def fig_06(region='hs3_2013', **kwargs):
	'''
	plot NAAPS data for HS3 2013 field campaign from climatology

	start July 2013 - Oct 1, 2013, maximum values for all species,
						and then just smoke and dust
	SAME AS FIG_05 but only for flight times
	'''
	from libnva import read_naapsaod as read_naaps
	from libnva import sub_region, draw_map
	from libtools import newdtg, shrink_ticks
	from libmeta import fields, regions
	import matplotlib.pyplot as plt
	from numpy import append, nan, meshgrid, isnan
	from libtools import rgbgen, epoch2iso, dtg2iso
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	from numpy.ma import masked_where

	#_pull out region data
	reg_meta = regions()[region]

	#_read in all flight track data for 2013
##	track = read_tracks() 
	track = read_raw_tracks()

	#_initialize draw space
	fig = plt.figure()
	ax = fig.add_subplot(111)

	#_plot around tropical atlantic
	map_args = {
		'corners' : reg_meta['corn'], 
        'grid'    : reg_meta['grid'], 
        'delta'   : reg_meta['delta'],
		'laby'    : [0,1,0,0],
		}
	
	#_get out contour levels
	m = draw_map(ax=ax, **map_args)
	m.plot(track.lon, track.lat, 'k-', linewidth=0.3, zorder=7)

	#_labet_
	ax.set_title('HS3 2013 Flight Tracks', size='small')
	shrink_ticks(ax)

	#_save
	pname = 'flight_tracks.png'.format(region)
	pname = os.path.join(DIR_PLOT, pname)
	print pname
	plt.savefig(pname)


def read_raw_tracks(**kwargs):
	
	from glob import glob
	import re
	from netCDF4 import Dataset as D
	from numpy import recarray, append, array, diff
	lat, lon = array([]), array([])

	for fcpl in glob('/data/wsessions/cpl/nc/OP_13*'):
		res = re.search('(\d{5})\w_(\d{6})_(\d{6})_(\d{6}).nc', fcpl)	
		dtg = res.group(2)
		if dtg < '130816':
			continue #_skip transit flights

		print fcpl

		#_get lat, lon, time
		cdf = D(fcpl, 'r')
		nt = len(cdf.dimensions['time'])
		lat = append(lat, cdf.variables['Lat'][:])
		lon = append(lon, cdf.variables['Lon'][:])

	latdx = abs(diff(lat)) > 20 
	londx = abs(diff(lon)) > 20 
	idx = append(latdx[:,None], londx[:, None], 1).any(1) 
	idx = append(idx, True)
	dtype = [('lat','f8'), ('lon','f8'), ('epoch','f8')]
	data = recarray((len(lat),), dtype)
	data.lat[:] = lat
	data.lon[:] = lon
	data.lat[idx] = None
	data.lon[idx] = None
	return data


def read_tracks(fname='/data/wsessions/hs3/hs3_2013_flighttracks.dat',**kwargs):
    from numpy import loadtxt, recarray, append, diff
    import matplotlib.pyplot as plt
    from libtools import dtg2epoch

    track = loadtxt(fname, delimiter=',')
    dtype = [('epoch', 'f8'), ('lat', 'f8'), ('lon', 'f8')]
    data = recarray((track.shape[0],), dtype=dtype)

    data.epoch = track[:,0]
    data.lat = track[:,1]
    data.lon = track[:,2]

    #_add nones to break up flights
    widx = data.lon > -100
    eidx = data.epoch > dtg2epoch('2013081000')

    idx = widx * eidx

    final = recarray((idx.sum(),), dtype=dtype)
    final.epoch = data.epoch[idx]
    final.lat = data.lat[idx]
    final.lon = data.lon[idx]

    idx = diff(final.epoch) > 10
#   idx = append(idx, [False])
    idx = append([False], idx)
    final.epoch[idx] = None
    final.lat[idx] = None
    final.lon[idx] = None

    return final


if __name__ == '__main__':
	run_main(**namelist)
