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
	'region'	: 'hs3_2013b',	
	}


################################################################################
################################################################################
################################################################################


def run_main(vnames=None, **kwargs):
	''' plot all figures '''

	#_plot spectral signature of habits, (currently kaolinite, gypsum)
##	fig_00(**kwargs)

	#_plot uncertainty by variable
	if 0:
	  for fname in kwargs.get('fnames'):
		for current_var in vnames:
			namelist = kwargs.copy()
			label = 'ERR_{0}'.format(current_var)
			label = re.sub('_$', '', label)
			label = label.replace('.','_')
			namelist.update({'out_label' : label})
			namelist.update({'fsim' : fname})
			[fig_01(v, **namelist) for v in ['ref','z']]

	  for sv in vnames:
		label = 'ERR_total'
		namelist.update({'out_label' : label})
		fig_01(sv, files=fnames, **namelist)


##	#_plot total uncertainty

	#_plot surf_temp variance
#	fig_02(**kwargs)

	#_plot surf_temp whisker plots
#	fig_03(**kwargs)

	#_plot matrices
#	fig_04(**kwargs)

	#_PLOT NAAPS SMOKE/DUST
	#_do for each flight then for entire period
	from flight_namelists import experiments
	from hs3_utils import Flight_segment as F
	from libtools import epoch2dtg, newdtg
	dtgs = []
	if 1:
 	  for dtg, values in experiments['HS3'].iteritems():

		#_read in flight file to get start and end times
		dbg(dtg)
		flight = F(dtg=dtg)
		dtg0 = '{0}00'.format(epoch2dtg(flight.SHIS_epoch.min())[:8])
		dtg1 = '{0}18'.format(epoch2dtg(flight.SHIS_epoch.max())[:8])
		dtg1 = newdtg(dtg1, 6) 

		dtg = dtg0
		while dtg <= dtg1:
			dtgs.append(dtg)
			dtg = newdtg(dtg, 6)

		#_individual flights
	#	fig_05(dtg0=dtg0, dtg1=dtg1, **kwargs)
	#_all times (USE FIG_06)
##	fig_05(**kwargs)

	#_plot only ones during flight times
	dtgs = list(set(dtgs))
	fig_06(dtgs, **kwargs)

	#_plot habit comparison of KAOL and GYP
##	fig_07(**kwargs)

	#_spectral comparison of KAOL and GYP
##	fig_08(**kwargs)	

	#_dashboard light
##	fig_09(**kwargs)



################################################################################
################################################################################
################################################################################


def fig_00(pname='delta_bt.png', dir_plot='.', dv=-27, **kwargs):
	'''
	Spectral signatures of habits by wavenumber

	I believe the input files I used to make these were done manually,
	so if you want to add things like "cloud" or other habits, you'll
	need to rerun the test cases in std_tropical with some new inputs.

	x == wavenumber
	y == change in radiance from clear sky?
	'''
	import matplotlib.pyplot as plt
	import os
	from lblrtm_utils import read_lbldis_out as read_lbldis
	from libtools import shrink_ticks
	from libgeo import planck_inv
	from lblrtm_utils import read_microwindows_header as mw
	from numpy import abs

	dir_in = '/data/wsessions/LBL-RTM_simulated/std_tropical'
	habits = {	#'quartz'	: 'blue', 
				'kaolinite' : 'red', 
				'gypsum'	: 'green', 
				'cloud'		: 'blue',
			}

	fig, axis = plt.subplots(1)

	lines = []
	names = []
	#_read in all data and plot
	for h, c in habits.iteritems():

		#_gen filename
		fname = os.path.join(dir_in, '{0}.cdf'.format(h))
		fn_in = os.path.join(dir_in, 'lbldis_input.{0}'.format(h))
	
		#_read in output file
		data = read_lbldis(fname, fname_input=fn_in)	

		#_convert to brightness temperatures
		wavs = data[0].wavenumber
		bt0 = planck_inv(data[0].radiances, wavs*100, domain='wavenumber') 
		bt1 = planck_inv(data[1].radiances, wavs*100, domain='wavenumber') 
		d_bt = bt0 - bt1

		line = axis.plot(wavs, d_bt, c=c, linewidth=0.4,
			label='changes by {0}'.format(h), zorder=1)
	
	#_plot clear sky bt
	atwn = axis.twinx()
	line = atwn.plot(wavs, bt0, 'k-', linewidth=0.2, label='clear sky')

	#_plot microwindow numbers
	chans = mw()[abs(dv+1)]
	for n, x in chans:
		v = (n + x) / 2.
		axis.plot([v, v], [0, 20], 'k--', linewidth=0.2)

	axis.set_xlim(600, 1300)
	axis.set_ylim(0,20)

	axis.set_ylabel("$\Delta B_{\\nu}T$ $(K)$", size='small')
	axis.set_xlabel("wavenumber $(cm^{-1})$", size='small')
	atwn.set_ylabel("$B_{\\nu}T$ $(K)$", size='small', rotation=-90)

	axis.legend(loc=2, frameon=False, fontsize=8)
	atwn.legend(loc=1, frameon=False, fontsize=8)

	#_save
	[shrink_ticks(a) for a in [axis, atwn]]
	pname = os.path.join(dir_plot, pname)
	dbg(pname)
	plt.savefig(pname)


#def plot_simulated_retrieval_by_var(var, out_label='', plot_var='tau',
def fig_01(var, out_label='', plot_var='tau',
	fsim=None, dir_hs3=os.path.join(os.environ['PRODUCTS'],'hs3'),
	dir_plot='.', files=[], **kwargs):
	'''
	produce some plots for this god damn stuff

	out_label	str,	experiment is not used so that it can
						denote specific input files as opposed 
						to output.  in_label should have been used.
	'''
	import matplotlib.pyplot as plt
	from libtools import mkdir_p, shrink_ticks
	from numpy import array, append
	import re
	from oe_sim import read_simulated_retrieval

	if len(files) == 0:
		print type(fsim)
		res = re.search('SIMULATED.(.*?).nc', fsim)
		sim_experiment = res.group(1)

		#_read in retrieval data
		fname = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.{1}.retrieved'.format(
			sim_experiment, out_label)
		fname = os.path.join(dir_hs3, fname)
		try:
			truth, retr, uncrt = read_simulated_retrieval(fname, **kwargs)
		except IOError:
			dbg(('Did you complete this retrieval yet?', out_label))
			os._exit(23)
	else:
		for i, f in enumerate(files):
			res = re.search('SIMULATED.(.*?).nc', f)
			sim_experiment = res.group(1)
			f = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.{1}.retrieved'.format(
				sim_experiment, out_label)
			f = os.path.join(dir_hs3, f)
			try:
				if i == 0:	#_initialize
					truth, retr, uncrt = read_simulated_retrieval(f, **kwargs)
				else:		#_append
					t, r, u = read_simulated_retrieval(f, **kwargs)
					for sv in truth:
						truth[sv] = append(truth[sv], t[sv]) 
					for sv in retr:
						retr[sv] = append(retr[sv], r[sv]) 
					for sv in uncrt:
						uncrt[sv] = append(uncrt[sv], u[sv]) 
			except IOError:
				dbg(('Did you complete this retrieval yet?', f))
				os._exit(23)

		res = re.search('SIMULATED.(.*?).nc', fsim)
		sim_experiment = res.group(1).split('.')[0]
		out_label = 'summary'

	#_what are we working with?
	state_vars = retr.keys() # truth.keys()
	dbg('');dbg(state_vars)

	#_kick back if not available
	if plot_var not in state_vars:
		dbg(('Variable not retrieved:', plot_var))
		return

	#_color dictionary
	order = ['low', 'mid', 'high', 'very_high']
	ranges = {
		'low'	: {
			'color'	: 'blue',
			'tau'	: (-9999., 0.5),
			'ref'	: (0.5, 0.675),
		},
		'mid'	: {
			'color'	: 'green',
			'tau'	: (0.5, 1.0),
			'ref'	: (0.625, 0.75),
		},
		'high'	: {
			'color'	: 'orange',
			'tau'	: (1.0, 1.5),
			'ref'	: (0.75, 0.875),
		},
		'very_high'	: {
			'color'	: 'red',
			'tau'	: (1.5, 9999),
			'ref'	: (0.875, 1.0),
		},
		}

	#_output directory for plots
	ppath = '{1}.{0}'.format(out_label, sim_experiment)
	ppath = os.path.join(dir_plot, ppath)
	if not os.path.exists(ppath):
		mkdir_p(ppath)

	#_get potential settings for value
	values = list(set(truth[var]))
	nval = len(values)
	values.sort()

	''' make loop over pages, 6 per page '''
	nrow = 6
	npag = nval / nrow #_do we need this? 
	page = 0

	#_produce for each state variable
	for sv in state_vars:

		for j, value in enumerate(values):
			dbg((sv, var, value))

			#_initialize new figure
			if not j % nrow:
				page += 1
				fig = plt.figure()
				k = 0
				fig.suptitle('Fractional Uncertainty')

			#_get indices of limit
			idx = truth[var] == value

			#_get residual
			res = truth[sv][idx] - retr[sv][idx]
	
			#_initialize figure and axes
			axis = fig.add_subplot(nrow, 1, k+1)
	
			ymax, xmax, xmin = 0, 0, 0
			nbin = 20
	
			#_plot uncertainty by truth
			fract = uncrt[sv][idx] / truth[sv][idx]
			x_truth, y_fract = [], []

			#_pull out values across all over vars, at truth value
			for value_truth in set(truth[sv]):
				x_truth.append(value_truth)
				y_fract.append(fract[truth[sv][idx] == value_truth])	

			#_pull out number of retrievals going into these
			y_bar = [len(y) for y in y_fract]

			##	axis.scatter(truth[sv][idx], fract, marker='x', color=col)
			axis.boxplot(y_fract, positions=x_truth, widths=0.05)
			axis.set_xticks(list(set(truth[sv])))

			abar = axis.twinx()
			bar_width = (max(x_truth) - min(x_truth)) / float(len(x_truth))\
				 - 0.05
			
			abar.bar(x_truth-bar_width/2., y_bar, bar_width, 
				alpha=0.24, color='brown')	
			#_get how high we can go and still be 95 % certain of the AOD
		##	m95 = truth[sv][certainty_garbage(res)]
			
			#_various things for uncertainty
		#	axis.text(truth[sv].max()-0.2, 1.6, 'n={0}'.format(len(y_fract))
			axis.set_ylim(bottom=0, top=2)
			axis.set_xlim(truth[sv].min()-0.1, truth[sv].max()+0.1)
			axis.set_ylabel('{0} = {1}'.format(var, value), size='xx-small')
			axis.set_ylabel('{0} = {1}'.format(var, value), size='x-small')
##			axis.set_xlabel('{0}, max_95 {1:5.2f}'.format(sv,m95), size='x-small')
			axis.grid(False)
			
			#_shrink...ticks...
			shrink_ticks(axis)
			shrink_ticks(abar)

			#_save when page full	
			
			if k == (nrow-1) or j == nval-1: 
				#_create image name
				pname = 'hist_uncert.{0}.{1}.{2}.by_{3}.p{4:02d}.png'.format(
					sim_experiment, out_label, sv, var, page)
				pname = os.path.join(ppath, pname)
				plt.savefig(pname)
				plt.close()		
		
				dbg(pname)

			k += 1


def fig_01_check(var, out_label='', plot_var='tau', fsim=None, dir_lblrtm=None,
	dir_hs3=os.path.join(os.environ['PRODUCTS'],'hs3'),  **kwargs):
	'''
	Uncertainty by optical depth
	
	First does plot of all median uncertainties?  Or by specific ref/z

	Plot uncertainty broken down by ref/z with box/whisker plots.

	Plot them in total.

	plot_simulated_by_var
	produce some plots for this god damn stuff

	out_label	str,	experiment is not used so that it can
						denote specific input files as opposed 
						to output.  in_label should have been used.
	'''
	import matplotlib.pyplot as plt
	from libtools import mkdir_p, shrink_ticks
	from numpy import array, append

	dbg('running...')
	#_read in retrieval data
	if type(fsim) != list:
		res = re.search('SIMULATED.(.*?).nc', fsim)
		sim_experiment = res.group(1)
		fname = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.{1}.retrieved'.format(
			sim_experiment, out_label)
		fname = os.path.join(dir_hs3, fname)

		try:
			truth, retr, uncrt = read_simulated_retrieval(fname, **kwargs)
		except IOError:
			dbg(('Did you complete this retrieval yet?', out_label))
			os._exit(23)

	else:
		for fname in fsim:
			res = re.search('SIMULATED.(.*?).nc', fname)
			sim_experiment = res.group(1)
			fname = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.{1}.retrieved'.format(
					sim_experiment, out_label)
			fname = os.path.join(dir_hs3, fname)
			t, r, u = read_simulated_retrieval(fname, **kwargs)
			try:
				[append(truth[sv], t[sv]) for sv in t.keys()]
				[append(retr[sv], r[sv]) for sv in r.keys()]
				[append(uncrt[sv], u[sv]) for sv in u.keys()]
			except UnboundLocalError:
				truth, retr, uncrt = t, r, u

		sim_experiment = 'NOISE'

	#_what are we working with?
	state_vars = retr.keys() # truth.keys()

	#_kick back if not available
	if plot_var not in state_vars:
		dbg(('Variable not retrieved:', plot_var))
		return

	#_output directory for plots
	ppath = 'PLOTS_{0}'.format(out_label)
	if not os.path.exists(ppath):
		mkdir_p(ppath)

	#_get potential settings for value
	values = list(set(truth[var]))
	nval = len(values)
	values.sort()

	nrow = 6
	npag = nval / nrow #_do we need this? 
	page = 0

	#_produce for each state variable
	for sv in state_vars:

		for j, value in enumerate(values):
			dbg((sv, var, value))

			#_initialize new figure
			if not j % nrow:
				page += 1
				fig = plt.figure()
				k = 0

			#_get indices of limit
			idx = truth[var] == value

			#_get residual
			res = truth[sv][idx] - retr[sv][idx]
	
			#_initialize figure and axes
			axis = fig.add_subplot(nrow, 1, k+1)
	
			ymax, xmax, xmin = 0, 0, 0
			nbin = 20
	
			#_plot uncertainty by truth
			fract = uncrt[sv][idx] / truth[sv][idx]
			x_truth, y_fract = [], []
			for value_truth in set(truth[sv]):
				x_truth.append(value_truth)
				y_fract.append(fract[truth[sv][idx] == value_truth])	
	
			##	axis.scatter(truth[sv][idx], fract, marker='x', color=col)
			axis.boxplot(y_fract, positions=x_truth, widths=0.05)
	
			axis.set_xticks(list(set(truth[sv])))
	
			#_various things for uncertainty
			axis.set_ylim(bottom=0, top=2)
			axis.set_xlim(truth[sv].min()-0.1, truth[sv].max()+0.1)
			axis.set_ylabel('{0} = {1}'.format(var, value), size='x-small')
			axis.grid(False)
			
			#_shrink...ticks...
			shrink_ticks(axis)

			#_save when page full	
			if k == (nrow-1) or j == nval-1: 
				#_create image name
				pname = 'hist_uncert.{0}.{1}.{2}.by_{3}.p{4:02d}.png'.format(
					sim_experiment, out_label, sv, var, page)
				pname = os.path.join(ppath, pname)
				plt.savefig(pname)
				plt.close()		
		
				dbg(pname)

			k += 1


def fig_02(sv='tau', fnames=[], delta_ts=[0], dir_plot='.', dir_hs3='.',
	by_var=False, threshold_aod=0.1, **kwargs):
	'''
	plot histograms of all the surf temperature tests (alpha it up)

	delta_ts	array of departures from the true surface temp
	fnames		list,		names of input files (SIMULATED)
	'''
	from oe_sim import read_simulated_retrieval
	import matplotlib.pyplot as plt
	from numpy import array, append, linspace
	from libtools import shrink_ticks, colors_spag

	fig, axes = plt.subplots(2,1)
	nshift = len(delta_ts)

	axis, ax0 = axes

	#_loop over shift in T 
	labels = []
	for j, shift_temp in enumerate(delta_ts): 
		#_ignore non-shifted
		if shift_temp == 0:
			continue

		#_define out_label
		out_label = 'SFC_TEMP_{0:4.2f}'.format(shift_temp)

		#_loop over input files
		for i, fname in enumerate(fnames):
			lab = re.search('SIMULATED.(.*?).nc$', fname).group(1)
			lab = '{0}.{1}'.format(lab, out_label)
			fname = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.retrieved'.format(lab)
			fname = os.path.join(dir_hs3, fname)

			if not i: #_read in data
				truth, retr, uncrt = read_simulated_retrieval(fname, **kwargs)
			else: #_append data together
				t, r, u = read_simulated_retrieval(fname, **kwargs)
				for sv in truth:
					truth[sv] = append(truth[sv], t[sv])
				for sv in retr:
					retr[sv] = append(retr[sv], r[sv])
				for sv in uncrt:
					uncrt[sv] = append(uncrt[sv], u[sv]) 

		#_calculate residuals
		res = array(truth[sv] - retr[sv])

		#_setup plot
		args = {
			'bins' : 100, 
			'alpha' : 0.5, 
			'range' : (-1, 1),
			'stacked' : True,
			'edgecolor' : None,
			} 
			
		#_Make and image separating out heights
		#_make whisker plot

		#_collect data and probabilities 
		if by_var != False:
			vars = list(set(truth[by_var]))
			vars.sort()
			nvar = len(vars)
	
			#_only do this crap first time through
			for k, v in enumerate(vars):
				idx = truth[by_var] == v 
				#init data or append
				if j == 0 and k == 0:
					data = res[idx].reshape(idx.sum(),1)
				else:
					data = append(data, res[idx].reshape(idx.sum(),1), axis=1)

				#_calculate percentage within thresh?
				P = (abs(res[idx]) <= threshold_aod).sum() / float(idx.sum())
				P *= 100
				fmt = '$\Delta$T={0:>5.2f}, {4}={1:>4.2f},{2:>6.1f}%' +\
					'$\pm${3:>5.2f}'
				labels.append(fmt.format(shift_temp,v,P,threshold_aod,by_var))

		#_collect data and probabilities for total
		else:
			if j == 0:
				#_only do this crap first time through
				k = 0
				data = res.reshape(res.size,1)
			else:
				data = append(data, res.reshape(res.size,1), axis=1)
	
			#_calculate percentage within thresh?
			P = (abs(res) <= threshold_aod).sum() / float(res.size)
			P *= 100
			fmt = '$\Delta$T={0:>5.2f},{1:>6.1f}%$\pm${2:>5.2f}'
			labels.append(fmt.format(shift_temp, P, threshold_aod)) 
				
	args.update({ 'label' : labels })
	(n, bins, patches) = axis.hist(data, **args)

	#_whisker plot	
	if not by_var:
		x_truth = set(truth[sv])
		print data.shape, len(x_truth)
		ax0.boxplot(data, positions=x_truth, widths=0.05)

	#_add legend
	rem = 1 if (j+k) % 26 else 0
	ncol = (j+k) / 26 + rem 
	lgd = axis.legend(fontsize=8, loc=2, bbox_to_anchor=(1.05,1.),
		borderaxespad=0., frameon=False, ncol=ncol)
	#_highlight within TOLERANCE level (90% within 0.01?) 

	#_particulars
	axis.set_xlabel('Residuals (True OD - S-HIS OE Retrieval)', size='small')
	axis.set_xlim(-1, 1)
	axis.set_xticks(linspace(-1,1,21))
	shrink_ticks(axis)

	#_savefig
	tf = type(by_var) == str
	pname = 'surf_temp_tolerance{0}.png'.format(tf*'_{0}'.format(by_var))
	pname = os.path.join(dir_plot, pname)
	dbg(pname)
	plt.savefig(pname, bbox_extra_artists=(lgd,), bbox_inches='tight')


def fig_03(sv='tau', fnames=[], delta_ts=[0], dir_plot='.', dir_hs3='.',
	by_var=False, threshold_aod=0.1, **kwargs):
	'''
	plot box/whisker things for changes in surf temp

	delta_ts	array of departures from the true surface temp
	fnames		list,		names of input files (SIMULATED)
	'''
	from oe_sim import read_simulated_retrieval
	import matplotlib.pyplot as plt
	from numpy import array, append, linspace, zeros
	from libtools import shrink_ticks, colors_spag
	from numpy.ma import masked_where

	fig, axis = plt.subplots(len(delta_ts)-1)
	nshift = len(delta_ts)

	#_loop over shift in T 
	labels = []
	delta_ts = delta_ts[delta_ts != 0]
	for j, shift_temp in enumerate(delta_ts): 
		#_ignore non-shifted
		if shift_temp == 0:
			continue

		#_define out_label
		out_label = 'SFC_TEMP_{0:4.2f}'.format(shift_temp)

		#_loop over input files
		for i, fname in enumerate(fnames):
			lab = re.search('SIMULATED.(.*?).nc$', fname).group(1)
			lab = '{0}.{1}'.format(lab, out_label)
			fname = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.retrieved'.format(lab)
			fname = os.path.join(dir_hs3, fname)

			if not i: #_read in data
				truth, retr, uncrt = read_simulated_retrieval(fname, **kwargs)
			else: #_append data together
				t, r, u = read_simulated_retrieval(fname, **kwargs)
				for sv in truth:
					truth[sv] = append(truth[sv], t[sv])
				for sv in retr:
					retr[sv] = append(retr[sv], r[sv])
				for sv in uncrt:
					uncrt[sv] = append(uncrt[sv], u[sv]) 

		#_calculate residuals
		res = abs(array(truth[sv] - retr[sv]))
		dbg(('TRUTH', retr[sv].max(), retr[sv].min()))

		#_fill in residual array
		x_truth = list(set(truth[sv]))
		x_truth.sort()
		nx = len(x_truth)

		data = zeros((len(truth[sv]) / nx, nx)) - 9999.
		for m, value in enumerate(x_truth): 
			idx = truth[sv] == value	
			data[:,m] = array(res[idx])

		#_setup plot
		args = {
			'bins' : 100, 
			'alpha' : 0.5, 
			'range' : (-1, 1),
			'stacked' : True,
			'edgecolor' : None,
			} 
			
		#_Make and image separating out heights
		#_make whisker plot
			
		#_whisker plot	
		data = array(data)
		data = masked_where(data > 1000, data)
	
		#_make fractional values
		data = data / array(x_truth)
		axis[j].boxplot(data, positions=x_truth, widths=0.05)

		#_particulars
#		axis.set_xlabel('Residuals (True OD - S-HIS OE Retrieval)', size='small')
		axis[j].set_xlim(0, x_truth[-1] + x_truth[0])
		axis[j].set_ylim(0,3)
		axis[j].xaxis.set_visible(False)

		axis[j].set_ylabel('{0:5.2f} K'.format(shift_temp), size='xx-small')

		axis[j].text(x_truth[-4],1, 'N={0:d}'.format(data[:,0].size),
			size='xx-small')

	#_set some stuff for last plot
	else:
		axis[j].xaxis.set_visible(True)
		axis[j].set_xticks(x_truth)
		axis[j].set_xlabel('Optical Depth', size='xx-small')

	fig.suptitle('Fractional Residual (SHISoe - True OD) for Simulated Data\n' +			'when surface temp arbitrarily borked $\pm$5 K', size='small')
	#_clean up	
	[shrink_ticks(ax, 5) for ax in axis]

	#_savefig
	tf = type(by_var) == str
	pname = 'surf_temp_boxplot{0}.png'.format(tf*'_{0}'.format(by_var))
	pname = os.path.join(dir_plot, pname)
	dbg(pname)
	plt.savefig(pname)#, bbox_extra_artists=(lgd,), bbox_inches='tight')


def fig_04(aods=[.24,.83,1.42,2.0], ref=1.1, z=3., **kwargs):
	''' plot matrix of delta whatever '''
	from oe_sim import namelist, read_sensitivity_qsub
	import matplotlib.pyplot as plt
	from numpy import zeros
	from libtools import combination, shrink_ticks
	from matplotlib.cm import RdBu_r
	
	#_what we're looping over
	arbitrary_fudge = namelist['arbitrary_fudge']

	#_setup draw area
	ncol = 2
	extra = 1 if len(aods) % ncol else 0 
	nrow = len(aods) / ncol + extra
	fig, axes = plt.subplots(nrow, ncol)

	#_Set
	order = arbitrary_fudge.keys()
	nx = len(arbitrary_fudge[order[0]])
	ny = len(arbitrary_fudge[order[1]])

	#_for simplicity if not lucidity 
	xvar, yvar = order
	xvals, yvals = arbitrary_fudge[xvar], arbitrary_fudge[yvar]

	#_read in all retrievals
	data = read_sensitivity_qsub(**namelist)

	#_loop over, fill in plt_data
	for j, tau in enumerate(aods):
		plt_data = zeros((nx*ny)) - 9999.

		#_initialize
		d0 = getattr(data, xvar) == 0
		d1 = getattr(data, yvar) == 0
		d2 = data.tau == tau
		bidx = (getattr(data, xvar) == 0) * (getattr(data, yvar) == 0) *	\
			(data.tau == tau)
		aod_base = data.retr_tau[bidx]

		#_sanity check
		if bidx.sum() != 1:
			raise RuntimeError, 'NO BASE CASE'

		#_loop over x and y and fill in
		combos = combination((arbitrary_fudge[xvar], arbitrary_fudge[yvar]))
		for i, vals in enumerate(combos):
			x, y = vals
			xidx = getattr(data, xvar) == x 
			yidx = getattr(data, yvar) == y
			tidx = data.tau == tau
			idx = xidx * yidx * tidx 
		
			if idx.sum() != 1:
				raise RuntimeError, 'NO BASE CASE'
				
			#_fill in plot data
			print aod_base 
			plt_data[i] = aod_base - data.retr_tau[idx] 
		
		k0 = j / ncol
		k1 = j % ncol
		plt_data = plt_data.reshape(nx,ny)
		ax = axes[k0, k1]
		cs = ax.pcolor(xvals, yvals, plt_data, cmap=RdBu_r)
		if k0:
			ax.set_xlabel(xvar)
		if not k1:
			ax.set_ylabel(yvar)
		ax.set_xlim(xvals[0], xvals[-1])
		ax.set_xticks(xvals)
		ax.set_yticks(yvals)
		
		ax.text(0,0,'AOD={0:5.1f}'.format(tau))

		shrink_ticks(ax)
	plt.colorbar(cs, shrink=0.9)
	plt.savefig('/home/wsessions/test.png')


def fig_06(dtgs, dir_naaps=DIR_NRL,
	species=['sulfate_aod','seasalt_aod', 'dust_aod','smoke_aod', 'total_aod'],
	region='hs3_2013', **kwargs):
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

	#_generate filename for intermediate step
	xname = 'naaps_max_{{0}}_flights.dat'.format(region)
	xname = os.path.join(DIR_PROD, 'NRL', xname)

	#_initialize collection dictionary
	max_aod = {}
	[max_aod.update({ v : None }) for v in species]

	#_read in all flight track data for 2013
	track = read_tracks() 

	#_see if files created
	if not all([os.path.exists(xname.format(s)) for s in species]):
	
		#_loop over and read dtgs until finished
		fmt = '{0}/NAAPSAOD/{{1}}/{{0}}_aod'.format(DIR_NRL)
		for i, dtg in enumerate(dtgs):
	
			#_generate filename
			fname = fmt.format(dtg, dtg[:6])
			dbg((fname, os.path.exists(fname)))
	
			#_if file does not exist, continue
			if not os.path.exists(fname):
				dbg(('missing {0}'.format(fname)))
				continue

			#_read in data
			data = read_naaps(fname)

			#_loop over, pull out data
			for s in species:
				sidx = data.variable == s
	
				#_pull out current max_aod data
				tmp = max_aod[s]
	
				#_subset by region
				read = sub_region(data, region=region)
				if i == 0: 
					max_aod[s] = read.values[sidx][0][:,:,None]
				else:
					max_aod[s] = append(tmp, read.values[sidx][0][:,:,None], 2)
					max_aod[s] = max_aod[s].max(axis=2)[:,:,None]

		#_write max values for period
		lats = read.values[0].lat
		lons = read.values[0].lon
		for s in species:
			sxname = xname.format(s)
			write_max(max_aod[s], lats, lons, sxname) 

	else:
		#_if data already computed, read in
		for s in species:
			sxname = xname.format(s)
			max_aod[s], lats, lons = read_max(sxname)	

	#_get coordinates
	lon2d, lat2d = meshgrid(lons, lats)

	#_initialize draw space
	nspec = len(species)
	if nspec < 5:
		fig, ax = plt.subplots(nspec, 1)
	else:
		fig = plt.figure()
		ax0 = plt.subplot2grid((3,2), loc=(0,0))#, colspan=2)
		ax1 = plt.subplot2grid((3,2), loc=(0,1))#, colspan=2)
		ax2 = plt.subplot2grid((3,2), loc=(1,0))#, colspan=2)
		ax3 = plt.subplot2grid((3,2), loc=(1,1))#, colspan=2)
		ax4 = plt.subplot2grid((3,2), loc=(2,0), colspan=2)
		ax = [ax0, ax1, ax2, ax3, ax4]

	#_plot around tropical atlantic
	map_args = {
		'corners' : reg_meta['corn'], 
        'grid'    : reg_meta['grid'], 
        'delta'   : reg_meta['delta'],
		'laby'    : [0,1,0,0],
		}
	
	#_get out contour levels
	levs = fields()['aod']['levs']
	for i, s in enumerate(species):

		#_get max for gridpoints
		xaod = max_aod[s].squeeze()
		m = draw_map(ax=ax[i], **map_args)
		m.plot(track.lon, track.lat, 'k-', linewidth=0.3, zorder=7)
		CS = m.contourf(lon2d, lat2d, xaod, levels=levs, 
			colors=rgbgen('aod',nens=None), zorder=1)

		#_labet_
		ax[i].set_ylabel(s.upper().replace('_',' '), size='x-small')
	
	fig.suptitle('NAAPS Max Aerosol Optical Depths during HS3 2013',
		size='small')

	fig.subplots_adjust(bottom=0.1)
	cax = fig.add_axes([0.3, 0.05, 0.42, 0.02])
	fig.colorbar(CS, cax=cax, orientation='horizontal')
	shrink_ticks(cax)

	#_save
	pname = 'naaps_{0}_flights.png'.format(region)
	pname = os.path.join(DIR_PLOT, pname)
	dbg(pname)
	plt.savefig(pname)


def fig_05(dtg0='2013080100', dtg1='2013093018', dir_naaps=DIR_NRL,
	dt=6, species=['sulfate_aod','seasalt_aod', 'dust_aod','smoke_aod',
	'total_aod'], region='hs3_2013', **kwargs):
	'''
	plot NAAPS data for HS3 2013 field campaign from climatology

	start July 2013 - Oct 1, 2013, maximum values for all species,
						and then just smoke and dust
	'''
	from libnva import read_naapsaod as read_naaps
	from libnva import sub_region, draw_map
	from libtools import newdtg, shrink_ticks
	from libmeta import fields, regions
	import matplotlib.pyplot as plt
	from numpy import append, meshgrid
	from libtools import rgbgen, epoch2iso, dtg2iso
	from mpl_toolkits.axes_grid1 import make_axes_locatable

	#_pull out region data
	reg_meta = regions()[region]

	#_generate filename for intermediate step
	xname = 'naaps_max_{{0}}_{0}_{1}_{2}.dat'.format(dtg0, dtg1, region)
	xname = os.path.join(DIR_PROD, 'NRL', xname)

	#_initialize collection dictionary
	max_aod = {}
	[max_aod.update({ v : None }) for v in species]
	#_see if files created
	if not all([os.path.exists(xname.format(s)) for s in species]):
	
		#_loop over and read dtgs until finished
		fmt = '{0}/NAAPSAOD/{{1}}/{{0}}_aod'.format(DIR_NRL)
		dtg = dtg0
		while dtg <= dtg1:
			
			#_generate filename
			fname = fmt.format(dtg, dtg[:6])
			dbg((fname, os.path.exists(fname)))
	
			#_if file does not exist, continue
			if not os.path.exists(fname):
				dbg(('missing {0}'.format(fname)))
				continue

			#_read in data
			data = read_naaps(fname)

			#_loop over, pull out data
			for s in species:
				sidx = data.variable == s
	
				#_pull out current max_aod data
				tmp = max_aod[s]
	
				#_subset by region
				read = sub_region(data, region=region)
				if dtg == dtg0:
					max_aod[s] = read.values[sidx][0][:,:,None]
				else:
					max_aod[s] = append(tmp, read.values[sidx][0][:,:,None], 2)
					max_aod[s] = max_aod[s].max(axis=2)[:,:,None]

			#_increment dtg
			dtg = newdtg(dtg, dt)

		#_write max values for period
		lats = read.values[0].lat
		lons = read.values[0].lon
		for s in species:
			sxname = xname.format(s)
			write_max(max_aod[s], lats, lons, sxname) 

	else:
		#_if data already computed, read in
		for s in species:
			sxname = xname.format(s)
			max_aod[s], lats, lons = read_max(sxname)	

	#_get coordinates
	lon2d, lat2d = meshgrid(lons, lats)

	#_initialize draw space
	nspec = len(species)
	if nspec < 5:
		fig, ax = plt.subplots(nspec, 1)
	else:
		fig = plt.figure()
		ax0 = plt.subplot2grid((3,2), loc=(0,0))#, colspan=2)
		ax1 = plt.subplot2grid((3,2), loc=(0,1))#, colspan=2)
		ax2 = plt.subplot2grid((3,2), loc=(1,0))#, colspan=2)
		ax3 = plt.subplot2grid((3,2), loc=(1,1))#, colspan=2)
		ax4 = plt.subplot2grid((3,2), loc=(2,0), colspan=2)
		ax = [ax0, ax1, ax2, ax3, ax4]

	#_plot around tropical atlantic
	map_args = {
		'corners' : reg_meta['corn'], 
        'grid'    : reg_meta['grid'], 
        'delta'   : reg_meta['delta'],
		'laby'    : [0,1,0,0],
		}
	
	#_get out contour levels
	levs = fields()['aod']['levs']
	for i, s in enumerate(species):
	#	if not i:
	#		ax[i].set_title('NAAPS Max Aerosol Optical Depths for HS3 2013',
	#			size='small')

		#_get max for gridpoints
		xaod = max_aod[s].squeeze()

		m = draw_map(ax=ax[i], **map_args)
		CS = m.contourf(lon2d, lat2d, xaod, levels=levs, 
			colors=rgbgen('aod',nens=None))

		#_labet_
		ax[i].set_ylabel(s.upper().replace('_',' '), size='x-small')
	
	dtg = (dtg2iso(dtg0), dtg2iso(dtg1))
	fig.suptitle('NAAPS Max Aerosol Optical Depths for {0} to {1}'.format(*dtg),
		size='small')

	'''
	divider = make_axes_locatable(ax[-1])
	cax = divider.append_axes('bottom', size='5%', pad=.05)
	'''
	fig.subplots_adjust(bottom=0.1)
	cax = fig.add_axes([0.3, 0.05, 0.42, 0.02])
	fig.colorbar(CS, cax=cax, orientation='horizontal')
	shrink_ticks(cax)

	#_save
	pname = 'naaps_{0}_{1}-{2}.png'.format(region, dtg0, dtg1)
	pname = os.path.join(DIR_PLOT, pname)
	dbg(pname)
	plt.savefig(pname)


def fig_07(mw=-27, dir_hs3=DIR_SHIS, vars=['tau','ref','z'], var_to_plot='tau',
	dir_plot='.', **kwargs):
	'''
	Plot comparison of retrieval results of mismatched 
	pseudo-observations and mineral type.

	Do only for one mw option at a time.

	'''
	from oe_sim import read_simulated_retrieval
	import matplotlib.pyplot as plt
	from libtools import unique, shrink_ticks

	#_stuff
	v = '_'.join(vars)
	mw = abs(mw)

	#_initialize artist object
	fig, ax = plt.subplots(4)
	for i, mw in enumerate([26, 27]):
		i0 = i * 2 
		i1 = i0 + 1
	
		#_setup crap 
		fmt = 'KAOLINITE_SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.{1}.retrieved'
		types = ['KAL', 'GYP']
		out_label0 = '{0}on{1}_{2:d}'.format(types[0], types[1], mw)	
		out_label1 = '{1}on{0}_{2:d}'.format(types[0], types[1], mw)	
		fname0 = os.path.join(dir_hs3, fmt.format(v, out_label0))
		fname1 = os.path.join(dir_hs3, fmt.format(v, out_label1))
	
	
		name = {'KAL' : 'Kaolinite', 'GYP' : 'Gypsum' }
		label0 = '{0} SSP, {1} layer'.format(name[types[0]], name[types[1]])
		label1 = '{0} SSP, {1} layer'.format(name[types[1]], name[types[0]])
		lab0 = 'mw={0}'.format(mw)
		lab1 = 'mw={0}'.format(mw)
		
		#_read in simulated retrievals
		t0, r0, u0 = read_simulated_retrieval(fname0, **kwargs)
		t1, r1, u1 = read_simulated_retrieval(fname1, **kwargs)
	
		#_get list of all possible real optical depths
		absc = unique(t0[var_to_plot])	
	
		#_get sizes
		ntau = len(unique(t0['tau']))
		nref = len(unique(t0['ref']))
		nzee = len(unique(t0['z']))
	
		#_plot the fractional differences
		tau0 = t0[var_to_plot].reshape(ntau, nref*nzee)
		tau1 = t1[var_to_plot].reshape(ntau, nref*nzee)
		ret0 = r0[var_to_plot].reshape(ntau, nref*nzee)
		ret1 = r1[var_to_plot].reshape(ntau, nref*nzee)
		dif0 = (ret0 - tau0) / tau0 
		dif1 = (ret1 - tau1) / tau1
	
	
		cs0 = ax[i0].boxplot(dif0.T, positions=absc, widths=0.05)
		cs1 = ax[i1].boxplot(dif1.T, positions=absc, widths=0.05)
	
		#_make labels on right
		ax[i0].set_ylabel(lab0, size='x-small')
		ax[i1].set_ylabel(lab1, size='x-small')
		ax[i0].text(absc[-6], 4, label0, size='xx-small') 
		ax[i1].text(absc[-6], 4, label1, size='xx-small') 
	
	[x.set_xlim(absc[0]-0.05, absc[-1]+0.05) for x in ax]
	[x.set_ylim(-1, 5) for x in ax]
	[shrink_ticks(x) for x in ax]
	ax[-1].set_xlabel('True AOD', size='x-small')
	ax[0].set_title('Fractional Difference (mw_option={0})'.format(mw),
		size='small') 
	fig.text(0.04, 0.5, '(SHIS-AOD - True) / True',
		va='center', size='x-small', rotation='vertical')
					
	#_save figure
	pname = 'pseudo.{0}.{1}.png'.format('-'.join(types), v)
#	pname = 'pseudo.{0}.{1}.{2}.png'.format('-'.join(types), v, mw)
	pname = os.path.join(dir_plot, pname)
	dbg(pname)
	plt.savefig(pname)


def fig_08(mw=-27, dir_hs3=DIR_SHIS, vars=['tau','ref','z'], var_to_plot='tau',
	dir_lbl=os.path.join(DIR_PROD, 'LBL-RTM_simulated'), dir_plot='.', **kwargs):
	'''
	Plot comparison of retrieval results of mismatched 
	pseudo-observations and mineral type.

	Do only for one mw option at a time.

	'''
	from oe_sim import read_simulated_retrieval
	import matplotlib.pyplot as plt
	from libtools import unique, shrink_ticks
	from hs3_utils import Flight_segment as F
	from libgeo import planck_inv
	from numpy import array, append
	from lblrtm_utils import read_lbldis_out, microwindow_average

	#_stuff
	v = '_'.join(vars)
	mw = abs(mw)

	#_read in simulated flight file
	fmt = 'lbldis_output.SIMULATED.{0:04d}.{1}.{2}.final.cdf'
	fmt_seg = '{1}_SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.nc'
	types = ['GYP', 'KAL']
	name = {'KAL' : 'Kaolinite', 'GYP' : 'Gypsum' }

	fig, ax = plt.subplots(4)

	for ii, typ in enumerate(types):	#_loop over mineral layer, then mw
		#_	
		layer = name[typ]
		#-really just trying to get this shit done, huh?
		ssp = 'GYP' if layer == 'Kaolinite' else 'KAL'
		label = '{0} SSP, {1} layer'.format(name[ssp], name[typ])
		
		#_read in sim retrieval output
		#_open flight file
		file_seg = fmt_seg.format(v, name[typ].upper())
		file_seg = os.path.join(dir_hs3, file_seg)
		flight = F(file_seg=file_seg)
		w_obs = flight.SHIS_wavenumber
	
		l26 = '{0}on{1}_{2:d}'.format(ssp, typ, 26)	
		l27 = '{0}on{1}_{2:d}'.format(ssp, typ, 27)	

		#_read in all radiances and append difference in BvT
		#_ planck_inv(lblr / 1e3, n_clr*100, domain='wavenumber')
		for i, fov in enumerate(flight):

			#_get field of views radiances
			r_obs = fov.SHIS_radiances
			r26_obs, w26 = microwindow_average(r_obs, w_obs, -26)
			r27_obs, w27 = microwindow_average(r_obs, w_obs, -27)
			t26_obs = planck_inv(r26_obs, w26*100, domain='wavenumber')
			t27_obs = planck_inv(r27_obs, w27*100, domain='wavenumber')

			#_generate retrieval output file names for both mw options
		#	'lbldis_output.SIMULATED.1319.tau_ref_z.KALonGYP_26.final.cdf'
			f26 = fmt.format(i, v, l26) 
			f27 = fmt.format(i, v, l27)
			f26 = os.path.join(dir_lbl, 'test_fov', f26) 			
			f27 = os.path.join(dir_lbl, 'test_fov', f27) 			

			#_read in radiances from retrievals
			f26 = read_lbldis_out(f26) 
			f27 = read_lbldis_out(f27) 

			#_get radiances and brightness temperatures
			r26 = f26.radiances
			r27 = f27.radiances
			t26 = planck_inv(r26, 100*f26.wavenumber, domain='wavenumber')
			t27 = planck_inv(r27, 100*f27.wavenumber, domain='wavenumber')

			#_get difference between truth and retrieval
			d26 = t26 - t26_obs
			d27 = t27 - t27_obs

			#_add to master array
			try:
				diff_26 = append(diff_26, d26[:,None], 1)
				diff_27 = append(diff_27, d27[:,None], 1)
			except UnboundLocalError:
				diff_26 = d26[:,None]
				diff_27 = d27[:,None]

		#_set indices
		i26 = ii % 2	#_26 plot
		i27 = i26 + 2	#_27 plot
		width_26 = (w26.max() - w26.min()) / 50
		width_27 = (w27.max() - w27.min()) / 50

			
		cs0 = ax[i26].boxplot(diff_26.T, positions=w26, widths=width_27)
		cs1 = ax[i27].boxplot(diff_27.T, positions=w27, widths=width_27)
	
		#_make labels on right
		ax[i26].set_ylabel('mw=26', size='x-small')
		ax[i27].set_ylabel('mw=27', size='x-small')
	#	ax[i26].text(805., -7.5, label, size='xx-small')
		ax[i26].text(w27[0], -6.6, label, size='xx-small') 
		ax[i27].text(w27[0], -20, label, size='xx-small') 

		#_fix xticks
##		ax[i26].set_xticks(w26, rotation=-30)
##		ax[i27].set_xticks(w27, rotation=-30)
	#	ax[i26].set_xlim(w26[0]-15, w26[-1]+15)	
		ax[i26].set_xlim(w27[0]-15, w27[-1]+15)	
		ax[i27].set_xlim(w27[0]-15, w27[-1]+15)	
		ax[i26].set_ylim(-10, 10)
		ax[i27].set_ylim(-30, 30)


	ax[0].set_title('Habit Sensitivity Test (Kaolinite and Gypsum)',
		size='small')
	[x.plot([-9999,9999], [0,0], 'k--', linewidth=0.4) for x in ax]
#	axtop = ax[0].twiny()
##	axtop.set_xlim(w26[0]-15, w26[-1]+15)	
##	axtop.set_xticks(w26)
#	axtop.set_xticks(w27)
#	axtop.set_xlim(w27[0]-15, w27[-1]+15)	
	ax[0].xaxis.set_visible(False)
	ax[1].xaxis.set_visible(False)
	ax[2].xaxis.set_visible(False)
	ax[3].set_xlabel('wavenumber ($cm^{-1}$)', size='small')

	[plt.setp(x.xaxis.get_majorticklabels(), rotation=70) for x in ax]
#	plt.setp(axtop.xaxis.get_majorticklabels(), rotation=70)
	[shrink_ticks(x) for x in ax]
#	shrink_ticks(axtop)
##	ax[-1].set_xlabel('True AOD', size='x-small')
##	ax[0].set_title('Fractional Difference (mw_option={0})'.format(mw),
##		size='small') 
	fig.text(0.04, 0.5, '$\Delta B_{\\nu}T$ (SHIS - True, $K$)',
		va='center', size='x-small', rotation='vertical')
					
	#_save figure
	pname = 'pseudo.BT.{0}.{1}.png'.format('-'.join(types), v)
	pname = os.path.join(dir_plot, pname)
	dbg(pname)
	plt.savefig(pname)


def fig_09(idx=[240,400,520], dtg='20130824175600', dz_cpl=25.,
	dir_plot='.', **kwargs):
	'''
	Image to show influence of dust layer on S-HIS
	Top Panel, three spectral IR lines at three FOVS
	Bot Panel, CPL with lines at appropriate panels.
	'''
	from hs3_utils import Flight_segment as F
	import matplotlib.pyplot as plt
	from libgeo import planck_inv
	from libtools import rgbgen
	from numpy import linspace, arange
	from libtools import epoch2iso, shrink_ticks

	#_read in flight segment file
	flight = F(dtg=dtg, **kwargs)

	#_init draw object
	fig, (ax_cpl, ax_rad) = plt.subplots(2)

	#_plot CPL
	cpl_plt, cpl_x, cpl_y = flight.map_cpl(dz_cpl=dz_cpl, **kwargs)
	cpl_nx = len(cpl_x)
	ax_cpl.pcolormesh(cpl_x, cpl_y, cpl_plt.T, vmin=-4e-7, vmax=1e-4)
	ax_cpl.set_title('CPL BACKSCAT 532_EXT', size='xx-small')
	ax_cpl.set_xlim(ax_cpl.xaxis.get_data_interval())
	ax_cpl.set_xticks(arange(0, cpl_nx, cpl_nx/5))
	ax_cpl.set_xticklabels([epoch2iso(ttt)[-8:] for
		ttt in flight.CPL_epoch[::cpl_nx/5]])
	ax_cpl.set_ylim(0, 6000)
	ax_cpl.set_yticks(linspace(0, 6000, 11))
##	ax_cpl.set_yticks(linspace(0, cpl_y[-1], 11))

	#_get number of lines
	ny = len(idx)
#	colors = rgbgen('spag', ny)
	colors = ['red', '#33ccff', '#33cc33']
	labels = ['clear', 'cloud + dust', 'dust']

	#_pull out x/y for rads
	x = flight.SHIS_wavenumber
	for ii, i in enumerate(idx):
		#_pull out radiances, convert to brightness temp
		r = flight[i].SHIS_radiances 
		y = planck_inv(r, x*100, domain='wavenumber')

		#_plot
		arg = {'zorder' : ii+2, 'color' : colors[ii], 'label' : labels[ii] }
		p = ax_rad.plot(x, y, linewidth=0.3, **arg) 

		#_plot location on cpl plot
		ax_cpl.plot([i,i], [0, cpl_y[-1]], linewidth=2., **arg) 

	#_
	ax_cpl.legend(loc=2, fontsize=8) # frameon=False, fontsize=8)
	
	#_limit plot area
	x0 = ax_rad.xaxis.get_data_interval()
	ax_rad.set_ylim(285, 298)
#	ax_rad.set_ylim(270, 300)
	ax_rad.set_xlim(750, 1100)
	ax_rad.grid(True)
	ax_rad.set_xlabel('wavenumber ($cm^{-1}$)', size='xx-small')
	ax_rad.set_ylabel('brightness temperature (K)', size='xx-small')

	#_
	[shrink_ticks(ax) for ax in [ax_cpl, ax_rad]]

	#_save image
	pname = 'spectral_shift_from_dust.png'
	pname = os.path.join(dir_plot, pname)
	dbg(pname)
	plt.savefig(pname)
	plt.close()


def write_max(data, lat, lon, fname, **kwargs):
	''' write maximum data for period '''
	from numpy import savetxt
	latname = fname.replace('max', 'lat')
	lonname = fname.replace('max', 'lon')
	savetxt(fname, data, delimiter=',')
	savetxt(latname, lat, delimiter=',')
	savetxt(lonname, lon, delimiter=',')
	

def read_max(fname):
	''' read in max aod values for period '''
	from numpy import loadtxt
	latname = fname.replace('max', 'lat')
	lonname = fname.replace('max', 'lon')
	
	data = loadtxt(fname, delimiter=',')
	lat = loadtxt(latname, delimiter=',')
	lon = loadtxt(lonname, delimiter=',')

	return data, lat, lon


def plot_simulated_histogram(out_label='', plot_var='tau',
	dir_hs3=os.path.join(os.environ['PRODUCTS'],'hs3'), **kwargs):
	'''
	produce some plots for this god damn stuff

	out_label	str,	experiment is not used so that it can
						denote specific input files as opposed 
						to output.  in_label should have been used.
	'''
	import matplotlib.pyplot as plt
	from libtools import mkdir_p, shrink_ticks
	from numpy import array, append, zeros

	#_read in retrieval data
	fname = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.retrieved'.format(out_label)
	fname = os.path.join(dir_hs3, fname)
	try:
		truth, retr, uncrt = read_simulated_retrieval(fname, **kwargs)
	except IOError:
		dbg(('Did you complete this retrieval yet?', out_label))
		os._exit(23)

	#_what are we working with?
	state_vars = retr.keys() # truth.keys()
	dbg(state_vars)

	#_kick back if not available
	if plot_var not in state_vars:
		dbg(('Variable not retrieved:', plot_var))
		return

	#_color dictionary
	order = ['low', 'mid', 'high', 'very_high']
	ranges = {
		'low'	: {
			'color'	: 'blue',
			'tau'	: (-9999., 0.5),
			'ref'	: (0.5, 0.675),
		},
		'mid'	: {
			'color'	: 'green',
			'tau'	: (0.5, 1.0),
			'ref'	: (0.625, 0.75),
		},
		'high'	: {
			'color'	: 'orange',
			'tau'	: (1.0, 1.5),
			'ref'	: (0.75, 0.875),
		},
		'very_high'	: {
			'color'	: 'red',
			'tau'	: (1.5, 2.0),
			'ref'	: (0.875, 1.0),
		},
		}

	#_output directory for plots
	ppath = 'PLOTS_{0}'.format(out_label)
	if not os.path.exists(ppath):
		mkdir_p(ppath)

	#_produce for each state variable
	for sv in state_vars:
		#_get residual
		res = truth[sv] - retr[sv]

		#_initialize figure and axes
		fig = plt.figure()
		axes = []
		figtext = []

		ymax, xmax, xmin = 0, 0, 0
		nbin = 100
		axes.append(plt.subplot2grid((1,1), (0,0)))

		#_plot histogram of residuals
		(n, bins, patches) = axes[0].hist(res, bins=nbin,
				facecolor='b', edgecolor='b')
		xmax = max((xmax, bins.max()))
		xmin = min((xmin, bins.min()))
		ymax = max((ymax, n.max()))
			
		axes[0].plot([0,0],[0,n.max()], 'k--')

		axes[0].set_xlabel('diff in tau, {0}'.format(sv), size='x-small')
		axes[0].set_title('retrieval residuals'.format(out_label), 
							size='xx-small')

		axes[0].set_ylabel('residuals, truth - oe_retrieval', size='x-small')

		#_max symmetric
		xbnd = max((abs(xmin), abs(xmax)))

		#_shrink...ticks...
		[ax.set_xlim(-xbnd, xbnd) for ax in axes]
		[ax.set_ylim(top=100) for ax in axes]
		[shrink_ticks(ax) for ax in axes]
	
		#_create image name
		pname = 'hist_{0}_{1}.png'.format(out_label, sv)
		pname = os.path.join(ppath, pname)
		plt.savefig(pname)
		plt.close()		
		
		dbg(pname)	


def plot_simulated_retrieval(out_label='', plot_var='tau', fsim=None,
	dir_hs3=os.path.join(os.environ['PRODUCTS'],'hs3'), **kwargs):
	'''
	produce some plots for this god damn stuff

	out_label	str,	experiment is not used so that it can
						denote specific input files as opposed 
						to output.  in_label should have been used.
	'''
	import matplotlib.pyplot as plt
	from libtools import mkdir_p, shrink_ticks
	from numpy import array, append

	res = re.search('SIMULATED.(.*?).nc', fsim)
	sim_experiment = res.group(1)

	#_read in retrieval data
	fname = 'SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.{1}.retrieved'.format(
		sim_experiment, out_label)
	fname = os.path.join(dir_hs3, fname)
	try:
		truth, retr, uncrt = read_simulated_retrieval(fname, **kwargs)
	except IOError:
		dbg(('Did you complete this retrieval yet?', out_label))
		os._exit(23)

	#_what are we working with?
	state_vars = retr.keys() # truth.keys()
	dbg(state_vars)

	#_kick back if not available
	if plot_var not in state_vars:
		dbg(('Variable not retrieved:', plot_var))
		return

	#_color dictionary
	order = ['low', 'mid', 'high', 'very_high']
	ranges = {
		'low'	: {
			'color'	: 'blue',
			'tau'	: (-9999., 0.5),
			'ref'	: (0.5, 0.675),
		},
		'mid'	: {
			'color'	: 'green',
			'tau'	: (0.5, 1.0),
			'ref'	: (0.625, 0.75),
		},
		'high'	: {
			'color'	: 'orange',
			'tau'	: (1.0, 1.5),
			'ref'	: (0.75, 0.875),
		},
		'very_high'	: {
			'color'	: 'red',
			'tau'	: (1.5, 9999),
			'ref'	: (0.875, 1.0),
		},
		}

	#_output directory for plots
	ppath = 'PLOTS_{0}'.format(out_label)
	if not os.path.exists(ppath):
		mkdir_p(ppath)

	#_produce for each state variable
	for sv in state_vars:
		#_get residual
		res = truth[sv] - retr[sv]

		#_get indices for ranges
	 	for label, opts in ranges.iteritems():
			min_val, max_val = opts[sv]
		#	idx0 = truth[sv] >= min_val 
		#	idx1 = truth[sv] < max_val
		#	idx = append(idx0[:,None], idx1[:,None], axis=1).all(axis=1)
			idx = (truth[sv] >= min_val) * (truth[sv] < max_val)
			ranges[label]['idx'] = idx

		#_initialize figure and axes
		fig = plt.figure()
		axes = []
		axis = plt.subplot2grid((2,4), (1,0), colspan=len(ranges.keys()))
		figtext = []

		ymax, xmax, xmin = 0, 0, 0
		nbin = 20
		for i, label in enumerate(order): #ranges):
			min_val, max_val = ranges[label][sv]

			idx = ranges[label]['idx']
			col = ranges[label]['color']
			if idx.sum() == 0: 
				continue	

			axes.append(plt.subplot2grid((2,4), (0,i)))

			#_plot histogram of residuals
			(n, bins, patches) = axes[i].hist(res[idx], bins=nbin,
											facecolor=col, edgecolor=col)
			xmax = max((xmax, bins.max()))
			xmin = min((xmin, bins.min()))
			ymax = max((ymax, n.max()))
			
			#_add label
			figtext.append('{0}\nmin {1:5.1f}\nmax {2:5.1f}'.format(label,
								min_val, max_val))

			axes[i].set_xlabel('difference in {0}'.format(sv), size='x-small')
			
			#_plot uncertainty by truth
			fract = uncrt[sv][idx] / truth[sv][idx]
			x_truth, y_fract = [], []
			for value in set(truth[sv][idx]):
				x_truth.append(value)
				y_fract.append(fract[truth[sv][idx] == value])	

		##	axis.scatter(truth[sv][idx], fract, marker='x', color=col)
			axis.boxplot(y_fract, positions=x_truth, widths=0.05)
		##	axis.plot(truth[sv][idx], fract, marker='x', color=col)
		##	axis.scatter(truth[sv][idx], uncrt[sv][idx], marker='x', color=col)

		axis.set_xticks(list(set(truth[sv])))
		axes[0].set_title('{0} (experiment label)'.format(out_label), 
							size='xx-small')

		#_various things for uncertainty
		axis.set_ylim(bottom=0)
		axis.set_xlim(truth[sv].min()-0.1, truth[sv].max()+0.1)
		axis.set_ylabel('posterior uncercertainty', size='x-small')
		axis.set_xlabel('{0}'.format(sv), size='x-small')
		axis.grid(False)
		
		axes[0].set_ylabel('residuals, truth - oe_retrieval', size='x-small')

		#_max symmetric
		xbnd = max((abs(xmin), abs(xmax)))
		for i, label in enumerate(order):	
			axes[i].text(0, ymax-ymax/5., figtext[i], size='xx-small')

		#_shrink...ticks...
		[ax.set_xlim(-xbnd, xbnd) for ax in axes]
		[ax.set_ylim(top=ymax) for ax in axes]
		[shrink_ticks(ax) for ax in axes]
		shrink_ticks(axis)
	
		#_create image name
		pname = 'hist_uncert.{0}.{1}.{2}.png'.format(sim_experiment,
			out_label, sv)
		pname = os.path.join(ppath, pname)
		plt.savefig(pname)
		plt.close()		
		
		dbg(pname)


def plot_flight_tracks(fname='/data/wsessions/hs3/hs3_2013_flighttracks.dat',
	**kwargs):
	''' plot flight tracks for particular dates '''
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
#	idx = append(idx, [False])
	idx = append([False], idx)
	print idx.size, final.size
	final.epoch[idx] = None 
	final.lat[idx] = None 
	final.lon[idx] = None 

	print (diff(final.epoch) > 10).sum(), 'TESTING'

	plt.plot(final.lon, final.lat)
	plt.savefig('/home/wsessions/tracks.png')
	plt.close()
	
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
#	idx = append(idx, [False])
	idx = append([False], idx)
	final.epoch[idx] = None 
	final.lat[idx] = None 
	final.lon[idx] = None 

	return final


def read_flight_tracks(dtg0='1900010100', dtg1='3001010100', **kwargs):
	''' read in SHIS flight track data from netdf files '''
	from glob import glob
	from libtools import merge
	from hs3_utils import read_cpl_flighttrack

	#_get filenames desired
	dir_cplnc = '/data/wsessions/cpl/nc'
	files = glob("{0}/OP_13*nc".format(dir_cplnc))

	for i, file in enumerate(files):
		if not i:
			data = read_cpl_flighttrack(file)

		else:
			data = merge((data, read_cpl_flighttrack(file)))

	return data


def dbg(msg, l=1, err=False):
	''' 
	if global debug is set to true, be more verbose 
	msg : str, Message to be printed
	l   : int, Debug level of message.  Set higher for lower level 
	        messages.  As debug increases, noisiness should also.
	'''
	import inspect

	if hasattr(msg, '__iter__'):
		msg = ' '.join([str(m) for m in msg])
	elif type(msg) != str:
		msg = str(msg)

	if DEBUG >= l:
		curf    = inspect.currentframe()
		calf    = inspect.getouterframes(curf,2)
		file, line, method = calf[1][1:4]
		file     = '.'.join(file.split('/')[-1].split('.')[:-1])
		scream  = '[%s.%s.%i] %s' % (file,method,line,msg)

		if not err:
			print scream
		else:
			raise RuntimeError, scream


################################################################################
#_SHELL_########################################################################
################################################################################


if __name__ == '__main__':
	run_main(**namelist)
	
