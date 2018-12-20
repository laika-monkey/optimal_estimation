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
	'delta_ts'		: arange(8) - 5, 
	'by_var'		: 'ref', #_z, ref, tau

	#_threshold for optical depth variance
	'threshold_aod'	: 0.1, 

	}


################################################################################
################################################################################
################################################################################


def run_main(vnames=None, **kwargs):
	''' plot all figures '''

#	fig_00(**kwargs)

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
		namelist.update({out_label : label})
		fig_01(sv, files=fnames, **namelist)


##	#_plot total uncertainty

	#_plot surf_temp variance
	fig_02(**kwargs)


################################################################################
################################################################################
################################################################################


def fig_00(pname='delta_bt.png', dir_plot='.', dv=-26, **kwargs):
	'''
	Spectral signatures of habits by wavenumber

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
				'gypsum'	: 'green' }

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

	fig, axis = plt.subplots(1)
	nshift = len(delta_ts)

	#_loop over shift in T 
	for j, shift_temp in enumerate(delta_ts): 
		#_ignore non-shifted
		if shift_temp == 0:
			continue

		#_define out_label
	##	kwargs.update({'out_label' : 'SFC_TEMP_{0:4.2f}'.format(shift_temp)})
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
		res = truth[sv] - retr[sv]

		#_setup plot
		args = {
			'bins' : 100, 
			'alpha' : 0.5, 
			'range' : (-1, 1)
			} 
		
		#_Make and image separating out heights
		#_make whisker plot

		#_make histogram
		if by_var != False:
			vars = list(set(truth[by_var]))
			nvar = len(vars)
	
			#_only do this crap first time through
			if j == 0:
				cols = colors_spag(nshift*nvar)
	
			for k, v in enumerate(vars):
				idx = truth[by_var] == v 

				#_calculate percentage within thresh?
				P = (abs(res[idx]) <= threshold_aod).sum() / float(idx.sum())
				P *= 100
				fmt = '$\Delta$T={0:>5.2f}, {4}={1:>4.2f},{2:>6.1f}%' +\
					'$\pm${3:>5.2f}'
				args.update({
					'facecolor' : cols[j*nvar+k],
					'edgecolor' : cols[j*nvar+k],
					'label' : fmt.format(shift_temp, v, P, threshold_aod, by_var
					)}) 
				(n, bins, patches) = axis.hist(res[idx], **args)

		#_make histogram
		else:
			if j == 0:
				#_only do this crap first time through
				cols = colors_spag(nshift)
				k = 0
	
			#_calculate percentage within thresh?
			P = (abs(res) <= threshold_aod).sum() / float(res.size)
			P *= 100
			fmt = '$\Delta$T={0:>5.2f},{1:>6.1f}%$\pm${2:>5.2f}'
			args.update({
					'facecolor' : cols[j],
					'edgecolor' : cols[j],
					'label' : fmt.format(shift_temp, P, threshold_aod) }) 
			(n, bins, patches) = axis.hist(res, **args)

	#_add legend
	rem = 1 if (j+k) % 26 else 0
	ncol = (j+k) / 26 + rem 
	dbg((ncol, rem, j+k))
	lgd = axis.legend(fontsize=8, loc=2, bbox_to_anchor=(1.05,1.),
		borderaxespad=0., frameon=False)
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
	from numpy import array, append

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
	
