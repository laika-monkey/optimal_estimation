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
if 'DISPLAY' not in os.environ:
    import matplotlib
    matplotlib.use('Agg')

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
	'fractional' : False,	#_get residuals for real cases as fraction of CPLOD

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
	'by_var'		: 'tau', #_False, z, ref, tau

	#_threshold for optical depth variance
	'threshold_aod'	: 0.1, 

	#_which naaps species to plot
	'species'	: ['dust_aod', 'smoke_aod', 'total_aod'],
	'region'	: 'hs3_2013b',	

	#_probably going to end up looping these later
	'out_label' : 'NVACLIMO1_M2_DV26_CLD_CTP-NAAPS_3KM', 
    'out_labels': [ 
				#	'NVACLIMO1_tau-ref_M2_DV26_CLD_CTP-NAAPS_3KM',
				#	'NVACLIMO1_tau-ref_M2_DV27_CLD_CTP-NAAPS_3KM',
					'NVACLIMO1_tau-ref_M2_DV26_CLD_NOCTP_3KM',
					'NVACLIMO1_tau-ref_M2_DV27_CLD_NOCTP_3KM',
                #	'NVACLIMO1_tau-ref_M2_DV26_NOCLD_CTP-NAAPS_3KM',
				#	'NVACLIMO1_tau-ref_M2_DV27_CLD_CTP-NAAPS_3KM',
				#	'NVACLIMO1_tau-ref_M2_DV27_NOCLD_CTP-NAAPS_3KM',
				#	'NVACLIMO1_tau-ref_M2_DV27_NOCLD_NOCTP_3KM',
				#	'NVACLIMO1_tau-ref_M2_DV27b_CLD_CTP-NAAPS_3KM',
				#	'NVACLIMO1_tau-ref_M2_DV27b_NOCLD_CTP-NAAPS_3KM'
					],
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
	if 0:
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
#	fig_06(dtgs, **kwargs)

	#_plot habit comparison of KAOL and GYP
##	fig_07(**kwargs)

	#_spectral comparison of KAOL and GYP
##	fig_08(**kwargs)	

	#_dashboard light (test case of three FOVs for jeff)
#	fig_09(**kwargs)

	#_KALonKAL, GYPonGYP (fig-07 redone)
#	fig_10(**kwargs)
#	fig_11(**kwargs)

	#_BOTHonKAL, BOTHonGYP, by AOD, spectrally
#	fig_14(**kwargs)
#	fig_15(**kwargs)

	#_compare cloud retrievals between CPL and S-HISoe
	if 0:
 	  for dtg, values in experiments['HS3'].iteritems():
		for out_label in kwargs.get('out_labels'):
			kwargs.update({'out_label' : out_label})
			dbg((out_label, dtg))
			kwargs.update(values)
			fig_12(dtg, **kwargs)

	#_plot histogram comparing real test case residuals
#	fig_13(**kwargs)

	#_do box/whisker for entire compaign
	fig_CAMPAIGN(**kwargs)

	#_write table of probabilities of borked sfc t getting right answer
#	fig_16(**kwargs)

##	fig_20(**kwargs)



################################################################################
################################################################################
################################################################################


def fig_00(pname='delta_bt_{0}.png', dir_plot='.', dv=-27, **kwargs):
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

	dir_in = '/odyssey/isis/users/wsessions/LBL-RTM_simulated/std_tropical'
	habits = {
		#		'quartz'	: 'purple', 
				'kaolinite'			: 'red', 
#				'kaolinite_thick'	: 'orange', 
#				'gypsum'			: 'green', 
#				'gypsum_thick'		: 'cyan', 
#				'gypsum_small'		: 'blue', 
#				'cloud'		: 'blue',
#				'sand'		: 'brown',
#				'dust'		: 'yellow',
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
	lin0 = atwn.plot(wavs, bt1, 'k--', linewidth=0.2, label='dusty sky')

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
	pname = os.path.join(dir_plot, pname.format('-'.join(habits.keys())))
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
			axis.set_ylim(bottom=0, top=2)
			axis.set_xlim(truth[sv].min()-0.1, truth[sv].max()+0.1)
			axis.set_ylabel('{0} = {1}'.format(var, value), size='xx-small')
			axis.set_ylabel('{0} = {1}'.format(var, value), size='x-small')
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


def fig_16(sv='tau', fnames=[], delta_ts=[0], dir_plot='.', dir_hs3='.',
	by_var=False, threshold_aod=0.1, tname='table_3.txt', **kwargs):
	'''
	Make table of probability of being within threshold value of AOD
	for arbitrarily wrong surface temperatures

	delta_ts	array of departures from the true surface temp
	fnames		list,		names of input files (SIMULATED)
	'''
	from oe_delete import read_simulated_retrieval
#	from oe_sim import read_simulated_retrieval
	import matplotlib.pyplot as plt
	from numpy import array, append, linspace, recarray
	from libtools import shrink_ticks, colors_spag

	fig, axes = plt.subplots(2,1)
	nshift = len(delta_ts)

	axis, ax0 = axes

	#_open output table
	f = open(tname, 'w')

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

					#_initialize recarray
					dtype = [('P', 'f4'), ('dT', 'f4'), ('OD', 'f4')]
				#	probs = recarray((nvar * len(delta_ts)), dtype)
					probs = recarray((0,), dtype)
				else:
					data = append(data, res[idx].reshape(idx.sum(),1), axis=1)

				#_calculate percentage within thresh?
			#	P = (abs(res[idx]) <= threshold_aod).sum() / float(idx.sum())
				#_within 10%
				P = (abs(res[idx]) <= v*.15).sum() / float(idx.sum())
				P *= 100
				fmt = '$\Delta$T={0:>5.2f}, {4}={1:>4.2f},{2:>6.1f}%' +\
					'$\pm${3:>5.2f}'
				labels.append(fmt.format(shift_temp,v,P,threshold_aod,by_var))

				#_add information into recarray
				probs.resize(probs.size + 1)
				probs[-1] = (P, shift_temp, v)

		#_collect data and probabilities for total
		else:
			raise RuntimeError, 'Not implemented.'
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

		#_generate output in form of nvar columns x n_dT rows
		tmp = probs[probs.dT == shift_temp]
		tmp = tmp[tmp.OD.argsort()]

		#_write header	
		if j == 0:
			f.write(' '*6 + ','.join(['{0:>6.1f}'.format(v) for v in tmp.OD]))	
			f.write('\n')

		out = '{0:>6.0f}'.format(shift_temp)
		out += ','.join(['{0:>6.1f}'.format(v) for v in tmp.P])
		f.write(out + '\n')

	dbg(tname)
	f.close()

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
	if nspec == 3:
		fig = plt.figure()
		ax0 = plt.subplot2grid((2,4), loc=(0,0), colspan=2)
		ax1 = plt.subplot2grid((2,4), loc=(0,2), colspan=2)
		ax2 = plt.subplot2grid((2,4), loc=(1,1), colspan=2)
		ax = [ax0, ax1, ax2]
	elif nspec < 5:
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
		if nspec == 3 and i == 0:
			map_args = {
			'corners' : reg_meta['corn'], 
		    'grid'    : reg_meta['grid'], 
		    'delta'   : reg_meta['delta'],
			'laby'    : [0,0,0,0],
			}
		else:
			map_args = {
			'corners' : reg_meta['corn'], 
		    'grid'    : reg_meta['grid'], 
		    'delta'   : reg_meta['delta'],
			'laby'    : [0,1,0,0],
			}
			


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
	if nspec == 3: plt.tight_layout()
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


def fig_15(dir_hs3=DIR_SHIS, vars=['tau','ref','z'], var_to_plot='tau',
	dir_lbl=os.path.join(DIR_PROD, 'LBL-RTM_simulated'), dir_plot='.',**kwargs):
	'''
	same as fig_08 (differences in Bt by MW channel), but with same 
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

	#_read in simulated flight file
	fmt = 'lbldis_output.SIMULATED.{0:04d}.{1}.{2}_tau-ref.final.cdf'
	fmt_seg = '{1}_SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.nc'
	types = ['GYP', 'KAL']
	name = {'KAL' : 'Kaolinite', 'GYP' : 'Gypsum' }

	fig, ax = plt.subplots(4)

	for ii, typ in enumerate(types):	#_loop over mineral layer, then mw
		#_	
		layer = name[typ]

		#-really just trying to get this shit done, huh?
		label = '{0} layer'.format(name[typ])
		
		#_read in sim retrieval output
		#_open flight file
		file_seg = fmt_seg.format(v, name[typ].upper())
		file_seg = os.path.join(dir_hs3, file_seg)
		flight = F(file_seg=file_seg)
		w_obs = flight.SHIS_wavenumber
	
		l26 = 'BOTHon{0}_{1:d}'.format(typ, 26)	
		l27 = 'BOTHon{0}_{1:d}'.format(typ, 27)	

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
		ax[i26].text(w27[0], -.6, label, size='xx-small') 
		ax[i27].text(w27[0], -.6, label, size='xx-small') 

		#_fix xticks
		ax[i26].set_xlim(w27[0]-15, w27[-1]+15)	
		ax[i27].set_xlim(w27[0]-15, w27[-1]+15)	
		ax[i26].set_ylim(-1, 1)
		ax[i27].set_ylim(-1, 1)

	ax[0].set_title('Discrimination Test (Kaolinite and Gypsum)',
		size='small')
	[x.plot([-9999,9999], [0,0], 'k--', linewidth=0.4) for x in ax]
	ax[0].xaxis.set_visible(False)
	ax[1].xaxis.set_visible(False)
	ax[2].xaxis.set_visible(False)
	ax[3].set_xlabel('wavenumber ($cm^{-1}$)', size='small')

	[plt.setp(x.xaxis.get_majorticklabels(), rotation=70) for x in ax]
	[shrink_ticks(x) for x in ax]
	fig.text(0.04, 0.5, '$\Delta B_{\\nu}T$ (SHIS - True, $K$)',
		va='center', size='x-small', rotation='vertical')
					
	#_save figure
	pname = 'pseudo.both0.BT.{0}.{1}.png'.format('-'.join(types), v)
	pname = os.path.join(dir_plot, pname)
	dbg(pname)
	plt.savefig(pname)



def fig_14(dir_hs3=DIR_SHIS, vars=['tau','ref','z'], var_to_plot='tau',
	dir_plot='.', **kwargs):
	'''
	Plot comparison of retrieval results of mismatched 
	pseudo-observations and mineral type.
	BOTHonKAL
	BOTHonGYP

	'''
	from oe_sim import read_simulated_retrieval
	import matplotlib.pyplot as plt
	from libtools import unique, shrink_ticks
	import re

	#_stuff
	v = '_'.join(vars)

	#_initialize artist object
	fig, ax = plt.subplots(4)
	for i, mw in enumerate([26, 27]):
		i0 = i * 2 
		i1 = i0 + 1
	
		#_setup crap 
		fmt = 'KAOLINITE_SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.{1}_tau-ref.retrieved'
#		fmt = 'GYPSUM_SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.{1}_tau-ref.retrieved'
		types = ['KAL', 'GYP']
		out_label0 = 'BOTHon{0}_{1:d}'.format(types[0], mw)	
		out_label1 = 'BOTHon{0}_{1:d}'.format(types[1], mw)	
		fname0 = os.path.join(dir_hs3, fmt.format(v, out_label0))
		fname1 = os.path.join(dir_hs3, fmt.format(v, out_label1))
	
	
		name = {'KAL' : 'Kaolinite', 'GYP' : 'Gypsum' }
		label0 = '{0} layer'.format(name[types[0]])
		label1 = '{0} layer'.format(name[types[1]])
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

		for idx, ssp in enumerate(r0.keys()):
			if re.search(name[types[0]].lower(), ssp):
				idx0 = ssp #idx	
			elif re.search(name[types[1]].lower(), ssp):
				idx1 = ssp #idx	
	
		#_plot the fractional differences
		tau0 = t0[var_to_plot].reshape(ntau, nref*nzee)
		tau1 = t1[var_to_plot].reshape(ntau, nref*nzee)
		ret0 = r0[idx0][var_to_plot].reshape(ntau, nref*nzee)
		ret1 = r1[idx1][var_to_plot].reshape(ntau, nref*nzee)
		dif0 = (ret0 - tau0) / tau0 
		dif1 = (ret1 - tau1) / tau1
	
		cs0 = ax[i0].boxplot(dif0.T, positions=absc, widths=0.05)
		cs1 = ax[i1].boxplot(dif1.T, positions=absc, widths=0.05)
	
		#_make labels on right
		ax[i0].set_ylabel(lab0, size='x-small')
		ax[i1].set_ylabel(lab1, size='x-small')
		ax[i0].text(absc[-6], .3, label0, size='xx-small') 
		ax[i1].text(absc[-6], .3, label1, size='xx-small') 
	
	[x.set_xlim(absc[0]-0.05, absc[-1]+0.05) for x in ax]
	[x.set_ylim(-.5,.5) for x in ax]
	[shrink_ticks(x) for x in ax]
	ax[-1].set_xlabel('True AOD', size='x-small')
	ax[0].set_title('Fractional Difference'.format(mw),
		size='small') 
	fig.text(0.04, 0.5, '(SHIS-AOD - True) / True',
		va='center', size='x-small', rotation='vertical')
					
	#_save figure
	pname = 'pseudo.both.{0}.{1}.png'.format('-'.join(types), v)
	pname = os.path.join(dir_plot, pname)
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


def fig_10(mw=-27, dir_hs3=DIR_SHIS, vars=['tau','ref','z'], var_to_plot='tau',
	dir_plot='.', **kwargs):
	'''
	fig 07, but with KAL on KAL and GYP on GYP
	'''
	from oe_sim import read_simulated_retrieval
	import matplotlib.pyplot as plt
	from libtools import unique, shrink_ticks

	#_stuff
	v = '_'.join(vars)
	mw = abs(mw)

	#_initialize artist object
	fig, ax = plt.subplots(4)
	name = {'KAL' : 'Kaolinite', 'GYP' : 'Gypsum' }
	for i, mw in enumerate([26, 27]):
		i0 = i * 2 
		i1 = i0 + 1
	
		#_setup crap 
		fmt = '{2}_SHIS.CPL.GDAS.COLLOC.SIMULATED.{0}.{1}.retrieved'
		types = ['KAL', 'GYP']
		out_label0 = '{0}on{0}_{1:d}'.format(types[0], mw)	
		out_label1 = '{0}on{0}_{1:d}'.format(types[1], mw)	
		fname0 = os.path.join(dir_hs3, fmt.format(v, out_label0,
			name[types[0]].upper()))
		fname1 = os.path.join(dir_hs3, fmt.format(v, out_label1,
			name[types[1]].upper()))
	
		label0 = '{0} SSP, {0} layer'.format(name[types[0]])
		label1 = '{0} SSP, {0} layer'.format(name[types[1]])
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
		ax[i0].text(absc[-6], .2, label0, size='xx-small') 
		ax[i1].text(absc[-6], .2, label1, size='xx-small') 
	
	[x.set_xlim(absc[0]-0.05, absc[-1]+0.05) for x in ax]
	[x.set_ylim(-.1, .1) for x in ax]
	[shrink_ticks(x) for x in ax]
	ax[-1].set_xlabel('True AOD', size='x-small')
	ax[0].set_title('Fractional Difference'.format(mw),
		size='small') 
	fig.text(0.04, 0.5, '(SHIS-AOD - True) / True',
		va='center', size='x-small', rotation='vertical')
					
	#_save figure
	pname = 'pseudo.sames.{0}.{1}.png'.format('-'.join(types), v)
	pname = os.path.join(dir_plot, pname)
	dbg(pname)
	plt.savefig(pname)


def fig_11(mw=-27, dir_hs3=DIR_SHIS, vars=['tau','ref','z'], var_to_plot='tau',
	dir_lbl=os.path.join(DIR_PROD, 'LBL-RTM_simulated'), dir_plot='.',**kwargs):
	'''
	same as fig_08 (differences in Bt by MW channel), but with same 
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
		label = '{0} SSP, {1} layer'.format(name[typ], name[typ])
		
		#_read in sim retrieval output
		#_open flight file
		file_seg = fmt_seg.format(v, name[typ].upper())
		file_seg = os.path.join(dir_hs3, file_seg)
		flight = F(file_seg=file_seg)
		w_obs = flight.SHIS_wavenumber
	
		l26 = '{0}on{0}_{1:d}'.format(typ, 26)	
		l27 = '{0}on{0}_{1:d}'.format(typ, 27)	

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
		ax[i26].text(w27[0], -.6, label, size='xx-small') 
		ax[i27].text(w27[0], -.6, label, size='xx-small') 

		#_fix xticks
		ax[i26].set_xlim(w27[0]-15, w27[-1]+15)	
		ax[i27].set_xlim(w27[0]-15, w27[-1]+15)	
		ax[i26].set_ylim(-1, 1)
		ax[i27].set_ylim(-1, 1)

	ax[0].set_title('Habit Sensitivity Test (Kaolinite and Gypsum)',
		size='small')
	[x.plot([-9999,9999], [0,0], 'k--', linewidth=0.4) for x in ax]
	ax[0].xaxis.set_visible(False)
	ax[1].xaxis.set_visible(False)
	ax[2].xaxis.set_visible(False)
	ax[3].set_xlabel('wavenumber ($cm^{-1}$)', size='small')

	[plt.setp(x.xaxis.get_majorticklabels(), rotation=70) for x in ax]
	[shrink_ticks(x) for x in ax]
	fig.text(0.04, 0.5, '$\Delta B_{\\nu}T$ (SHIS - True, $K$)',
		va='center', size='x-small', rotation='vertical')
					
	#_save figure
	pname = 'pseudo.sames.BT.{0}.{1}.png'.format('-'.join(types), v)
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
	pname = 'pseudo.sames.BT.{0}.{1}.png'.format('-'.join(types), v)
	pname = os.path.join(dir_plot, pname)
	dbg(pname)
	plt.savefig(pname)


##def fig_09(idx=[122,242,400,520], dtg='20130824175600', dz_cpl=25.,
def fig_09(idx=[242,400,520], dtg='20130824175600', dz_cpl=25.,
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
##	colors = ['magenta', 'red', '#33ccff', '#33cc33']
	colors = ['red', '#33ccff', '#33cc33']
	labels = ['clear', 'cloud + dust', 'dust']
##	labels = ['clear', 'clear', 'cloud + dust', 'dust']

	#_create list to dump
	dump_dict = {} 
	time_labels = []

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

		#_squirrel away
		time_labels.append(epoch2iso(flight.CPL_epoch[i]))
		dump_dict[epoch2iso(flight.CPL_epoch[i])] = y 

	#_plot differences between dust/cloud 
	if 0:
		r_clr = flight[idx[0]].SHIS_radiances
		r_cld = flight[idx[1]].SHIS_radiances
		r_dst = flight[idx[2]].SHIS_radiances

		y_clr = planck_inv(r_clr, x*100, domain='wavenumber')
		y_cld = planck_inv(r_cld, x*100, domain='wavenumber')
		y_dst = planck_inv(r_dst, x*100, domain='wavenumber')

	#	y_diff = y_dst - y_cld 
		y_diff = y_dst - y_clr
	
		from numpy import abs

		ax_rad_twin = ax_rad.twinx()
		ax_rad_twin.plot(x, abs(y_diff), linewidth=0.3, color='black')

		ax_rad_twin.set_ylim(0,8)

	#_dump BT
	with open('jeff_bt.txt', 'w') as f:
		tmp = ['wavenumber']
		tmp.extend(time_labels)
		f.write(','.join(tmp) + '\n')
		
		#_loop over each value
		for i, v in enumerate(flight.SHIS_wavenumber[:]):
			tmp2 = [dump_dict[l][i] for l in time_labels]
			numb = ['{0:10.5f}'.format(n) for n in tmp2]
			line = ['{0:10.5f}'.format(v)]
			line.extend(numb)
			f.write(','.join(line) + '\n')
 
		

	#_
	ax_cpl.legend(loc=2, fontsize=8) # frameon=False, fontsize=8)
	
	#_limit plot area
	x0 = ax_rad.xaxis.get_data_interval()
	ax_rad.set_ylim(285, 298)
#	ax_rad.set_ylim(270, 300)
#	ax_rad.set_xlim(600, 1300)
	ax_rad.set_xlim(750, 1300)
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


def fig_12(dtg, fidx=None, out_label=None, cloud_threshold=1., 
	dir_plot='.', dz_cpl=25., thresh=1.2, **kwargs):
	''' plot comparison between cpl and shis cloud retrievals '''
	import matplotlib.pyplot as plt
	from hs3_utils import Flight_segment as F
	from libcmap import rgbcmap #_calipso
	from numpy import meshgrid, append, tile, linspace, array
	from libtools import epoch2iso, shrink_ticks

	#_read in flight segment file
	flight = F(dtg=dtg, **kwargs)
	if fidx is None:
		fidx = range(len(flight))
	fidx = array(fidx)
	
	#_get COD data (cod.COD_CPL, cod.COD_SHIS)
	cod = read_cod(flight, dtg, fidx=fidx, out_label=out_label, **kwargs)

	#_find where it's the various combos (shis only, cpl only, both, neither)
	mask_c = cod.COD_CPL > thresh
	mask_s = cod.COD_SHIS > thresh
	both = mask_c * mask_s
	neit = append((mask_c == False)[:,None],(mask_s == False)[:,None], 1).all(1)
	cpll = append((mask_s == False)[:,None], (mask_c == True)[:,None], 1).all(1)
	shis = append((mask_s == True)[:,None], (mask_c == False)[:,None], 1).all(1)
	real = len(fidx)
	total = neit.sum() + cpll.sum() + shis.sum() + both.sum()

	#_initialize draw objects
	fig = plt.figure(figsize=(8,2))
	ax = fig.add_subplot(111)
#	fig, ax = plt.subplots(1)

	#_plot up backscatter
	cpl_plt, cpl_x, cpl_y = flight.map_cpl(dz_cpl=dz_cpl, **kwargs)

	#_generate cpl colormap
	calipso = rgbcmap('calipso')

	#_plot cpl data
	cpl_max = cpl_plt.shape[1]
	cpl_nx, cpl_ny = cpl_plt.shape
	cb = ax.pcolormesh(cpl_x, cpl_y, cpl_plt.T, vmin=-4e-7,
        vmax=1e-4, cmap=calipso, zorder=0)

	#_plot where there is CPL or SHIS
	y0 = tile(cpl_y[-1]-14000, cpl_nx)[shis]
	y1 = tile(cpl_y[-1]-14100, cpl_nx)[cpll]
	y2 = tile(cpl_y[-1]-14200, cpl_nx)[both]
##	y2 = tile(cpl_y[-1]-14200, cpl_nx)[cpll*shis]
	ax.scatter(fidx[shis], y0, marker='x', color='red', s=0.3)
	ax.scatter(fidx[cpll], y1, marker='x', color='yellow', s=0.3)
	ax.scatter(fidx[both], y2, marker='x', color='magenta', s=0.3)
##	ax.scatter(fidx[cpll*shis], y2, marker='x', color='magenta', s=0.3)

	#_setup ticks
	xlim = ax.xaxis.get_data_interval()
	ax.set_xlim(xlim)
#	ax.set_ylim(0, cpl_y[-1]) 
#	ax.set_yticks(linspace(0, cpl_y[-1], 11))
	ax.set_ylim(0, 8000)
	ax.set_yticks(linspace(0, 8000, 6))
	try:
		ax.set_xticks(arange(0, cpl_nx, cpl_nx/5))
		ax.set_xticklabels([epoch2iso(ttt)[-8:] for 
			ttt in flight.CPL_epoch[::cpl_nx/5]])
	except ZeroDivisionError:
		dbg('not sure why this error is happening')

	title = 'CPL_only={0:d}, SHIS_only={1:d}, BOTH={2:d}, NEITHER={3:d},'
	title +=  ' {4:d}/{5:d}'
	title = title.format(cpll.sum(), shis.sum(), both.sum(), neit.sum(),
		total, real)
	ax.set_ylabel('HEIGHTS TOO LOW ~1KM', size='xx-small')
	ax.set_title(title, size='xx-small')

	shrink_ticks(ax)

	#_save figure
	pname = 'cloud-mask_{0}.{1}_{2}.png'.format(flight.dtg0, flight.dtg1,
		out_label) 
	pname = os.path.join(DIR_PLOT, out_label, pname)
	dbg(pname)
	plt.savefig(pname) 


def fig_20(out_labels=['NVACLIMO1_tau-ref_M2_DV26_CLD_CTP-NAAPS_3KM'],
	residuals=None, dir_plot='.', out_label='.', **kwargs):
	'''
	fname	str,	path to output file of optimal estimation
	
	do full campaign comparison between different options

	box whisker version (still implementing)	
	'''
	import matplotlib.pyplot as plt
	from numpy import append, vstack, empty, linspace, tile 
	from hs3_utils import Flight_segment as F
	from libtools import epoch2iso, shrink_ticks
	from numpy import random
	from libtools import shrink_ticks
	import locale
##	import matplotlib
	from oe_real import get_all_residuals


	locale.setlocale(locale.LC_ALL, 'en_US')

	#_initialize artist objects
##	fig, ax = plt.subplots(1)
	fig = plt.figure()
	ax = fig.add_subplot(111)

	#_setup histogram bins
	xmin, xmax = -1, 1
	bins = linspace(xmin, xmax, 101)
	maxy = 0	
	for i, out_label in enumerate(out_labels):
 
		#_initialize residual
		residuals = get_all_residuals(out_label, **kwargs)
		nstr = locale.format("%d", residuals.size, grouping=True)

		#_create histogram
		tmp = ax.hist(residuals, bins=bins, alpha=0.5, label=out_label)
		maxy = max(max(tmp[0]), maxy)

	#_place a line at 0
	ax.plot([0,0], [0,9999], 'k--', linewidth=0.4, zorder=10)
	ax.set_ylim(0, round(maxy, 2) + 50)
	box = ax.get_position()

	#_shrink current axes by ten percent at bottom
	ax.set_position([	box.x0,		box.y0 + box.height * 0.25,
						box.width,	box.height * 0.75])
	
##	ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=8)
	ax.set_xlabel('S-HIS AOT minus CPL AOT', size='xx-small')
	ax.set_ylabel('Fields of View', size='xx-small')
	ax.set_xticks(linspace(xmin, xmax, 11))
	ax.text(-0.9, 200, 'n={0}'.format(nstr), fontsize=6)
##	ax.text(-0.9, 200, 'n={0}\n$\sigma_blue={1}\n\sigma_green={2}$'.format(nstr, total_ok.std(), total_nk.std()), fontsize=6)

	#_shrink	
	shrink_ticks(ax, 12)

	#_name file and save
	pname = 'histogram.shisoe-cpl.summary.{0}.png'.format(len(out_labels))	
	pname = 'histogram.shisoe-cpl.summary.{0}.png'.format('.'.join(out_labels))	
	pname = os.path.join(dir_plot, pname)
	plt.savefig(pname)
	dbg(pname)


def fig_13(out_labels=['NVACLIMO1_tau-ref_M2_DV26_CLD_CTP-NAAPS_3KM'],
	residuals=None, dir_plot='.', out_label='.', **kwargs):
	'''
	fname	str,	path to output file of optimal estimation
	
	do full campaign comparison between different options
	
	out_label kwargs is pulled out to not interfere with get_all_res	
	'''
	import matplotlib.pyplot as plt
	from numpy import append, vstack, empty, linspace, tile, recarray 
	from hs3_utils import Flight_segment as F
	from libtools import epoch2iso, shrink_ticks
	from numpy import random
	from libtools import shrink_ticks
	import locale
	import matplotlib
	from oe_real import get_all_residuals

	locale.setlocale(locale.LC_ALL, 'en_US')

	#_initialize artist objects
	fig, ax = plt.subplots(1)
	fig.suptitle('Residual Optical Depths when Aerosol Layer found in CPL',
		size='small')
	#_setup histogram bins
	xmin, xmax = -1, 1
	bins = linspace(xmin, xmax, 101)
	maxy = 0	
	for i, out_label in enumerate(out_labels):
 
		#_initialize residual
		residuals = get_all_residuals(out_label,**kwargs)
		nstr = locale.format("%d", residuals.size, grouping=True)

		#_create histogram
		tmp = ax.hist(residuals, bins=bins, alpha=0.5, label=out_label)
		maxy = max(max(tmp[0]), maxy)

	#_place a line at 0
	maxy = round(maxy, -1) + round(maxy/9, -1) 
	ax.plot([0,0], [0,9999], 'k--', linewidth=0.4, zorder=10)
	ax.set_ylim(0, maxy)

	box = ax.get_position()

	#_shrink current axes by ten percent at bottom
	ax.set_position([	box.x0,		box.y0 + box.height * 0.25,
						box.width,	box.height * 0.75])
	
	ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=8)
	ax.set_xlabel('S-HIS AOT minus CPL AOT', size='xx-small')
	ax.set_ylabel('Fields of View', size='xx-small')
	ax.set_xticks(linspace(xmin, xmax, 21))
	ax.text(-0.9, 200, 'n={0}'.format(nstr), fontsize=6)
##	ax.text(-0.9, 200, 'n={0}\n$\sigma_blue={1}\n\sigma_green={2}$'.format(nstr, total_ok.std(), total_nk.std()), fontsize=6)

	#_shrink	
	shrink_ticks(ax)

	#_name file and save
	pname = 'histogram.shisoe-cpl.summary.{0}.png'.format(len(out_labels))	
	pname = 'histogram.shisoe-cpl.summary.{0}.png'.format('.'.join(out_labels))	
	pname = os.path.join(dir_plot, pname)
	plt.savefig(pname)
	dbg(pname)

def match(d, comparison_sensor, diff_type='aerosol'):
		''' derp
		comparison sensor		CPL or MODIS

		diff type,		aerosol or cloud. Sorts through 
						the similar habits and what not
		'''
		from numpy import array

		shis = d[d.sensor == 'S-HIS']
		comp = d[d.sensor == comparison_sensor]

		e_shis = array(list(set(shis.epoch)))
		e_comp = array(list(set(comp.epoch)))

		epoch = e_shis if e_shis.size < e_comp.size else e_comp

		pairs = []
		for e in epoch:
			sidx = shis.epoch == e
			cidx = comp.epoch == e

			s = shis[shis.epoch == e]
			c = comp[comp.epoch == e]

			if c.size == 0:
				raise RuntimeError, 'NOTHING TO COMPARE TO'

			if diff_type == 'aerosol':
				#_sum up all the aerosol types
				shis_tau = 0.
				for layer in s:
					if layer.ssp in ['mie_quartz', 'mie_kaolinite','mie_gypsum']:
						shis_tau += layer.tau
				
				comp_tau = 0.
				for layer in c:
					print layer
					if layer.ssp == 'aerosol':
						comp_tau += layer.tau

			elif diff_type == 'cloud':
				#_sum up all the aerosol types
				shis_tau = 0.
				for layer in s:
					if layer.ssp in ['mie_wat', 'mie_ice']: 
						shis_tau += layer.tau
				
				comp_tau = 0.
				for layer in c:
					if layer.ssp == 'cloud':
						comp_tau += layer.tau
			
			if comp_tau == 0 or comp_tau == -9999:
				continue
	
			print 'NON-ZERO', shis_tau, comp_tau		
			pairs.append((shis_tau - comp_tau) / comp_tau)	

		return pairs


def fig_CAMPAIGN2(out_labels=['NVACLIMO1_tau-ref_M2_DV26_CLD_CTP-NAAPS_3KM'],
	bins=[[0,0.15], [0.15, 0.3], [0.3, 0.45], [0.45, 0.6], [0.6, 0.8], [0.8,1.],
		[1, 1.25], [1.25, 1.5], [1.5,1.75], [1.75,2.0], [2,3.]],
	residuals=None, dir_plot='.', out_label='.', comparison='CPL', **kwargs):
	'''
	fname	str,	path to output file of optimal estimation
	
	do full campaign comparison between different options
	
	out_label kwargs is pulled out to not interfere with get_all_res	
	'''


	import matplotlib.pyplot as plt
	from numpy import append, vstack, empty, linspace, tile, recarray 
	from numpy import random, mean
	from hs3_utils import Flight_segment as F
	from flight_namelists import experiments
	from libtools import epoch2iso, shrink_ticks
	from libtools import shrink_ticks
	import locale
	import matplotlib
	from oe_real_hs3 import get_all_residuals, read_real_retrieval
	from libtools import merge


	bins = linspace(0,3,31)
	bins = zip(bins[:-1], bins[1:])

	segments = experiments['HS3']

	locale.setlocale(locale.LC_ALL, 'en_US')

	#_initialize artist objects
	np = len(out_labels)
	fig, ax = plt.subplots(np)
	fig.suptitle('Residual Optical Depths when Aerosol Layer',
		size='small')

	#_initialize plot data, fill with fields to be passed to plot.boxplot
	y = [None] * len(bins)

	#_get mean values within bins
	x = mean(bins, 1)

	#_read all AOD data into recarray for out_labels
	for i, out_label in enumerate(out_labels):
		
		#_list of the number of obs in each bin
		n = []	
	
		#_read in data
		if comparison == 'MODIS':
			modis = True
			cpl = False
		elif comparison == 'CPL':
			modis = False
			cpl = True

		data = read_retrieval_campaign(modis=modis, cpl=cpl, **kwargs)
		shis = data[data.sensor == 'S-HIS']
		
		#_loop over bins, build data
		for b, (bmin, bmax) in enumerate(bins):

			
			#_shis/cpl/modis recarrays should have same size and order
			if comparison == 'MODIS':
				idx = (modis.aod >= bmin) * (modis.aod < bmax)
				idx = (modis.aod > 0) * idx
				y[b] = (shis.aod[idx] - modis.aod[idx]) / modis.aod[idx]	

			elif comparison == 'CPL':
				idx = (cpl.aod >= bmin) * (cpl.aod < bmax)
				idx = (cpl.aod > 0) * idx
				y[b] = (shis.aod[idx] - cpl.aod[idx]) / cpl.aod[idx]

			#_record number of obs used
			n.append(idx.sum())
	
		#_select plot axis to use
		axis = ax[i]

		#_make box plot
		axis.plot([0, len(bins)+2], [0, 0]) 
		axis.boxplot(y, positions=x, widths=0.1)
		axis.set_xlim(0, bins[-1][-1])
		axis.set_ylim(-5, 5)
		axis.text(1.5, 1.5, out_label, size='xx-small')

		#_include the number of retrievals along top
		axis_top = axis.twiny()
		axis_top.set_xticks(x)
		axis_top.set_xticklabels(n)

		#_shrink	
		shrink_ticks(axis)
		shrink_ticks(axis_top)

#	box = ax.get_position()
#
#	#_shrink current axes by ten percent at bottom
#	ax.set_position([	box.x0,		box.y0 + box.height * 0.25,
#						box.width,	box.height * 0.75])
#	
#	ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=8)
#	ax.set_xlabel('S-HIS AOT minus CPL AOT', size='xx-small')
#	ax.set_ylabel('Fields of View', size='xx-small')
#	ax.set_xticks(linspace(xmin, xmax, 21))
#	ax.text(-0.9, 200, 'n={0}'.format(nstr), fontsize=6)

	#_name file and save
	pname = 'histogram.shisoe-cpl-modis.{1}.{0}.png'.format(
		'.'.join(out_labels), comparison)	
	pname = os.path.join(dir_plot, pname)
	plt.savefig(pname)
	dbg(pname)


def fig_CAMPAIGN(out_labels=['NVACLIMO1_tau-ref_M2_DV26_CLD_CTP-NAAPS_3KM'],
	bins=[[0,0.15], [0.15, 0.3], [0.3, 0.45], [0.45, 0.6], [0.6, 0.8], [0.8,1.],
		[1, 1.25], [1.25, 1.5], [1.5,1.75], [1.75,2.0], [2,3.]],
	residuals=None, dir_plot='.', out_label='.', comparison='CPL', **kwargs):
	'''
	fname	str,	path to output file of optimal estimation
	
	do full campaign comparison between different options
	
	out_label kwargs is pulled out to not interfere with get_all_res	
	'''


	import matplotlib.pyplot as plt
	from numpy import append, vstack, empty, linspace, tile, recarray 
	from numpy import random, mean
	from hs3_utils import Flight_segment as F
	from flight_namelists import experiments
	from libtools import epoch2iso, shrink_ticks
	from libtools import shrink_ticks
	import locale
	import matplotlib
	from oe_real_hs3 import get_all_residuals, read_real_retrieval
	from libtools import merge

	bins = linspace(0,3,31)
	bins = zip(bins[:-1], bins[1:])

	segments = experiments['HS3']

	def rmod(flight, **kwargs):
		''' read in modis '''
		from numpy.ma import masked_where
		from netCDF4 import Dataset
	
		#_read in flight data
		fmod = flight.FLIGHT_file.replace('CPL.GDAS','MODIS')
		with Dataset(fmod, 'r') as cdf_mod:

			#_mask missing
			aod_mod = cdf_mod.variables['aod'][:]
			dis_mod = cdf_mod.variables['distance'][:]
			aod_mod = masked_where(aod_mod <= 0, aod_mod)
			aod_mod = masked_where(dis_mod > 3e4, aod_mod)

			dtype = [('aod','f4'), ('epoch', 'f8')]
			modis = recarray((aod_mod.size), dtype)
			modis.aod[:] = aod_mod[:]
			modis.epoch[:] = flight.SHIS_epoch[:] 
		##	modis.epoch[:] = cdf_mod.variables['epoch'][:]	

		return modis

	def rshis(d, flight, **kwargs):
		''' clean up retrieval read '''
		dtype = [('aod', 'f4'), ('epoch', 'f8')] 

		if d.size != flight.size:
			raise RuntimeError, 'What {0:d} {1:d}'.format(d.size, flight.size)

		dtype = [('aod','f4'), ('epoch', 'f8')]
		shis = recarray((d.size), dtype)
		shis.aod[:] = d.tau	
		shis.epoch[:] = flight.SHIS_epoch[:]	

		#_Maybe try to work in cloud optical depth?
		# or work that into generate_retrieved or whatever it is called
		return shis

	def rcpl(d, flight, **kwargs):
		''' clean up read of CPL '''
		dtype = [('aod', 'f4'), ('epoch', 'f8'), ('cod', 'f4')]

		#_put CPL data into arrays  
		cpl_typ = vstack(d.CPL_type)
		cpl_tau = vstack(d.CPL_tau)

		#_trim down to number of layers possibly defined in type
		cpl_tau = cpl_tau[:,:5]
		cpl_typ = cpl_typ[:,::2]

		#_find fovs with any aerosol
		idx_pbl = (cpl_typ == 1).any(1) #_index of PBL type
		idx_ele = (cpl_typ == 2).any(1) #_index of AEROSOL type
		idx_aer = append(idx_pbl[:,None], idx_ele[:,None], axis=1).any(1)
	##	idx_cld = (cpl_typ == 3).any(1) * idx_aer   #_aerosol w/ cloud
	##	idx_clr = (cpl_typ != 3).all(1) * idx_aer   #_aerosol w/o cloud

		#_get column totals of elevated and pbl aersol when not missing
		col_pbl = ((cpl_typ == 1) * (cpl_tau)).sum(1) #max(1) #sum(1)
		col_ele = ((cpl_typ == 2) * (cpl_tau)).sum(1) #max(1) #sum(1)
		cpl_aod = col_ele + col_pbl
		
		cpl = recarray((cpl_aod.size,), dtype)
		cpl.aod[:] = cpl_aod[:]
		cpl.epoch[:] = flight.SHIS_epoch[:]

		return cpl
		
	locale.setlocale(locale.LC_ALL, 'en_US')

	#_initialize artist objects
	np = len(out_labels)
	fig, ax = plt.subplots(np)
	fig.suptitle('Residual Optical Depths when Aerosol Layer',
		size='small')

	#_initialize plot data, fill with fields to be passed to plot.boxplot
	y = [None] * len(bins)

	#_get mean values within bins
	x = mean(bins, 1)

	#_read all AOD data into recarray for out_labels
	for i, out_label in enumerate(out_labels):
		
		#_list of the number of obs in each bin
		n = []	
	
		for dtg, values in segments.iteritems():
			flight = F(dtg=dtg, **kwargs)

			#_generate input file name
			fname = flight.FLIGHT_file.replace('.nc',
				'.{0}.retrieved'.format(out_label))

			#_read in data
			data = read_real_retrieval(fname, **kwargs)
		
			if data.size == 0:
				dbg('NO DATA FOR {0}'.format(dtg))
				continue

			try:
				#_read in SHIS aod
				shis = merge((shis, rshis(data, flight))) 

				#_read in CPL aod
				cpl = merge((cpl, rcpl(data, flight)))

				#_read in MODIS aod
				modis = merge((modis, rmod(flight, **kwargs)))
 
			except UnboundLocalError:
				#_read in SHIS aod
				shis = rshis(data, flight)

				#_read in CPL aod
				cpl = rcpl(data, flight)

				#_read in MODIS aod
				modis = rmod(flight, **kwargs)
				

		#_loop over bins, build data
		for b, (bmin, bmax) in enumerate(bins):
			
			#_shis/cpl/modis recarrays should have same size and order
			if comparison == 'MODIS':
				idx = (modis.aod >= bmin) * (modis.aod < bmax)
				idx = (modis.aod > 0) * idx
				y[b] = (shis.aod[idx] - modis.aod[idx]) / modis.aod[idx]	

			elif comparison == 'CPL':
				idx = (cpl.aod >= bmin) * (cpl.aod < bmax)
				idx = (cpl.aod > 0) * idx
				y[b] = (shis.aod[idx] - cpl.aod[idx]) / cpl.aod[idx]

			#_record number of obs used
			n.append(idx.sum())
	
		#_select plot axis to use
		axis = ax[i]

		#_make box plot
		axis.plot([0, len(bins)+2], [0, 0]) 
		axis.boxplot(y, positions=x, widths=0.1)
		axis.set_xlim(0, bins[-1][-1])
		axis.set_ylim(-5, 5)
		axis.text(1.5, 1.5, out_label, size='xx-small')

		#_include the number of retrievals along top
		axis_top = axis.twiny()
		axis_top.set_xticks(x)
		axis_top.set_xticklabels(n)

		#_shrink	
		shrink_ticks(axis)
		shrink_ticks(axis_top)

#	box = ax.get_position()
#
#	#_shrink current axes by ten percent at bottom
#	ax.set_position([	box.x0,		box.y0 + box.height * 0.25,
#						box.width,	box.height * 0.75])
#	
#	ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.1), fontsize=8)
#	ax.set_xlabel('S-HIS AOT minus CPL AOT', size='xx-small')
#	ax.set_ylabel('Fields of View', size='xx-small')
#	ax.set_xticks(linspace(xmin, xmax, 21))
#	ax.text(-0.9, 200, 'n={0}'.format(nstr), fontsize=6)

	#_name file and save
	pname = 'histogram.shisoe-cpl-modis.{1}.{0}.png'.format(
		'.'.join(out_labels), comparison)	
	pname = os.path.join(dir_plot, pname)
	plt.savefig(pname)
	dbg(pname)
		

def read_cod(flight, dtg, fidx=None, out_label=None, thresh=0.3, 
	dir_lbl=os.path.join(DIR_PROD, 'LBL-RTM_hs3'), experiment='HS3', **kwargs):
	''' 
	Read in cloud optical thicknesses for CPL and SHIS
	Plot 532 backscatter of CPL
	Have two scatter dots respresenting if CPL and SHIS identified cloud layers

	'''
	from lblrtm_utils import Lbldis_input
	import re
	from numpy import recarray

	#_get flight information
	d0, d1 = flight.dtg0, flight.dtg1

	#_build base date dir structure
	dir_dtg = '{0}_{1}.{2}_fov{{0:06d}}'.format(experiment, d0, d1)

	#_create regular expression for water cloud ssp file
	reg = re.compile('ssp_db.mie_wat.gamma_sigma_0p100')

	#_initialize recarray
	dtype = [('epoch', 'f8'), ('fidx', 'i4'), 
			('COD_SHIS', 'f4'), ('COD_CPL', 'f4')]

#	#_get full length of flight if no indices passed
#	if fidx is None:
#	#	from hs3_utils import Flight_segment as F
	nf = len(flight) if fidx is None else len(fidx)
	data = recarray((nf,), dtype)

	def _get_cpl_cod(flight, fidx):
		from numpy.ma import masked_where
		from numpy import array

		''' pull out field of views cloud optical depth from CPL '''
		fov = flight[f]

		#_create string of text containing CPL aod values for layers
		tau_532     = fov.CPL_tau_532[:]
		tau_type    = fov.CPL_tau_type[::2]
	#	tau_type    = fov.CPL_tau_type[:]
		tau_type_str = [None, 'PBL', 'AEROSOL', 'CLOUD']
		tau_532 = masked_where(tau_532 <= 0, tau_532)

	#	tau_typ = array(zip(tau_532, tau_type))[tau_532.mask == False]

		#_go down where not masked.... I guesss?
		cod = 0.
		lll = min(10-tau_532.mask.sum(), 5) #_found cases with more OD vals than types...
		#_HS3_20130829133102.20130829182526_fov000384/lbldis_input.NVACLIMO1_tau-ref_M2_DV26_CLD_CTP-NAAPS_3KM.fina
		for ii in xrange(lll):
		#for ii in xrange(10 - tau_532.mask.sum()):
			try:
				if tau_type_str[tau_type[ii]] == 'CLOUD':
					cod += tau_532[ii] 
			except IndexError:
				print 'FAILED ON THIS'
				print ii
				print tau_532
				print tau_type
				print fov.CPL_tau_532
				print fov.CPL_tau_type
				raise IndexError, 'FAILED ON THIS'
	
	#	for tau, typ in tau_typ:
	#		print typ
	#		if tau_type_str[int(typ)] == 'CLOUD':
	#			cod += tau
	
		return cod

	#_loop over fields of view
	for i, f in enumerate(fidx):
		#_pull out flight time
		epoch = flight[f].CPL_epoch

		#_read in lbldis_input for this file
		fname = 'lbldis_input.{0}.final'.format(out_label)	
		lbl_input = os.path.join(dir_lbl, dir_dtg.format(f), fname) 

		#_read in layers
		lbl = Lbldis_input(lbl_input)

		cod = 0.
		for j, layer in enumerate(lbl.layers):
			#_pull out current layers ssp
			ssp = lbl.ssp_db[layer['dbnum']]
			
			#_if it's the water cloud, include its OD in total
			if reg.search(ssp):
				cod += sum(layer['tau'])

		#_add values for cod to record array for SHIS and CPL
		data[i] = (epoch, f, cod, _get_cpl_cod(flight, f))	

	return data


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
	
