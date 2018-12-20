def plot_simulated_retrieval_by_var(var, out_label='', plot_var='tau',
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

	if len(files) == 0:
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

