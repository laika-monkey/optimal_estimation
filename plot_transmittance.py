#!/usr/bin/env python
################################################################################
#_script:	plot_transmittance.py
# purpose:	Recreates Petty's fig. 7.6 from Intro to Radiation 
# author:	Walter R Sessions
################################################################################

from numpy import array
import os

namelist = {
	#_name of transmission file
	'fname'	: os.path.join(os.environ['PRODUCTS'], 'zenith_2.trans'),

	#_range of plotted wavelength
##	'xlim'		: (0,50),
##	'xlim'		: (9,13),	#_longer window channel + o3

	#_labels for wavelength
#	'lab_wl'	: array([.3,.4,.5,.6,.8,1,1.2,1.5,2,2.5,3,
#					4,5,6,7,8,9,10, 12,15,20,25,40,50]),
	'lab_wl'	: array([  8.62589494,   8.73286176,   8.86132034,  10.39392994,
						10.69919221,  11.08893324,  11.18755943,  11.43053095,
						12.02718143,  12.19735317,  12.33121647]),
	#_save figure to file or show?
	'save'		: True,  
	'pname'		: 'TIR_O3',

	#_which increases left to right, wavelength or number
	#_DO NOT USE
	'abscisa'	: 'wavenumber',

	#_skip any constituents? (use lowercase total, h20, co2, o3, o2, ch4, n2o,co
	'skip'		: ['co','n2o','o2', 'ch4'],

	#_add planck functions in the background for terrestrial or solar radiation
	'add_blackbody'	: {
		'earth'	: {
			'T'		: 300.,
			'color'	: 'blue',
			},
		'solar'	: {
			'T'		: 5778.,
			'color'	: 'orange',
			},
		},
	}


#############################################################################_80
#_MAIN_#########################################################################
################################################################################


def run_main(xlim=(0,50), lab_wl=None, add_blackbody=False,abscisa='wavelength',
	save=False,	pname='zenith_transmittance', **kwargs):
	'''
	xlim	tup,	range of wavelengths to plot
	'''
	import matplotlib.pyplot as plt

	#_read data
	data = read_transmittance(**kwargs)

	#_pull out variables names from dtype and drop wlength
	varnames = list(data.dtype.names)
	del varnames[varnames.index('wl')]
	nvar = len(varnames)
	varnames.sort()	#_for now, works to throw 'total' at end

	#_initialize artist
	fig, ax = plt.subplots(nvar, 1)
	x = data.wl

	#_initialize lists of twinned axes
	aT = []

	#_loop over variables and plot
	for i, name in enumerate(varnames):
		#_create a text label for variable name
		varname_label, lab_mod = build_label(name)

		#_get transmittance values
		y = getattr(data, name)

		#_copy right and left
		aT.append(ax[i].twiny())

		#_make plot of transmittance and dummy plot to extend top axes
		ax[i].plot(x, y, linewidth=0.3)
		aT[i].plot(x,[-1]*x.size)

		#_add outgoing blackbody radiances
		for source, opts in add_blackbody.iteritems():
			planck_plot(x, ax[i], **opts)

		#_put x axis on log scale
		ax[i].set_xscale('log')
		aT[i].set_xscale('log')

		#_set yaxis labels
		yhgt = 0.6 if lab_mod else 0.775
		ax[i].text(xlim[1]-2, yhgt, varname_label, #name.upper(), 
				verticalalignment='center',
				horizontalalignment='right')
		ax[i].set_ylim(0, 1)
		ax[i].set_yticks([0,.5,1])

		#_make wavenumber grid on top for all plots, but only visible on top 
		ax[i].set_xlim(*xlim)
		aT[i].set_xlim(*xlim)

		#_turn on grid
		ax[i].grid(True)
		aT[i].grid(True)

		#_use only labels within range
		lab_wl = lab_wl[(lab_wl >= xlim[0]) * (lab_wl <= xlim[1])]

		#_set location of wavel/n xticks
		ax[i].set_xticks(lab_wl)
		aT[i].set_xticks(lab_wl)
		ax[i].set_xticklabels([])
		aT[i].set_xticklabels([])

		#_make top x axis visible if 0
		if not i:
			#_calculate labels and put into format
			lab_wn = ['{0:5.0f}'.format(n*1e6) \
					for n in wavelength2number(lab_wl)]

			#_set locations and converted labels
			aT[i].set_xticklabels(lab_wn, rotation=90)
			aT[i].set_xlabel('wavenumber ($cm^{-1}$)', size='x-small')

		#_make bottom x axis visible if nvar-1
		elif i == nvar-1:
			#_trim up values if they're greater than 1
			labels = []
			for label in lab_wl:
				if label > 3:
					labels.append('{0:3.1f}'.format(label))
				else:
					labels.append('{0:3.1f}'.format(label))

			#_set location and labels of xticks
			ax[i].set_xticks(lab_wl)
			ax[i].set_xticklabels(labels, rotation=90)
			
			#_set dimensional label
			ax[i].set_xlabel('wavelength (${\mu}m$)', size='x-small')

		#_set ylabel around middle
		if i == nvar/2:
			ax[i].set_ylabel('transmittance', size='medium')

	'''
	if abscisa == 'wavenumber':
		[axis.invert_xaxis() for axis in fig.axes]
	'''

	#_reduce size of all ticks
	[shrink_ticks(axis, 5) for axis in fig.axes]
	
	#_rotate labels
	if save:
		pname = '{0:s}.png'.format(pname)
		print 'saving', pname
		plt.savefig(pname)
	else:
		plt.show()


def shrink_ticks(ax, size=6):
	''' reduce font size of ticks '''
	for t in ax.yaxis.get_major_ticks():
		t.label1.set_fontsize(size)
	
		#_if twinned
		if hasattr(t, 'label2'):
			t.label2.set_fontsize(size)

	for t in ax.xaxis.get_major_ticks():
		t.label1.set_fontsize(size)

		#_if twinned
		if hasattr(t, 'label2'):
			t.label2.set_fontsize(size)


''' pass l in microns, returns in cm-1 '''
def wavelength2number(l): return (l*100)**-1 


def read_transmittance(fname='/Users/wsessions/Data/zenith_2.trans',
	skip=False, **kwargs):
	'''
	Read in space delimited transmission file used for Petty fig. 7.6

	arguments
	fname	str,	path to input file

	returns
	numpy recarray with column headers as attributes
	'''

	from numpy import loadtxt, recarray

	#_get variable names
	with open(fname,'r') as f:
		variables = f.readline().split()

	#_remove variables if not desired
	try:
		for varname in skip:
			del variables[variables.index(varname)]
	except:
		pass

	#_read in data
	data = loadtxt(fname, skiprows=1)

	#_build dtype
	dtype = [(name, 'f8') for name in variables]
	out = recarray((data.shape[0],), dtype) 

	#_load into recarray
	for i, name in enumerate(variables):
		out.__getattribute__(name)[:] = data[:,i]

	return out 


def planck_plot(wl, axis, T=300, color='k', **kwargs):
	''' add normalized outgoing radiances to plot area '''
	from numpy import append, exp
	def planck(t, wl):
		C = 2.99792458e8  #_Speed of Light        [m/s]
		H = 6.6260957e-34 #_Planck's Constant     [m2kg/s]
		K = 1.381e-23     #_Boltzmann's Constant  [J/K]
		e = exp(H*C/K/wl/t)
		return 2*H*C*C/wl**5/(e-1.)

	#_create normalized values
	bt = planck(T, wl*1e-6)
	bt = bt / bt.max()

	#_add last point as first to fill
	x = append(0, wl)
	x = append(x, 0)
	y = append(0, bt)
	y = append(y, 0)

	#_plot fill region
	axis.fill(x, y, color, alpha=0.1)


def build_label(name):
	'''
	make prettier label names
	
	returns
	str,	name but with all numbers dropped to subscripts
	bool,	True if subscript added
	'''
	lab_mod = False

	varname_label = []
	for c in name:
		if c.isdigit():
			varname_label.append('$_{0:s}$'.format(c))
			lab_mod = True
		else:
			varname_label.append(c.upper())
			
	return ''.join(varname_label), lab_mod


if __name__ == '__main__':
	import os, sys
	if len(sys.argv) != 3:
		print 'usage: ./plot.py <start_wavelength> <end_wavelength>'
		os._exit(0)

	namelist['xlim'] = tuple([float(a) for a in sys.argv[1:]])
	run_main(**namelist)
