#!/usr/bin/env python
#############################################################################_80
# Odds and ends for vimont's class. Individual problems should be handled      #
#	externally to this in a script.					       #
################################################################################

import os, sys

kwargs = {


	'show'	: True,

	#_AERONET sites
	'sites'	: ['capove'], #_crozet, baengn, kanpur
	
	}

'''
def read_grib(fname, **kwargs):
	import pygrib as pyg
	
	var = '10 metre U wind component'
	grb = pyg.open(fname)
	
	#_pull out list of object pointing to whatever
	grb.select(name=var)
	
	#_
	grb[0].keys()	#_attributes
	grb[0].values 	#_data
'''

def read_exp2(aod_thresh=0.5, **kwargs):
	from libtools import epoch2dtg, dtg2epoch, subset, unique
	from numpy import vstack
	
	#_read icap data
	icap = read_forecasts(**kwargs)
	
	#_read GDAS data at same times
	dust = subset(icap, variable='dust_aod', fhr=24)
	
	#_drop missing values
	iiii = vstack(dust.values).mean(axis=1)
	dust = dust[iiii>0.5]
	
	#_get gdas
	dtgs = epoch2dtg(dust.epoch)
	gdas = read_GDASGRIB(dtgs=dtgs, **kwargs)

	#_limit to when times are same
	dust = subset(dust, epoch=dtg2epoch(gdas.dtg))
	
	return dust, gdas
	
def exp2(vars=['geopot','uwind'], data=None, thresh=0.8,
	**kwargs):
	from numpy import vstack, arange, linspace, sqrt, meshgrid
	import matplotlib.pyplot as plt
	
	if data == None:
		dust, gdas = read_exp2(**kwargs)
	else:
		dust, gdas = data
	
	#_stack and reshape g??
	lat = linspace(-90, 90, 181)
	lon	= linspace(-180, 180, 360)
	lon2d, lat2d = meshgrid(lon, lat)

	#_stack
	dust 	= vstack(dust.values)
	nt 		= dust.shape[0]
	nx, ny 	= 360, 181 #_no time
	idx		= arange(nt)
	
	#_remove mean
	dust -= dust.mean(axis=0)
	
	for var in vars:
		#_init figure
		fig = plt.figure()
		nax = 321
		
		#_remove temporal mean
		if var == 'uwind':
			u = vstack(gdas.uwind).reshape(nt, ny, nx)
			v = vstack(gdas.vwind).reshape(nt, ny, nx)
			g = sqrt(u**2 + v**2)
		elif var == 'vwind':
			continue
		else:
			g = vstack(gdas.__getattribute__(var)).reshape(nt, ny, nx)
			
		#_get anomaly
		g -= g.mean(axis=0)
	
		#_do icap first
		fcst = dust.mean(axis=1)
		reg_map = regression_map(fcst, g, **kwargs)
		
		#_get significance
		pidx = idx[fcst >= thresh]
		nidx = idx[fcst < thresh]
		H, diff = composite(g, pidx, nidx, **kwargs)
		
		plot_reg(reg_map, H, lat2d, lon2d, label='Test', figure=fig, nax=nax)
		
		#_loop over model
		models = ['NAAPS','GEOS-5','MACC','MASINGAR','NGAC']
		for nmod in xrange(5):
			nax += 1

			fcst 	= dust[:,nmod]
			reg_map = regression_map(fcst, g, **kwargs)
			
			#_get significance
			pidx = idx[fcst >= thresh]
			nidx = idx[fcst < thresh]
			H, diff = composite(g, pidx, nidx, **kwargs)

			plot_reg(reg_map, H, lat2d, lon2d, label=var.upper() + ' ' 
				+ models[nmod], figure=fig, nax=nax)
			
		plt.show()

def eof(d, full_matrices=False, centered=False, **kwargs):
	'''
	given a dataset, returns eigenvectors, values, pcs and percent var 
	assumes data is a two dimensional matrix in sample and state space
	d	np.ndarray,		ndimensional array with sample space in the 0th axis
	
	_OUTPUT
	eof	np.ndarray,		EOFS will be in ROWS (NxN)
	lam	np.ndarray,		eigenvalues will be in a vector (N)
	pcs	np.ndarray,		Principal components will be in ROWS (NxM)
	per	np.ndarray,		Percent variance will also be in a vector
	'''
	from numpy import diag, matrix, dot, sqrt, product, zeros, array, ones
	from numpy.linalg import eig, svd
	
	multi = True if d.ndim != 2 else False
	if multi:	#_resize to 2d matrix
		orig_shape = d.shape[:]
		d = d.reshape( orig_shape[0], product(orig_shape[1:]) )

	m,n	= d.shape
	mm 	= max(m, n)
	nn 	= min(m, n)
	d 	= matrix(d)
	
	if not centered: d -= (1./n)*matrix(ones((m, n)))		#_remove mean

	#_apply eof analysis (s is the sqrt of eigenvalues)
	l, s, r = svd(d) #,full_matrices=False)		#_using SVD
	tmp		= zeros( mm )						#_create full array eigenvalues
	tmp[:nn]= s**2 / (m-1.)						#_fill with what we have
	eigL	= tmp[:n]
	eigV	= r.T								#_right matrix col eigenvectors
	pca		= dot(l, diag(tmp)[:m,:n])			#_left matrix col holds the pcs

##	#_using EIG (HASNT BEEN UPDATED FOR ROTATION)
##	#_get covariance matrix of data
##	covar	= np.dot( d.T, d ) / (nt-1.)
##	eigL, eigV 	= eig( covar )				#_get eigValues and eigVectors
##	idx 	= np.argsort(eigL)[::-1]		#_get eigenvalues in order
##	eigL 	= eigL[idx]						#
##	eigV 	= eigV[:,idx]					#_and eigenvectors
##	pca 	= np.dot( dataN, eigV )			#_calculate principal comps

	#_calc variance explained by each EOF
	percent = eigL / eigL.sum() * 100.

	eof = dot(eigV, diag(sqrt(eigL))).T		#_get eofs and put vectors in rows
	pcs = (pca / sqrt(eigL)).T				#_calculate principle components
	eof, pcs = array(eof), array(pcs)		#_just in case

	#_return spatial component back to original shape
	if multi: eof = array(eof[:m]).reshape(orig_shape)

	return eof, eigL, pcs, percent
	

def eof_anal(data, show=False, label='', **kwargs):
	'''
	eof/pca to synthetic dataset
	'''
	import numpy as np
	import matplotlib.pyplot as plt
	from numpy.random import randn
	import matplotlib.cm as cm
	from numpy.linalg import svd, eig
	from numpy import dot, product, arange
	import libtools as lt
	
	#_if this is a 3d dataset, blorgh
	if data.ndim > 2:
		shape = data.shape
		data  = data.reshape( shape[0], product( shape[1:]) )
		multi = True
	else:
		multi = False	
		
	#_create fake data
	nt, nx	= data.shape
	x, t	= arange(nx), arange(nt)
	
	#_get covariance matrix of data
	covar	= dot( data.T, data ) / (nt-1.)

	eigV, eigL, pca, percent = eof(data)

	#_create fake data
	nt, nx	= data.shape
	x, t	= np.arange(nx), np.arange(nt)
	
	#_apply eof analysis (s is the sqrt of eigenvalues)
	#_ using SVG
##	l, s, r = svd(dataN, full_matrices=False)
##	eigV	= r.T.copy()
##	eigL	= s**2 / (nt-1.)
##	pca		= np.dot(l, np.diag(eigL))

##	#_using EIG
##	eigL, eigV 	= eig( covar )				#_get eigValues and eigVectors
##	idx 	= np.argsort(eigL)[::-1]		#_get eigenvalues in order
##	eigL 	= eigL[idx]						#
##	eigV 	= eigV[:,idx]					#_and eigenvectors
##	pca 	= np.dot( dataN, eigV )			#_calculate principal comps

	#_calc variance explained by each EOF
	percent = eigL / eigL.sum() * 100.
	lamberr = eigL * np.sqrt(2./(nt-1))

	#_initialize draw object
	fig = plt.figure( )
	arg = { 
		'extend': 'both', 
		'cmap' 	: cm.RdBu_r, 
		'levels': np.linspace(-1,1,11) }
	shape = (5,2)

	#_make unimportant eigenvectors plot lightly
	lw = np.ones((percent.size)) * 0.5
	lw[percent < 10.] = 0.2
	
	if not multi:
		#_create hovmuller diagram of synthetic data
		ax0 = plt.subplot2grid( shape, (0,0), rowspan=2 )
		ap0 = ax0.contourf( x, t, data, **arg )
		ax0.set_title( label ); ax0.set_ylabel('t')
		plt.colorbar( ap0 )
	
		#_create hovmuller diagram of synthetic data
		ax1 = plt.subplot2grid( shape, (0,1), rowspan=2 )
		ap1 = ax1.contourf( x, t, dataN, **arg )
		ax1.set_title( 'with white noise' )
		plt.colorbar( ap1 )
	
		#_plot the eigenspectrum with error bars (100 indep samples)
		ax2 = plt.subplot2grid( shape, (2,0), colspan=2 )
	##	ax2.errorbar( x+0.01, percent, yerr=lamberr, fmt='o', markersize=1.3)
		ax2.errorbar( x+0.01, eigL, yerr=lamberr, fmt='o', markersize=1.3)
		ax2.set_xlim(-1,nx); ##ax2.set_ylim(0,40)
		ax2.set_title('eigenspectrum scree plot')
	
		#_plot the EOFs scaled by the sqrt of the eigenvalues
		eof0 = dot( eigV[:,0], np.sqrt(eigL[0]) )
		eof1 = dot( eigV[:,1], np.sqrt(eigL[1]) )
		eof2 = dot( eigV[:,2], np.sqrt(eigL[2]) )
		ax3 = plt.subplot2grid( shape, (3,0), colspan=2 )
		ax3.plot( x, eof0, '-k', linewidth=lw[0] )
		ax3.plot( x, eof1, '-r', linewidth=lw[1] )
		ax3.plot( x, eof2, '-b', linewidth=lw[2] )
		ax3.set_title('first three eofs : black=1,red=2,blue=3')
		ax3.set_xlim(0,len(x)-1)

		#_plot first three pcs
		ax4 = plt.subplot2grid( shape, (4,0), colspan=2 )
		pc0 = pca[:,0] / np.sqrt(eigL[0])
		pc1 = pca[:,1] / np.sqrt(eigL[1])
		pc2 = pca[:,2] / np.sqrt(eigL[2])
		ax4.plot( t, pc0, '-k', linewidth=lw[0] )
		ax4.plot( t, pc1, '-r', linewidth=lw[1] )
		ax4.plot( t, pc2, '-b', linewidth=lw[2] )
		ax4.set_title('first three pcs : black=1,red=2,blue=3')
	##	ax4.set_ylim(-4,4)
		ax4.set_xlim(0,len(t)-1)
		
	elif multi:
		pass
	
	for ax in fig.axes: lt.shrink_ticks(ax,size=8)
	
	plt.tight_layout()
	if show: plt.show()
	else: plt.savefig('problem_'+label+'_.png')

def composite(data, idxp, idxn, alphaT=.95, alphaF=.90, **kwargs):
	''' returns boolean array of where null is rejected '''
	import numpy as np
	import matplotlib.pyplot as plt
	import libtools as lt
	import libnva as ln
	from scipy.stats import f
	from scipy.stats import t
	from numpy import nan, tile, sqrt, round, ones, append
	
	Np, Nn = len(idxp), len(idxn)
	ny, nx = 181, 360
	
	#_find differences between means
	diff = data[idxp,:,:].mean(axis=0) - data[idxn,:,:].mean(axis=0)
	sigp = data[idxp].std(axis=0) # (x,y)
	sign = data[idxn].std(axis=0) # (x,y)
	
	'''
	#  3.1:  Identify points with equal variance
	#  1.  90% significance level
	#  2.  H0:  variance is the same
	#      HA:  Variance is different
	%  3.  Statistic:  2-tailed F-test
	%  4.  Critical region:  finv(0.05,Np-1,Nn-1) < F < finv(0.95,Np-1,Nn-1)
	%  5.  Evaluate the statistic (below)
	''' 
	#_The equality of variance evaluation is required to use the appropriate
	# T test calculation.
	#_calc f test crap (f.isf is ONE TAILED and for UPPER prob, opposite
	# of MLAB)
	F = sigp**2 / sign**2 #_variance of positive cases over neg (x,y)
	flw = f.isf(0.95, Np-1, Nn-1) #_get upper level (flt)

	#_find where cases variance is within F limits, and cases where it's diff
	varIdxSame = flw<F
	varIdxDiff = F>=flw

	'''
	%  3.2  Calculate t-test with pooled std. dev.

	%  Set up t-test with equal variance.
	%  1.  95% significance level
	%  2.  H0:  SST during El Nino events is different than La Nina events
	%      HA:  SST is no different during the two events
	%  3.  Statistic:  2-tailed t-test with N1+N2-2 DOFs
	%  4.  Critical region:  tlow < T < thi
	%  5.  Evaluate the statistic ...
	'''
	#_empty array for T test	
	T = tile(nan, (ny,nx))
	
	#_calculate pooled stdev. is std weighted by pop size
	#_all t tests with unequal sample size
	#_t test with same variance [t = (u1-u2)/s12/sqrt(n1^-1+n2^-1)]
	#_in locations with the same variance
	spool = sqrt(((Np-1)*sigp**2 + (Nn-1)*sign**2) / (Np+Nn-2))
##	T[varIdxSame] = diff[varIdxSame] / spool[varIdxSame] / sqrt(1./Np+1./Nn)
					
	#_t test with diff variance [t= (u1-u2)/sqrt(s1^2/n1 + s2^2/n2)]
	# diff var means no pooled std
##	T[varIdxDiff] = diff[varIdxDiff] / \
##				sqrt( sigp[varIdxDiff]**2/Np + sign[varIdxDiff]**2/Nn )
	print Np, Nn
#	T = diff / spool / sqrt(1./Np+1./Nn)
	T = diff / spool / sqrt(1./Nn)
	
	# DOF, unequal var, DOF = (s1^2/n1 +  s2^2/n2) /
	#					(s1^2/n1)^2 / (n1-1) + (s2^2/n2)^2 / (n2-1)
	# DOF, equal var, diff sample = n1+n2-2
	DOF = ( sigp**2/Np + sign**2/Nn )**2 / \
		((sigp**2/Np)**2/(Np-1.) + (sign**2/Nn)**2/(Nn-1.))
	DOF[varIdxSame] = Np + Nn - 2
	DOF = round(DOF)
	
	#_fill in H with 0's where H0 is accepted and 1's where it is rejected
	thresh_hgh = (1-alphaT)/2
	thi = t.isf( thresh_hgh, DOF )

	#_H0 accepted => 0 [SSTs same], H0 rejected => 1 [SSTs diff]
	#_find where T is within significance range, set those values to 0
##	idx = append( (tlw<T),(T<thi) ).reshape(2,ny,nx).all(axis=0) 
	H	= ones((ny, nx)) 
	H[thi < T] = 0

	return H, diff

def read_GDASGRIB(dir_in='/'.join((os.environ['PRODUCTS'], 'gdas')),
	dtgs=None, level=1000, vars=['geopot', 'uwind', 'vwind'], **kwargs):
	'''
	fname	string,		file path
	
	read files gdas1.PGrbF##.YYMMDD.HHz
	'''
	import pygrib as pyg
	from glob import glob
	import re
	from numpy import recarray, ndarray, meshgrid, linspace
	
	re_dtg = re.compile('\.(\d{6})\.(\d{2})z$')

	trans = {
		'geopot'	: 'Geopotential height anomaly',
		'uwind'		: 'U component of wind',
		'vwind'		: 'V component of wind',
	}

	dtype	= [('dtg', 'a10'), ('lat', ndarray), ('lon', ndarray)]
	dtype.extend([(v, ndarray) for v in vars])
	data 	= recarray((0,), dtype=dtype)

	lat = linspace(-90, 90, 181)
	lon	= linspace(-180, 180, 360)
	lon2d, lat2d = meshgrid(lon, lat)

	for fname in glob(dir_in + '/gdas*'):
		tmp = re_dtg.search(fname)
		dtg = '20' + tmp.group(1) + tmp.group(2)
		if dtgs != None and dtg not in dtgs: continue	#_skip if we care
		 
		#_expand recarray
		data.resize(data.size + 1)
		data.dtg[-1] = dtg
		data.lat[-1] = lat2d
		data.lon[-1] = lon2d
		
		#_open file
		print(fname)
		grb = pyg.open(fname)
	
		#_loop over vars
		for var in vars:
			gvar = trans[var]	#_get gdas name
			locations = grb.select(name=gvar)
			
			#_find proper level
			for loc in locations:
				if loc.level == level: data[-1].__setattr__(var, loc.values)
	
		grb.close()
		
	return data

def regression_map(fcst, met, **kwargs):
	''' Dan's regression map example '''
	import numpy as np
	from numpy import dot, matrix
	
	#_get dimension sizes
	nt, ny, nx = met.shape
	
	#_reshape into 2d array
	skt = met.copy().reshape( nt, ny*nx ) 
	
	#_MAP OF COVARIANCE?? ( anom(x)*anom(y)/ N-1 )
	#_means are already removed, get the covariance of the n34 with each column
	# of skt-monthly.  The variance of n34 should be 1, so we don't have to div
	map_reg = dot(fcst, skt) / fcst.var() / (nt-1.)
	
	#_MAP OF CORRELATION COEF ( cov(x,y) / sig(x)sig(y) )
	#_corr coeff requires extra step. divide each col by its std then get cov
##	map_cor = n34 * np.matrix( skt / skt.std(axis=1).reshape(nt,1) ) / (nt-1.)
##	map_cor = (n34 / n34.std()) * (skt / skt.std(axis=1).reshape(nt,1)) / (nt-1.)
##	print n34.shape, skt.shape, map_cor.shape, map_reg.shape
##	print 'Why matrix require division...of... oi... '
	
	map_reg = map_reg.reshape( ny, nx )
##	map_cor = map_cor.reshape( ny, nx )
	
##	return map_reg, map_cor
	return map_reg
	
def plot_reg(regres, H, lat, lon, figure=None, nax=111, label='', **kwargs):
	from libnva import draw_map
	import matplotlib.pyplot as plt
	from matplotlib.cm import RdBu_r
	from numpy import linspace
	
	if figure == None: 
		figure 	= plt.figure()
		ax		= figure.add_subplot(111)
	else:
		ax = figure.add_subplot(nax)
	
	#_draw and plot
	m 	= draw_map(corners=[-90,90,-180,180], laby=[1,0,0,0])
	C5 	= m.contourf( lon, lat, regres, cmap=RdBu_r, latlon=True, extend='both')
##		latlon=True, levels=np.linspace(-.5,.5,11), extend='both')

	ax.set_title('Regression Map:' + label, size='x-small')
	ax.set_xlabel('Basis: CAPOVERDE AOD', size='x-small')
	
##	C2 = m.contour( lon, lat, H, levels=[-.5,.5], latlon=True, 
##		 linestyles=['--','-'], colors=['magenta','magenta'] )
	
	br3 = plt.colorbar(C5, orientation='vertical', shrink=0.8)
	[ t.set_fontsize(6) for t in br3.ax.get_yticklabels() ]
	
def read_forecasts( path_fcst=os.environ['PRODUCTS']+'/icap/', label='2012FL',
	sites=None, **kwargs ):
	'''
	Reads merged forecast files produced by ICAP.merge_data.py
	These are collacted with AERONET sites	
	'''
	from glob import glob
	import libtools as lt
	import libnva   as ln
	import libclass as cl
	import numpy    as np
	from multiprocessing import Queue, Process

	path_fcst = '/'.join((path_fcst, 'seasonal_statistics-'+label, 'points'))
	if sites == None:
		files = glob( path_fcst + '/*-*-'+label+'.nc' ) #_AFTER
	else: 
		files = []
		for site in sites:
			files.extend(glob(path_fcst+'/*-'+site+'-'+label+'.nc'))

	records = False
	print( path_fcst )
	dtype = [       ('values', cl.var),
	                ('model', 'a16'),
	                ('epoch', 'i8'),
	                ('fhr', 'i4'),
	                ('code', 'a10'),
	                ('lat', 'f4'),
	                ('lon', 'f4'),
	                ('variable', 'a20'),
	                ('ensemble', 'i4'),
	                ('dimname', tuple),
	                ('dimsize', tuple) ]
	records = np.recarray( (0,), dtype=dtype )

	groups = lt.setup_groups( files, **kwargs )
	for group in groups:
	        l = len(group)
	        q = [None]*l
	        t = [None]*l
	        for i in xrange(l):
	                file = group[i]
	                q[i] = Queue()
	                args = (file,)
	                kwargs.update({'pipe':q[i]})
	                t[i] = Process( target=ln.read_period_points,
	                        args=args, kwargs=kwargs )
	                t[i].start()

	        for i in xrange(l):
	                tmp, attrv = q[i].get()
	                for rec in tmp:
	                        rec.values = cl.var( rec.values,
	                                attrn=['member'], attrv=attrv )
	                records = lt.merge([records,tmp])

	ln.check_member_order( records )

	return records

def parse_aeronet( file, **kwargs ):
	'''reads netcdf aeronet file for one site'''
	from netCDF4 import Dataset
        import numpy as np
        import libtools as lt
        print(file)

        ncdf    = Dataset( file, mode='r', format='NETCDF3_CLASSIC' )
        nt      = len( ncdf.dimensions['time'] )
        modes   = [str(a) for a in ncdf.modes.split(',')]

        #_construct recarray datatyp
        values = {}
        for mode in modes:
                tmp = ncdf.variables[mode][:]
                tmp = np.ma.masked_where( tmp == -9999., tmp )
                mask = tmp.mask == False
                values.update({ mode : tmp[mask] })

        nt = mask.tolist().count(True)
        dtype = [('epoch','i8'),('code', 'a6'),('lat','f4'),('lon','f4')]
        [ dtype.append((mode,'f4')) for mode in modes ]
	records         = np.recarray((nt,), dtype=dtype)
        records.dtg_start = ncdf.dtg_start
        records.dtg_end = ncdf.dtg_end
        records.code    = [ncdf.code] * nt
        records.modes   = modes
        records.lat     = [ncdf.lat] * nt
        records.lon     = [ncdf.lon] * nt

        #_some of the files didn't save the time correctly
        # this corrects that.
        epochs = ncdf.variables['time'][:].copy()
        if (epochs[1:] - epochs[:-1]).max() == 1:
                epoch_start = np.min( epochs )
                epochs -= epoch_start
                epochs *= 21600
                epochs += epoch_start
        epochs = epochs[mask]
        records.epoch = epochs[:]

        records.total   = values['total']
        records.coarse  = values['coarse']
        records.fine    = values['fine']
        return records

def read_aeronet( path_aero=os.environ['PRODUCTS'] + '/aeronet', **kwargs):
        ''' read in all aeronet sites in path_aero '''
        from glob import glob
        import libtools as lt
        files = glob( path_aero + '/*-*-*nc' )

        records = False
        for file in files:
                if records != False:
                        tmp = parse_aeronet( file )
                        records = lt.merge([records, tmp])
                else:
                        records = parse_aeronet( file )

        return records

#############################################################################_80
#############################################################################_80
#############################################################################_80

def dbg( msg ):
        ''' 
        if global debug is set to true, be more verbose 
        msg     : str, Message to be printed
        l       : int, Debug level of message.  Set higher for lower level 
                        messages.  As debug increases, noisiness should also.
        '''
        import libtools as lt
        msg = lt.to_string( msg )
        if hasattr( msg, '__iter__'): msg = ' '.join( msg )

        import inspect
        curf = inspect.currentframe()
        calf = inspect.getouterframes( curf, 2 )
        file, line, method = calf[1][1:4]
        file = file.split('/')[-1]
        print '[%s.%s.%i] %s' % ( file, method, line, msg )
