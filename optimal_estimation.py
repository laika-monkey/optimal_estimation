#!/usr/bin/env python
#############################################################################_80
# file:		optimal_estimation.py	(sensor agnostic)
# author:	Walter R. Sessions, 2014.10.04	(overhaul 2016.04.20.blazeit
# purpose:	Attempts to use optimal estimation techniques taken from Clive
#			Rodgers's Inverse Methods for Atmospheric Sounding, Vol II to
#			retrieve profile information from hyperspectral infrared imagers
#			such as S-HIS and AIRS. Code also includes cases to run tests on
#			pseudo-obs.
#
#			Much of this code is intended to be split out into separate 
#			files, such as the namelist and hs3_2013 dictionaries. All
#			in here just for testing.
################################################################################
################################################################################


from subprocess import call
import numpy as np
from numpy import matrix
from numpy.linalg import inv, svd
from scipy.linalg import sqrtm
from scipy.io import netcdf
import os
import math
import time
import sys


'''
To turn off uncertainties
	'model'			set 'uncertainty' to 0
	'instrument'	set 'uncertainty' to 0
'''

# DIR_LOG  = os.path.expanduser('~/qsub_logs/')
DIR_LOG  = os.environ['LOG']
DIR_PROD = os.environ['PRODUCTS']
DIR_LBL  = 'LBL-RTM' #_generally want to use just LBL-RTM
DIR_TMP  = '/data/wsessions/TMP'
DEBUG    = 1


################################################################################
#_retrieval_####################################################################
################################################################################


def optimal_estimation(obs_rads, obs_wave, Sr,
	#_Which sensor are we using for radiances?

	surf_temp=-1,
	dir_lblrtm='.', 
	ref_wn=900, 
	L_curve_index=4, 
	L_curve_size=151, 
	form='rodgers',
	max_iter=25,
	out_type=dict,
	pipe=None,
	apriori={},
	posterior_dump=None,
	bias_correction=False,
	**kwargs):
	'''
	phys_info
	radiances		Array containing observations to be matched
	wavenumber		Array containing associated wavenumber dimension

	form    str,    carissimo or rodgers, defines the form of OE
                    as following either:

	                CARISSIMO, 2005
	                The physical retrieval methodology for IASI: the d-IASI code

	                RODGERS, 2000
	                Inverse Methods for Atmospheric Sounding
		
					TURNER, 2014
					Information Content and Uncertainties in Thermodynamic
					Profiles and Liquid Cloud Properties Retrieved from the
					Ground-Based Atmospheric Emitted Radiance Interferometer
					(AERI)
	  
                    Both including the usage of Gamma found in Turner, 2014
	L_curve_index	int,	range over which to do L-curve calculations

	Sr		matrix,	instrument error covariance matrix

	ret_vars	list,	contains strings for var to actually do a retrieval upon

	out_type	type, dict will output a dictionary by state var,
						defaults to matrix
	surf_temp	float,	Surface temperature value to use. Set to -1 to 
						use LBL-ATM surface temperature, set to -2 to use
						GDAS (not implemented yet), or >0 to be arbitrary
	dynamic_CTP		bool||list,	if false, use static apriori cloud top pressure
								otherwise use what's available based on
								priority defined by list. If nothing available,
								will default to static.
	eosterior_dump	str,		for qsub write. dumps the dictionary with 
								posterior uncertainties to a pickle. Arg is
								name of file	
	apriori		list(dicts),	definition of initial layers

	bias_correction	dict,		apply a spectrally dependent bias correction
								to the measurement state matrix, or pass as False
	persistence		dict,		If not false, use previous fov values as a priori
								as long as they are within a certain percent of
								value/time
	'''
	from lblrtm_utils import microwindow_average, get_surf_temp
	from lblrtm_utils import check_reff_v_ssp as check_ssp
	from numpy import append, array, arange, c_, diag, eye, matrix, power, sqrt
	from numpy import zeros, diff, tile, trace, linspace, concatenate, ndarray
	from scipy.linalg import sqrtm
	from libtools import strictly_increasing
	from copy import deepcopy
	from hs3_utils import get_naaps_ctp
	from pickle import dump
	from libgeo import planck_inv

	#_for testing
	import matplotlib.pyplot as plt
	from libgeo import planck, p2z
	from tape7 import tape7

	############################################################################
	#_LINEAR_LEAST_SQUARES_#####################################################
	############################################################################

	#_pull out sensitivity dictionary
	sensitivity = kwargs.get('sensitivity')
	uncertainty = kwargs.get('uncertainty') 
	state_vars = kwargs.get('state_vars', ['tau'])
	uncert_vars = kwargs.get('uncert_vars', ['tau'])
	out_label = kwargs.get('out_label', '')

	nstate = len(state_vars)
	nlayer = len(apriori)

	#_get height of observer from LBL-RTM to limit convergence
	t7 = tape7(os.path.join(dir_lblrtm, 'TAPE7'))
	max_z = t7.H1 - sensitivity['z']*2

	#_define an a priori cloud layer dictionary. 
	cld = deepcopy(apriori) #_first guess

	'''

	NOTHING ABOUT LOCATION IS REQUIRED FOR OE, MOVE THIS INTO WRAPPER

	#_use static cloud top height or use DR SHIS guess
	#_convert pressure to AGL... sorta...
	for layer in cld:
		# MAKE THIS A METHOD ELSEWHERE 
		if dynamic_ATP:
			for ctp_source in dynamic_ATP:
				if ctp_source == 'NAAPS':
					arg = (epoch, latitude, longitude) 
					CTP, CBP = get_naaps_ctp(*arg, **kwargs)
					CTZ = p2z(CTP) / 1e3 
					CBZ = p2z(CBP) / 1e3 

					#_if apprioriate layer and value NAAPS, use
					if layer['type'] == 'aerosol' and	((type(CTZ) == list or \
														type(CTZ)==ndarray) and\
															max(CTZ) > 0) or   \
														(type(CTZ) == float and\
															CTZ > 0): 
				#	if layer['type'] == 'aerosol' and CTZ > 0:
						layer['z_top'] = max(CTZ) #CTZ
						layer['z'] = max(CBZ) #CTZ-0.1
						dbg(('using just the top layer', max(CTZ), max(CBZ)))
					else:
						continue

				elif ctp_source == 'SHIS':
					CTZ = p2z(fov.SHIS_CTP) / 1e3

					#_if apprioriate layer and value present, use SHIS DR 
					if layer['type'] == 'aerosol' and CTZ > 0:
						layer['z_top'] = CTZ
						layer['z'] = CTZ-0.1
					else:
						continue
					
				#_set, move on (break condition only met when not NAAPS/SHIS	
				break			

	for layer in cld:
		if dynamic_CTP:
			for ctp_source in dynamic_CTP:
				if ctp_source == 'SHIS':
					CTZ = p2z(fov.SHIS_CTP) / 1e3

					#_if apprioriate layer and value present, use SHIS DR 
					if layer['type'] == 'cloud' and CTZ > 0:
						layer['z_top'] = CTZ
						layer['z'] = CTZ-0.1

						#_check if ice layer
						obs_rads = fov.SHIS_radiances
			
						#_get brightness temperatures for a few windows	
						r, w = microwindow_average(obs_rads, obs_wave, -26)
						bt = planck_inv(r,w*100,domain='wavenumber')[-3:].mean()

						if bt < 273.15: #_maybe lower...
							#_make this not dumb later
							ice = '/data/wsessions/lbldis_inputs/' \
								+ 'ssp_db.mie_ice.gamma_sigma_0p100'
							kwargs['ssp_db_files'].append(ice)
							layer['dbnum'] = len(kwargs['ssp_db_files'])
							
					else:
						continue
					
				#_set, move on	
				break			
	'''

	#_If using full obs spectrum, convolve lbldis to retrieval
	# If using microwindows, average retrieval to channels
	dv       = kwargs.get('dv')

	#_DV less than zero tell LBL-DIS to use select channels instead of full spec
	if dv < 0:
		#_when using microwindows, average observation data
		y, obs_wave = microwindow_average(obs_rads, obs_wave, dv)

		#_if the window chanels are averaging a shift of more than 5 K from
		# surface temperature, adjust cloud apriori
		Bt = planck_inv(y, obs_wave*100, domain='wavenumber')
		Bt = Bt.mean()

		for layer in cld:
		#	print "TYPE", layer['type']
		#	print "SURF", surf_temp
		#	print "BTTT", Bt
			if layer['type'] == 'cloud' and surf_temp-Bt > 4.5:
				layer['tau'] = 4.0
	else:
		#_use full observation spectrum (NOT RECOMMENDED FOR OE)
		y = obs_rads[:]

	#_update cloud into kwargs and add analysis sfc temperature
	kwargs.update({'clddef' : cld, 'surf_temp' : surf_temp})

	#_MAKE THE NEXT THREE OPERATIONS NOT NECESSARY
	#_and convert to matrix
	y = matrix(y)

	#_make sure matrix is [m,1]
	if y.shape[1] != 1:
		y = y.T

	#_maybe I don't know
	y *= 1e3

	#_generate a clear sky run based on surf temp guess (for conv plots)
	kwargs['r_clr'], kwargs['n_clr'] = clear_sky_radiance(dir_lblrtm, **kwargs)

	#_get values for uncertainty jacobians
	dbg('launching uncertainty calculations')
	Ki, Fx, wave_lbl = jacobian_forker(dir_lblrtm,
		jacobian_vars=uncert_vars, **kwargs)
	dbg('finishing uncertainty calculations')

	#_if bias correction being applied, put in correct order and type
	if bias_correction:
		bias_correction = _sort_bias(wave_lbl, bias_correction, dv)
	
	#_a priori and definition of the vector of the first guess
	# v and vg in Carismo
	xi = matrix([[a[v] for v in state_vars] for a in apriori]).flatten().T
	xa = matrix([[a[v] for v in state_vars] for a in apriori]).flatten().T

	#_get sizes of measurement and state space
	m = len(y)
	n = len(xi)

	'''
	#_NOT GAMMA DIST, but a regularization factor. (see pveg about range) 
	# See Tikhonov-Miller / Phillips-Twomey / Constrained linear invers. method
	#_CARISSIMO, 2.2.3
	gamma = 10.**linspace(L_curve_index, L_curve_index, L_curve_size) 
	gamma = tile(eye(n), (L_curve_size, 1, 1)) * gamma[:,None,None]

	#_See Carissimo 2005, 2.2 Inversion Methodology
	N			= len(y)			#_size of DOF
	chi2alpha	= N + 2*sqrt(2.*N)	# condition to be <= in eq(3) 
									#_assumes 2sig tolerance interval
	chi22n		= array([1e4, 1e3])	# used in conditional escape below 
	chi2steps	= 1000. 
	'''

	#_(v-vg).T*L*(v-vg), L = B.T = apriori covar
	#_(xi-xa).T*L*(xi-xa)
	############################################################################
	#_optimal_estimation_#######################################################
	############################################################################

	#_initialize collection of all iterations
	all_xi	= xi.copy()		#_keep a running list of state iterations
	all_J	= []			#_list of cost function values
	all_di2	= []			#_DELETE WRS all convergence criteria

	#_arbitrary gamma values for testing
	gamma_arb = [1000, 300, 100, 30, 10, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1]
	
	#_i => OD, Reff (2xMeas)
	#_define covariance matrices
	Sm = matrix(diag(Fx.A.squeeze() * uncertainty['model']))
	Sr = Sr * uncertainty['instrument']
	Sa = array([[uncertainty[unc]	for unc in state_vars] \
									for i in range(nlayer)])
	Sa = matrix(diag(Sa.flatten()))
##	Sa = matrix(diag([[uncertainty[unc]	\
##		for unc in state_vars]			\
##		for i in range(nlayer)]))
##	Sa = matrix(diag(array([uncertainty[unc] for unc in state_vars])))
									#_uncertainty covar of apriori
##	St = sensitivity['surf_temp']	#_include these even though not used?
##	St = uncertainty['surf_temp']	#_include these even though not used?
##	Sh = sensitivity or sensitivity ['z']
	#_s-1 precision/determination of information in prior
	#_total covariance  (dF/dstate * d_unc_state => variance in F by state)
	# uncertainty of SHIS comes from sensor
	# uncertainty of model appears... arbitrary? Do I just trust docs?
	#_assumes uncorrelated error.  Used to normalize ellipse?
	#_Se should include all error that would prevent the Forward model
	# from preproducing y from x.
	Se = Sr + Sm

	#_hrm..... noooooo
	for lidx in range(nlayer):
 	  for jv in uncert_vars:

														#_skip vars w/ uncert in Sa
		if	jv in state_vars or uncertainty[jv] == 0 or	\
			(lidx != 0 and (jv[:4] == 'surf' or jv == 'rel_hum')):	#_only do surface vars once
			continue

		#_initialize jacobian
		Ktmp = matrix(Ki[lidx][jv]) * uncertainty[jv]

		#_covar for
		Stmp = Ktmp * sensitivity[jv] * Ktmp.T
		Stmp = matrix(diag(diag(Stmp)))

		Se = Se + Stmp

	''' look at eigen decomposition of covariance matrices '''
	'''
	LOOK BELOW FOR OTHER QUOTE and WHY IS SENSITIVITY TO TAU/REF NOT HERE

	CHECK DIM V2 TO SEE IF HE'S INCLUDES THE DESIRED ONES IN Sm and THE NON
	IN Se
	'''
	for oe in xrange(max_iter):
		#_update cloud for jacobians to current guess (x_n), updates at end 
		[[layer.update({sv : xi[i+j*nstate, 0]})	\
			for i, sv in enumerate(state_vars)]	\
			for j, layer in enumerate(kwargs['clddef'])]

		#_Fx == Forward Model Obs at Current Guess
		#_calculate jacobian around first guess, get K(jac), F(x), dim of F(x)  
		K, Fx, wave_lbl = jacobian_forker(dir_lblrtm,
			jacobian_vars=state_vars, **kwargs)

		#_plot current iteration DELETE::WRS
		plot(y.copy(), obs_wave, Fx.copy(), wave_lbl, label='OE', **kwargs)

		#_jacobians (build matrix from jacobian dictionary)
		for lidx in range(nlayer):
			
			#_current state of K is K[<layer_number>]{<state_var>}
			# desired state of K is K[n_state, n_measurement]
			for i, sv in enumerate(state_vars):
				#_pull out this layers jacobian for this state_var
				tmp = K[lidx][sv] if K[lidx][sv].ndim == 2 \
								else K[lidx][sv][:,None]

				#_if first state variable and first layer, initialize
				if not lidx and not i:
					Ki = tmp
				else:	
					Ki = append(Ki, tmp, 1)

		Ki = matrix(Ki)
##		Ki = array([[K[lidx][v].flatten() \
##			for v in state_vars] \
##			for lidx in range(nlayer)])
##		Ki = matrix([[K[lidx][v] \
##			for v in state_vars] \
##			for lidx in range(nlayer)]).T
#_#		#_s-1 precision/determination of information in prior
#_#		#_total covariance  (dF/dstate * d_unc_state => variance in F by state)
#_#		# uncertainty of SHIS comes from sensor
#_#		# uncertainty of model appears... arbitrary? Do I just trust docs?
#_#		#_assumes uncorrelated error.  Used to normalize ellipse?
#_#		#_Se should include all error that would prevent the Forward model
#_#		# from preproducing y from x.
#_#		Se = Sr
#_#		for jv in jacobian_vars:
#_#			#_uncertainty of state vars in Sa
#_#			if jv in state_vars or uncertainty[jv] == 0:
#_#				continue
#_#
#_#			#_covar for
#_#			Stmp = K[jv] * sensitivity[jv] * K[jv].T
#_#			Se = Se + Stmp 
##		Se = Sr + Sm # + Kt.T*St*Kt	# + Kh*Sh*Kh.T
##		Se = Sr + Sm # + Kt*St*Kt.T	# + Kh*Sh*Kh.T	#_orig?
##		sum([K[sv] * uncertainty[sv] * K[sv].T for sv in state_vars]) 

		#_choose between formulation of Carissimo (2005) or Rodgers (2000)
		if form == 'carissimo':
			pass
			'''
			#_change of variables to make the problem dimless, easier to
			# manipulate and more efficient.
			#_CARISSIMO ET AL. (2005) - The physical retrieval methodology for
			#							IASI: the delta-IASI code
			Yn	= (y - Fx) + Ki*(xi - xa)#_current x minus first guess minus xa?
	
			#_the inverse of the variance of all errors * dF/dstate
			Ji	= matrix(inv(sqrtm(Se)) * Ki)	#_CARISSIMO eq.18~
			G	= Ji*sqrtm(Sa)			#_Gn = Jn*sqrt(B), B=a priori covar mtrx
			U,D2,V = svd(G.T*G)
			V	= U						#_matrix being factorized is symmetric
			z	= inv(sqrtm(Se))*Yn		#_same reference

			#_new radiance values
			# (SHIS minus LBL-DIS rads + dF/dx(first guess minux a priori)
	
			#_HANSEN (1992) - 	Analysis of discrete ill-posed problems by means
			#					of the L-curve.
			curv = []
			for g in range(0): #range(L_curve_size):	#_can this be indexed
	
				#_CARISSIMO eq.20: Solve for u == inv(Sa)*x_n+1
				# u_n+1 = V inv(gI+D2)V.T G.T z
				u = V * inv(gamma[g] + D2) * V.T * G.T * z
	
				#_first and second derivative of u
				u1 = -V * inv(gamma[g]*eye(n) + D2) * V.T*u
				u2 = -2. * inv(gamma[g]*eye(n) + D2) * V.T*u1

##				Should all be dimensionless
##				a(gamma)= u.Tu
##				a(g)'	= u.T'u + u.Tu'
##				a(g)"	= u.T"u + 2u.T'u' + u.Tu" 	
##	
##				b(gamma)= (Gu-z).T(Gu-z)
##				b(g)'	= (Gu').T(Gu-z) + (Gu-z).T(Gu')
##				b(g)"	= (Gu").T(Gu-z) + 2(Gu').T(Gu') + (Gu-z).T(Gu")

				a1 = u1.T*u + u.T*u1				#_first der of param a
				a2 = u2.T*u + 2.*u1.T*u1 + u.T*u2	#_first der of param a
				b1 = (G*u1).T*(G*u-z) + (G*u-z).T*(G*u1)
				b2 = (G*u2).T*(G*u-z) + 2.*(G*u1).T*(G*u1) + (G*u-z).T*(G*u2)
	
				#_Curvature of the L-curve may be then expressed in terms of gamma
				# CARISSIMO eq.23
				cgamma	= abs(a1*b2 - a2*b1) / power((a1**2 + b1**2), 1.5)
				curv	= append(curv, cgamma)
##				this is currently a mess.  dimensions wrong? 
	
			#_search for the best gamar regularization param for the current
			# iter  Gamma is a ratio of prior information to obs, so g > 1
			#
			# should an array of gamma 1000, 300, 100, 30... show up, it's from 
			# TURNER 2014
			gammamax = 1 # curv.max()

			#_build giant matrices, calculate all at once.
			#_optimal estimation rearranged according to Carissimo 2005
			# eq.21
			u = V * inv(gammamax*eye(n) + D2) * V.T * G.T * z
	
			#_??? state vector for current iteration, w/ dims
			# Sa*u == x_n+1  (is this square root correct and is it just a step)
			#_carissimo above eq.20, u_n+1 definition
			xoe = sqrtm(Sa) * u + xa 

			#_approximate the current state vector
			xoe[0] = round(xoe[0].real, 3)
			xoe[1] = round(xoe[1].real, 2)

			#_chisquared test stuff
			chi2current	= (y-Fx).T * inv(Se) * (y-Fx)
			chi2n		= append(chi2n, chi2current)
			chi2steps	= (chi2n[-2] - chi2n[-1]) / chi2n[-2]
	
			print chi2current.shape, chi2n[-1], N, n, gam
			#_set exit conditions
			#_CHARISSIMO (3), 2.2.1 Iterization
			# (R - Fv).T*S^-1*(R-Fv) <= chi2alpha
			if chi2n[-1] < chi2alpha and gam == 1:		#_and di2 < 0.01:
	##		if chi2current < chi2alpha and gam == 1:	#_and di2 < 0.01:
				out_status = 1	#_good
				print 'GAMMA STATIC'
				break
			elif chi2steps <= 0.1:
				out_status = -1	#_bad, oscillation or slow convergeance
				break
			elif chi2n[-1] >= chi2n[-2]:
				out_status = -3	#_bad, residuals increasing
				break
			else:
				#_still converging... this is superfluous.
				xi = xoe.copy()
				continue
			'''

		#_RODGER_###############################################################
		elif 'rodgers':
			try:
				gam = gamma_arb[oe]	# curv.max()
			except IndexError:
				gam = 1.

			#_RODGERS eq.5.9
			#_n-form
			xoe = xa + inv(gam*inv(Sa) + Ki.T*inv(Se)*Ki) * Ki.T \
					* inv(Se) * (y-Fx-bias_correction + Ki*(xi-xa)) 
		##	xoe = xa + inv(gam*inv(Sa) + Ki.T*inv(Se)*Ki) * Ki.T \
		##			* inv(Se) * (y-Fx + Ki*(xi-xa)) 

			#_m-form
			#	= xa + Sa * Ki.T * inv(Ki*Sa*Ki.T + Se) * (y-Fx + Ki*(xi-xa))

			#_ensure output is within bounds
			xoe = check_state_boundaries(xoe, **kwargs)

			#_add to previously computed retrievals
			all_xi = concatenate((all_xi, xoe), axis=1)

			#_TURNER 2014 eq.3-5, calculating the posterior error covar
			B	= gam*inv(Sa) + Ki.T*inv(Se)*Ki
			A	= inv(B) * Ki.T*inv(Se)*Ki		#_averaging kernel	
			Sxi = inv(B) * (gam**2 * inv(Sa) + Ki.T*inv(Se)*Ki) * inv(B)
			degf= trace(A)

			#_RODGERS eq.5.29
			# di2 = (x_i - x_i+1)' * S^-1 * (x_i - x_i+1) << n
			di2 = (all_xi[:,oe]-all_xi[:,oe+1]).T * inv(Sxi) \
				* (all_xi[:,oe]-all_xi[:,oe+1])
			all_di2.append(di2)	#_DELETE WRS for testing

			#_update initial point
			xi = xoe.copy()
			dbg(('CURRENT STATE','\n',xi))
	
			#_MAP / Minimum cost function.  Find CF
			J = (xi-xa).T*inv(Sa)*(xi-xa) + (Fx-y).T*inv(Se)*(Fx-y)
			all_J.append(J)

##			crit_val = (xi - xoe).T*inv(Sa)*(xi - xoe) 
##			if crit_val < n/100:
			dbg(('TESTING_CONV', oe, di2, J, di2 < 0.01, J.A[0] < 0.01))
			dbg(('TESTING_TYPE', type(J.A[0]), J.A[0].shape))
			dbg(('TESTING_____', bias_correction, Fx.shape))
			try:
				dbg(bias_correction.shape)
			except:
				pass

			if di2 < 0.001 and gam == 1:
##			if J.A[0] < 0.01 and gam == 1:	
				#_convergeance
				break

		elif 'jacobs':
			#_SEE SECTION 7 FOR JACOBS TREATMENT USING ADJOINT
			pass

	#_update final iteration and rerun with no cloud and converged cloud
	[[layer.update({sv : xi[i+j*nstate,0]})			\
		for i, sv in enumerate(state_vars)]	\
		for j, layer in enumerate(kwargs['clddef'])]
	lbldis_final(dir_lblrtm, **kwargs)

	#_make sure plot is complete DELETE WRS
	plot(y.copy(), obs_wave, Fx.copy(), wave_lbl,
		label='OE_LAST', last=True, **kwargs)
	dbg(os.getcwd())
	try:
		plot_close(**kwargs)
	except IOError:
		dbg('closing image failed')

	#_prepare output type	
	if out_type == dict:

		#_initialize output dictionary
		out = [{}] * nlayer		#_not really used for anything right now...
		xi = array(xi)
	
		#_create dictionary for retrieval uncertainty
		posterior = {} # [] * nlayer
		posterior_std = sqrt(diag(Sxi))
	
		#_build output jacobian dictionary keyed with state_var order
		fname_unc = 'uncertainties_LOOP-FOR-LAYERS_{0}.txt'.format(out_label)
		with open(os.path.join(dir_lblrtm, fname_unc), 'w') as f:
			for l, cld_opt in enumerate(cld): #range(nlayer):
				key = cld_opt['ssp_db']
				posterior.update({key : {}})
 
				for i, v in enumerate(state_vars):
					idx = i + nstate*l
					out[i].update({v : xi[idx, 0]})
												#_calling function will
												# have access to apriori
												# to know the ssp order
					posterior[key].update({v : posterior_std[idx]})
	
					f.write('{0:10s} {1:10.5f}\n'.format(v, posterior_std[idx]))

	#	dbg(('WALTER_11', posterior))
	#	dbg(('WALTER_12', posterior_std))
	#	dbg(('WALTER_13', xi.shape))
	#	dbg(('WALTER_14', xi))
	#	dbg(('WALTER_15', state_vars))
	#	dbg(('WALTER_16', out))

		#_dump posterior
		if posterior_dump is not None:
			pkname = os.path.join(dir_lblrtm, posterior_dump)
			dump(posterior, open(pkname, 'wb'))
			
		if pipe is None:
			return out, all_xi, posterior

		else:
			pipe.send((out, all_xi, posterior))
			pipe.close()

	#_don't use this... it's functionally broken
	else:
		raise RuntimeError, 'don\'t use this'
		return xi	


def _sort_bias(wave_lbl, bias_correction, dv):
	''' make sure bias correction is in right order '''
	from numpy import argsort, array, matrix
	from libtools import strictly_decreasing

	#_pull out values to work with
	wave_bias = bias_correction[dv].keys()
	bias = array(bias_correction[dv].values())

	#_put in ascending order
	bias = bias[argsort(wave_bias)]
		
	#_if order disagrees, reverse it
	if strictly_decreasing(wave_lbl):
		bias = bias[::-1]

	#_put in used variable
	bias_correction = matrix(bias)

	dbg('BIAS_CORRECTION:')
	dbg(wave_lbl)
	dbg(bias_correction)
	return bias_correction.T


def jacobians(dir_lblrtm, sensitivity={}, jacobian_vars=['tau', 'ref'],
	clddef=[], output_array=True, pipe=None, apriori_rads=None, nsend=1,
	clean=False, lidx=0, sensor='AIRS', **kwargs):
	'''
	Builds an MxN matrix for the jacobian describing the relationship 
	of the forward model and the state variables linearized around
	certain values.  

	sensitivity		dict,	dictionary of the sensitivities of various
							state variables
	jacobian_vars	list,	strings of keys to sensitivity.dict above
							for the variables that will make up the jacobian

	clddef			list,	list for cloud layers, see write_lbldis_params
							gets updated here to that method.
	output_array	bool,	loosey goosey way of doing this.  When true,
							output Jacobian is in a single array with the
							column order following jacobian_vars.  Otherwise,
							output as 1D arrays of dictionary labeled by 
							jacobian_vars
	lidx			int,	layer index within clddef to perturb

	Each run is performed depending on the presence of the label in the
	jacobian_vars list.  Do not try to loop over that and be slick about it.
	There are too many variable specific things that will come into play
	down the road to do that.

	The Jacobian will be linearized around the cloud values passed in
	clddef, then shifted according to the corresponding change values
	found in sensitivity
	'''
	from lblrtm_utils import write_lbldis_parameters as lbldis_write
	from numpy import append, matrix, zeros
	from copy import deepcopy
	import time
	from tape7 import tape7
	from convolution import sensors

	#_load appropriate convolution
	convolve = sensors[sensor]

	dbg(dir_lblrtm)
	dbg(jacobian_vars)
	
	#_get processing ID
	pid = os.getpid()

	#_create a shortcut
	path = os.path.join

	#_move to working directory
	dir_pwd = os.getcwd()
	os.chdir(dir_lblrtm)

	#_check to make sure tau values are in list
	for cloud in clddef:
		if type(cloud['tau']) != list:
			cloud['tau'] = [cloud['tau']]

	#_initialize dictionary
	jacobian = {}

	#############################
	#_optical_depth_sensitivity_#
	#############################
	#_always do optical depth, but keep tabbed for symmetry
	if 'tau' in jacobian_vars:	

		#_make a copy to futz with here
		cloud = deepcopy(clddef) 

		#_add second instance of lbldis run to tau field in clddef
	##  2015.04.21 WRS removed. Let write_lbldis split this?
	##	also with the inclusion of lidx, I don't want every tau arb inc'd 
	##	dOD = sensitivity['tau'] / len(cloud)	#_split inc among layers(1)
	##	for cld in cloud:
	##		cld['tau'].append(cld['tau'][0] + dOD) #_make sure this ok 
		dOD = sensitivity['tau']
		for l, cld in enumerate(cloud):
			cld['tau'].append(cld['tau'][0] + dOD * (lidx == l))

		#_create a name for this run based upon state variables tested
		label = 'tau_{0:5.3f}_{1:07d}'.format(cloud[lidx]['tau'][-1], pid)
		fname = path(dir_lblrtm, 'lbldis.input.' + label)
		oname = 'lbldis.out.{0}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kwargs)

		if call(['lbldis', fname, '0', oname]):
			dbg('lbldis failed to run {0}'.format(os.getcwd))
			raise RuntimeError, 'lbldis failed to run\n{0}\n{1}'.format(fname,
																		oname)

		ncdf = netcdf.netcdf_file(oname+'.cdf','r')
		rads = ncdf.variables['radiance'][:] #_[nchan, ninstances]
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave_tmp = wave.copy()
			wave, tmp0 = convolve(wave_tmp, rads[:,0])
			wave, tmp1 = convolve(wave_tmp, rads[:,1])
		##	wave, tmp0 = convolve2shis(wave_tmp, rads[:,0])
		##	wave, tmp1 = convolve2shis(wave_tmp, rads[:,1])
	
			n_rads = zeros((wave.size, 2))
			n_rads[:,0] = tmp0.copy().squeeze()
			n_rads[:,1] = tmp1.copy().squeeze()
			rads = n_rads.copy()

		#_find first column of Jacobian (dMeasurement / dState)
		try:
			jacobian['tau'] = ((rads[:,1] - rads[:,0]) / dOD)
		except IndexError:
			jacobian['tau'] = ((apriori_rads - rads[:,0]) / dOD)

		#_keep a copy of the apriori rads
		apriori_rads = rads[:,0].copy()

		if clean:
			[os.unlink(ff) for ff in [fname, '{0}.cdf'.format(oname)]]

		if pipe is not None:
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			sys.exit(0)	#_don't exit? Need apriori below?
	
	################################
	#_effective_radius_sensitivity_#
	################################
	elif 'ref' in jacobian_vars: 
		
		#_make a copy to futz with here
		cloud = deepcopy(clddef)
		cloud[lidx]['ref'] = cloud[lidx]['ref'] + sensitivity['ref']

		label = 'ref_{0:5.3f}_{1:07d}'.format(cloud[lidx]['ref'], pid)
		fname = path(dir_lblrtm, 'lbldis.input.{0}'.format(label))
		oname = 'lbldis.out.{0}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kwargs)
 
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run', fname
		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve(wave, rads)
		##	wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_this version assumes that rads is from one instance
		tmp = (rads - apriori_rads) / sensitivity['ref']
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T
		jacobian['ref'] = tmp 

		if clean:
			[os.unlink(ff) for ff in [fname, '{0}.cdf'.format(oname)]]

		if pipe is not None:
			t= [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)

	##############################
	#_cloud_position_sensitivity_#
	##############################
	elif 'z' in jacobian_vars: 

		#_make a copy to futz with here
		cloud = deepcopy(clddef)
		cloud[lidx]['z'] = cloud[lidx]['z'] + sensitivity['z']
		cloud[lidx]['z_top'] = cloud[lidx]['z_top'] + sensitivity['z']

		label = 'z_{0:5.3f}_{1:07d}'.format(cloud[lidx]['z'], pid)
		fname = path(dir_lblrtm, 'lbldis.input.{0}'.format(label))
		oname = 'lbldis.out.{0:s}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kwargs)
 
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run'
		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve(wave, rads)
		##	wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_find third column of jacobians
		#_this version assumes that rads is from on instance
		tmp = (rads - apriori_rads) / sensitivity['z']
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T

		jacobian['z'] = tmp 

		if clean:
			[os.unlink(ff) for ff in [fname, '{0}.cdf'.format(oname)]]

		if pipe is not None:
			t = [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)
	
	###########################
	#_cloud_depth_sensitivity_#
	###########################
	elif 'thickness' in jacobian_vars: 

		#_make a copy to futz with here
		cloud = deepcopy(clddef)
		thickness = cloud[lidx]['thickness'] + sensitivity['thickness']
		cloud[lidx]['z_top'] = cloud[lidx]['z'] + thickness 

		label = 'thickness_{0:5.3f}_{1:07d}'.format(thickness, pid)
		fname = path(dir_lblrtm, 'lbldis.input.{0}'.format(label))
		oname = 'lbldis.out.{0:s}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kwargs)
 
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run'
		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve(wave, rads)
		##	wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_this version assumes that rads is from on instance
		tmp = (rads - apriori_rads) / sensitivity['thickness']
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T

		jacobian['thickness'] = tmp 

		if clean:
			[os.unlink(ff) for ff in [fname, '{0}.cdf'.format(oname)]]

		if pipe is not None:
			t = [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)

	###################################
	#_surface_temperature_sensitivity_#
	###################################
	elif 'surf_temp'in jacobian_vars:
		if lidx:
			sys.exit(0)

		#_make a copy to futz with here
		kw = kwargs.copy()
		kw['surf_temp'] = kwargs['surf_temp'] + sensitivity['surf_temp']
		cloud = deepcopy(clddef)

		label = 'surf_temp_{0:5.3f}_{1:07d}'.format(kw['surf_temp'], pid)
		fname = path(dir_lblrtm, 'lbldis.input.' + label)
		oname = 'lbldis.out.{0:s}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kw)
 
		#call(['lbldis', fname, '0', oname])
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run'
		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve(wave, rads)
		##	wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_find third column of jacobians
		#_this version assumes that rads is from on instance
		tmp = (rads - apriori_rads) / sensitivity['surf_temp']
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T

		jacobian['surf_temp'] = tmp 

		if clean:
			[os.unlink(ff) for ff in [fname, '{0:s}.cdf'.format(oname)]]

		if pipe is not None:
			t = [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)
		pass

	##################################
	#_surface_emissivity_sensitivity_# 
	##################################
	elif 'surf_emis' in jacobian_vars:
		if lidx:
			sys.exit(0)
		#_make a copy to futz with here
		kw = kwargs.copy()
		kw['surf_emis'] = clddef[0]['surf_emis'][:]
		for k, e in enumerate(kw['surf_emis'][1]):
			kw['surf_emis'][1][k] = min((e + sensitivity['surf_emis'], 1))
		cloud = deepcopy(clddef)

		label = 'surf_emis_{0:5.3f}_{1:07d}'.format(kw['surf_emis'][1][0], pid)
		fname = path(dir_lblrtm, 'lbldis.input.' + label)
		oname = 'lbldis.out.{0:s}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kw)
 
		#call(['lbldis', fname, '0', oname])
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run'
		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve(wave, rads)
		##	wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_find third column of jacobians
		#_this version assumes that rads is from on instance
		tmp = (rads - apriori_rads) / sensitivity['surf_emis']
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T

		jacobian['surf_emis'] = tmp 

		if clean:
			[os.unlink(ff) for ff in [fname, '{0:s}.cdf'.format(oname)]]

		if pipe is not None:
			t = [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)

	###############################
	#_cloud_thickness_sensitivity_#
	###############################
	#_WRS WHY IS THIS IN HERE TWICE
	elif 'thickness' in jacobian_vars: 
		dbg("MAKE SURE THERE IS A WAY FOR THICKNESS TO BE SHRUNK")

		#_make a copy to futz with here
		cloud = deepcopy(clddef)
		for cld in cloud:
			#_need to pass current layer base and top, get returned new base/top
			dz_orig = cld['z_top'] - cld['z']
			cld['z'], cld['z_top'], dz = increment_thickness(dir_lblrtm, **cld)

		#_change in thickness
		dz -= dz_orig

		label = 'thickness_{0:5.3f}_{1:07d}'.format(dz, pid)
		fname = path(dir_lblrtm, 'lbldis.input.' + label)
		oname = 'lbldis.out.{0:s}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kwargs)
 
		#call(['lbldis', fname, '0', oname])
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run'
		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve(wave, rads)
		##	wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_find third column of jacobians
		#_this version assumes that rads is from on instance
		tmp = (rads - apriori_rads) / dz 
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T

		jacobian['thickness'] = tmp 

		if pipe is not None:
			t = [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)

	################################
	#_effective_radius_sensitivity_#
	################################
	elif 'ssp_dist' in jacobian_vars: 
		if lidx:
			sys.exit(0)
		
		#_make a copy to futz with here
		cloud = deepcopy(clddef)
		for cld in cloud:
			cld['ref'] = cld['ref'] + sensitivity['ref']

		label = 'ref_{0:5.3f}_{1:07d}'.format(cloud[0]['ref'], pid)
		fname = path(dir_lblrtm, 'lbldis.input.{0}'.format(label))
		oname = 'lbldis.out.{0}'.format(label)
		lbldis_write(dir_lblrtm, filename=fname, clddef=cloud, **kwargs)
 
		#call(['lbldis', fname, '0', oname])
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run'
		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve(wave, rads)
		##	wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_this version assumes that rads is from one instance
		tmp = (rads - apriori_rads) / sensitivity['ref']
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T
		jacobian['ref'] = tmp 

		if clean:
			[os.unlink(ff) for ff in [fname, '{0:s}.cdf'.format(oname)]]

		if pipe is not None:
			t= [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)

	############################################################################
	#_state_mods_lblrtm_########################################################
	#_water_vapor_profile_#
	#######################
	elif 'rel_hum' in jacobian_vars:
		dbg('rel_hum {0}'.format(lidx))
		if lidx:
			dbg('exiting')
			sys.exit(0)

		from libtools import mkdir_p
		from lblrtm_utils import link_tape3, read_tape5, write_tape5
		from shutil import rmtree

		#_pull out shift
		delta = sensitivity['rel_hum']

		#_construct labels
		label = 'rel_hum_{0:5.3f}_{1:07d}'.format(delta, pid)

		#_create new LBLRTM run directory 
		dir_lblrtm_n = '{0}.{1}'.format(dir_lblrtm, label)

		#_read in old tape5
		prof, pres, t5kw = read_tape5(os.path.join(dir_lblrtm, 'TAPE5'))

		#_directory crap
		cwd = os.getcwd()
		mkdir_p(dir_lblrtm_n)
		os.chdir(dir_lblrtm_n)
		dbg("Changing directoris {0} {1}".format(cwd, dir_lblrtm_n))

		#_shift new tape5 and write
		varr = '{0:s}_{1:s}'.format(t5kw['profile_source'], 'relative_humidity')
		setattr(prof, varr, getattr(prof, varr) + delta)
		write_tape5(prof, pres, **t5kw)

		try:
			link_tape3(**kwargs)
		except OSError:
			dbg(('WARNING: tape3 link failure', dir_lblrtm, dir_lblrtm_n))

		#_run LBLRTM and move back into run directory
		dbg('running lblrtm')
		call(['lblrtm'])
		dbg('completed lblrtm')	
	
		#_lbldis inputs
		cloud = deepcopy(clddef)
		fname = path(dir_lblrtm_n, 'lbldis.input.' + label)
		oname = 'lbldis.out.{0:s}'.format(label)
		lbldis_write(dir_lblrtm_n, filename=fname, clddef=cloud, **kwargs)
 
	#	call(['lbldis', fname, '0', oname])
		dbg('calling lbldis')
		if call(['lbldis', fname, '0', oname]):
			raise RuntimeError, 'lbldis failed to run'
		dbg('lbldis_complete')

		ncdf = netcdf.netcdf_file(oname+'.cdf', 'r')
		rads = ncdf.variables['radiance'][:].squeeze()
		wave = ncdf.variables['wnum'][:]

		#_if not using microwindows, convolve to sensor
		if kwargs.get('dv') > 0:
			wave, tmp = convolve(wave, rads)
		##	wave, tmp = convolve2shis(wave, rads)
			rads = tmp.copy().squeeze()

		#_find third column of jacobians
		#_this version assumes that rads is from on instance
		tmp = (rads - apriori_rads) / delta 
		if tmp.ndim > 1 and tmp.shape[0] < tmp.shape[1]:
			tmp = tmp.T
		else:
			tmp = matrix(tmp).T

		jacobian['rel_hum'] = tmp 
	
		if clean:
			dbg('Think about how to clean up LBLRTM sections')
			[os.unlink(ff) for ff in [fname, '{0}.cdf'.format(oname)]]

		#_WRS 2015.09.14 moves this from below lblrtm call... unsure why
		# it was up there.
		os.chdir(cwd)

		#_delete tmp rel_hum directory
		if clean:
			rmtree(dir_lblrtm_n)
			
		if pipe is not None:
			t = [vv[1].shape for vv in jacobian.iteritems()]
			pipe.send((jacobian, matrix(apriori_rads).T, wave))
			pipe.close()
			os._exit(0)
	
	#_move back into calling directory
	os.chdir(dir_pwd)

	'''
	#_put jacobian in  
	if output_array and type(jacobian_vars) == list:
		#_loop over jacobian_vars and put into a 2D array
		for i, sv in enumerate(jacobian_vars):
			if not i:
				tmp = jacobian[sv][:,None]
			else:
				tmp = append(tmp, jacobian[sv][:,None], axis=1)

		#_update return value
		jacobian = tmp
	'''
	return jacobian, matrix(apriori_rads).T, wave


def clear_sky_radiance(dir_lblrtm, clddef=None, lbldis_clear='lbldis.clear.cdf',
	sensor='AIRS', **kwargs):
	'''
	dummy run to drop aerosol / cloud layers from retrieval

	Produces a file with the expected name of lbldis_output.clear.cdf
	that is commonly used in other places as a clear sky bounding for plots.
	
	lbldis_clear	str,	name of clear file output
	'''
	from scipy.io import netcdf
	from lblrtm_utils import write_lbldis_parameters as lbldis_write
	from numpy import zeros
	from convolution import sensors

	convolve = sensors[sensor]

	cwd = os.getcwd()
	res = os.chdir(dir_lblrtm)
##	dbg(('TES', dir_lblrtm, cwd, 'res', res)); exit()

	#_get a label from lbldis clear
	label = '.'.join(lbldis_clear.split('.')[1:-1])

	#_create a name for this run based upon state variables tested
	fname = os.path.join(dir_lblrtm, 'lbldis_input.' + label)
	oname = 'lbldis_output.{0}'.format(label)
	lbldis_write(dir_lblrtm, filename=fname, clddef=[], **kwargs)

	if call(['lbldis', fname, '0', oname]):
		raise RuntimeError, 'lbldis failed to run'
	ncdf = netcdf.netcdf_file(oname+'.cdf','r')
	rads = ncdf.variables['radiance'][:] #_[nchan, ninstances]
	wave = ncdf.variables['wnum'][:]

	#_if not using microwindows, convolve to sensor
	if kwargs.get('dv') > 0:
		wave_tmp = wave.copy()
		wave, tmp0 = convolve(wave_tmp, rads[:,0])
	##	wave, tmp0 = convolve2shis(wave_tmp, rads[:,0])

		n_rads = zeros((wave.size, 2))
		n_rads[:,0] = tmp0.copy().squeeze()
		rads = n_rads.copy()

	os.chdir(cwd)
	return rads, wave 


def wavenum_mismatch(obs_wave, rad_wave, rads):
	'''
	I introduced some crappiness when the mw channels are in 
	dead sensor zones, so they end up being dropped, which messes
	up some matrix sizes.  So I should really just limit which
	MW channels we can use for these and now go nuts here.  Hold
	on.
	'''
	pass


def lbldis_final(dir_lblrtm, clddef=None, clean=True,
	lbldis_output='lbldis_output.final.cdf', **kwargs):
	'''
	dummy run to create final retrieval cdf file

	Produces an file of the lbldis_output.final.cdf that is
	expected elsewhere.  This should probably be updated to 
	include out_label.
	'''
	from lblrtm_utils import write_lbldis_parameters as lbldis_write

	pid = os.getpid()
	cwd = os.getcwd()
	os.chdir(dir_lblrtm)

	#_extract label
	label = '.'.join(lbldis_output.split('.')[1:-1])

	#_create a name for this run based upon state variables tested
	fname = os.path.join(dir_lblrtm, 'lbldis_input.{0}'.format(label))
	oname = 'lbldis_output.{0}'.format(label)
	lbldis_write(dir_lblrtm, filename=fname, clddef=clddef, **kwargs)

	if call(['lbldis', fname, '0', oname]):
		raise RuntimeError, 'lbldis failed to run'

	os.chdir(cwd)


def jacobian_forker(dir_lblrtm, jacobian_vars=['tau', 'ref'],  **kwargs):
	'''
	Treads each column of jacobian matrix
	dir_lblrtm		str,	location of parent lblrtm output for LBL-DIS
	jacobian_vars	list,	list of state variables for columns of jacobian

	returns
	K		dict,	jacobians		
	'''
	from libtools import setup_groups
	from multiprocessing import Process, Queue, Pipe
	from numpy import append
	dbg(jacobian_vars)

	#_make list of state vars without tau
	state = jacobian_vars[:]
	try:
		del state[state.index('tau')]
	except ValueError:
		pass

	#_Ki{sv}[lidx] || Ki[lidx]{sv}
	nlayer = len(kwargs.get('clddef'))

	#_start process for each column of jacobian
	args	= (dir_lblrtm,)
	nv		= len(state) 
	thread	= [None]*nv
	pipe	= [None]*nv
	pope	= [None]*nv

	#_initialize output for change, a priori values, and wavenumbers
	K    = [{} for i in range(nlayer)]
	rads = [{} for i in range(nlayer)]
	wave = [{} for i in range(nlayer)]

	#_need to get apriori_rads from first tau case, launch that first
	#_dict, matrix, array
	for lidx in range(nlayer):
		dbg('Getting tau a apriori radiances: nlayer {0}'.format(nlayer))
		K[lidx], rads[lidx]['tau'], wave[lidx]['tau'] =	\
			jacobians(dir_lblrtm, lidx=lidx, jacobian_vars=['tau'], **kwargs)

	apriori = rads[0]['tau'].copy()

	#_loop over each layer
	for lidx in range(nlayer):

		#_setup fork groups
		groups = setup_groups(state, **kwargs)
		for group in groups:

		  #_launch threads for non-Tau cases
 		  for i, sv in enumerate(group):

			#_skip first layer's tau
			#if (lidx and sv[:4] == 'surf'):
			if (not lidx and sv == 'tau')	\
			or (lidx and (sv[:4] == 'surf' or sv == 'rel_hum')):
				continue

			#_open communication pipe
			pipe[i], pope[i] = Pipe()

			#_launch subprocess
			kwargs.update({	'pipe'			: pope[i],
							'lidx'			: lidx, 
							'jacobian_vars'	: [sv], 
							'apriori_rads'	: apriori.A.squeeze() })
			thread[i] = Process(target=jacobians, args=args, kwargs=kwargs)
			thread[i].start()
	
		  #_collect output
 		  for i, sv in enumerate(group):

			#_skip first layer's tau
			if (not lidx and sv == 'tau') \
			or (lidx and (sv[:4] == 'surf' or sv == 'rel_hum')):
				continue

			#_get data
			Ktmp, rads[lidx][sv], wave[lidx][sv] = pipe[i].recv()
			K[lidx].update(Ktmp)
	
			#_wait for it to complete
			thread[i].join()

	#_set return vals... to just the last one in the list?  seems bad.
	wavenum = wave[0]['tau']
	apriori = rads[0]['tau']

	#_put in desired format
	return K, apriori, wavenum
	

################################################################################
#_tools_########################################################################
################################################################################


def check_settings(model=None, uncert_vars=[], uncertainty={},
	surf_temp=None, **kwargs):
	''' 
	a not even nearly complete dumping ground to make sure
	that the settings set in namelists make sense.
	'''

	#_if model is not False, then profile shifts cannot be performed
	if model and 'rel_hum' in uncert_vars:
		raise ValueError, 'Cannot test uncertainty on profile jacobians'
	elif not model and 'rel_hum' not in uncert_vars:
		raise ValueError, 'You can test this, but let\'s not!'

	#_don't allow for flat uncertainty to be added
	if uncertainty['model'] != 0:
		raise ValueError, "What the hell are you doing?"

	#_check that surface temp matches things
	if model != 0 and surf_temp != 'LBL-ATM':
		raise ValueError, 'Use LBL-ATM with standard model profiles'
	elif model == 0 and surf_temp == 'LBL-ATM':
		raise ValueError, 'User GDAS or ECMWF with user defined profile'


def check_state_boundaries(x, state_vars, sensitivity, style='max', 
	clddef=[], **kwargs):
	'''
	Looks through current state values against permitted
	values in the SSP file and other constraints such
	as observation height.
	'''
	from lblrtm_utils import check_reff_v_ssp as check_ssp
	from tape7 import tape7

	nlayer = len(clddef)
	nstate = len(state_vars)

	#_do check for each layer
	for lidx in range(nlayer):

		#_check optical depth
		if 'tau' in state_vars:
			tidx = state_vars.index('tau') + lidx * nstate
			x[tidx,-1] = 0.0001 if x[tidx,-1] < 0 else x[tidx,-1]
	
		#_check effective radius
		if 'ref' in state_vars:
			ridx = state_vars.index('ref') + lidx * nstate	
 
			x[ridx,-1], style = check_ssp(x[ridx,-1],
				delta=sensitivity['ref'],
				dbnum=clddef[lidx]['dbnum'], **kwargs) 

			dbg(('SSP_CRASH_CHECK_00', x[ridx,-1], style))
			if style == 'max':
				x[ridx,-1] -= sensitivity['ref']*2.
			dbg(('SSP_CRASH_CHECK_01', x[ridx,-1], style))
	
		#_check height of layer top and bottom
		if 'z' in state_vars:
			max_z = tape7(os.path.join(dir_lblrtm, 
					'TAPE7')).H1-sensitivity['z']*1.5
	
			zidx = state_vars.index('z') + lidx * nstate
			x[zidx,-1] = 0.01 if x[zidx,-1] < 0 else x[zidx,-1]
			x[zidx,-1] = max_z if x[zidx,-1] > max_z else x[zidx,-1]
	
		#_make sure bimodal ssp dist fraction is not outside [0,1]
		if 'ssp_dist' in state_vars:
			tidx = state_vars.index('ssp_dist')
			if x[tidx,-1] < 0:
				x[tidx,-1] = 0
			elif x[tidx,-1] > 1:
				x[tidx,-1] = 1

	return x


def plot_close(pname=None, fig=None, **kwargs):
	''' call to plot figure in case not plotted elsewhere '''
	import matplotlib.pyplot as plt
	import os
	if not os.path.exists(pname):
		dbg(pname)
		fig.tight_layout()
		plt.savefig(pname)


def plot(r_obs, n_obs, r_lbl, n_lbl, fig=None, clddef={}, surf_temp=0, dv=0,
    label='', last=False, r_clr=None, n_clr=None, sup=None,
    ax=None, **kwargs):
    '''
    create plot of r values in brightness temperature space
    To have multiple on a single plot, pass axes through ax
    '''
    import os
    import matplotlib.pyplot as plt
    from libgeo import planck_inv
    from numpy import array, tile, linspace
    from libtools import shrink_ticks
    import os

    #_check if clear sky passed
    clear = False if r_clr is None else True
    first = True if len(ax.lines) == 0 else False

    #_kill matrices
    r_obs = r_obs if r_obs.ndim == 1 else array(r_obs).squeeze()
    r_lbl = r_lbl if r_lbl.ndim == 1 else array(r_lbl).squeeze()

    try:
        tau = clddef[0]['tau'][0]
    except:
        tau = clddef[0]['tau']

    ref = clddef[0]['ref']
    zkm = clddef[0]['z']
    ztp = clddef[0]['z_top']
    thk = ztp - zkm
    tit = label + ' '
    tit += '{0:4.2f} tau, {1:4.2f} ref, {2:4.2f} z, {3:4.2f} thk '.format(tau,
                                                             ref, zkm, thk)
    tit += '(red = obs, green = lbl)'

    #_check if in brightness temperature space
    t       = tile(surf_temp, n_obs.size)
    t_obs   = planck_inv(r_obs/1e3, n_obs*100, domain='wavenumber')
    t_lbl   = planck_inv(r_lbl/1e3, n_lbl*100, domain='wavenumber')

    #_line settings
    arg = { 'linewidth' : 0.3 }
    if dv < 0:
        arg.update({'marker' : 'x' })

    #_clear and obs should be plotted first, then skipped for alpha
    if first:
        #_plot obs in red
        ax.plot(n_obs, t_obs, color='r', **arg)

        #_if clear sky radiances are passed
        if clear:
            #_output from lbldis is many flavored if tau vs taus given
            if r_clr.ndim != 1:
                r_clr = array(r_clr).squeeze()

            #_plot clear sky radiances in blue
            t_clr = planck_inv(r_clr/1e3, n_clr*100, domain='wavenumber')
            ax.plot(n_clr, t_clr, color='#CCCCFF', **arg)

    #_plot this iteration's radiances
    ax.plot(n_lbl, t_lbl, color='g', **arg)

    #_if last, setup axes and the like
    if last:
        #_loop over lines, skip clear sky and obs, set alphas
        iskip = 2 if clear else 1
        nline = len(ax.lines) - iskip
        alphas = linspace(0.9, 1, nline)

        #_set supertitle if not set
        if fig is not None and sup is not None:
            fig.suptitle('{0}, n={1:d}'.format(sup, nline) , size='xx-small')

        ax.set_ylabel('n={0:d}'.format(nline))

        dbg(("LAST", nline))

        for i, alpha in enumerate(alphas):
            if i == nline - 1:
                continue 

            #_fade old lines
            plt.setp(ax.lines[i+iskip], 'alpha', alpha)
            plt.setp(ax.lines[i+iskip], 'linestyle', '-.')
            plt.setp(ax.lines[i+iskip], 'marker', None)

        #_shrink axis tick marks
        shrink_ticks(ax)
        ax.grid(True)
        ax.set_xlim(n_lbl[0]-20, n_lbl[-1]+20)
        ax.set_title(tit, size='xx-small')

        #_do something where the list of values tried is somewhere


def plot_delete(r_obs, n_obs, r_lbl, n_lbl, fig=None, clddef={}, surf_temp=0, dv=0,
	pname='default.png', label='', last=False, r_clr=None, n_clr=None, sup=None,
	**kwargs):
	'''
	create plot of r values in brightness temperature space
	To have multiple on a single plot, pass axes through ax
	'''
	import os
	import matplotlib.pyplot as plt
	from libgeo import planck_inv
	from numpy import array, tile
	from libtools import shrink_ticks
	import os

	#_check if clear sky passed
	clear = False if r_clr is None else True

	#_kill matrices
	r_obs = r_obs if r_obs.ndim == 1 else array(r_obs).squeeze()
	r_lbl = r_lbl if r_lbl.ndim == 1 else array(r_lbl).squeeze()

	try:
		tau = clddef[0]['tau'][0]
	except:
		tau = clddef[0]['tau']

	ref = clddef[0]['ref']
	zkm = clddef[0]['z']
	ztp = clddef[0]['z_top']
	thk = ztp - zkm 
	tit = label + ' '
	tit += '{0:4.2f} tau, {1:4.2f} ref, {2:4.2f} z, {3:4.2f} thk '.format(tau,
															 ref, zkm, thk)
	tit += '(red = obs, green = lbl)'

	#_check if in brightness temperature space
	t		= tile(surf_temp, n_obs.size)
	t_obs	= planck_inv(r_obs/1e3, n_obs*100, domain='wavenumber')	
	t_lbl	= planck_inv(r_lbl/1e3, n_lbl*100, domain='wavenumber')	

	if clear:
		if r_clr.ndim != 1:
			r_clr = array(r_clr).squeeze()

		t_clr = planck_inv(r_clr/1e3, n_clr*100, domain='wavenumber')

	#_find which axes are empty
	i = 0
	for i, ax in enumerate(fig.axes):
		#_save final plot for solution
		if last:
			ax = fig.axes[-1]

			#_set supertitle if not set
			if fig is not None and sup is not None:
				fig.suptitle(sup, size='xx-small')

			break	
		elif not len(ax.lines) and i != len(fig.axes)-1: #_don't fill last
			#_plot empty, get out of loop
			break
		else:
			#_line not empty, continue
			continue
	else:
		
		return 

	#_plot new brightness temperatures and difference between last
	if dv > 0:
		arg = { 'linewidth' : 0.3 }
		ax.plot(n_obs, t_obs, color='r', **arg)
		ax.plot(n_lbl, t_lbl, color='g', **arg)
		if clear:
			ax.plot(n_clr, t_clr, color='#CCCCFF', **arg)

	elif dv < 0: #_only plot for microwindows
		arg = { 'linewidth' : 0.3, 'marker' : 'x' }
		ax.plot(n_obs, t_obs, color='r', **arg)
		ax.plot(n_lbl, t_lbl, color='g', **arg)
		ax.scatter(n_obs, t_obs, color='r', **arg)
		ax.scatter(n_lbl, t_lbl, color='g', **arg)
		if clear:
			ax.plot(n_clr, t_clr, color='#CCCCFF', **arg)

##	ax.plot(n_obs, t, color='k', **arg)
	#_if going to do differences, match them up
		
	#_set labels
	shrink_ticks(ax) 
	ax.grid(True)
	ax.set_xlim(n_lbl[0]-20, n_lbl[-1]+20)
	ax.set_title(tit, size='xx-small')



def increment_thickness(dir_lbl, z=None, z_top=None, **kwargs):
	'''
	Expands passed layer the smallest geometric amount based
	on atmosphere in lblrtm directory

	dir_lbl	str,	path to lblrtm output
	z		flt,	layer base
	z_top	flt,	layer top
	'''
	from tape7 import tape7
	from numpy import append, diff

	atm = tape7(os.path.join(dir_lbl, 'TAPE7'))
	
	zmin = atm.zlevel1.min()
	zmax = atm.zlevel2.max()

	#_get max index
	imax = atm.zlevel1.size - 1

	#_get thickness and midlevel points
	depths = append(diff(atm.zlevel1), diff(atm.zlevel2)[-1])	
	midlev = atm.zlevel1 + depths/2

	#_get current level indices, idx0 will correspond to depth/midlev, idx1 not
	idx = range(atm.zlevel1.size)
	idx0 = abs(z - atm.zlevel1).argmin()
	idx1 = abs(z_top - atm.zlevel2).argmin()
	idx0 = idx0 if z > atm.zlevel1[idx0] else idx0-1
	idx1 = idx1 if z_top <= atm.zlevel2[idx1] else idx1+1

	#_figure out if we're going up or down
	#_First, are we at the edge already?
	if idx0 == 0 and idx1 == imax: 
		return z, z_top, atm.zlevel2[-1] - atm.zlevel1[0]
	if idx0 == 0:
		direction = 'up'
	elif idx1 == imax:
		direction = 'down'
	else: #_expand in shortest direction
		direction = 'up' if depths[idx0-1] > depths[idx1+1] else 'down'

	#_expand
	if direction == 'up':
		idx1 += 1
		z_top = midlev[idx1] #atm.zlevel2[idx1+1]
		z = midlev[idx0]
	elif direction == 'down':
		idx0 -= 1
		z_top = midlev[idx1]
		z = midlev[idx0] #atm.zlevel1[idx0-1]

	#_get new geometric depth
	dz = atm.zlevel2[idx1] - atm.zlevel1[idx0]

	return z, z_top, dz


def add_ssps(apriori=None, **kwargs): 
	''' pull out ssp files form apriori and put them in kwargs '''

	#_setup path to ssp databases
	dir_ssp = os.path.join(os.environ['PRODUCTS'], 'lbldis_inputs')

	#_pull out all ssp files from layers
	ssps = [layer['ssp_db'] for layer in apriori]

	#_add ssp file to list of all (add a check later to prevent dupes)
	[kwargs['ssp_db_files'].append(os.path.join(dir_ssp, f)) for f in ssps]

	#_add database number associated with ssp
	[layer.update({'dbnum' : ssps.index(layer['ssp_db'])}) for layer in apriori]


def dbg(msg, l=1, err=False):
	''' 
	if global debug is set to true, be more verbose 
	msg : str, Message to be printed
	l   : int, Debug level of message.  Set higher for lower level 
	        messages.  As debug increases, noisiness should also.
	'''
	import inspect
	import re

	if hasattr(msg, '__iter__'):
		msg = ' '.join([str(m) for m in msg])
	elif type(msg) != str:
		msg = str(msg)

	if DEBUG >= l:
		curf    = inspect.currentframe()
		calf    = inspect.getouterframes(curf,2)
		file, line, method = calf[1][1:4]
		file     = '.'.join(file.split('/')[-1].split('.')[:-1])
		scream  = '[%s.%s.%i] %s' % (file,method,line,msg) \
			if not re.search('^\s*$', msg) else ' ' 

		if not err:
			print scream
		else:
			raise RuntimeError, scream


################################################################################
#_SHELL_########################################################################
################################################################################


if __name__ == '__main__':
	dbg('Use wrapper scripts such as oe.real|sim.py')
