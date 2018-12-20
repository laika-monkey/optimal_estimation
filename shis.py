def convolve2shis(v_mono, r_mono, zflag=0, v_rolloff=None, v_laser=15799.0,
	**kwargs): 
	'''
	function [nu,spec] = aeri(v_mono,r_mono,zflag,pflag,v_laser,v_rolloff)

	apply AERI01 (SGP CART central fac.) instrument function to a monochromatic 
	calculation.  FFOV self apodization not included.

	Input variables:
		v_mono are the input monochromatic wavenumbers
		r_mono are the input monochromatic radiances in units of W/(m^2 ster cm^-1)

	zflag determines if the apodized interferogram is truncated or zero filled.
	if zflag == 1, the interferogram points past MOPD are retained, but zero 
	filled and the output spectrum has a finer point spacing.  
	if zflag == 0 (the default), the interferogram is trucated past MOPD 
	and the output spectrum has the minimum wavenumber point spacing.

	v_laser is the metrology laser wavenumber (typically 15799.0 1/cm)

	v_rolloff is a two element array defining where the boundaries
	(in 1/cm) of where the rolloff function is applied.  See
	rolloff.m for info.

	Output variables:
		nu are the wavenumbers of the apodized spectrum
		spec are the apodized radiances [W/(m^2 ster cm^-1)]

	External functions used:
		rolloff.m (appended)

	DCT, CIMSS/SSEC/UW-Madison
	'''	
	from numpy.fft import fft
	from numpy import array, matrix, matrixlib, append, zeros, real, linspace
	from scipy.interpolate import splev, splrep, interp1d

	#_meh, get it going. WRS 
	r_mono = array(r_mono).squeeze()
	return_matrix = True
##	if type(r_mono) == matrixlib.defmatrix.matrix:
##		return_matrix = True
##		r_mono = array(r_mono).squeeze()
##	else:
##		return_matrix = False	

	#_check radiance values
	if r_mono.max() > 100.:
		print 'MAKE SURE RADIANCE VALUES IN W/m2/cm!  Values appear high!'	

	#_DEFAULT VALUE FOR ROLL-OFF WIDTH IS 5 PERCENT OF BANDWIDTH
	if v_rolloff is None:
		v_rolloff = [None, None]
		v_rolloff[0] =  0.05*(v_mono.max()-v_mono.min()) + v_mono.min()
		v_rolloff[1] = -0.05*(v_mono.max()-v_mono.min()) + v_mono.max()

	#_APPLY COSINE ROLL-OFFS TO TAKE INPUT SPECTRUM SMOOTHLY TO ZERO AT ENDS
	r_mono = rolloff(v_mono, r_mono, v_rolloff)

	# CREATE A WAVENUMBER SCALE WHICH GOES FROM 0 CM^-1 TO THE NYQUIST
	# WAVENUMBER (V_NYQUIST = V_LASER/2/REDUCTION_FACTOR) WITH [(2^N) +1] 
	# POINTS.  (CHOOSE N SUCH THAT THE WAVENUMBER POINT SPACING IS MUCH 
	# SMALLER THAN THE FINAL APODIZED WAVENUMBER POINT SPACING).
	red_factor = 1		#_samples every nth fringe
	nsamples = 2**14	#_number of points in the (single sided) interferogram 
						# from 0 to MOPD cm
	n = 20
	npts = 2**n + 1
	wnbr = linspace(0, v_laser/2./red_factor, npts)

	#_interpolate to new wavenumber scale
	rad = zeros((npts,))
	index = append(	(wnbr >= min(v_mono))[:,None],
					(wnbr <= max(v_mono))[:,None], axis=1).all(axis=1)
	
	#_interpolateeee
	f1d = splrep(v_mono, r_mono, k=3)
##	f1d = interp1d(v_mono, r_mono, kind='cubic')
##	rad = splev(wnbr[index], f1d)
	rad[index] = splev(wnbr[index], f1d)

	# NOW, CREATE A SPECTRUM WHICH ALSO CONTAINS THE NEGATIVE WAVENUMBER
	# COMPONENTS.  TO DO THIS, GRAB THE 2ND THROUGH 2^N POINTS (GOING FROM 
	# DV CM^-1 TO (V_NYQUIST-DV) CM^-1 AND APPEND THEM (IN DESCENDING INDEX 
	# ORDER) TO THE ORIGINAL (2^N)+1 POINTS, (GOING FROM 0 TO V_NYQUIST CM^-1).
	# THE RESULTING VECTOR (WHICH WILL BE FFT'D) HAS [(2^N)+1]+[(2^N)-1]=
	# 2^(N+1) POINTS AND IS REAL AND EVEN.  THE IMPORTANT POINT (NO PUN 
	# INTENDED) HERE IS THAT THE 0 CM^-1 AND V_NYQUIST CM^-1 POINTS ARE UNIQUE 
	# AND DO NOT GET REPEATED.  (THE FACT THE NEGATIVE WAVENUMBER POINTS ARE 
	# PLACED TO THE RIGHT OF THE POSITIVE WAVENUMBER POINTS IS JUST MATLAB'S 
	# WAY OF INTERPRETING THE DATA, JUST AS IT ASSOCIATES THE FIRST POINT OF 
	# AN FFT WITH THE ZERO OPTICAL PATH DIFFERENCE POINT).

	rad2 = append(rad[:npts], rad[npts-1:1:-1])
	npts2 = 2**(n+1)

	# TAKE THE COMPLEX FFT THE SPECTRUM
	yfft = fft(rad2)

	# FOR THE APODIZED SPECTRUM WITH THE MINUMUM NUMBER OF POINTS:
	# TRUNCATE THE INTERFEROGRAM TO THE DESIRED NUMBER OF POINTS AND FFT IT.  
	# FOR THE TRUNCATED INTERFEROGRAM, EXTRACT THE FIRST [(2^M)+1] POINTS 
	# (ASSOCIATED WITH THE 0 CM TO MAX OPTICAL PATH DELAY CM POINTS) AND THE 
	# LAST [(2^M)-1] POINTS OF THE ORIGINAL INTERFEROGRAM.  FFT'ING THIS YIELDS 
	# AN APODIZED SPECTRUM WHICH CONTAINS THE POSITIVE AND NEGATIVE WAVENUMBER 
	# COMPONENTS WITH [2^(M+1)] POINTS.  PULL OUT THE FIRST [(2^(M)+1)] POINTS 
	# OF THE SPECTRUM TO GET THE POINTS ASSOCIATED WITH THE POSITIVE WAVENUMBERS, 
	# AND CREATE A NEW WAVENUMBER SCALE GOING FROM 0 CM^-1 TO V_NYQUIST CM^-1 
	# WITH [(2^M)+1] POINTS.

	if not zflag: 
		yfft1 = zeros((2*nsamples,))
		yfft1[:nsamples] = yfft[:nsamples]
		yfft1[nsamples+1:2*nsamples] = yfft[npts2-nsamples+1:npts2]
		spec = fft(yfft1)/npts2
		spec = real(spec[:nsamples+1])
		nu = linspace(0, v_laser/2/red_factor, nsamples+1)

	# FOR THE APODIZED SPECTRUM WITH INTERPOLATION BETWEEN MINIMUM INFORMATION
	# POINTS USING A ZERO PADDED INTERFEROGRAM: JUST FOLLOW THE SAME
	# PROCEDURE AS ABOVE, BUT DO NOT REDUCE THE NUMBER OF POINTS IN THE 
	# TRUNCATED INTERFEROGRAM BUT SET THE VALUES TO ZERO IF THEY EXCEED THE 
	# MAXIMUM OPTICAL PATH DIFFERENCE POINT.
	
	elif zflag:
		yfft2 = yfft.copy()
		yfft2[nsamples+1:npts2-nsamples+1] \
			= yfft2[nsamples+1:npts2-nsamples+1]*0	#_wrs ????
		spec = fft(yfft2)/npts2
		spec = real(spec[:npts])
		nu = wnbr.copy()

	# EXTRACT PORTION OF SPECTRUM WITHIN ROLL-OFF BOUNDARIES
	idx = append(	(nu >= v_rolloff[0])[:,None],
					(nu <= v_rolloff[1])[:,None],
					axis=1	).all(axis=1)

	nu = nu[idx]
	spec = spec[idx]

	#_kludge WRS
	if return_matrix:
		spec = matrix(spec)

		#_ok flip it
		ooo, iii = spec.shape
		if ooo < iii:
			spec = spec.T

	return nu, spec #_wavenumbers, radiances 


def rolloff(vin, rin, vb=None):
	'''
	function rout = rolloff(vin,rin,vb)

	Apply cosine roll-offs to edges of an input spectrum

	Inputs:
	[vin,rin]: input wavenumbers,radiances
	vb       : [2x1] boundaries where roll-off function is applied
             DEFAULT values are computed using 5% of the input
             bandwidth.
  
	Outputs:
	 rout: modified radiances

	The input spectrum for v < v_rolloff(1) and for v > v_rolloff(2)
	is replaced with the appropriate cosine function.

	DCT 2-12-99
	
	_Make this smarter about matrices versus arrays, and check
	for shape and dimensional order.
	'''
	from numpy import array, float32, where, sin, cos, pi

	#_DEFAULT VALUE FOR ROLL-OFF WIDTH IS 5 PERCENT OF BANDWIDTH
	if vb is None: 
		vb[0] =  0.05*(max(vin)-min(vin)) + min(vin)
		vb[1] = -0.05*(max(vin)-min(vin)) + max(vin)

	rout = array(rin.copy()).squeeze()

	#_for lower wavenumber end
	ind = where(vin < vb[0])[0]
	tmp = ind-min(ind)
	tmp = tmp/max(tmp)
	a =  rin[max(ind)] 
	b =  (sin(tmp*pi-pi/2, dtype=float)*0.5+0.5)
	rout[ind] = rin[max(ind)] * (float32(sin(tmp*pi-pi/2))*0.5+0.5)

	#_for upper wavenumber end
	ind = where(vin > vb[1])[0]
	tmp = ind-min(ind)
	tmp = tmp/max(tmp)
	rout[ind] = rin[min(ind)-1] * (float32(cos(tmp*pi))*0.5+0.5)

	return rout
