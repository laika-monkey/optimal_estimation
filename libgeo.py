#_constants
C 	= 2.99792458e8	#_Speed of Light		[m/s]
H	= 6.6260957e-34	#_Planck's Constant		[m2kg/s]
K	= 1.381e-23		#_Boltzmann's Constant	[J/K]
Rd	= 289.058		#_gas constant, dry		[J/kg/K]
Rv	= 461.91		#_gas constant, moist	[J/kg/K]
R	= 8.3144621		#_universal gas costant	[J/kg/mol]
N	= 6.0221413e23	#_Avogadro Constant		[molecules/mole]
Re	= 6.378e6		#_radius of earth (avg)	[m]
G	= 9.81			#_sfc gravitational acc [N/kg]

from numpy import pi as PI
from libtools import dbg

#_standard temperature conversions
def F2C(T): return (T - 32) * 5. / 9.
def C2F(T): return T * 9./5. + 32		#_double, subtract a tenth, add 32
def C2K(T): return T + 273.15
def F2K(T): return C2K(F2C(T))
def F2R(T): raise RuntimeError, 'fuck you'


#_standard speed conversions
def ms2kn(spd): return spd*1.94384
def kn2ms(spd): return spd/1.94384
def ms2mph(spd): return spd*2.23694
def mph2ms(spd): return spd/2.23694
def kph2mph(spd): return spd*.621371
def mph2kph(spd): return spd/.621371
def kn2mph(spd): return spd*1.15078
def mph2kn(spd): return spd/1.15078


#_radiation, wavenumber functions are not ANGULAR wavenumbers
def wavenumber2length(k): return 1./k/100	#_K IS IN CM-1
def wavelength2frequency(l): return C/l
def wavenumber2frequency(k): return wavelength2frequency(wavenumber2length(k))

#_convert to sigma coords
def pres2sig(p, psfc=1013.15, ptop=100, **kwargs): return (p-ptop)/(psfc-ptop)
def sig2pres(sig, psfc=1013.15, ptop=100, **kwargs): return sig*(psfc-ptop)+ptop

class molecules(object):
	def __init__(self):
		self.retard = 'strength'
		molecules = {
			'ozone' 	: {
				'molecular_mass' : .0479982, 	},		#_kg/mol}
			'benzene' 	: {
				'molecular_mass' : .07811,		},
		}		

		for molecule, attributes in molecules.items():
			self.__setattr__(molecule,attributes)

			
############################################################################_80
#_FUNCTIONS_###################################################################
###############################################################################


def planck(T, dom, rayleigh_jeans=False, domain='wavelength'):
	'''
    PASS SI UNITS, RETURNS NON-SI
    dom,            float,      wavelength of information {m|m-1|s}
    T,              float,      temperature of blackbody
    rayleigh-jeans, bool,       for u-wave, can use R-J approx
    domain          string,     what is the spectral domain
                                {wavenumber|length|frequency}
    
    returns 
    radiance,       float|ndarray,    W/sr/m/{um|cm-1|s}
    '''
	from numpy import exp

	#_calculate the emission as a function of wavelength and temp
	if domain == 'wavelength':
		if not rayleigh_jeans:
		    e = exp(H*C/K/dom/T)
		    B = 2*H*C*C/dom**5/(e-1.)
		else:
		    B = 2*C*K*T/dom**4			
		B *= 1e-6 						#_W/m2/sr/um
	elif domain == 'frequency':
		e = exp(H*dom/K/T)
		B = 2*H*dom**3/C/C/(e-1.)		#_W/m2/sr/s
	elif domain == 'wavenumber':
		e = exp(H*C*dom/K/T)
		B = 2*H*C*C*dom**3/(e-1)		
		B *= 100 						#_W/m2/sr/cm
	else:
		raise RuntimeError,	'invalid domain space ' + domain

	return B
	

def planck_inv(B, doms, domain='wavelength'):
	'''
	PASS IN NON-SI (per um, cm-1)
	B,              float(N),       radiances (W/m2/sr/{um|s|cm-1})
	lamb,           float(N),       wavelength (m)
	rayleigh-jeans, bool,           for u-wave, can use R-J approx
	domain,         string,         'wavelength' or 'frequency'
	
	returns
	Brightness Temperature, float|ndarray,	K
	'''
	from numpy import array, log
	dom = doms.copy()

	if B.ndim != 1:
		B = array(B).squeeze()

	#_for now
	if B.flatten().size != dom.flatten().size:
		raise RuntimeError, 'sizes must match'
	dom[dom == 0] = 1e-34

	if (B < 1e-50).sum():
		#_very small or negative radiance values being truncated'
		B[B < 1e-50] = 1e-50

	#_inverse planck function
	if domain == 'wavelength':
		T = H*C/dom/K/log(1.+2.*H*C*C/B/1.e6/dom**5.)		
	elif domain == 'frequency':
		T = H*dom/K/log(1.+2.*H*dom**3./C**2./B)
	elif domain == 'wavenumber':
	##	if dom.max() < 1e5:
	##		dbg('WARNING: Wavenumbers may not be in proper units (m)')	
		T = H*C*dom/K/log(1.+2.*H*C*C*dom**3./(B/1.e2))
	else:
		raise RuntimeError, 'invalid domain space ' + domain

	return T


def beer_lamber(tau, mu=0):
	'''
	t=exp(-tau)
	tau		float,	optical depth of layer
	mu		float,	zenith angle through layer in radians from vertical
	returns:
	transmittance	float,	
	'''
	from numpy import exp, cos
	return exp(-tau/cos(mu))

	
def size_parameter(r,l):
	'''
	returns size parameter for given particle radii and wavelengths
	r,	ndarray(flt),	radii of particle
	l,	ndarray(flt),	wavelengths of radiation
	
	x << 1,		rayleigh regime
	x ~ 1-100	mie regime, techincally requires spherical particle
	x > 100		geometric ray tracing
	'''
	from numpy import array, pi
	return 2.*pi*(array(r)/array(l))

	
def size_parameter_plot(lrng=(1e-7,1e-1),rrng=(1e-10,1e-1),logscale=True):
	'''
	produces size parameter contour plt
	'''
	from numpy import linspace, meshgrid
	import matplotlib.pyplot as plt
	from matplotlib import cm, ticker
	import matplotlib

	l = linspace(lrng[0],lrng[1],1e6)	#_wavelengths
	r = linspace(rrng[0],rrng[1],100)	#_radii

	#_kills memory
	l2, r2 = meshgrid(l,r)		#_make two dimensional arrays
	x = size_parameter(r2,l2)	#_calculate 2d size parameter

	#_initialize plotting area
	ax = plt.figure().add_subplot(111)
	
	#_contour levels
	levels = [0,.002,.2,2000]	#_Geometric Optics, Mie, Rayleight, Neg[::-1]
	cn = ax.contourf(l,r,x,locator=ticker.LogLocator(),cmap=cm.jet)
	cx = ax.contour(l,r,x,levels=levels,linestyles='--')
	
	#_plot colorbar and inline labels for noted regimes
	plt.colorbar(cn)
	plt.clabel(cx,use_clabeltext=True)
	
	#_set axes limits and style
	ax.set_xlim(lrng); ax.set_ylim(rrng)
	ax.set_xscale('log'); ax.set_yscale('log')
	
	if logscale:
		loc = [	1e-7, 1e-6,  1e-5,   1e-4, 1e-3, 1e-2,  1e-1]
		lab = ['0.1um','1um','10um','100um','1mm','1cm','10cm']
		ax.set_xticks(loc); ax.set_xticklabels(lab)
	
		loc = [	 1e-9,  1e-8,   1e-7, 1e-6,  1e-5,   1e-4, 1e-3, 1e-2]
		lab = [ '1nm','10nm','0.1um','1um','10um','100um','1mm','1cm']
		ax.set_yticks(loc); ax.set_yticklabels(lab)
	
	ax.set_xlabel('Wavelength'); ax.set_ylabel('Radius')
	ax.set_title('Petty Figure 12.1: Size Parameter')
	plt.savefig('figure_12.1.png')

	
def rh2wv(rh,t,p):
	'''
	convert from relative humidity to vapor mass mixing ratio
	rh	float,	relative humidity, percent
	T	float,	temperature, K
	P	float,	pressure, Pa
	
	returns:
	w	float,	mass mixing ratio in g/kg
	'''
	es = e_s(t)
##	ws = es * Rd / Rv / (p-es)
	ws = .622 * es / p
	return rh * ws / 100 * 1000	#_convert to g/kg

	
def e_s(t):
	'''
	calculate the approximate saturation vapor pressure at T
	T	float,	temperature, K
	'''
	from numpy import exp
	A = 2.53e11	#_Pa
	B = 5420. 	#_K
	return A*exp(-B/t)


def p2z(p):
	''' estimated height of pressure level in hPa, returned in m '''
	#return 1013.25*(1-2.25577e-7*z)**5.25588
	try:
		return ((p*100/101325.)**(1./5.25588) - 1) / -2.25577e-5
	except:
		return -9999.


def Z2z(Z):
	''' convert geopotential height to altitude '''
	return Re / (Re/Z-1)
##	return Z*9.80/G


def ppmv2kg(ppmv_gas, solute='ozone', **kwargs):
    '''
    convert mass mixing ratio (kg/kg) to parts per million, volume
    ppmV = 1 V_solute / 1e6 V_solution

    r   float,  part per million, volume

    1e-6 * ppmv_gas * (Rair/Rgas)
	1e-3 * ppmv(MOL) * (MW(mol) / MW(Dry Air))
    returns
    ppmV
    '''
    mol = molecules().__getattribute__(solute)['molecular_mass']
    return 1e-6 * ppmv_gas * (Rd / R * mol)


def kg2ppmv(m,p,t,solute='ozone',**kwargs):
	'''
	convert mass mixing ratio (kg/kg) to parts per million, volume
	ppmV = 1 V_solute / 1e6 V_solution
	
	m	float,	mass mixing ratio	[kg/kg]
	p	float,	air pressure		[Pa]
	t	float,	air temperature		[K]
	
	returns
	ppmV
	'''
	M	= molecules()
	n	= m / M.__getattribute__(solute)['molecular_mass']	#_[kg]
	vs	= n*R*t/p
	vd	= Rd*t/p
	return vs/vd


###_START USING CLASSES MORE, GOOBER.  everything should have a plot option
### with an axes keyword/default to whatever
'''this section cribbed from geopy module'''
class Distance(object):
    """
    Base for :class:`.great_circle` and :class:`.vincenty`.
    """

    def __init__(self, *args, **kwargs):
        kilometers = kwargs.pop('kilometers', 0)
        if len(args) == 1:
            # if we only get one argument we assume
            # it's a known distance instead of
            # calculating it first
            kilometers += args[0]
        elif len(args) > 1:
            for a, b in zip(args): #util.pairwise(args):
                kilometers += self.measure(a, b)

        kilometers += units.kilometers(**kwargs)
        self.__kilometers = kilometers

    def __add__(self, other):
        if isinstance(other, Distance):
            return self.__class__(self.kilometers + other.kilometers)
        else:
            raise TypeError(
                "Distance instance must be added with Distance instance."
            )

    def __neg__(self):
        return self.__class__(-self.kilometers)

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        return self.__class__(self.kilometers * other)

    def __div__(self, other):
        if isinstance(other, Distance):
            return self.kilometers / other.kilometers
        else:
            return self.__class__(self.kilometers / other)

    __truediv__ = __div__

    def __abs__(self):
        return self.__class__(abs(self.kilometers))

    def __nonzero__(self):
        return bool(self.kilometers)

    __bool__ = __nonzero__

    def measure(self, a, b):
        raise NotImplementedError()

    def __repr__(self): # pragma: no cover
        return 'Distance(%s)' % self.kilometers

    def __str__(self): # pragma: no cover
        return '%s km' % self.__kilometers

    def __cmp__(self, other):
        if isinstance(other, Distance):
            return cmp(self.kilometers, other.kilometers)
        else:
            return cmp(self.kilometers, other)

    @property
    def kilometers(self): # pylint: disable=C0111
        return self.__kilometers

    @property
    def km(self): # pylint: disable=C0111
        return self.kilometers

    @property
    def meters(self): # pylint: disable=C0111
        return units.meters(kilometers=self.kilometers)

    @property
    def m(self): # pylint: disable=C0111
        return self.meters

    @property
    def miles(self): # pylint: disable=C0111
        return units.miles(kilometers=self.kilometers)

    @property
    def mi(self): # pylint: disable=C0111
        return self.miles

    @property
    def feet(self): # pylint: disable=C0111
        return units.feet(kilometers=self.kilometers)

    @property
    def ft(self): # pylint: disable=C0111
        return self.feet

    @property
    def nautical(self): # pylint: disable=C0111
        return units.nautical(kilometers=self.kilometers)

    @property
    def nm(self): # pylint: disable=C0111
        return self.nautical


class Point(object):
    def __new__(cls, latitude=None, longitude=None, altitude=None, **kwargs):
        latitude = float(latitude or 0.0)
        if abs(latitude) > 90:
            latitude = ((latitude + 90) % 180) - 90
	
        longitude = float(longitude or 0.0)
        if abs(longitude) > 180:
            longitude = ((longitude + 180) % 360) - 180

        altitude = float(altitude or 0.0)

        self = super(Point, cls).__new__(cls)
        self.latitude = latitude
        self.altitude = altitude
        self.longitude = longitude
        self._items = [self.latitude, self.longitude, self.altitude]
        return self 

    def __getitem__(self, index):
        return self._items[index]

    def _iter__(self):
        return iter((self.latitude, self.longitude, self.altitude))

    def __repr__(self):
        return "Point(%r, %r, %r)" % tuple(self._items)


def great_circle(a, b, **kwargs):
    """
    Use spherical geometry to calculate the surface distance between two
    geodesic points. This formula can be written many different ways,
    including just the use of the spherical law of cosines or the haversine
    formula.

    Set which radius of the earth to use by specifying a 'radius' keyword
    argument. It must be in kilometers. The default is to use the module
    constant `EARTH_RADIUS`, which uses the average great-circle radius.

    Example::

        >>> from geopy.distance import great_circle
        >>> newport_ri = (41.49008, -71.312796)
        >>> cleveland_oh = (41.499498, -81.695391)
        >>> great_circle(newport_ri, cleveland_oh).miles
        537.1485284062816

    """
    from math import radians, sin, cos, tan, pi, sqrt, atan2, asin

    a, b = Point(*a), Point(*b)
    lat1, lng1 = radians(a.latitude), radians(a.longitude)
    lat2, lng2 = radians(b.latitude), radians(b.longitude)

    sin_lat1, cos_lat1 = sin(lat1), cos(lat1)
    sin_lat2, cos_lat2 = sin(lat2), cos(lat2)

    delta_lng = lng2 - lng1
    cos_delta_lng, sin_delta_lng = cos(delta_lng), sin(delta_lng)

    d = atan2(sqrt((cos_lat2 * sin_delta_lng) ** 2 +
                   (cos_lat1 * sin_lat2 -
                    sin_lat1 * cos_lat2 * cos_delta_lng) ** 2),
              sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_delta_lng)

    return Re * d


################################################################################
#_CELESTIAL_BODIES_#############################################################
################################################################################

''' change to include option for local time '''
def sunset(elevation=0, latitude=43.0667, longitude=-89.4, **kwargs):
	''' calculate the next sunset occurance in epoch UTC time '''
	import ephem
	from calendar import timegm

	#_initialize a body object to set and calc current position
	sun = ephem.Sun()
	sun.compute()

	#_set an observer location
	obs = ephem.Observer()
	obs.lat = str(latitude) #_takes radians, strings converted from deg
	obs.lon = str(longitude)
	obs.elev = elevation

	utc_tup = obs.next_setting(sun).tuple()
	return timegm(utc_tup)
