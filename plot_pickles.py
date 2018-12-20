#!/usr/bin/env python
# author	Walter Sessions
# purpose	plot up rgba values of data counts
# date		2015.03.30


def main():
	from mpl_toolkits.basemap import Basemap
	import matplotlib.pyplot as plt
	from numpy import array, zeros

	#_read in data and map coords
	d_warm, d_cool, d_invs = calc_fract()
	x2, y2 = read_coords()

	#_get data shape
	n0, n1 = d_warm.shape
	
	#_initialize map
	m = Basemap(projection='npstere',boundinglat=66,lon_0=0,resolution='l')
	m.drawcoastlines(linewidth=1.0)

	#_put into a tuple with each fraction representing rgb values 0-1
	rgb = array(zip(d_warm.flatten(), d_invs.flatten(), d_cool.flatten()))
	rgb = rgb.reshape(n0*n1,3)
	
	mesh = m.pcolor(x2, y2, zeros((n0,n1))) 
##	mesh = m.imshow(rgb, interpolation='none')
	mesh.set_alpha(1.0)
	mesh.set_facecolors(rgb)

	plt.savefig('test.png')
	plt.show()


def read_coords(region='arctic', **kwargs):
	from pickle import load

	#_read in map coordinates
	x = load(open("{0:s}_x2.p".format(region), 'rb'))
	y = load(open("{0:s}_y2.p".format(region), 'rb'))

	return x, y


def calc_fract():
	'''
	read in the three files and calculate fraction
	of total occurances
	'''
	from pickle import load

	fmt = 'lapse_arcticaug_n_sfc_bin1{0:s}.p'

	#_read in data
	d_warm = load(open(fmt.format('_greater'), 'rb'))
	d_cool = load(open(fmt.format('_less'), 'rb'))
	d_invs = load(open(fmt.format(''), 'rb'))

	#_read in total
	d_totl = load(open('lapse_arcticaug_flux_valid_count_cloudy.p', 'rb'))

##	#_sanity check
##	if not (d_warm + d_cool + d_invs == d_totl).all():
##		raise RuntimeError, 'count not the same as reported total'

	return d_warm, d_cool, d_invs 


if __name__ == '__main__':
	main()
