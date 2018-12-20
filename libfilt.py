# collection of dumb bad stupid filters

def onethreeone(a, win=5):
	'''
	scipy.signal should have something like this in it that works right.
	I can't be bothered to look today
	'''
	from numpy import max, mean, array, tile
	a = array(a)
	p = win-2
	m = win-3
	tmp = tile(a, (3,1)).T.flatten()
	tmp = array([mean(tmp[max([i-m, 0]):i+p]) for i in range(a.size*3)[1::p]])

	return tmp
