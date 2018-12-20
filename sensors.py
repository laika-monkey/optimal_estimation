'''
specifications for LBL-RTM and LBL-DIS computations by sensor

'''

airs = {
	#_only the window we care about
	'v1'	: 649.,
	'v2'	: 1614., 
	'dv'	: 0.1,		#_not great around 10um
}

full = {
	'v1' : 1000.,
	'v2' : 2800.,
	'dv' : 0.1,
}

ir_window = {
	'v1' : 600.,
	'v2' : 1430.,
	'dv' : 0.1,
}

nir_to_co2 = {
	'v1' : 600.,
	'v2' : 2500.,
	'dv' : 0.1,
}

shis = {
	#_get to caring about solar angles
}
