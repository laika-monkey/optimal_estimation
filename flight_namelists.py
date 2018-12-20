#!/usr/bin/env python
#############################################################################_80
# file:		hs3_namelists.py
# author:	Walter R. Sessions, 2015.06.01
# purpose:	Central file for dictionaries containing segment information
#			during the HS3 campaign being used for my dust studies.
#			Intended for use initially within optimal_estimation.py	
################################################################################
################################################################################


experiments = {
	'HS?' : { #_tester
		'20130820200000' :	{	#_just make dtg somewhere in range
			'desc'	: 'Edge ELE 0.2 aerosol optical dept event over Atlantic',
			'fidx'	: range(5),
			},
		},

	#_What the hell are these? Seems to be some transitional cases?
	'HS?' : { #_plotlabel SEPT_2015 
		'20130820200000' :	{	#_just make dtg somewhere in range
			'desc'	: 'Edge ELE 0.2 aerosol optical dept event over Atlantic',
			'fidx'	: range(1370),
			},
		'20130821000000' : {
			'desc'	: 'ELE 0.2 AOD event north of SA. CTT below dust layer.',
			'fidx'	: range(100),
			},
		'20130821020000' : {
			'desc'	: 'ELE 0.2 AOD event north of SA. CTT below dust layer.',
			'fidx'	: range(1450),
			},
		'20130821060000' : {
			'desc'	: 'PBL aerosol activity, otherwise clear',
			'fidx'	: range(435),
			},
		'20130824170000' : {
			'desc'	: 'PBL aerosol activity transition',
			'fidx'	: range(570),
			},
		'20130824190000' : {
			'desc'	: 'ELE aerosol event west of CAPOVE. 0.2-3 AOD.',
			'fidx'	: range(350),
			},
		'20130824230000' : {
			'desc'	: 'LELE/PBL, AOD 0.5',
			'fidx' : range(928),
			},
		'20130825060000' : {
			'desc'	: 'ELE/PBL 0.3 AOD event.',
			'fidx'	: range(560),
			},
		'20130830060000' : {
			'desc'	: 'Elevated layer',
			'fidx'	: range(500),
			},
		},

	'HS3_USE-FOR-CLOUD-ATP-TEST-2016.02'	: {
		'20130821020000' : {
			'desc'	: 'Elevated layer',
			},

		'20130824190000' : {
			'desc'	: 'PBL aerosol activity transition',
			},

		'20130830060000' : {
			'desc'	: 'PBL aerosol activity transition',
			},	},

	#_2016.02 USE THIS ONE
	'HS3_short'	: {
		'20130824170000' : {
			'desc'	: 'PBL aerosol activity transition',
		#	'fidx'	: [], 
			},},

		
	'HS3'	: {
		'20130824170000' : {
			'desc'	: 'PBL aerosol activity transition',
		#	'fidx'	: [], 
			},

		'2013082106000' : {
			'desc'	: 'Long elevated layer',
		#	'fidx'	: [228],
			},
		'20130820200000' : {
			'desc'	: 'PBL aerosol activity transition',
		#	'fidx'	: [],
			},
		'20130821020000' : {
			'desc'	: 'Elevated layer',
		#	'fidx'	: [],
			},
		'20130824190000' : {
			'desc'	: 'PBL aerosol activity transition',
		#	'fidx'	: [], 
			},
		'20130824220000' : {
			'desc'	: 'PBL aerosol activity transition',
		#	'fidx'	: [], 
			},
		'20130829150000' : {
			'desc'	: 'PBL aerosol activity transition',
		#	'fidx'	: [], 
			},
		'20130829190000' : {
			'desc'	: 'PBL aerosol activity transition',
		#	'fidx'	: [], 
			},
		'20130830060000' : {
			'desc'	: 'PBL aerosol activity transition',
		#	'fidx'	: [], 
			},	},

	#_segments of interest for TC4
#	'''
#	Taking notes while skimming through cpl site, pre-colloc
#	09JUL2007 has layers off Baja CA
#	19JUL2007 SO MUCH DUST OFF COSTA RICA
#	25JUL2007 more dust, but lots of high clouds. Good test of separation
#	29JUL2007 lower clouds and this dust layers
###	09AUG2007 dust in Gulfo de Mexico before landfall
#	
#	'''
	'TC4' : {
		'20070719160000'	: {
			'desc'	: 'dust after land transit, low, GOM',
			'fidx'	: range(440, 661),
			},
		'20070725160000'	: {
			'desc'	: 'dust off Central Am.',
			'fidx'	: range(330, 441),
		##	'fidx'	: range(810, 1091),
			},
		},


	#_segments of interest for THORPEX over Pacific theater
#	'''
#	22FEB03	Dust and other assundries, low, clouds
#	24FEB03 " " " " "
#	01MAR03 JUST DO ALL OF THESE AND FIGURE OUT IF SAME CRAP WORKS
#	Probably seasalt?  Layer low lwow lowowo
#	03MAR03
#	11MAR03	High layer of aerosols
#	12MAR03	lower layers again
#	13MAR03	Lower layers, goes over the Big Island, less low clouds	
#	'''
	'PTOST' : {
		'DTG'	: {
			'desc'	: 'of segment',
			'fidx'	: range(0),
			},
		},	


	}	#_end dict.segments_####################################################
