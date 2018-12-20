#!/usr/bin/env python
################################################################################
# script	run_flexpart.py
# author	walter raymond sessions
# date		2015 april 1
# purpose	wrapper script to run flexpart wrf
################################################################################


from os import path, environ


def main(**kwargs):
	pass


def write_available():
	pass


def dtg2str(dtg):
	''' ya know '''
	return '{0:s} {1:s}'.format(dtg[:6], dtg[6:])


def write_input(
	#_PATHNAMES options
	outdir_fmt=path.join(environ['PRODUCTS'],'flexwrf','{0:s}.{1:s}.{2:d}'),
	indir_wrf=path.join(environ['PRODUCTS'], 'wrf'),
	available=path.join(environ['PRODUCTS'], 'wrf', 'AVAILABLE'),

	#_COMMAND
	forward=True,	#_direction of simulation, bool 
	dtg0='20150328000000',	#_start time
	dtg1='20150329000000',	#_end time
	interval=3600,	#_output interval, s
	time_avg=3600,	#
	sampling=180,	#_sampling interval, s
	split=99999999,	#_time between spliting particles, s
	sync=180,		#_synchronisation interval, s
	ctl=10.,		#_factor time step must be smaller than tl
	ifine=10,		#_decrease of time step for vertical motion
	iout=1,			#_{conc, mixrat, both, plum traj, 1+4}, 1-5
	ipout=0,		#_particle dump, 0=no, 1=@output int, 2=@end
	lsubgrid=0,		#_subgrid terrain effect paramiterization
	lconvection=0,	#_convection, 3=yes, 0=no
	dt_conv=3600.,	#_time interval to call convection, seconds
	lagespectra=1,	#_age spectra: 1 yes, 0 no
	ipin=0,			#_restart w/ dumped particle data: 1 yes, 0 no
	iflux=0,		#_calculate fluxes
	ioutputforeachrel=0,	#_create output for each location
	mdomainfill=0,	#_domain-filling traj, 0=n,1=y,2=strat_o3 
	ind_source=1,	#_1=mass unit , 2=mass mixing ratio unit
	ind_receptor=2,	#_1=mass unit , 2=mass mixing ratio unit
	nested_output=0,#_nested output, 0 no, 1 yes
	linit_cond=0,	#_ICOND4BW RUNS: 0=NO,1=MASSUNIT,2=MASSMIXRAT
	turb_option=1,	#_0=no; 1=diagnosed as flexpart_ecmwf; 2 and 3=tke.
	cbl=0,			#_0=no, 1=yes. works if TURB_OPTION=1
	sfc_option=1,	#_0=def comp of u*, hflux, pblh, 1=from wrf
	wind_option=1,	#_0=snapshot winds,1=mean winds,
					# 2=snapshot eta-dot,-1=w based on divergence
    time_option=1,	#_1=corr of time validity for time-average wind
					# 0=no need
	outgrid_coord=1,#_0=wrf grid(meters), 1=regular lat/lon grid
	release_coord=1,#_0=wrf grid(meters), 1=regular lat/lon grid
	iouttype=2,		#_0=default binary, 
					# 1=ascii (for particle dump only),
					# 2=netcdf
	ntimerec=3,		#_time frames per output file, only used for netcdf
	verbose=0,		#_VERBOSE MODE,0=minimum, 100=maximum

	#_OUTGRID opts
	llon=-120.,		#_lower left lon
	llat=40.,		#_lower left lat
	nx=24,			#_grid points in x
	ny=28,			#_grid points in y
	outgriddef=0,	#_0=grid dist, 1=upperright corner coord
	dxlon=0.125,	#_delta x
	dylon=0.125,	#_delta x
	zlevs=[100.,500.,1000.,2000.,5000.,7000.,9000.,12000.,20000.],
	nllon=-120.,	#_lower left lon
	nllat=40.,		#_lower left lat
	nnx=24,			#_grid points in x
	nny=28,			#_grid points in y
	noutgriddef=0,	#_0=grid dist, 1=upperright corner coord
	ndxlon=0.125,	#_delta x
	ndylon=0.125,	#_delta x

	#_RECEPTOR opts (used??)
	nreceptor=0,

	#_EVENTUALLY WORK OUT SPECIES
	#_RELEASES opts
	npec=1,			#_total number of species emitted
	emitvar=0,		#_1=emit variation
	link0=1,		#_index of species in file SPECIES
	link1=2,		#_index of species in file SPECIES
	releases={
		'test' : { 
			'llon'	: -120., 'ulon' : -119., 'llat' : 40., 'ulat' : 41., 
			'z_top' : 150., 'z_bot' : 100.,	#_vert
			'npart' : 5e5,					#_number of particles
			'spec_mass' : [1e3,],
			},
		},
	dtg0r=None,		#_beginning time of release
	dtg1r=None,		#_ending time of release
	zkind=1,		#_1 m agl, 2 m asl, 3 pressure
	**kwargs):
	'''	generate input file for flexpart run '''
	import os

	pid = os.getpid()

	#_open file for writing
	fname = 'flexpart.input.{0:07d}'.format(pid)
	f = open(fname, 'w')
	w = f.write

	#_write PATHNAMES section 
	w('=====================FORMER PATHNAMES FILE===================\n')
	w(outdir_fmt.format(dtg0, dtg1, pid) + '\n')
	w(indir_wrf + '\n')
	w(available + '\n')
	w('=============================================================\n')

	#_write COMMAND section
	w('=====================FORMER COMMAND FILE=====================\n')
	w('{0:d}\n'.format(1 if forward else -1))
	w('{0:s}\n'.format(dtg2str(dtg0)))
	w('{0:s}\n'.format(dtg2str(dtg1)))	
	w('{0:d}\n'.format(interval))	
	w('{0:d}\n'.format(time_avg))
	w('{0:d}\n'.format(sampling))		
	w('{0:d}\n'.format(split))				
	w('{0:d}\n'.format(sync))					
	w('{0:f}\n'.format(ctl))
	w('{0:d}\n'.format(ifine))
	w('{0:d}\n'.format(iout))
	w('{0:d}\n'.format(ipout))
	w('{0:d}\n'.format(lsubgrid))
	w('{0:d}\n'.format(lconvection))
	w('{0:f}\n'.format(dt_conv))
	w('{0:d}\n'.format(lagespectra))
	w('{0:d}\n'.format(ipin))
	w('{0:d}\n'.format(iflux))
	w('{0:d}\n'.format(ioutputforeachrel))
	w('{0:d}\n'.format(mdomainfill))
	w('{0:d}\n'.format(ind_source))
	w('{0:d}\n'.format(ind_receptor))
	w('{0:d}\n'.format(nested_output))
	w('{0:d}\n'.format(linit_cond))
	w('{0:d}\n'.format(turb_option))
	w('{0:d}\n'.format(cbl))
	w('{0:d}\n'.format(sfc_option))
	w('{0:d}\n'.format(wind_option))
	w('{0:d}\n'.format(time_option))
	w('{0:d}\n'.format(outgrid_coord))
	w('{0:d}\n'.format(release_coord))
	w('{0:d}\n'.format(iouttype))
	w('{0:d}\n'.format(ntimerec))
	w('{0:d}\n'.format(verbose))

	#_write AGESCLASSES section
	w('=====================FORMER AGECLASESS FILE==================\n')
	w('{0:d}\n'.format(2))
	w('{0:d}\n'.format(7200))
	w('{0:d}\n'.format(999999))

	#_write OUTGRID section
	w('=====================FORMER OUTGRID FILE=====================\n')
	w('{0:7.2f}\n'.format(llon))
	w('{0:7.2f}\n'.format(llat))
	w('{0:d}\n'.format(nx))
	w('{0:d}\n'.format(ny))
	w('{0:d}\n'.format(outgriddef))
	w('{0:8.4f}\n'.format(dxlon))
	w('{0:8.4f}\n'.format(dylon))
	w('{0:d}\n'.format(len(zlevs)))
	[w('{0:f}\n'.format(z)) for z in zlevs]

	#_write OUTGRID_NEST	
	w('================OUTGRID_NEST==========================\n')
	w('{0:7.2f}\n'.format(nllon))
	w('{0:7.2f}\n'.format(nllat))
	w('{0:d}\n'.format(nnx))
	w('{0:d}\n'.format(nny))
	w('{0:d}\n'.format(noutgriddef))
	w('{0:8.4f}\n'.format(ndxlon))
	w('{0:8.4f}\n'.format(ndylon))

	#_write RECEPTOR section (backward) 
	w('=====================FORMER RECEPTOR FILE====================\n')
	w('{0:d}\n'.format(nreceptor))

	print 'MAKE A SPECIES DICTIONARY'
	w('=====================FORMER SPECIES FILE=====================\n')
	line = \
    '''2	NUMTABLE        number of variable properties. The following lines are fixed format
XXXX|NAME    |decaytime |wetscava  |wetsb|drydif|dryhenry|drya|partrho  |parmean|partsig|dryvelo|weight |
    AIRTRACER     -999.9   -9.9E-09         -9.9                 -9.9E09                   -9.99   29.00
    Cs-137        -999.9    1.0E-04  0.80   -9.9                  2.5E03  6.0E-7  3.0E-1   -9.99   -9.99\n'''
	w(line)

	#_write RELEASES section
	w('=====================FORMER RELEEASES FILE===================\n')
##	w('{0:d}\n'.format(nspec))
	w('{0:d}\n'.format(2))
	w('{0:d}\n'.format(emitvar))
	w('{0:d}\n'.format(link0))
	w('{0:d}\n'.format(link1))
	w('{0:d}\n'.format(len(releases)))
	w('{0:s}\n'.format(dtg2str(dtg0 if dtg0r is None else dtg0r)))
	w('{0:s}\n'.format(dtg2str(dtg1 if dtg1r is None else dtg1r)))
	for name, opts in releases.iteritems():
		w('{0:9.4f}\n'.format(opts['llon']))	
		w('{0:9.4f}\n'.format(opts['llat']))
		w('{0:9.4f}\n'.format(opts['ulon']))	
		w('{0:9.4f}\n'.format(opts['ulat']))
		w('{0:d}\n'.format(zkind))
		w('{0:10.3f}\n'.format(opts['z_bot']))
		w('{0:10.3f}\n'.format(opts['z_top']))
		[w('{0:10e}\n'.format(sm)) for sm in opts['spec_mass']]
		w('{0:s}\n'.format(name))

	#_close file
	f.close()
	print fname
	return fname


if __name__ == '__main__':
	#_internal or external namelist?
	main(**namelist)
