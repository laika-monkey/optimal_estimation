#############################################################################_80
#_NAME		: libnet.py						       #
#_AUTHOR	: Walter R. Sessions, 8 Mar 2013			       #
#_USAGE		: As library only					       #
#_PURPOSE	: No clue.  Just threw some download functions in here.        #
#############################################################################_80
debug = 3 
################################################################################
#_PROTOCOLS_####################################################################
################################################################################

def http( url ):
	''' downloads through http '''
	import urllib2
	dbg( url )
	file = url.split('/')[-1]
	u = urllib2.urlopen(url)
	f = open(file, 'wb')
	meta = u.info()
	dbg( meta, l=2 )
	f_size = int( meta.getheaders("Content-Length")[0] )
	dbg( "Downloading: %s Bytes: %s" % (file, f_size), l=2 )

	fs_dl = 0
	blk_sz = 8192
	while True:
		buf = u.read(blk_sz)
		if not buf: break

		fs_dl += len(buf)
		f.write(buf)
		status = r"%10d  [%3.2f%%]" % (fs_dl, fs_dl * 100. / f_size)
		status = status + chr(8)*( len(status)+1 )
		print status,

	f.close()

def ftp( url, file, username='anonymous', password='none', dir='' ):
	''' downloads file though ftp '''
	from ftplib import FTP
	dbg( 'ftp://'+url+dir+'/'+file )

	#_Open handle
	ftp_handle = FTP( url, username, password )
###	ftp_handle.login()

	#_Move to directory containing file
	ftp_handle.cwd( dir )
	
	#_Open local binary and stuff it in
	local_file = open( file, 'wb' )
	ftp_handle.retrbinary( 'RETR ' + file, local_file.write, 8*1024 )
	local_file.close()

	ftp_handle.close()

################################################################################
#_MODEL_SPECIFIC_MODULES_#######################################################
################################################################################

def opendap_example( pipe, variable, file ):
	'''
	GEOS5 uses an OPENDAP server that is easier handled on its own

	dtg	: str*10, Expected 00.  GEOS5 is initialized at 22z
	pipe	: pydap.open_url() object, connected to GSFC
	variable: string of nrl variable name
	file	: species output file name
	'''
	import numpy as np
	from netCDF4 import Dataset
        #_Conversion between GSFC variables and local ones
        var_nrl2gsfc = {'dust_aod' 		: 'duexttau',
       	                'sulfate_aod' 		: 'suexttau', 
       	                'blackcarbon_aod' 	: 'bcexttau',
       	                'organiccarbon_aod' 	: 'ocexttau',
       	                'seasalt_aod'		: 'ssexttau'	}

	#_Get diminsional values
	spec = var_nrl2gsfc[ variable + '_aod' ]
	nt, ny, nx = pipe[spec][spec].shape

        #_setup netcdf
	#_Get GSFC name for species and create ncdf variable
        ncdf = Dataset( file, 'w', format='NETCDF3_CLASSIC' )
        dbg( file )
	
        #_Create dimensions
        ncdf.createDimension( 'lon', nx )
        ncdf.createDimension( 'lat', ny )
        ncdf.createDimension( 'time', nt )
	
	#_Add dimensional variables to file
        v = ncdf.createVariable('lats','f4',('lat'))
	v[:] = pipe.lat[:]
     	v = ncdf.createVariable('lons','f4',('lon'))
       	v[:] = pipe.lon[:]
       	v = ncdf.createVariable('times','f8',('time'))
       	v[:] = pipe.time[:]

       	#_Download variables
	v = ncdf.createVariable( variable, 'f4', ('time','lat','lon') )

	#_Download data to file
	v[:] = pipe[spec][spec][:]

	#_Close file and make readable 
	ncdf.close()

def dbg( msg, l=1 ):
        ''' if global debug is set to true, be more verbose '''
        import inspect
        if debug >= l:
                curf = inspect.currentframe()
                calf = inspect.getouterframes( curf, 2 )
                file, method = calf[1][1], calf[1][3]
                print '[%s.%s] %s' % ( file, method, msg )

################################################################################
#_COMMAND_LINE_METHOD_##########################################################
################################################################################

if __name__ == '__main__':
	pass
