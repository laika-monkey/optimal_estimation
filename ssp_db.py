import os
import numpy as np

# one-liner helper function 
def _read_array_line(f, expected_dtype=np.float, expected_sep=' '):
    return np.fromstring(f.readline(), 
                            dtype = expected_dtype, 
                            sep = expected_sep)

def get_ssp_dtype(num_phase_angles):
    """
    get_ssp_dtype, given a number of phase angles.
    the dtype is a compound record type for numpy.ndarray.
    """
    return np.dtype( [ ('s_ext', np.float),
                       ('s_sca', np.float),
                       ('s_abs', np.float),
                       ('w0', np.float),
                       ('g', np.float),
                       ('Q_ext', np.float),
                       ('Q_abs', np.float),
                       ('Q_sca', np.float),
                       ('vol', np.float),
                       ('proj_area', np.float),
                       ('phase_func', (np.float, num_phase_angles)) ] )
   
def write_ssp_db(filename, data, comments):
    """
    write_ssp_db(filename, data, comments)

    Write a single scattering properties database to disk, 
    using the text format.
    Inputs:
    filename: string containing filename of new database file.
       raises IOError if this file already exists.
    data: database data structure. Should be a dictionary with 
       the following keys (or, equivalently, an ssp_db object):
       wavenum: list or array containing N wavenumbers [1/cm]
       wavelen: list or array containing N wavelengths [um]
       r_eff: list or array containing M effective radii [um]
       angles: list or array containing P phase function angles [deg]
       ssp: dictionary or ndarray (with compound data type), 
           containing the following database fields, each with 
           [N,M] elements:
           Q_ext: extinction efficiency
           Q_sca: scattering efficiency
           Q_abs: absorption efficiency
           s_ext: extinction cross section per particle [um^2]
           s_sca: scattering cross section per particle [um^2]
           s_abs: absorption cross section per particle [um^2]
           proj_area: projected area per particle [um^2]
           vol: volume per particle [um^3]
           g: asymmetry parameter
           w0: single scattering albedo
           phase_func: array of P phase function values
       
    """
    
    f = open(filename, 'w')

    f.write(comments[0] + '\n')
    f.write(comments[1] + '\n')

    nlines = data['wavenum'].shape[0] * data['r_eff'].shape[0]
    nangles = data['angles'].shape[0]

    f.write('{0:<7d} Number of data lines\n'.format(nlines))
    f.write('{0:<7d} Number of angles in the phase function\n'.format(nangles))

    f.write('Phase function angles:\n')

    f.write(" ".join( ['{0:12.5e}'.format(p) for p in data['angles']] ))
    f.write('\n')

    hdrs = ['wave [um]', 'wave [cm-1]', 'reff [um]', 'ext [um^2]', 
            'scat [um^2]', 'abs [um^2]', 'w0', 'g', 'Qext', 'Qabs', 
            'Qsca', 'Vol [um^3]', 'Proj_area [um^2]', 
            '...and phase function values...']
    f.write(" ".join( ['{0:>12s}'.format(s) for s in hdrs] ))
    f.write('\n')

    ordered_flds_set1 = ['s_ext', 's_sca', 's_abs']
    ordered_flds_set2 = ['w0', 'g', 'Q_ext', 'Q_abs', 'Q_sca']
    ordered_flds_set3 = ['vol', 'proj_area']

    # to account for this as a compound dtype, or python dictionary,
    # we use the shape of the g ssp variable; Note this assumes if the
    # variables are in a python dictionary that all are the same shape.
    for n in np.ndindex(data['ssp']['g'].shape):
        f.write('{0:12.4f}'.format(data['wavelen'][n[0]]) + " ")
        f.write('{0:12.4f}'.format(data['wavenum'][n[0]]) + " ")
        f.write('{0:12.4f}'.format(data['r_eff'][n[1]]) + " ")
        f.write(" ".join( ['{0:12.4f}'.format(data['ssp'][fld][n]) 
                           for fld in ordered_flds_set1] ))
        f.write(" ".join( ['{0:12.6f}'.format(data['ssp'][fld][n]) 
                           for fld in ordered_flds_set2] ))
        f.write(" ".join( ['{0:13.6e}'.format(data['ssp'][fld][n]) 
                           for fld in ordered_flds_set3] ))
        f.write(" ".join( ['{0:12.5e}'.format(p) 
                           for p in data['ssp']['phase_func'][n]] ))
        f.write('\n')

    f.close()

def read_ssp_db(filename):
    """
    Read an ssp_db file, and populate a python dictionary with the 
    various datasets. Each dataset is a numpy ndarray. The main dataset 
    (the single scattering properties) is a compound dtype.
    
    Contents: see write_ssp_db.
    """
    f = open(filename, 'r')

    dat = {}
    # First 5 lines are comments. stripping endline characters;
    # not sure if that's the best thing to do
    dat['comments'] = []
    for n in range(2):
        dat['comments'].append(f.readline()[0:-1])

    # get the array layout from the comments - not sure if that 
    # is the standard ICD or not, either.
    num_lines = int(f.readline().split()[0])
    num_angles = int(f.readline().split()[0])
    # dummy comment 'phase function angles:' - dont need to save it
    f.readline()

    nums_per_row = num_angles + 13

    # phase function angles from header
    if num_angles > 0:
        dat['angles'] = _read_array_line(f)
    else:
        dat['angles'] = np.zeros(0)
        # this should be a blank line, if there are no angles
        f.readline()

    # last header line - not saving this either
    f.readline()

    raw_arr = np.empty((nums_per_row, num_lines))
    for n in range(num_lines):
        raw_arr[:,n] = _read_array_line(f)

    f.close()

    # attempt to figure out the 2D
    # do a crude check first
    num_wave = sum(raw_arr[2,:] == raw_arr[2,0])
    num_r_eff = sum(raw_arr[0,:] == raw_arr[0,0])
    if num_lines != (num_r_eff * num_wave):
        raise ValueError, 'Database does not appear to have [R,W] elements'

    # Find unique r_eff and wave, by subindexing each by a 
    # fixed value of the other
    sorter = raw_arr[2,:].argsort()
    dat['wavenum'] = raw_arr[1,sorter[0:num_wave]]
    dat['wavenum'].sort()
    dat['wavelen'] = raw_arr[0,sorter[0:num_wave]]
    dat['wavelen'].sort()
    dat['wavelen'][:] = dat['wavelen'][::-1]
    sorter = raw_arr[0,:].argsort()
    dat['r_eff'] = raw_arr[2,sorter[0:num_r_eff]]
    dat['r_eff'].sort()

    new_dtype = get_ssp_dtype(num_angles)
    
    raw_arr_trim = raw_arr[3:,:]
    # this reshape assumes the file is organized with the wavenum
    # fixed in blocks while the r_eff cycles through all possible 
    # values: line1 (wnum_1, re_1), line2 (wnum_1, re_2) ...
    dat['ssp'] = \
        raw_arr_trim.T.ravel().view(new_dtype).reshape((num_wave, num_r_eff))

    return dat


def _reform_ssp_db(d):

    # reforms the data into a set of 2-D arrays indexed 
    # by R_eff and wavenumber - This is assuming the ssp_db has 
    # [N_r x N_w] lines, which may not be true (?)

    # This should not be needed if the reshape works (see last line 
    # of the read function, above.); just keeping this around in case 
    # the method is needed in the future, but at the moment this is 
    # not used.
    
    # do a crude check first
    num_wave = sum(d['r_eff'] == d['r_eff'][0])
    num_r_eff = sum(d['wavelen'] == d['wavelen'][0])
    if d['num_lines'] != (num_r_eff * num_wave):
        raise ValueError, 'Database does not appear to have [R,W] elements'

    # Find unique r_eff and wave, by subindexing each by a 
    # fixed value of the other
    sorter = d['r_eff'].argsort()
    wavenum = d['wavenum'][sorter][0:num_wave]
    wavenum.sort()
    wavelen = d['wavelen'][sorter][0:num_wave]
    wavelen.sort()
    sorter = d['wavenum'].argsort()
    r_eff = d['r_eff'][sorter][0:num_r_eff]
    r_eff.sort()
    
    new_dtype = get_ssp_dtype(d['num_angles'])
    arr = np.zeros( (num_wave, num_r_eff), dtype=new_dtype )

    twod_fields = d.keys()
    twod_fields.remove('wavenum')
    twod_fields.remove('wavelen')
    twod_fields.remove('r_eff')
    twod_fields.remove('num_angles')
    twod_fields.remove('num_lines')
    twod_fields.remove('comments')
    twod_fields.remove('angles')

    for w in range(num_wave):
        for r in range(num_r_eff):
            try:
                k = np.argwhere( 
                    np.logical_and(d['r_eff'] == r_eff[r], 
                                      d['wavenum'] == wavenum[w]) )
            except:
                import pdb; pdb.set_trace()
                print 'igh'
            if len(k) != 1:
                import pdb; pdb.set_trace()
                raise ValueError, 'array indexing got lost'
            for fld in twod_fields:
                try:
                    arr[fld][w,r] = d[fld][:,k].flatten()
                except (IndexError, TypeError, ValueError), theError:
                    import pdb; pdb.set_trace()
                    print 'you suck'
                    print str(theError)

    return {'wavenum':wavenum, 'wavelen':wavelen, 'r_eff':r_eff, 
            'angles':d['angles'].copy(), 'dat':arr}

class ssp_db(dict):
    """
    ssp_db class

    This is a subclass of the python builtin dictionary.
    Mainly this is currently used just to have a nicer __str__ for 
    the ssp_db in order to print it nicely.

    Use:

    db_obj = ssp_db(filename='ssp_db.mie_water')

    This will load the ssp_db text file.

    optional keywords:
    data=None : create a new ssp_db object from an already loaded ssp_db 
        that is contained in a python dictionary or compound ndarray.
        This would be used instead of a filename load.
    copydata=True : If the data keyword is used, this controls whether 
        the input data is copied.
    comments=None : set to some strings for some optional comments 
        (can be a 2-element list/tuple)
    """
    def __init__(self, filename=None, data=None, copydata=True, comments=None):

        if filename is not None:
            d = read_ssp_db(filename)
            for k,v in d.iteritems():
                self[k] = v
            self['filename'] = filename
        elif data is not None:
            # predefine list of fields we need to copy (if we just copied 
            # all, we might get other things we don't want)
            klist = ['wavenum', 'wavelen', 'r_eff', 'angles']
            ssp_klist = ['Q_abs', 'Q_sca', 'Q_ext', 'w0', 'g', 
                        's_abs', 's_sca', 's_ext', 'phase_func', 
                        'proj_area', 'vol']
            # use np.array, for the copy keyword (can just send the copydata
            # flag directly), and to make sure self will contain ndarrays.
            # this would imply self will contain a copy, even if copydata 
            # is false, if the data input had python lists anywhere.
            for k in klist:
                self[k] = np.array(data[k], copy=copydata)

            # copy ssp as either dict or ndarray (having compound dtype)
            if isinstance(data['ssp'],dict):
                self['ssp'] = {}
                for k in ssp_klist:
                    self['ssp'][k] = np.array( data['ssp'][k], 
                                                  copy=copydata )
            else:
                self['ssp'] = np.array( data['ssp'], copy=copydata )

            self['comments'] = ['','']
            if comments is not None:
                self['comments'][0] = comments[0]
                if len(comments) > 1:
                    self['comments'][1] = comments[1]

        else:
            raise TypeError, 'ssp_db class requires a filename or ' + \
                'data structure as input'


    def __str__(self):
        output_line = \
            'ssp database:\n' + \
            'number of wnum:  {3:6d} ({4:f} - {5:f} [1/cm])\n' + \
            'number of radii: {0:6d} ({1:f} - {2:f} [um])\n' + \
            'number of angles: {6:6d}\n' + \
            'comment 1: {7:s}\n' + \
            'comment 2: {8:s}'
        return output_line.format(self['r_eff'].shape[0], 
                                  self['r_eff'][0], self['r_eff'][-1], 
                                  self['wavenum'].shape[0], 
                                  self['wavenum'][0], self['wavenum'][-1], 
                                  self['angles'].shape[0], 
                                  self['comments'][0], 
                                  self['comments'][1])
