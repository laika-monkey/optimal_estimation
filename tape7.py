import numpy as np

class tape7:
    """
    TAPE7 reader for LBLRTM output.

    Usage: t7 = lblrtm_utils.tape7(tape7_filename)
    If the filename is omitted, then the default is to attempt to read 
    TAPE7 in the current directory.

    Note that this assumes a certain TAPE7 format, which should be the 
    format created by the TAPE5 writer inside this module. A general 
    TAPE7 reader would be problematic since I don't think there is 
    any metadata in the file that describes the format.

    Output is a tape7 class, which has the various data items 
    from the TAPE7 as data members, as follows:

    header fields:
      nlayer: number of layers (scalar int)
      nmol: number of molecular species (scalar int)
      desc: string description of the standard atmosphere used
      H1: start of slant path (from TAPE5 input)
      H2: end of slant path
      ang: angle of slant path
    data fields: (all floating point numy arrays):
      pbar: average pressure within layer
      tbar: average temperature within layer
      zlevel1: altitude of layer bottom
      zlevel2: altitude of layer top
      tlevel1: temperature at layer bottom
      tlevel2: temperature at layer top
      plevel1: pressure at layer bottom
      plevel2: pressure at layer top
      q: water vapor mass mixing ratio [kg wv/kg dry air]
      mol_vmr: volume mixing ratios, expressed with respect 
        to dry air. This array is [nmol+1, nlayer], with the 
        last row containing the vmr of all remaining skipped 
        molecular species. So, the sum of mol[1:,1] will be 1, 
        since the first row is always water vapor.
    """

    def __init__(self, filename='TAPE7'):
        self.filename = filename
        self._read_header()
        self._create_arrays()
        self._read_data()

    def _read_header(self):

        f = open(self.filename, 'r')
        hdrline = f.readline()
        self.comment = hdrline[1:-1]
        hdrline = f.readline()
        f.close()

        self.nlayer = int(hdrline[2:5])
        self.nmol = int(hdrline[7:10])
        self.desc = hdrline[20:36]
        self.H1 = float(hdrline[40:48])
        self.H2 = float(hdrline[52:60])
        self.ang = float(hdrline[65:73])

    def _create_arrays(self):
        self.pbar = np.zeros(self.nlayer)
        self.tbar = np.zeros(self.nlayer)
        self.zlevel1 = np.zeros(self.nlayer)
        self.zlevel2 = np.zeros(self.nlayer)
        self.tlevel1 = np.zeros(self.nlayer)
        self.tlevel2 = np.zeros(self.nlayer)
        self.plevel1 = np.zeros(self.nlayer)
        self.plevel2 = np.zeros(self.nlayer)
        self.q = np.zeros(self.nlayer)
        self.mol_vmr = np.zeros((self.nmol+1, self.nlayer))
        self.mol_den = np.zeros((self.nmol+1, self.nlayer))

    def _read_data(self):

        # Since the output is formatted fixed width, 
        # we can pull out the numbers with hard-coded character 
        # numbers instead of doing something scanf-like; this is 
        # a little nicer since python by default does not have a 
        # scanf type function.
        # Note that since I don't have an ICD or file definition 
        # for the TAPE7 , the field boundaries might be incorrect 
        # since I am just eyeballing them.

        # first, need to figure out how many mol lines to expect.
        # the 8th value is always printed as the total number of 
        # molecules for all the skipped species. 8 numbers are 
        # always printed even if nmol < 7, with zero values padding.
        # if nmol > 8, then additional lines are printed (8 
        # numbers per line). For example, if nmol is 18, then 3 lines
        # will be written - 8, 8, 3. (remember the first line has only
        # 7 mol values since position 8 is the remainder.)

        min_nmol = max([self.nmol, 7])
        fract_num_mol_lines = (min_nmol+1)/8.0;
        num_mol_lines = int(np.ceil(fract_num_mol_lines))
        num_mol_per_line = np.zeros(num_mol_lines,dtype='int') + 8
        if np.ceil((min_nmol+1)/8.0) != np.floor((min_nmol+1)/8.0):
            num_mol_per_line[-1] = 8*(fract_num_mol_lines - 
                                      np.floor(fract_num_mol_lines))

        # ok, now start reading the file.
        f = open(self.filename)

        fline = f.readline()
        fline = f.readline()

        fline = f.readline()

        self.pbar[0] = float(fline[0:11])
        self.tbar[0] = float(fline[15:25])
        self.zlevel1[0] = float(fline[41:48])
        self.plevel1[0] = float(fline[48:56])
        self.tlevel1[0] = float(fline[56:63])
        self.zlevel2[0] = float(fline[63:70])
        self.plevel2[0] = float(fline[70:78])
        self.tlevel2[0] = float(fline[78:85])

        moldata = np.zeros( min_nmol + 1 )
        p = 0
        for n in range(num_mol_lines):
            fline = f.readline()
            tmpdata = [float(fline[c*15:c*15+15]) 
                       for c in range(num_mol_per_line[n])]
            moldata[p:p+num_mol_per_line[n]] = tmpdata
            p = p + num_mol_per_line[n]

        # calculation to convert the molecular column density [cm^-2]
        # back into VMR with respect to dry air;
        dry_air_col = sum(moldata[1:])
        if self.nmol <= 7:
            self.mol_den[0:self.nmol,0] = moldata[0:self.nmol]
        else:
            self.mol_den[0:7,0] = moldata[0:7]
            self.mol_den[7:self.nmol,0] = moldata[8:]
        self.mol_den[self.nmol,0] = moldata[7]
        self.mol_vmr[:,0] = self.mol_den[:,0] / dry_air_col

        # Now convert water vapor VMR into mass mixing ratio Q [kg/kg].
        R_d_over_R_v = 287.0 / 461.0
        self.q[0] = R_d_over_R_v * moldata[0] / dry_air_col

        for v in range(1,self.nlayer):

            fline = f.readline()

            self.pbar[v] = float(fline[0:11])
            self.tbar[v] = float(fline[15:25])
            self.zlevel2[v] = float(fline[63:70])
            self.plevel2[v] = float(fline[70:78])
            self.tlevel2[v] = float(fline[78:85])

            p = 0
            for n in range(num_mol_lines):
                fline = f.readline()
                tmpdata = [float(fline[c*15:c*15+15]) 
                    for c in range(num_mol_per_line[n]) ]
                moldata[p:p+num_mol_per_line[n]] = tmpdata
                p = p + num_mol_per_line[n]
                
            dry_air_col = sum(moldata[1:])
            if self.nmol <= 7:
                self.mol_den[0:self.nmol,v] = moldata[0:self.nmol]
            else:
                self.mol_den[0:7,v] = moldata[0:7]
                self.mol_den[7:self.nmol,v]  = moldata[8:]
            self.mol_den[self.nmol,v] = moldata[7]
            self.mol_vmr[:,v] = self.mol_den[:,v] / dry_air_col
        
            R_d_over_R_v = 287.0 / 461.0
            self.q[v] = R_d_over_R_v * moldata[0] / dry_air_col

        f.close()

        self.zlevel1[1:] = self.zlevel2[0:-1]
        self.plevel1[1:] = self.plevel2[0:-1]
        self.tlevel1[1:] = self.tlevel2[0:-1]

    def __repr__(self):
        return \
            'filename : {0:s}\n'.format(self.filename) + \
            'comment  : {0:s}\n'.format(self.comment) + \
            'nlayer   : {0:d}\n'.format(self.nlayer) + \
            'desc     : {0:s}\n'.format(self.desc) + \
            'H1       : {0:g}\n'.format(self.H1) + \
            'H2       : {0:g}\n'.format(self.H2) + \
            'ang      : {0:g}\n'.format(self.ang)

