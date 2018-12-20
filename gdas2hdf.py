#!/usr/bin/env python

# This program takes as input an NCEP GDAS file in GRIB format and produces
# an output file in HDF format. The output file contains a subset of the
# information contained in the GRIB file. The subset includes only some of
# the GDAS variables (see below) and is also resticted to locations taken
# from the latitude and longitude arrays as determined from the match file
# and variable names specified on the command line. The following data is
# written to the HDF output file:
#
#     - Latitude
#     - Longitude
#           (These simply contain the N latitude / longitude pairs pulled
#            out of the match file.
#     - Pressure
#     - Surface_Temperature
#           (These also have dimension N: one value per input latitude /
#            longitude.)
#     - Geopotential_Height
#     - Temperature
#     - Wind_U_Component
#     - Wind_V_Component
#           (These are arrays with dimensions Nx26. The first dimension
#           indicates the corresponding latitude / longitutde; the second
#           dimension indicates the pressure level.)
#     - Realtive_Humidity
#           (Like the above only with dimensions Nx21.)

# import pygrib
import Nio
import numpy
import os
from pyhdf.SD import SD, SDC
import sys

# process command line arguments
if len(sys.argv) == 1:
    print ('usage: gdas2hdf.py <output_file> <gdas_file> ' +
                              '<input_file> <lat_var_name> <lng_var_name> ' +
                              '[<output_sds_prefix>]')
    sys.exit(1)
out_path = sys.argv[1]
gdas_path = sys.argv[2]
match_path = sys.argv[3]
lat_var_name = sys.argv[4]
lng_var_name = sys.argv[5]
output_sds_prefix = sys.argv[6] if len(sys.argv) == 7 else ''

# pull the latitudes and longitudes out of the match file
# (note that we expect the dimensions of the given latitude and longitude
# variables to be 1xN)
match = SD(match_path, SDC.READ)
sds = match.select(lat_var_name)
lats = sds[:].flatten()
sds.endaccess()
sds = match.select(lng_var_name)
lngs = sds[:].flatten()
sds.endaccess()
match.end()
if lats.shape != lngs.shape:
    print ('ERROR: size mismatch between latitude array ' + str(lats.shape) +
           ' and longitude array ' + str(lngs.shape))
    sys.exit(1)

# determine the GRIB indices we're interested in
# (the GDAS grid's corners are (90N,0) at indices (0,0) and
# (90S, 1W) at (180, 359)
lat_indices = []
lng_indices = []
fill_indices = []
for i in range(len(lats)):
    lat = lats[i]
    lng = lngs[i]
    bad_lat = False
    bad_lng = False
    if (lat < -90.0) or (lat > 90.0):
        bad_lat = True
    if (lng < -180.0) or (lng > 180.0):
        bad_lng = True
    if bad_lat or bad_lng:
        if not (bad_lat and bad_lng):
            print ('WARNING: invalid latitude / longitude pair: ' +
                   str((lat, lng)))
        fill_indices.append(i)
        lat_indices.append(0)
        lng_indices.append(0)
    else:
        lat_indices.append(90 - int(round(lat)))
        index = int(round(lng))
        if index < 0:
            index += 360
        lng_indices.append(index)

# helper function for checking GRIB array sizes
def check_gdas_array_size(var_name, expected_shape):
    if gdas.variables[var_name].shape != expected_shape:
        print ('ERROR: unexpected shape for GDAS variable ' + var_name + ': ' +
               str(gdas.variables[var_name].shape))
        os.remove(out_path)
        sys.exit(1)

# open input and output files
gdas = Nio.open_file(gdas_path, format='grib')
# gdas = pygrib.open(gdas_path)
out = SD(out_path, SDC.WRITE | SDC.CREATE)

# deal with the single-dimension arrays (everything but the profiles)
def write_variable(hdf_var, dim_name, units, data_type, data):
    sds = out.create(output_sds_prefix + hdf_var, data_type, len(data))
    sds.units = units
    sds.dim(0).setname(dim_name)
    sds[:] = data
    sds.endaccess()
write_variable('Latitude', 'Location', 'deg', SDC.FLOAT64, numpy.round(lats))
write_variable('Longitude', 'Location', 'deg', SDC.FLOAT64, numpy.round(lngs))
check_gdas_array_size('TMP_3_SFC_10', (181, 360))
data = gdas.variables['TMP_3_SFC_10'].get_value()[lat_indices,lng_indices]
data[fill_indices] = -8888.0
write_variable('Surface_Temperature', 'Location', 'K', SDC.FLOAT32, data)
pressure_levels = (1000.0, 975.0, 950.0, 925.0, 900.0,
                   850.0, 800.0, 750.0, 700.0, 650.0,
                   600.0, 550.0, 500.0, 450.0, 400.0,
                   350.0, 300.0, 250.0, 200.0, 150.0,
                   100.0, 70.0, 50.0, 30.0, 20.0,
                   10.0)
write_variable('Pressure', 'Level', 'hPa', SDC.FLOAT32, pressure_levels)

# deal with the profiles
def transfer_profile(grib_var, hdf_var, units, min_level, max_level,
                     lat_indices, lng_indices):
    sds = out.create(output_sds_prefix + hdf_var, SDC.FLOAT32,
                     (len(lat_indices), 26))
    sds.units = units
    sds.dim(0).setname('Location')
    sds.dim(1).setname('Level')
    check_gdas_array_size(grib_var, (max_level - min_level, 181, 360))
    data = numpy.empty((len(lat_indices), 26), numpy.float32)
    data[:,:min_level] = 0.0
    data[:,min_level:max_level] = (gdas.variables[grib_var].get_value()
                                       [::-1,lat_indices,lng_indices]
                                       .transpose())
    data[:,max_level:] = 0.0
    data[fill_indices,:] = -8888.0
    sds[:] = data
    sds.endaccess()
profiles = (('HGT_3_ISBL_10', 'Geopotential_Height', 'gpm', 0, 26),
            ('TMP_3_ISBL_10', 'Temperature', 'K', 0, 26),
            ('R_H_3_ISBL_10', 'Relative_Humidity', '%', 0, 21),
            ('U_GRD_3_ISBL_10', 'Wind_u_Component', 'm/s', 0, 26),
            ('V_GRD_3_ISBL_10', 'Wind_v_Component', 'm/s', 0, 26),
            ('O3MR_3_ISBL_10', 'Ozone_Mixing_Ratio', 'kg/kg', 20, 26))
for grib_name, hdf_name, units, min_level, max_level in profiles:
    transfer_profile(grib_name, hdf_name, units, min_level, max_level,
                     lat_indices, lng_indices)

# close up shop
gdas.close()
out.end()
