import Nio
import os
import subprocess

def make_filename(time_info):
    return 'gdas1.PGrbF00.%02d%02d%02d.%02dz' % (time_info.tm_year % 100,
                                                 time_info.tm_mon,
                                                 time_info.tm_mday,
                                                 time_info.tm_hour)

# for the given struct_time (as returned by time.localtime() for
# example), fetch the corresponding GDAS file
def fetch(time_info):
    filename = make_filename(time_info)
    if os.path.exists('gdas/' + filename):
        return
    dirname = '%04d_%02d_%02d_%03d' % (time_info.tm_year, time_info.tm_mon,
                                       time_info.tm_mday, time_info.tm_yday)
    url = ('ftp://ftp.ssec.wisc.edu/pub/eosdb/ancillary/' + dirname + '/' +
           filename)
    subprocess.check_call(['wget', '--directory-prefix=gdas', url])

def open(time_info):
    fetch(time_info)
    filename = make_filename(time_info)
    return Nio.open_file('gdas/' + filename, format='grib')
