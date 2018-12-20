import numpy
import os
import subprocess
import time

fill_value = 999.0

# for the given buoy ID and year, fetch the corresponding data
# file from the NDBC site
def fetch(buoy_id, year):
    filename = str(buoy_id) + 'h' + str(year) + '.txt'
    if os.path.exists('ndbc/' + filename):
        return
    filename += '.gz'
    url = 'http://www.ndbc.noaa.gov/data/historical/stdmet/' + filename
    subprocess.check_call(['wget', '--directory-prefix=ndbc', url])
    subprocess.check_call(['gunzip', 'ndbc/' + filename])

class BuoyData:

    def __init__(self, buoy_id, year):
        filename = 'ndbc/' + str(buoy_id) + 'h' + str(year) + '.txt'
        fetch(buoy_id, year)
        self.water_temp = numpy.empty((366 * 24,), numpy.float64)
        self.water_temp[:] = fill_value
        buoy_file = open(filename, 'r')
        while True:
            line = buoy_file.readline()
            if not line:
                break
            if line[0] == '#':
                continue
            line = line.split()
            ts = time.strptime(' '.join(line[:5]), '%Y %m %d %H %M')
            if ts.tm_min != 0:
                raise Exception('unexpected measurement at ' +
                                time.asctime(ts))
            slot = 24 * (ts.tm_yday - 1) + ts.tm_hour
            self.water_temp[slot] = float(line[14]) + 273.15
