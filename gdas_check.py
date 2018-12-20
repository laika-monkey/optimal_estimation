#!/usr/bin/env python

import gdas
import ndbc

import matplotlib.pyplot as plt
import Nio
import sys
import time

if len(sys.argv) != 6:
    print 'usage: ' + sys.argv[0] + ' <year> <month> <buoy_id> <latitude> <longitude>'
    sys.exit(1)
year = int(sys.argv[1])
month = int(sys.argv[2])
buoy_id = sys.argv[3]
lat = float(sys.argv[4])
lon = float(sys.argv[5])

lat_idx = 90 - int(round(lat))
lon_idx = int(round(lon))
if lon_idx < 0:
    lon_idx += 360

ts = time.strptime(str(year) + ' ' + str(month), '%Y %m')
time_str = time.strftime('%b %Y', ts)

gdas_temps = []
buoy_temps = []
buoy_data = ndbc.BuoyData(buoy_id, year)
while True:

    if ts.tm_mon != month:
        break

    gdas_file = gdas.open(ts)
    gdas_temps.append(gdas_file.variables['TMP_3_SFC_10'][lat_idx,lon_idx])
    gdas_file.close()

    buoy_idx = 24 * (ts.tm_yday - 1) + ts.tm_hour
    buoy_temp = buoy_data.water_temp[buoy_idx]
    if buoy_temp != ndbc.fill_value:
        buoy_temps.append(buoy_temp)
    else:
        buoy_temps.append(float('nan'))

    ts = time.localtime(time.mktime(ts) + 6 * 60 * 60)

plt.title('%s Lat: %.1f, Lon: %.1f' % (time_str, lat, lon))
plt.ylabel('Surface Temparature')
plt.plot(gdas_temps, label = 'GDAS')
plt.plot(buoy_temps, label = 'Buoy ' + str(buoy_id))
plt.legend()
plt.show()
