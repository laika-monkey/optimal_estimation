
"""Determine filename to use for SHIS/CPL collocation data

As decided via e-mail discussion 2013-02-05.
"""

from argparse import ArgumentParser
from datetime import datetime, timedelta
from pyhdf import SD

def main():

    parser = ArgumentParser(description=__doc__)
    parser.add_argument('colloc_file', help='HDF collocation output file')
    parser.add_argument('shis_file', help='NetCDF SHIS file')
    args = parser.parse_args()

    colloc_sd = SD.SD(args.colloc_file)
    shis_sd = SD.SD(args.shis_file)

    idx_sds = colloc_sd.select('SHIS_Index')
    idx_ini, idx_fin = idx_sds[0], idx_sds[-1]

    t_ini = shis_time(shis_sd, idx_ini)
    t_fin = shis_time(shis_sd, idx_fin)

    print('SHIS.CPL.COLLOC.{0}.{1}.hdf'.format(t_ini, t_fin))

def shis_time(shis_sd, idx):

    year, month, day, second = (
        shis_sd.select(v)[idx] for v in ['refTimeYear', 'refTimeMonth', 'refTimeDay',
                                         'refTimeSec'])
    t = datetime(int(year), int(month), int(day)) + timedelta(seconds=second)
    return t.strftime('%y%m%d%H%M%S')

main()

