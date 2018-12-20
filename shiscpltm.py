
"""Print the begin and end times of an SHIS file"""

from argparse import ArgumentParser
from datetime import datetime, timedelta
from pyhdf import SD

def main():

    parser = ArgumentParser(description=__doc__)
    parser.add_argument('shis_file', help='NetCDF SHIS file')
    args = parser.parse_args()

    shis_sd = SD.SD(args.shis_file)

    print shis_time(shis_sd, 0)
    print shis_time(shis_sd, -1)

def shis_time(shis_sd, idx):

    year, month, day, second = (
        shis_sd.select(v)[idx] for v in ['refTimeYear', 'refTimeMonth', 'refTimeDay',
                                         'refTimeSec'])
    return datetime(int(year), int(month), int(day)) + timedelta(seconds=second)

main()

