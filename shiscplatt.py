
"""Add shiscpl input file names to HDF output"""

from argparse import ArgumentParser
import os
from pyhdf import SD

def main():

    parser = ArgumentParser(description=__doc__)
    parser.add_argument('shis_file', help='NetCDF SHIS file')
    parser.add_argument('cpl_file', help='XDR CPL file')
    parser.add_argument('colloc_file', help='HDF collocation output file')
    args = parser.parse_args()

    colloc_sd = SD.SD(args.colloc_file, SD.SDC.WRITE)
    colloc_sd.SHIS_File = os.path.basename(args.shis_file)
    colloc_sd.CPL_File = os.path.basename(args.cpl_file)

main()

