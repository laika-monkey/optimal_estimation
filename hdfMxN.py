
"""hdfMxN - convert Fred's ASCII collocation output to HDF"""

from argparse import ArgumentParser
from itertools import count
import sys

import numpy as np
from pyhdf import SD

def main():

    parser = ArgumentParser(description=__doc__)
    parser.add_argument('input_file', help='ASCII collocation output file')
    parser.add_argument('-m', '--master-sds-names',
                        help='comma-separated list of master index SDS names')
    parser.add_argument('-s', '--slave-sds-names',
                        help='comma-separated list of slave index SDS names')
    parser.add_argument('-o', '--output-file',
                        help='name of output file (overrides name found in ASCII input)')
    args = parser.parse_args()

    # first pass through input file: determine dimension sizes needed
    # for output arrays
    num_masters = 0
    max_num_slaves = 0
    f = open(args.input_file)
    output_filename = f.readline()
    while True:

        # each record begins with ***; we're done when we reach EOF
        stars = f.readline()
        if not stars:
            break
        assert stars == '***\n'

        num_masters += 1

        # pull number of slave pixels for the current master out of the
        # record's first line
        master_line = f.readline()
        num_slaves = int(master_line.split()[-1])
        max_num_slaves = max(max_num_slaves, num_slaves)

        # skip all the per-slave-pixel lines
        for _ in range(num_slaves):
            slave_line = f.readline()

    # end here if there were no records
    if not num_masters:
        print 'No collocations; exiting'
        sys.exit(101)

    # use the number of tokens in the last master and slave lines read
    # in to determine how many indexes the master and slave each use
    # (could probably also get this from the the input filename)
    num_indexes = {
        'master': len(master_line.split()) - 1,
        'slave': len(slave_line.split()) - 1}

    # allocate output arrays
    master_index = np.zeros([num_indexes['master'], num_masters],
                            np.int32)
    slave_index = np.zeros([num_indexes['slave'], num_masters, max_num_slaves],
                           np.int32)
    weight = np.zeros([num_masters, max_num_slaves],
                      np.float32)

    # second pass through input file: fill-in output arrays
    f.seek(0)
    output_filename = f.readline().strip()
    for i in count():

        # stars for next record or EOF?
        stars = f.readline()
        if not stars:
            break
        assert stars == '***\n'

        # master indexes and number of slave pixels come from first line
        master_tokens = f.readline().split()
        master_index[:,i] = [int(t) for t in master_tokens[:-1]]

        # slave indexes and weights come from additional lines
        for j in range(int(master_tokens[-1])):
            slave_tokens = f.readline().split()
            slave_index[:,i,j] = [int(t) for t in slave_tokens[:-1]]
            weight[i,j] = float(slave_tokens[-1])

    # determine SDS names for output file; use names given by -m or -s;
    # default to names like Master_Index_1, ...
    names = {}
    for side in ['master', 'slave']:
        sds_names_str = side + '_sds_names'
        given_names = getattr(args, sds_names_str)
        if given_names is not None:
            names[side] = getattr(args, sds_names_str).split(',')
            if len(names[side]) != num_indexes[side]:
                print ('wrong number of %s SDS names given ' % side +
                           '(expected %s)' % num_indexes[side])
                sys.exit(1)
        else:
            names[side] = ['%s%s_Index_%s' % (side[0].upper(), side[1:], i + 1)
                               for i in range(num_indexes[side])]

    # write the output arrays to HDF. use compression
    if args.output_file:
        output_filename = args.output_file
    sd = SD.SD(output_filename, SD.SDC.WRITE | SD.SDC.CREATE | SD.SDC.TRUNC)
    for i in range(num_indexes['master']):
        sds = sd.create(names['master'][i], SD.SDC.INT32,
                        [num_masters])
        sds.setcompress(SD.SDC.COMP_DEFLATE, 6)
        sds[:] = master_index[i]
        sds.endaccess()
    for i in range(num_indexes['slave']):
        sds = sd.create(names['slave'][i], SD.SDC.INT32,
                        [num_masters, max_num_slaves])
        sds.setcompress(SD.SDC.COMP_DEFLATE, 6)
        sds[:] = slave_index[i]
        sds.endaccess()
    sds = sd.create('Weights', SD.SDC.FLOAT32, [num_masters, max_num_slaves])
    sds.setcompress(SD.SDC.COMP_DEFLATE, 6)
    sds[:] = weight
    sds.endaccess()
    sd.end()

if __name__ == '__main__':
   main()

