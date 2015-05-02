#!/bin/python

#
# Produce a list of target regions from a primer bed. Primer bed has the form:
#
# chr19 1206861 1206878 STK11_t1_2_FPrimer
# chr19 1206957 1206978 STK11_t1_2_RPrimer
# chr19 1206905 1206924 STK11_t1_3_FPrimer
# chr19 1206996 1207019 STK11_t1_3_RPrimer
#
# and output looks like:
#
# chr19 1206878 1206957 STK11_t1_2
# chr19 1206924 1206996 STK11_t1_3
#
# Usage:
#   python primers_to_target_regions.py < primer_locations.bed

import sys

def iterblocks(iter, n):
    while 1:
        rval = []
        for i in range(n):
            rval.append(iter.next())
        yield rval

block_count = 0
for block in iterblocks(sys.stdin, 2):
    # Split lines into fields.
    fwd_primer_fields = block[0].split('\t')
    rvs_primer_fields = block[1].split('\t')

    # Validate that block makes sense.
    if fwd_primer_fields[0] != rvs_primer_fields[0] or \
       fwd_primer_fields[3] != rvs_primer_fields[3]:
        print "Error: mismatched block:"
        print "\t", block[0],
        print "\t", block[1],
        break

    # Print target region for primers.
    print "\t".join([fwd_primer_fields[0], 
                     fwd_primer_fields[2], 
                     rvs_primer_fields[1], 
                     fwd_primer_fields[3]]),
