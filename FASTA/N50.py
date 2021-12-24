#!/usr/bin/env python

## calculate N50 from fasta file
## N50 = contig length such that half of the contigs are longer and 1/2 of contigs are shorter

import commands
import sys
import os
from itertools import groupby
import numpy

# from Bio import SeqIO

lengths = []

with open(sys.argv[1]) as fasta:
    ## parse each sequence by header: groupby(data, key)
    faiter = (x[1] for x in groupby(fasta, lambda line: line[0] == ">"))

    for record in faiter:
        ## join sequence lines
        seq = "".join(s.strip() for s in faiter.next())
        lengths.append(len(seq))

## sort contigs longest>shortest
all_len = sorted(lengths, reverse=True)
csum = numpy.cumsum(all_len)

print
"N: %d" % int(sum(lengths))
n2 = int(sum(lengths) / 2)

# get index for cumsum >= N/2
csumn2 = min(csum[csum >= n2])
ind = numpy.where(csum == csumn2)

n50 = all_len[ind[0]]
print
"N50: %s" % n50

## N90
nx90 = int(sum(lengths) * 0.90)

## index for csumsum >= 0.9*N
csumn90 = min(csum[csum >= nx90])
ind90 = numpy.where(csum == csumn90)

n90 = all_len[ind90[0]]
print
"N90: %s" % n90

## write lengths to file
with open('/mnt/Crucial/Frame_Pred/prodigal/trimmed_paired_SRR873595_combined_megahit_contigs_Min_1000.fa', 'w') as handle:
    handle.write('\n'.join(str(i) for i in lengths))