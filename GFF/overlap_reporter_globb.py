import gzip
import glob
import collections
import matplotlib.pyplot as plt
import sys
import matplotlib.ticker as mtick
from matplotlib.ticker import PercentFormatter
import seaborn as sns
import numpy as np


overlaps = collections.defaultdict(int)
cds_overlaps = collections.defaultdict(int)

count = 0
for gff_file in glob.glob('/home/nick/Desktop/Ensem/GFF/*.gff3.gz'):
    count +=1
    gff_in = gzip.open(gff_file, 'rt')
    prev_stop = 0
    cds_prev_stop = 0
    contig = ''
    for line in gff_in:  # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
        line_data = line.split()
        try:
            if line.startswith('\n'):  # Not to crash on empty lines in GFF
                continue
            elif 'ID=gene' in line_data[8]: # could be ID=gene in line[8]
                if line_data[0] != contig:
                    prev_stop = 0
                    cds_prev_stop = 0
                    contig = line_data[0]
                start = int(line_data[3])
                stop = int(line_data[4])
                if start <= prev_stop:
                    overlap = int(abs(start - prev_stop))
                    overlaps[overlap] += 1
                prev_stop = stop
            # elif 'CDS' in line_data[2]:
            #     if line_data[0] != contig:
            #         prev_stop = 0
            #         cds_prev_stop = 0
            #         contig = line_data[0]
            #     start = int(line_data[3])
            #     stop = int(line_data[4])
            #     if start <= cds_prev_stop:
            #         cds_overlap = int(abs(start - cds_prev_stop))
            #         cds_overlaps[cds_overlap] += 1
            #     cds_prev_stop = stop

        except IndexError:
            continue

    print(count)

    if count > 2000:
        break


#overlaps_tuple = sorted(overlaps.items())
temp = sorted(overlaps.items(), key=lambda x: int(x[0]))

overlaps = dict(temp)
#overlap_lengths = [tuple[0] for tuple in sorted(overlaps_tuple, key = lambda x:x[1])]
#
print("WW")
from statistics import median
#{k:median([v for v in overlaps[k]]) for k in overlaps.keys()}
#

# overlaps_tuple = sorted(overlaps_tuple, key=lambda  x:x[1])
#
# tmp = []
# for tuple in overlaps_tuple:
#     tmp.append(tuple[0])
# print(np.median(tmp))
# print(np.mean(tmp))
# print(max(tmp))

#sys.exit()
#print(np.)



ax = sns.displot(overlaps,binwidth=1, stat="probability")
#ax = sns.displot(cds_overlaps,binwidth=1, stat="probability")#,color='r')
ax = sns.ecdfplot(overlaps, stat="proportion")
#ax = sns.ecdfplot(cds_overlaps, stat="proportion",color='r')
# ax = sns.countplot(x = 'lengths', data = overlaps)
ax.set_title('Lengths of Gene Overlaps')
ax.set(ylabel='Proportion of Overlaps',xlabel='Overlap Lengths (nts)')
ax.set(ylim=(0, None))
ax.set(xlim=(0, 50))
plt.tight_layout()
ax.figure.savefig("Overlap_Lengths.pdf", bbox_inches='tight')
plt.xlim(0)

    # plt.bar(list(overlaps.keys()), overlaps.values(), color='b')
    # plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
#plt.show()
   # sys.exit()