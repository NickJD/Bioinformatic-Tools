# This script was designed to specifically work with Ensembl data format
import argparse
from collections import defaultdict

def gff_load(gff_in):
    #Will code in different versions for different types of GFF3 files (Prodigal,Ensembl etc)
    prev_stop = 0
    overlaps = defaultdict(int)
    for line in gff_in:  # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
        line_data = line.split()
        try:
            if line.startswith('\n'): # Not to crash on empty lines in GFF
                continue
            elif 'ID=gene' in line_data[8]:
                start = int(line_data[3])
                stop = int(line_data[4])
                if start <= prev_stop:
                    overlap = abs(start - prev_stop)
                    overlaps[overlap] +=1
                prev_stop = stop
        except IndexError:
            continue



def extractor(options):
    gff_in = open(options.gff,'r')
    gff_load(gff_in)





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-gff', action='store', dest='gff', help='GFF annotation file for the FASTA',
                        required=True)
    options = parser.parse_args()
    extractor(options)