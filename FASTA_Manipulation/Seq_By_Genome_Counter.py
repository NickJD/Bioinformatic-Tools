import argparse
from collections import defaultdict
import gzip


def fasta_load(fasta_in):
    first = True
    for line in fasta_in:
        line = line.strip()
        if line.startswith(';'):
            continue
        elif line.startswith('>') and not first:
            genomes[sequence_name] +=1
            sequence_name = line.split('>')[1].split('|')[0]
        elif line.startswith('>'):
            sequence_name = line.split('>')[1].split('|')[0]
        else:
            first = False
    genomes[sequence_name] +=1



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta_seq', action='store', dest='fasta_seq', required=True,
                        help='FASTA file')
    parser.add_argument('-o', '--output_file', action='store', dest='out_file', required=True,
                        help='Output file')

    options = parser.parse_args()
    genomes = defaultdict(int)
    try: # Detect whether fasta files are .gz or text and read accordingly
        fasta_in = gzip.open(options.fasta_seq,'rt')
        fasta_load(fasta_in)
    except:
        fasta_in = open(options.fasta_seq,'r')
        fasta_load(fasta_in)

    out = open(options.out_file, 'w', newline='\n', encoding='utf-8')
    for genome_id, sequence_count in sorted(genomes.items(), key=lambda item: item[1], reverse=True):
        out.write(genome_id+'\t'+str(sequence_count)+'\n')









