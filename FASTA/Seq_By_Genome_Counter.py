import argparse
from collections import defaultdict
import gzip


def fasta_load(fasta_in,string_match,string_count):
    first = True
    for line in fasta_in:
        line = line.strip()
        if line.startswith(';'):
            continue
        elif line.startswith('>') and not first:
            genomes[sequence_name] +=1
            sequence_name = line.split('>')[1].split('|')[0]
            if string_match:
                if string_match in line:
                    string_count[sequence_name] += 1
        elif line.startswith('>'):
            sequence_name = line.split('>')[1].split('|')[0]
            if string_match:
                if string_match in line:
                    string_count[sequence_name] += 1
        else:
            first = False
    genomes[sequence_name] +=1
    string_count[sequence_name] += 1
    return string_count



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta_seq', action='store', dest='fasta_seq', required=True,
                        help='FASTA file')
    parser.add_argument('-s', '--string_match', action='store', dest='string_match', required=False,
                        default='', help='Optional: Count number of sequences with "string" in ID')
    parser.add_argument('-o', '--output_file', action='store', dest='out_file', required=True,
                        help='Output file')

    options = parser.parse_args()
    string_count = defaultdict(int)
    genomes = defaultdict(int)
    try: # Detect whether fasta files are .gz or text and read accordingly
        fasta_in = gzip.open(options.fasta_seq,'rt')
        string_count = fasta_load(fasta_in,options.string_match,string_count)
    except:
        fasta_in = open(options.fasta_seq,'r')
        string_count = fasta_load(fasta_in,options.string_match,string_count)

    out = open(options.out_file, 'w', newline='\n', encoding='utf-8')
    out.write("Genome_ID\tNumber_Of_Seqs\t"+options.string_match+'\n')
    for genome_id, sequence_count in sorted(genomes.items(), key=lambda item: item[1], reverse=True):
        if string_count[genome_id] > 0:
            out.write(genome_id+'\t'+str(sequence_count)+'\t'+str(string_count[genome_id])+'\n')
        else:
            out.write(genome_id + '\t' + str(sequence_count) +'\n')








