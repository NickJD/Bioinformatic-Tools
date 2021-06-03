import argparse
from collections import  OrderedDict
import textwrap
import gzip

###################
gencode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
      'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}
############################


def translate_frame(sequence):
    translate = ''.join([gencode.get(sequence[3 * i:3 * i + 3], 'X') for i in range(len(sequence) // 3)])
    return translate



def revCompIterative(watson): #Gets Reverse Complement - May be used later
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                   'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M',
                   'M': 'K', 'V': 'B', 'B': 'V', 'H': 'D', 'D': 'H'}
    watson = watson.upper()
    watsonrev = watson[::-1]
    crick = ""
    for nt in watsonrev: # make dict to catch bad nts - if more than 1 output them to std error
        try:
            crick += complements[nt]
        except KeyError:
            crick += nt
    return crick

def fasta_load(fasta_in):
    first = True
    for line in fasta_in:
        line = line.strip()
        if line.startswith(';'):
            continue
        elif line.startswith('>') and not first:
            sequences.update({sequence_name: seq})
            seq = ''
            sequence_name = line
        elif line.startswith('>'):
            seq = ''
            sequence_name = line
        else:
            seq += str(line)
            first = False
    sequences.update({sequence_name: seq})


def write_fasta(sequence_id,amino_sequence):  # Some Lines commented out for BetaRun of ConStORF
    if options.gz == False:
        out_fasta = open(options.out_file, 'a', newline='\n', encoding='utf-8')
    elif options.gz == True:
        out_fasta = gzip.open(options.out_file +'.gz', 'at', newline='\n', encoding='utf-8')
    out_fasta.write(sequence_id+'\n')
    if options.line_wrap:
        wrapped = textwrap.wrap(amino_sequence, width=60)
        for wrap in wrapped:
            out_fasta.write(wrap + '\n')
    else:
        out_fasta.write(amino_sequence + '\n')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta_seq', action='store', dest='fasta_seq', required=True,
                        help='FASTA file')
    parser.add_argument('-o', '--output_file', action='store', dest='out_file', required=True,
                        help='Output file')
    parser.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')
    parser.add_argument('-lw', action="store", dest='line_wrap', default=False, type=eval, choices=[True, False],
                        help='Default - False: Line wrap FASTA sequence output at 60 chars')

    options = parser.parse_args()
    #### Blank output files
    if options.gz == False:
        out_fasta = open(options.out_file, 'w', newline='\n', encoding='utf-8')
    elif options.gz == True:
        out_fasta = gzip.open(options.out_file +'.gz', 'wt', newline='\n', encoding='utf-8')
    sequences = OrderedDict()
    try: # Detect whether fasta files are .gz or text and read accordingly
        fasta_in = gzip.open(options.fasta_seq,'rt')
        fasta_load(fasta_in)
    except:
        fasta_in = open(options.fasta_seq,'r')
        fasta_load(fasta_in)


    for sequence_id, sequence in sequences.items():
        amino_sequence = translate_frame(sequence[0:])
        write_fasta(sequence_id,amino_sequence)









