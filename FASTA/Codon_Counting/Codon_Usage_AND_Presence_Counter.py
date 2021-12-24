import argparse
import collections
import gzip
import re

## Currently only works for stop codons

def fasta_load(fasta_in):
    first = True
    dna_regions = collections.OrderedDict()
    for line in fasta_in:
        line = line.strip()
        if line.startswith('>') and first == False:  # Check if first seq in file
            dna_region_length = len(seq)
            dna_regions.update({dna_region_id: (seq, dna_region_length, list(), None)})
            seq = ''
            dna_region_id = line.split()[0].replace('>', '')
        elif line.startswith('>'):
            seq = ''
            dna_region_id = line.split()[0].replace('>', '')
        else:
            seq += str(line)
            first = False
    dna_region_length = len(seq)
    dna_regions.update({dna_region_id: (seq, dna_region_length, list(), None)})
    return dna_regions

def gff_load(gff_in,dna_regions):
    #Will code in different versions for different types of GFF3 files (Prodigal,Ensembl etc)
    for line in gff_in:  # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
        #temp
        #line = line.replace('NZ_','')
        line_data = line.split()
        if line.startswith('\n'): # Not to crash on empty lines in GFF
            continue
        elif 'ID=gene' in options.gene_ident:
            if line_data[0] in dna_regions and options. gene_ident in line_data[8]:
                pos = line_data[3] + line_data[6] + line_data[4]
                dna_regions[line_data[0]][2].append(pos) # This will add to list
        else:
            gene_types = options.gene_ident.split(',')
            try:
                if line_data[0] in dna_regions:
                    if any(gene_type in line_data[2] for gene_type in gene_types): # line[2] for normalrun

                        pos = line_data[3] + line_data[6] + line_data[4]
                        dna_regions[line_data[0]][2].append(pos) # This will add to list
            except IndexError:
                continue
    return dna_regions

def revCompIterative(watson): #Gets Reverse Complement
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

def codon_Counting(dna_regions,codons_Counted,codons):
    total_Codons_Count = 0
    total_Codons_Usage = 0
    for region, data in dna_regions.items():
        for codon in codons.split(','):  # Find all Stops in seq
            codons_Counted[codon+'_DNA'] += data[0].count(codon)
            sequence_rev = revCompIterative(data[0])
            codons_Counted[codon+'_DNA'] += sequence_rev.count(codon)

        total_Codons_Count = sum(codons_Counted.values())

        sequence_rev = revCompIterative(data[0])
        seq_length = len(sequence_rev)
        for position in data[2]:  # need to change this for other codons
            if '+' in position:
                end = int(position.split('+')[1])
                codon = data[0][end-3:end]
                if codon in codons:  # filter out unwanted codons
                    codons_Counted[codon+'_Usage'] += 1
                    total_Codons_Usage +=1
            if '-' in position:
                end_rev = int(position.split('-')[0])
                r_Stop = seq_length - end_rev
                codon = sequence_rev[r_Stop-2:r_Stop+1]
                if codon in codons: # filter out unwanted codons
                    codons_Counted[codon+'_Usage'] += 1
                    total_Codons_Usage +=1

    for codon in codons.split(','):  # Find all Stops in seq
        codons_Counted[codon + '_DNA_Percentage'] = (codons_Counted[codon + '_DNA'] / total_Codons_Count) * 100
    for codon in codons.split(','):  # Find all Stops in seq
        codons_Counted[codon + '_Usage_Percentage'] = (codons_Counted[codon + '_Usage'] / total_Codons_Usage) * 100

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', action='store', dest='fasta', required=True,
                        help='FASTA file for codon counting')
    parser.add_argument('-gff', action='store', dest='gff', help='GFF annotation file for FASTA file',
                        required=True)
    parser.add_argument('-c', '--codons', action='store', dest='codons', required=True,
                        default='', help='Which codonds to count - "TAG,TGA,TAA"')
    parser.add_argument('-gene_ident', action='store', dest='gene_ident', default='ID=gene',
                        help='Identifier used for extraction of "genic" regions "CDS,rRNA,tRNA": Default for Ensembl_Bacteria = "ID=gene"')

    parser.add_argument('-o', '--output_file', action='store', dest='out_file', required=False,
                        help='Output file')

    options = parser.parse_args()
    codons_Counted = collections.defaultdict(int)
    try:  # Detect whether fasta/gff files are .gz or text and read accordingly
        fasta_in = gzip.open(options.fasta, 'rt')
        dna_regions = fasta_load(fasta_in)
    except:
        fasta_in = open(options.fasta, 'r')
        dna_regions = fasta_load(fasta_in)
    try:
        gff_in = gzip.open(options.gff, 'rt')
        dna_regions = gff_load(gff_in, dna_regions)
    except:
        gff_in = open(options.gff, 'r')
        dna_regions = gff_load(gff_in, dna_regions)
    print(options.fasta)
    codon_Counting(dna_regions,codons_Counted,options.codons)

    if options.out_file:
        out = open(options.out_file, 'a', newline='\n', encoding='utf-8')
        out.write(options.fasta+'\n')
        for codon, count  in codons_Counted.items():
            print(codon+'\t'+format(count,'.2f'))
            out.write(codon+'\t'+format(count,'.2f')+'\n')
    else:
        for codon, count  in codons_Counted.items():
            print(codon+'\t'+format(count,'.2f'))








