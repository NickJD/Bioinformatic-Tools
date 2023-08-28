
import argparse
import textwrap




###### CURRENTLY ONLY WORKS ON SINGLE CONTIG FASTA

parser = argparse.ArgumentParser()
parser.add_argument('-dna', action="store", dest='dna', required=True,
                    help='Which DNA file to modify?')
parser.add_argument('-loci', action="store", dest='loci', required=True,
                    help='Which base pair positions to modify? - enter as "500|G,804|T,12044|A"')
parser.add_argument('-oname', action="store", dest='oname', required=False,
                    help='Default - Appends \'_Point_Mutated.fasta\' to end of input FASTA filename')
args = parser.parse_args()



def point_Mutation(dna,loci,oname):

    if oname:
        outfile = open(oname, 'w', newline='\n', encoding='utf-8')
    else:
        outfile = open(dna.split('.')[0] + '_Point_Mutation.fasta', 'w', newline='\n', encoding='utf-8')

    loci_dict = {}
    for i in loci.split(','):
        loci_dict.update({int(i.split('|')[0])-1:i.split('|')[1]}) # '-1' is to account for 1-base GFF system


    with open(dna, mode='r') as genome:
        first_line = genome.readline()
        genome_Seq = "".join(line.rstrip() for line in genome if not line.startswith('>'))

    # Replace multiple characters with different replacement characters
    mutated_genome_Seq = genome_Seq
    for index, replacement in loci_dict.items():
        mutated_genome_Seq = mutated_genome_Seq[:index] + loci_dict[index] + mutated_genome_Seq[index:]

    outfile.write(first_line)
    wrapped = textwrap.wrap(mutated_genome_Seq, width=60)
    for wrap in wrapped:
        outfile.write(wrap + '\n')


if __name__ == "__main__":
    point_Mutation(**vars(args))









