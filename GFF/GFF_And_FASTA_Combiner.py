import argparse
import textwrap
import glob
import pathlib


parser = argparse.ArgumentParser()
parser.add_argument('-dir', action="store", dest='directory', required=True,
                    help='directory containing matching FASTA and GFF files')
args = parser.parse_args()


def Combiner(directory):

    gff_list = list(pathlib.Path(directory).glob('*.gff'))
    gff_list.extend(pathlib.Path(directory).glob('*.gff3'))
    gff_list = list(map(str, gff_list))

    fasta_list = list(pathlib.Path(directory).glob('*.fasta'))
    fasta_list.extend(list(pathlib.Path(directory).glob('*.fna')))
    fasta_list.extend(list(pathlib.Path(directory).glob('*.fa')))
    fasta_list = list(map(str, fasta_list))

    for gff in gff_list:
        gff_genome = gff.split('/')[-1].split('.')[0]
        fasta = list(filter(lambda x: gff_genome in x, fasta_list))
        with open(gff.replace('.gff3','_Combined.gff3'), 'a') as out_file:
            with open(gff) as infile:
                out_file.write(infile.read())
            with open(gff) as infile:
                out_file.write('##FASTA\n')
            with open(fasta[0]) as infile:
                out_file.write(infile.read())
if __name__ == "__main__":
    Combiner(**vars(args))





