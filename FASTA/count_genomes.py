import collections
import gzip

def get_Genus(clustered):
    clustered_genus = clustered.split('|')[0]
    if '_' in clustered_genus[0]:  # Remove name error
        clustered_genus = clustered_genus.split('_')[1]
    else:
        clustered_genus = clustered_genus.split('_')[0]
    return str(clustered_genus).capitalize()

def get_Species(clustered):
    clustered_species = clustered.split('|')[0]
    if '_' in clustered_species[0]:  # Remove name error
        clustered_species = clustered_species.split('_')[1]
    else:
        clustered_species = clustered_species.split('_')[:2]
    return str('_'.join(clustered_species)).capitalize()


def fasta_read(fasta_in,genomes):
    for line in fasta_in:
        if line.startswith('>'):
            if '_' in line[1]:
                line = line.split('.')[0].replace('>_','')
            else:
                line = line.split('.')[0].replace('>','')
            if line not in genomes:
                genomes[line] +=1
                genus = get_Genus(line)
                genera[genus] +=1
                species = get_Species(line)
                specieses[species] +=1



genomes_in = '/home/nick/Desktop/Ensembl_CDS_Genes_Combined.fasta.gz'

genomes = collections.defaultdict(int)
genera = collections.defaultdict(int)
specieses = collections.defaultdict(int)

try:  # Detect whether fasta files are .gz or text and read accordingly
    fasta_in = gzip.open(genomes_in, 'rt')
    fasta_read(fasta_in,genomes)
except:
    fasta_in = open(genomes_in, 'r')
    fasta_read(fasta_in)

outfile_genera = open('counted_genera.txt','w')
outfile_species = open('counted_species.txt','w')

for key,value in genera.items():
    outfile_genera.write(key+','+str(value)+'\n')

for key,value in specieses.items():
    outfile_species.write(key+','+str(value)+'\n')


