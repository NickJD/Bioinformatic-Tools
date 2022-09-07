import gzip
import glob
import collections

import numpy as np


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

CDS_out = open('/home/nick/Nextcloud/General/Collab/Fiona/Corrected/Original_CDSs.fasta', 'w', newline='\n', encoding='utf-8')
StORF_out = open('/home/nick/Nextcloud/General/Collab/Fiona/Corrected/StORFs.fasta', 'w', newline='\n', encoding='utf-8')
Genome_INFO = collections.defaultdict(list)

Stats = []

for gff_file in glob.glob('/home/nick/Nextcloud/General/Collab/Fiona/Corrected/*StORF_Reporter_Combined.gff'):
    gff_in = open(gff_file, 'r')
    Genes = collections.defaultdict(list)
    Sequences = collections.defaultdict(str)
    Original,StORFs = 0,0
    FASTA = False
    for line in gff_in:  # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
        line_data = line.split()
        genome = gff_file.split('/Cleaned_')[1].split('_seq_StORF')[0]
        try:
            if line.startswith('\n'):  # Not to crash on empty lines in GFF
                continue
            elif 'PseudoCAP' in line_data[1] and 'CDS' in line_data[2]: # could be ID=gene in line[8]
                start = int(line_data[3])
                stop = int(line_data[4])
                pos = str(start) + '_' + str(stop)
                frame = line_data[6]
                contig = line_data[0]
                Genes[contig].append({pos:['PseudoCAP',start,stop,frame]})
                Original +=1
            elif 'StORF_Reporter' in line_data[1] and 'CDS' in line_data[2]:
                start = int(line_data[3])
                stop = int(line_data[4])
                pos = str(start) + '_' + str(stop)
                frame = line_data[6]
                contig = line_data[0]
                Genes[contig].append({pos:['StORF_Reporter',start,stop,frame]})
                StORFs +=1
        except IndexError:
            if line.startswith('>'):
                FASTA = True
                Current_Contig = line.strip().replace('>','')
            elif FASTA == True:
                Sequences[Current_Contig] = line.strip()

    proportion = (StORFs/Original)*100
    Stats.append(int(proportion))
    Genome_INFO[genome] = [Original,StORFs,int(proportion)]

    for contig, seqs in Genes.items():
        contig_sequence = Sequences[contig]
        reverse_contig_sequence = revCompIterative(contig_sequence)
        print(contig)
        contig_sequence_length = len(contig_sequence)
        for seq in seqs: # This needs to be rewritten
            print(seq)
            seq_id = list(seq.keys())
            seq_data = seq[seq_id[0]]
            start = seq_data[1]
            stop = seq_data[2]
            frame = seq_data[3]
            if '-' in frame:  # Reverse Compliment starts and stops adjusted
                r_start = contig_sequence_length - stop
                r_stop = contig_sequence_length - start
                if seq_data[0] == 'StORF_Reporter': # the remove first stop codon adjustment
                    current_seq = reverse_contig_sequence[r_start-3:r_stop+1]
                    StORF_out.write('>' + genome + '_' + contig + '_' + seq_data[
                        0] + '_' + str(start) + '_' + str(stop) + ';' + frame + '\n' + current_seq + '\n')
                else:
                    current_seq = reverse_contig_sequence[r_start:r_stop+1]
                    CDS_out.write('>'+genome+'_'+contig+'_'+seq_data[0]+'_'+str(start)+'_'+str(stop)+';'+frame+'\n'+current_seq+'\n')

            elif '+' in frame:
                if seq_data[0] == 'StORF_Reporter': # the remove first stop codon adjustment
                    current_seq = contig_sequence[start-4:stop]
                    StORF_out.write('>' + genome + '_' + contig + '_' + seq_data[
                        0] + '_' + str(start) + '_' + str(stop) + ';' + frame + '\n' + current_seq + '\n')
                else:
                    current_seq = contig_sequence[start-1:stop]
                    CDS_out.write('>' + genome + '_' + contig + '_' + seq_data[
                        0] + '_' + str(start) + '_' + str(stop) + ';' + frame + '\n' + current_seq + '\n')
            print(current_seq)




    print("Done")

print("Stats: " + str(np.average(Stats)))






