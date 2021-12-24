



fasta_in = open('/home/nick/Git/CoDing_Sequence_Frame_Prediction_From_DNA_Reads/classifier/Running/megahit_assembly_contigs_Min_1000_aligned_gff_unaligned_reads.fa','r')

IDs_in = open('/home/nick/Git/CoDing_Sequence_Frame_Prediction_From_DNA_Reads/classifier/Running/megahit_assembly_contigs_Min_1000_aligned_unaligned_reads_IDs','r')

out = open('/home/nick/Git/CoDing_Sequence_Frame_Prediction_From_DNA_Reads/classifier/Running/megahit_assembly_contigs_Min_1000_aligned_gff_BUT_not_CDS_aligned_reads.fa','w')


ID_List = {}

testing_count = 0
for line in IDs_in:
    line = line.strip()
    ID_List.update({line.strip():None})


print("Here")


first = True

for line in fasta_in:
    line = line.strip()
    if line.startswith('>') and first == False:  # Check if first seq in file
        if seq_id not in ID_List:
            out.write(seq_id+'\n'+seq+'\n')
        seq = ''
        seq_id = line
    elif line.startswith('>'):
        seq = ''
        seq_id = line
    else:
        seq += str(line)
        first = False

if seq_id not in ID_List:
    out.write(seq_id+seq)

print("DONE")













