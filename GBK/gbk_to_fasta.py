import sys
def extract_cds_proteins(genbank_file):
    protein_records = []
    current_protein_id = ""
    current_protein_sequence = ""
    reading_sequence = False

    with open(genbank_file, "r") as file:
        current_locus_tag = None
        for line in file:
            if line.startswith('ORIGIN'):
                break
            if line.endswith('"\n'):
                reading_sequence = False
            stripped_line = line.strip()
            if stripped_line.startswith("CDS"):
                if current_locus_tag and current_protein_sequence:
                    if 'complement' in line:
                        start_index = line.index('(') + 1
                        end_index = line.index(')')
                        numbers = line[start_index:end_index].split('..')
                        start_number = int(numbers[0])
                        end_number = int(numbers[1])
                    else:
                        numbers = line.split()[1].split('..')
                        start_number = int(numbers[0])
                        end_number = int(numbers[1])

                    description = current_locus_tag+'_'+str(start_number)+'_'+str(end_number)
                    protein_records.append((description, current_protein_sequence))
                current_locus_tag = ""
                current_protein_sequence = ""
                reading_sequence = False
            elif stripped_line.startswith("/locus_tag"):
                current_locus_tag = line.split('"')[1]
            elif stripped_line.startswith("/translation"):
                current_protein_sequence = ""
                current_protein_sequence += line.strip().replace('/translation="','')
                reading_sequence = True
            elif reading_sequence:
                current_protein_sequence += line.strip().replace('"','')

    if current_locus_tag and current_protein_sequence:
        description = current_locus_tag + '_' + str(start_number) + '_' + str(end_number)
        protein_records.append((description, current_protein_sequence))

    return protein_records

def write_fasta(records, output_file):
    with open(output_file, "w") as f:
        for description, sequence in records:
            f.write(f">{description}\n")
            f.write(f"{sequence}\n")

if __name__ == "__main__":
    def process_file(file_name):
        # Your file processing code here
        print("Processing file:", file_name)
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <file_name>")
    else:
        input_genbank_file = sys.argv[1]
        output_fasta_file = input_genbank_file.replace('.gbk','.fasta')

    proteins = extract_cds_proteins(input_genbank_file)
    write_fasta(proteins, output_fasta_file)
    print(f"Protein sequences extracted and saved to {output_fasta_file}")
