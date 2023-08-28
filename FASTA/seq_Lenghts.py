import os
import gzip
import math

# Define a function to parse a fasta file and return the number of sequences and the maximum sequence length
def parse_fasta_file(file_path):
    num_sequences = 0
    max_seq_len = 0
    with gzip.open(file_path, 'rt') as file:
        try:
            seq_len = 0
            for line in file:
                if line.startswith('>'):
                    num_sequences += 1
                    max_seq_len = max(max_seq_len, seq_len)
                    seq_len = 0
                else:
                    seq_len += len(line.strip())
            max_seq_len = max(max_seq_len, seq_len)
        except Exception as e:
            print('Exception: ' + e)
    return num_sequences, max_seq_len

# Define a function to calculate the median number of sequences per file
def median_sequences_per_file(directory_path):
    file_sequence_counts = []
    total_num_seqs = 0
    for file_name in os.listdir(directory_path):
        if file_name.endswith('singleline.fasta.gz'):
            file_path = os.path.join(directory_path, file_name)
            num_sequences, _ = parse_fasta_file(file_path)
            file_sequence_counts.append(num_sequences)
            total_num_seqs += num_sequences
    file_sequence_counts.sort()
    if len(file_sequence_counts) % 2 == 0:
        median_index = len(file_sequence_counts) // 2
        median_value = (file_sequence_counts[median_index] + file_sequence_counts[median_index - 1]) / 2
    else:
        median_index = len(file_sequence_counts) // 2
        median_value = file_sequence_counts[median_index]
    return median_value,total_num_seqs

# Define a function to calculate the median length and maximum length of all sequences
def median_and_max_seq_lengths(directory_path):
    all_seq_lengths = []
    for file_name in os.listdir(directory_path):
        try:
            if file_name.endswith('singleline.fasta.gz'):
                file_path = os.path.join(directory_path, file_name)
                with gzip.open(file_path, 'rt') as file:
                    seq_len = 0
                    for line in file:
                        if not line.startswith('>') and not line.startswith('#') and not line.startswith('\n'):
                            seq_len = len(line.strip())
                            all_seq_lengths.append(seq_len)
        except Exception as e:
            print('Exception: ' + e)

    all_seq_lengths.sort()
    if len(all_seq_lengths) % 2 == 0:
        median_index = len(all_seq_lengths) // 2
        median_value = (all_seq_lengths[median_index] + all_seq_lengths[median_index - 1]) / 2
    else:
        median_index = len(all_seq_lengths) // 2
        median_value = all_seq_lengths[median_index]
    seq_lengths_mean = sum(all_seq_lengths) / len(all_seq_lengths)
    seq_lengths_variance = sum([((x - seq_lengths_mean) ** 2) for x in all_seq_lengths]) / len(all_seq_lengths)
    seq_lengths_std = math.sqrt(seq_lengths_variance)
    max_seq_length = max(all_seq_lengths)
    return median_value, seq_lengths_std, max_seq_length, all_seq_lengths


directory_path = '/home/nick/Documents/StORF-Reporter/E-coli/Filtered/Filtered_UR_StORFs'
median_sequences, total_num_seqs  = median_sequences_per_file(directory_path)
median_length, std_length, max_length, all_seq_lengths = median_and_max_seq_lengths(directory_path)
print('total number of sequences : {}'.format(len(all_seq_lengths)))
print('Median number of sequences per file: {}'.format(median_sequences))
print('Median sequence length: {}'.format(median_length))
print('Standard deviation of sequence lengths: {}'.format(std_length))
print('Maximum sequence length: {}'.format(max_length))
