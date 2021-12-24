import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--seqs', default='', help='Which seqs to analyse?')
args = parser.parse_args()

def gc_count(dna):
    c = 0
    a = 0
    g = 0
    t = 0
    n = 0
    for i in dna:
        if "C" in i:
            c += 1
        elif "G" in i:
            g += 1
        elif "A" in i:
            a += 1
        elif "T" in i:
            t += 1
        elif "N" in i:
            n += 1
    gc_content = (g + c) * 100 / (a + t + g + c + n)
    n_per = n * 100 / (a + t + g + c + n)
    return gc_content


def seq_Lengths(seqs):
    all_gc = []
    seq = ""
    with open(seqs, 'r') as seqs:
        for line in seqs:
            line = line.replace("\n","")
            # if line.startswith('>'):
            #     seq = ""
            if not line.startswith('>'):
                seq += str(line)

    all_gc.append(gc_count(seq))

    print("Number of Seqs: " + str(len(all_gc)))
    print("Longest Seq: " + str(max(all_gc)))
    print("Shortest Seq: " + str(min(all_gc)))
    print("Median of Seqs: " + str(np.median(all_gc)))
    print("Mean of Seqs: " + format(np.mean(all_gc), '.2f'))
    print("Std: " + format(np.std(all_gc),'.2f'))
    print("25th Percentile " + format(np.percentile(all_gc, 25), '.2f'))
    print("75th Percentile " + format(np.percentile(all_gc,75),'.2f'))


if __name__ == "__main__":
    seq_Lengths(**vars(args))








