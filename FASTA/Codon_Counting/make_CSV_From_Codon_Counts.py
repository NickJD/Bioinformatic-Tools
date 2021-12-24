import argparse
import collections
import csv

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',  action='store', dest='input', required=True,
                        help='Input file for codon counting')
    parser.add_argument('-o', action='store', dest='out_file', required=False,
                        help='Output file')

    options = parser.parse_args()
s

    in_file = open(options.input, 'r')
    for line in in_file:
        line = line.strip()
        line = line.split('\t')
        try:
            countings[line[0]].append(line[1])
        except KeyError:
            continue
        except IndexError:
            continue

    out_file = open(options.out_file,'w')
    out_file = csv.writer(out_file)
    for group, data in countings.items():
        out_file.writerow(data)



