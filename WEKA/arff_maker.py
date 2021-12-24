import collections

test_data = open('/mnt/L_Data/Nextcloud/Nick/SpikeWork/colistin_all_genes.csv','r')

genomes = collections.defaultdict(list)

count = 0
First = True
for line in test_data:
    if First == False:
        line = line.split(',')
        data = line # pointless
        genome = data[0]+','+data[3]
        if data[5].startswith('1'): # gene hit
            if data[4] not in genomes[genome]: # only add 1 example of a gene for each genome
                genomes[genome].append(data[4])
        #count +=1
        # if len(genomes) == 5:
        #     break
    First = False

#the index of this list is important!! - it is used later
list_of_genes = list(set({ele for val in genomes.values() for ele in val}))

list_of_zeros = [0] * len(list_of_genes)

arff_data = collections.defaultdict(list)

for genome, genes in genomes.items():
    arff_data[genome] = [0] * len(list_of_genes)
    for gene in genes:
        idx = list_of_genes.index(gene)
        arff_data[genome][idx] = 1

print(','.join(map(str,list_of_genes)))
for genome, array in arff_data.items():
    array = ','.join(map(str, array))

    print(genome+','+array)


print("END")