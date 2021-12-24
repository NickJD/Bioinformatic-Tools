import collections

input = open('/home/nick/Desktop/Ensem/Biggest_Con-StORF_Only_Cluster_1326.fa','r')
output = open('/home/nick/Desktop/Ensem/Biggest_Con-StORF_Only_Cluster_1326_renamed.fa','w')

count = 0

genera = collections.defaultdict(int)

for line in input:
    if line.startswith('>'):
        line = line.split('_')
        genera[line[0]] +=1
        gen = line[0] + "_Con-StORF_" + str(genera[line[0]])
        output.write(gen+"\n")
    else:
        output.write(line.replace('*','-'))

print("W")



