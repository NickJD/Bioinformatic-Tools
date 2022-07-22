import collections
import sys

file_name = sys.argv[1:]


myfile = open(file_name[0],'r')

num_classified = 0

my_dict = collections.defaultdict(int)

for line in myfile:
    line_data = line.split('\t')
    if line_data[0] != "U":
        print(line)
        num_classified +=1
        my_dict[line_data[2]] +=1
print(num_classified)


print(my_dict)


tax_of_interest = my_dict[file_name[1]]

print(tax_of_interest)

per = 100 * (tax_of_interest/num_classified)

print(per)
