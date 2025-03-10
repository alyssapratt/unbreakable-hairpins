import sys
import numpy as np

def find_diffs(all_unbreakables, split_unbreakables):
    unsplit_unbreakables = np.setdiff1d(all_unbreakables, split_hairpins)
    return unsplit_unbreakables

def process_unbreakables(tabbed_file):
    unbreakable_hairpins=[]
    with open(tabbed_file, 'r') as F:
        for line in F:
            id,shuffles,seq=line.split('\t')
            unbreakable_hairpins.append(id)
    return unbreakable_hairpins

def process_split_list(list_file):
    split_hairpins = []
    with open(list_file, 'r') as F:
        for line in F:
            id = line.split()
            split_hairpins.append(id)
    return split_hairpins
##########
# MAIN   #
##########

usage = "Usage: " + sys.argv[0] + " <full unbreakable tabbed file> <split list of IDs>"
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

all_unbreakable_file = sys.argv[1]
all_split_hairpins_file = sys.argv[2]

unbreakable_hairpins = process_unbreakables(all_unbreakable_file)
split_hairpins = process_split_list(all_split_hairpins_file)

unsplit_unbreakables = find_diffs(unbreakable_hairpins, split_hairpins)

output_file = "unsplit_T50s.txt"
F = open(output_file, 'w')
for value in split_hairpins:
    F.write(str(value) + '\n')
