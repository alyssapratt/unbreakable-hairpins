import sys, re, math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess

#############
#Subroutines#
#############

def get_mfes(info_file):
    hairpin_rna_ids = []
    hairpin_lengths = []
    hairpin_mfes = []
    with open(info_file, 'r') as F:
        for line in F:
            rna_id,hp_length,mfe = line.strip().split('\t')
            if mfe is not None:
                hairpin_rna_ids.append(rna_id)
                hairpin_lengths.append(int(hp_length))
                hairpin_mfes.append(float(mfe))
        F.close()
    return hairpin_rna_ids,hairpin_lengths,hairpin_mfes

def get_mfes_control(unbreakable_rna_ids, control_info_file):
    hairpin_rna_ids = []
    hairpin_lengths = []
    hairpin_mfes = []
    with open(control_info_file, 'r') as F:
        for line in F:
            rna_id,hp_length,mfe = line.strip().split('\t')
            if rna_id not in unbreakable_rna_ids:
                hairpin_rna_ids.append(rna_id)
                hairpin_lengths.append(int(hp_length))
                hairpin_mfes.append(float(mfe))
        F.close()
    return hairpin_rna_ids,hairpin_lengths,hairpin_mfes

def get_hairpin_number(rna_id):
    return int(rna_id.split('H')[1].split('_')[0])

def make_filter_set(lineages_file):
    filter_ids = set()
    with open(lineages_file, 'r') as F:
        for line in F:
            id,accession,bpRNA_id,lineage,domain,url = line.strip().split('\t')
            filter_ids.add(int(id))
        F.close()
    return filter_ids

######
#Main#
######

usage = "usage: " + sys.argv[0] + " <unbreakable length/mfe info file> <control hairpins length/mfe info file> <lineages file for filter> <scatterplot filename>"
if len(sys.argv) != 5:
    print(usage)
    sys.exit()

unbreakable_file = sys.argv[1]
control_file = sys.argv[2]
lineages_file = sys.argv[3]
figure_file = sys.argv[4]
figure_file_no_type = figure_file[0:-4]

unbreakable_rna_ids,unbreakable_lengths,unbreakable_mfes = get_mfes(unbreakable_file)
control_rna_ids,control_lengths,control_mfes = get_mfes_control(unbreakable_rna_ids, control_file)
filter_set = make_filter_set(lineages_file)

unbreakable_hairpin_numbers = []
breakable_hairpin_numbers = []

for id in unbreakable_rna_ids:
    num_id = get_hairpin_number(id)
    if num_id in filter_set:
        unbreakable_hairpin_numbers.append(num_id)


for id in control_rna_ids:
    num_id = get_hairpin_number(id)
    if num_id in filter_set:
        breakable_hairpin_numbers.append(num_id)

my_bins = [1,2,3,4,5,6,7,8,9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,20, 21, 22, 23, 24, 25, 26, 27, 28, 29,30, 31, 32, 33, 34, 35, 36, 37, 38, 39,40, 41, 42, 43, 44, 45, 46, 47, 48, 49,50, 51, 52, 53, 54, 55, 56, 57, 58, 59,60, 61, 62, 63, 64, 65, 66, 67, 68, 69,70, 71, 72]

plt.figure()
plt.hist(breakable_hairpin_numbers,my_bins, density=True, color="#92c5de", alpha=0.5,label='Control Hairpins')
plt.hist(unbreakable_hairpin_numbers, my_bins, density=True, color="#ca0020",alpha=0.5,label='Unbreakable Hairpins')
plt.ylabel("Density")
plt.xlabel("bpRNA-1m Hairpin Number")
plt.legend(loc='upper left')
plt.savefig(figure_file)
