import sys, re, math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess

#############
#Subroutines#
#############

def get_stem_lengths(tabbed_struc_file):
    mod_stem_lengths = []
    orig_stem_lengths = []
    stem_length_differences = []
    with open(tabbed_struc_file, 'r') as F:
        for line in F:
            if not line.startswith("#"):
                rna_id,seq,mod_struc,orig_struc = line.strip().split('\t')
                if seq is not None:
                    orig_stem_length = min(orig_struc.count('('), orig_struc.count(')'))
                    mod_stem_length = min(mod_struc.count('('), mod_struc.count(')'))
                    mod_stem_lengths.append(mod_stem_length)
                    orig_stem_lengths.append(orig_stem_length)
                    stem_length_differences.append((mod_stem_length-orig_stem_length)/orig_stem_length)
        F.close()
    return mod_stem_lengths, orig_stem_lengths,stem_length_differences

######
#Main#
######

usage = "usage: " + sys.argv[0] + " <T1000 hairpin tabbed structure file> <control hairpin tabbed structure file> <scatterplot filename>"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

T1000_struc_file = sys.argv[1]
control_struc_file = sys.argv[2]
figure_file = sys.argv[3]

mod_T1000_stem_lengths, orig_T1000_stem_lengths,T1000_stem_change = get_stem_lengths(T1000_struc_file)
mod_control_stem_lengths, orig_control_stem_lengths,control_stem_change = get_stem_lengths(control_struc_file)



mplot = plt.figure()
vplot = plt.violinplot([control_stem_change,T1000_stem_change])
#plt.ylabel("Fractional change in stem length after random insertion")
colors =  ["#92c5de", "#ca0020"]
for patch,color in zip(vplot['bodies'], colors):
    patch.set_facecolor(color)

mplot.set_figwidth(3)
mplot.set_figheight(3)
#plt.scatter(orig_control_stem_lengths, mod_control_stem_lengths,s=5,color = "#ca0020")
#plt.scatter(orig_T1000_stem_lengths, mod_T1000_stem_lengths,s=5,color = "#92c5de")
#plt.xlabel('Stem Length of Hairpin before Deletion')
#plt.ylabel('Stem Length of Hairpin after Deletion')
plt.savefig(figure_file)
