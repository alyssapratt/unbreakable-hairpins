import sys, re, math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)

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

def filter_control_length(control_lengths,control_mfes,length_limit):
    filtered_mfes = []
    filtered_lengths=[]
    for i in range(len(control_lengths)):
        if control_lengths[i] <= length_limit:
            filtered_mfes.append(control_mfes[i])
            filtered_lengths.append(control_lengths[i])
    return filtered_mfes,filtered_lengths

######
#Main#
######

usage = "usage: " + sys.argv[0] + " <unbreakable length/mfe info file> <control hairpins length/mfe info file> <scatterplot filename>"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

unbreakable_file = sys.argv[1]
control_file = sys.argv[2]
figure_file = sys.argv[3]
figure_file_no_type = figure_file[0:-4]

unbreakable_rna_ids,unbreakable_lengths,unbreakable_mfes = get_mfes(unbreakable_file)
control_rna_ids,control_lengths,control_mfes = get_mfes_control(unbreakable_rna_ids, control_file)

unbreakable_mfes,unbreakable_lengths = filter_control_length(unbreakable_lengths,unbreakable_mfes,2000)
control_mfes,control_lengths = filter_control_length(control_lengths,control_mfes,2000)

limited_unbreakable_mfes,limited_unbreakable_lengths = filter_control_length(unbreakable_lengths,unbreakable_mfes,25)
limited_control_mfes,limited_control_lengths = filter_control_length(control_lengths,control_mfes,25)

plt.figure()
fig,ax1=plt.subplots(figsize=(8, 6))
ax1.scatter(control_lengths,control_mfes, alpha=0.5,s=5,color = "#92c5de",label='Control Hairpins')
ax1.scatter(unbreakable_lengths,unbreakable_mfes,s=5,alpha=0.5, color = "#ca0020",label='Unbreakable Hairpins')
plt.xlabel('Hairpin Length (nt)')
plt.ylabel('Predicted Hairpin Minimum Free Energy (kcal/mol)')
ax1.legend(loc='upper right')

left, bottom, width, height = [0.2, 0.2, 0.2, 0.2]
inset_ax=fig.add_axes([left, bottom, width, height])
mark_inset(ax1, inset_ax, loc1=2, loc2=1, fc="none", ec='0.5')
inset_ax.scatter(limited_control_lengths,limited_control_mfes, alpha=0.5,s=5,color = "#92c5de",label='Control Hairpins')
inset_ax.scatter(limited_unbreakable_lengths,limited_unbreakable_mfes,s=5,alpha=0.5, color = "#ca0020",label='Unbreakable Hairpins')





plt.savefig(figure_file)
