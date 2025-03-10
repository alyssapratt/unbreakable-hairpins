import sys, re, math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy import ndimage
from Bio import SeqIO
import subprocess

def run_RNAfold(RNA):
    output = subprocess.check_output("echo " + str(RNA) + " | /nfs0/Hendrix_Lab/bin/RNAfold --noPS", shell=True, universal_newlines=True)
    lines = output.split('\n')
    dot_bracket = lines[1].split(' ')[0]
    return dot_bracket

def create_RY_split_stats(fasta_file, ID_list):
    split_list = []
    left_Y_proportions = []
    right_R_proportions = []
    IDs = []
    for record in SeqIO.parse(fasta_file,'fasta'):
        rna_id = record.id
        seq = str(record.seq)
        dot_bracket = run_RNAfold(seq)
        if is_valid(seq,dot_bracket):
            RY_split,left_prop,right_prop = RY_split_score(seq, dot_bracket)
            split_list.append(RY_split)
            if rna_id not in ID_list:
                IDs.append(rna_id)
            left_Y_proportions.append(left_prop)
            right_R_proportions.append(right_prop)
    return split_list, left_Y_proportions, right_R_proportions, IDs

def RY_split_score(seq, dot_bracket):
    loop_match = re.search("\(\.*\)", dot_bracket)
    loop_start = loop_match.start()
    loop_end = loop_match.end()
    left = seq[0:loop_start]
    right = seq[loop_end:]
    left_Y_prop = (left.count("C") + left.count("U"))/float(len(left))
    left_R_prop = (left.count("A") + left.count("G"))/float(len(left))
    right_Y_prop = (right.count("C") + right.count("U"))/float(len(right))
    right_R_prop = (right.count("A") + right.count("G"))/float(len(right))
    # idea 1:
    score1 = (left_Y_prop - 0.5)**2 + (right_R_prop - 0.5)**2
    score2 = (left_R_prop - 0.5)**2 + (right_Y_prop - 0.5)**2
    #score = max(score1,score2)
    # idea 2:
    score = max(left_Y_prop,left_R_prop) + max(right_Y_prop,right_R_prop)
    #score2 = max(left_Y_prop+right_R_prop,left_R_prop+right_Y_prop)
    return (score,max(left_Y_prop,left_R_prop),max(right_Y_prop,right_R_prop))

def RY_axes_list(seq):
    if onlyACGU(seq):
        split_point = int(len(seq)/2)
        left = seq[0:split_point]
        right = seq[split_point:]
        left_Y_prop = (left.count("C") + left.count("U"))/float(len(left))
        left_R_prop = (left.count("A") + left.count("G"))/float(len(left))
        right_Y_prop = (right.count("C") + right.count("U"))/float(len(right))
        right_R_prop = (right.count("A") + right.count("G"))/float(len(right))
        return left_Y_prop, right_R_prop
    else:
        return None

def is_valid(seq, dot_bracket):
    loops = re.findall("\(\.*\)", dot_bracket)
    if onlyACGU(seq) and len(loops) ==1:
        return True
    return False

def onlyACGU(RNA):
    if re.search('[^ACGU]',str(RNA)):
        return False
    else:
        return True


##########
# MAIN   #
##########

usage = "Usage: " + sys.argv[0] + " <test fasta> <ctrl/background fasta> <scatterplot figure name>"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

test_fasta_file = sys.argv[1]
ctrl_fasta_file = sys.argv[2]
fig_name = sys.argv[3]

test_split_list,test_left_proportions,test_right_proportions,test_IDs = create_RY_split_stats(test_fasta_file, [])
ctrl_split_list,ctrl_left_proportions,ctrl_right_proportions,ctrl_IDs = create_RY_split_stats(ctrl_fasta_file, test_IDs)

avg_test_split_score = sum(test_split_list) / len(test_split_list)
avg_ctrl_split_score = sum(ctrl_split_list) / len(ctrl_split_list)

# create summary file:
output_file = "RY_split_list_T50.txt"
F = open(output_file, 'w')
F.write("Test: Number of only-ACGU sequences: " + str(len(test_split_list)) + '\n')
F.write("Test: Average score: " + str(avg_test_split_score) + '\n')
F.write("Ctrl: Number of only-ACGU sequences: " + str(len(ctrl_split_list)) + '\n')
F.write("Ctrl: Average score: " + str(avg_ctrl_split_score) + '\n')
for number in test_split_list:
    F.write(str(number) + '\n')

# plot figures
plt.figure()
plt.scatter(ctrl_left_proportions, ctrl_right_proportions, color="#2c7bb6")
plt.scatter(test_left_proportions, test_right_proportions, color="#d7191c")
plt.xlabel("Maximum left proportion")
plt.ylabel("Maximum right proportion")
plt.savefig(fig_name)

plt.figure()
plt.hist(test_split_list, normed=1, bins=100,cumulative=False,fill=False,histtype='step',label='RY Split Scores of T50 Hairpins')
plt.hist(ctrl_split_list,normed=1,bins=100,cumulative=False,fill=False,histtype='step',label='RY Split Scores of All Hairpins')
plt.xlabel('Purine/Pyrimidine Split Score')
plt.ylabel('Distribution')
plt.legend(loc='upper right')
plt.title('Histogram of RY split scores')
plt.savefig('hist_RY_split_scores.pdf')

F.close()
