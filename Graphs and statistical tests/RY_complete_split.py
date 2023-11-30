import sys, re, math, subprocess
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from Bio import SeqIO

###############
# SUBROUTINES #
###############

def run_RNAfold(RNA):
    output = subprocess.check_output("echo " + str(RNA) + " | RNAfold --noPS", shell=True, universal_newlines=True)
    lines = output.split('\n')
    dot_bracket = lines[1].split(' ')[0]
    return dot_bracket

def is_valid(seq, dot_bracket):
    loops = re.findall("\(\.*\)", dot_bracket)
    if onlyACGU(seq) and len(loops) ==1:
        return True
    return False

def is_split(seq, dot_bracket):
    i = 0
    crosses =
    if seq.count("AC")


    loop_match = re.search("\(\.*\)", dot_bracket)
    loop_start = loop_match.start()
    loop_end = loop_match.end()
    stem_left = seq[0:loop_start]
    stem_right = seq[loop_end:]
    if stem_left.count("A") + stem_left.count("G") == len(stem_left):
        if stem_right.count("C") + stem_right.count("U") == len(stem_right):
            return True
        else:
            return False
    elif stem_left.count("C") + stem_left.count("U") == len(stem_left):
        if stem_right.count("A") + stem_right.count("G") == len(stem_right):
            return True
        else:
            return False
    else:
        return False

def onlyACGU(RNA):
    if re.search('[^ACGU]',str(RNA)):
        return False
    else:
        return True

def get_split_hairpins(ctrl_fasta_file):
    split_hairpins = []
    for record in SeqIO.parse(ctrl_fasta_file,'fasta'):
        rna_id = record.id
        seq = str(record.seq)
        dot_bracket = run_RNAfold(seq)
        if is_valid(seq, dot_bracket) and is_split(seq, dot_bracket):
            split_hairpins.append(rna_id)
    return split_hairpins

def get_unbreakable_hairpins(filtered_fasta_file):
    unbreakable_hairpins = []
    for record in SeqIO.parse(filtered_fasta_file,'fasta'):
        rna_id = record.id
        seq = str(record.seq)
        unbreakable_hairpins.append(rna_id)
    return unbreakable_hairpins


##########
# MAIN   #
##########

usage = "Usage: " + sys.argv[0] + " <test fasta> <ctrl/background fasta> <venn diagram figure name>"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

test_fasta_file = sys.argv[1]
ctrl_fasta_file = sys.argv[2]
fig_name = sys.argv[3]


# test_split_list,test_number_split, test_number_not_split = create_RY_split_stats(test_fasta_file)
# ctrl_split_list,ctrl_number_split,ctrl_number_not_split = create_RY_split_stats(ctrl_fasta_file)

split_hairpins = get_split_hairpins(ctrl_fasta_file)
unbreakable_hairpins = get_unbreakable_hairpins(test_fasta_file)

# create summary file:
output_file = "RY_perfect_split_list_T50.txt"
F = open(output_file, 'w')
F.write("Number of split hairpins: " + str(len(split_hairpins)) + '\n')
F.write("Number of unbreakable hairpins: " + str(len(unbreakable_hairpins)) + '\n')
F.write("Split Hairpins:\n")
for value in set(split_hairpins):
    F.write(str(value) + '\n')
F.write("Unbreakable Hairpins:\n")
for value in set(unbreakable_hairpins):
    F.write(str(value) + '\n')


# plot figure
plt.figure()
venn2([set(split_hairpins), set(unbreakable_hairpins)], ("Purine/Pyrimidine Split Hairpins", "Unbreakable Hairpins"))
plt.savefig(fig_name)

F.close()
