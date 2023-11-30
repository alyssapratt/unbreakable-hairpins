import sys, re, math, subprocess
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from Bio import SeqIO


def create_RY_split_dict(ctrl_fasta_file):
    RY_split_dict = {}
    split_hairpins = []
    for record in SeqIO.parse(ctrl_fasta_file,'fasta'):
        rna_id = record.id
        seq = str(record.seq)
        RY_split = perfect_RY_split(seq)
        if RY_split:
        RY_split_dict[rna_id] = {"Seq": seq, "Unbreakable": False, "Split": is_split}
        split_hairpins.append(rna_id)
    return RY_split_dict,split_hairpins

def amend_RY_split_dict(filtered_fasta_file, RY_split_dict):
    unbreakable_hairpins = []
    for record in SeqIO.parse(filtered_fasta_file,'fasta'):
        rna_id = record.id
        seq = str(record.seq)
        seq_dict = RY_split_dict[rna_id]
        if seq_dict != None:
            seq_dict["Unbreakable"] = True
        unbreakable_hairpins.append(rna_id)
    return RY_split_dict,unbreakable_hairpins

def RY_split_score(seq):
    if onlyACGU(seq):
        split_point = int(len(seq)/2)
        left_Y_prop = (seq.count("C",0,split_point) + seq.count("U",0,split_point))/float(split_point)
        right_R_prop = (seq.count("A",split_point) + seq.count("G",split_point))/float(split_point)
        return (left_Y_prop - 0.5)**2 + (right_R_prop - 0.5)**2
    else:
        return None

def perfect_RY_split(seq, dot_bracket):
    if dot_bracket.count('(') == 1 and dot_bracket.count(')') == 1 and onlyACGU(seq):
        loop_match = dot_bracket.search("\(.\)")
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
    else:
        return None

def onlyACGU(RNA):
    if re.search('[^ACGU]',str(RNA)):
        return False
    else:
        return True


##########
# MAIN   #
##########

usage = "Usage: " + sys.argv[0] + " <test fasta> <ctrl/background fasta>"
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

test_fasta_file = sys.argv[1]
ctrl_fasta_file = sys.argv[2]

RY_split_dict,split_hairpins = create_RY_split_dict(ctrl_fasta_file)
RY_split_dict,unbreakable_hairpins =  amend_RY_split_dict(test_fasta_file, RY_split_dict)


# plot figure
plt.figure()
venn2([set(split_hairpins), set(unbreakable_hairpins)], ("Purine/Pyrimidine Split Hairpins", "Unbreakable Hairpins"))
plt.savefig("RY_T50_venn_diagram.jpg")
