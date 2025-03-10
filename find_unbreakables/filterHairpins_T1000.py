#!/local/cluster/bin/python
import sys, os, re, subprocess
from Bio import SeqIO
from altschulEriksonDinuclShuffle_forpy3 import dinuclShuffle

#################
## SUBROUTINES ##
#################

def read_hairpin_fasta_file(hairpin_fasta_file):
    hairpins = {}
    sequences = SeqIO.parse(hairpin_fasta_file,'fasta')
    for record in sequences:
        hairpins[record.id] = record.seq
    return hairpins

def run_RNAfold(RNA):
    output = subprocess.check_output("echo " + str(RNA) + " | RNAfold", shell=True, universal_newlines=True)
    lines = output.split('\n')
    dot_bracket = lines[1].split(' ')[0]
    return dot_bracket

def compute_basepair_density(dot_bracket):
    num_bp = dot_bracket.count('(')
    length= len(dot_bracket)
    bp_density = float(num_bp)/length
    return bp_density

def compute_balance_score(dot_bracket):
    left_bracket_count = float(dot_bracket.count('('))
    right_bracket_count = float(dot_bracket.count(')'))
    if left_bracket_count != right_bracket_count:
        print("Unbalanced brackets.\n")
        sys.exit()
    if left_bracket_count == 0:
        return 0.0
    half_length = round(len(dot_bracket)//2)
    left_half = dot_bracket[0:int(half_length)]
    right_half = dot_bracket[int(half_length):]
    left_on_left = float(left_half.count('('))
    right_on_right = float(right_half.count(')'))
    balance_score = float((left_on_left/left_bracket_count) + (right_on_right/right_bracket_count))
    return balance_score

def onlyACGU(RNA):
    if re.search('[^ACGU]',str(RNA)):
        return False
    else:
        return True


##########
## MAIN ##
##########

usage = "Usage: " + sys.argv[0] + " <hairpin fasta file>"
if len(sys.argv) != 2:
    print(usage)
    sys.exit()

hairpin_fasta_file = sys.argv[1]
BPD_threshold = 0.25
balance_threshold = 1.86
max_shuffle_count = 1000

hairpins = read_hairpin_fasta_file(hairpin_fasta_file)

outputFile = "T1000_list.txt"
F = open(outputFile,'w')

for hpID in hairpins:
    # if no Ns in sequence, only A,C,G,U
    if onlyACGU(hairpins[hpID]):
        dot_bracket = run_RNAfold(hairpins[hpID])
        BPD = compute_basepair_density(dot_bracket)
        balance_score = compute_balance_score(dot_bracket)
        # IF it passes "hairpin thresholds" i.e. is a hairpin
        if balance_score >= balance_threshold and BPD >= BPD_threshold:
            DNA = str(hairpins[hpID]).replace('U','T')
            # re-initialize the shuffle count
            shuffle_count = 1
            shuffled_DNA = dinuclShuffle(DNA)
            shuffled_dot_bracket = run_RNAfold(shuffled_DNA)
            shuffled_BPD = compute_basepair_density(shuffled_dot_bracket)
            shuffled_balance_score = compute_balance_score(shuffled_dot_bracket)
            # while this shuffled RNA passes "hairpin thresholds" i.e. is a hairpin
            while shuffled_balance_score >= balance_threshold and shuffled_BPD >= BPD_threshold:
                # create new shuffle, increment shuffle count
                shuffled_DNA = dinuclShuffle(DNA)
                shuffle_count += 1
                # recompute the dotbracket, and BPD etc
                shuffled_dot_bracket = run_RNAfold(shuffled_DNA)
                shuffled_BPD = compute_basepair_density(shuffled_dot_bracket)
                shuffled_balance_score = compute_balance_score(shuffled_dot_bracket)
                # exit while-loop because this is no longer forms a hairpin
                if shuffled_balance_score <= balance_threshold or shuffled_BPD <= BPD_threshold:
                    break
                # exit of shuffled too many times
                if shuffle_count >= max_shuffle_count:
                    F.write("%s\t%d\t%s\n" % (hpID, shuffle_count, hairpins[hpID]))
                    break
            # at this point, at some point the shuffle broke the HP, or max_shuffle_count was reached
            #print hpID, BPD, balance_score, shuffle_count

F.close()
