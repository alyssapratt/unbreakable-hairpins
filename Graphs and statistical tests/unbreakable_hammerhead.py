import random, re, subprocess, sys
from altschulEriksonDinuclShuffle import dinuclShuffle

###############
# SUBROUTINES #
###############

def run_RNAfold(RNA):
    output = subprocess.check_output("echo " + str(RNA) + " | RNAfold --noPS", shell=True, universal_newlines=True)
    lines = output.split('\n')
    dot_bracket = lines[1].split(' ')[0]
    return dot_bracket

def complement(seq):
    complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
    bases = list(seq)
    bases = [complement[base] for base in bases]
    return ''.join(bases)

def reverse_complement(s):
    return complement(s[::-1])

def get_random_base():
    bases = ['A', 'U', 'G', 'C']
    base_index = random.randint(0,3)
    return bases[base_index]

def construct_hammerhead():
    A = ""
    B = ""
    C = ""
    X1 = ""
    X2 = ""
    X3 = ""
    X4 = ""
    X5 = ""
    for i in range(0,5):
        A += get_random_base()
        B += get_random_base()
        C += get_random_base()
        X2 += get_random_base()
        X4 += get_random_base()
    for i in range(0,2):
        X1 += get_random_base()
        X3 += get_random_base()
        X5 += get_random_base()
    A_comp = reverse_complement(A)
    B_comp = reverse_complement(B)
    C_comp = reverse_complement(C)
    return(A+X1+B+X2+B_comp+X3+C+X4+C_comp+X5+A_comp)

def is_unbreakable(RNA):
    DNA = RNA.replace('U','T')
    shuffled_DNA = dinuclShuffle(DNA)
    shuffled_dot_bracket = run_RNAfold(shuffled_DNA)
    hammerhead_match = re.search("\(\.*\(+\.*\)+\.*\(+\.*\)+\.*\)", shuffled_dot_bracket)
    shuffle_count = 1
    while hammerhead_match is not None and shuffle_count < max_shuffle_count:
        shuffled_DNA = dinuclShuffle(DNA)
        shuffled_dot_bracket = run_RNAfold(shuffled_DNA)
        hammerhead_match = re.search("\(\.*\(+\.*\)+\.*\(+\.*\)+\.*\)", shuffled_dot_bracket)
        shuffle_count += 1
        if shuffle_count >= max_shuffle_count:
            return True
    return False

##########
# MAIN   #
##########
max_shuffle_count = 10
outputFile = sys.argv[1]
F = open(outputFile,'w')

usage = "Usage: " + sys.argv[0] + " <hammerhead list file name>"
if len(sys.argv) != 2:
    print(usage)
    sys.exit()

for i in range(0,50000):
    rna_hammerhead = construct_hammerhead()
    if is_unbreakable(rna_hammerhead):
        F.write(rna_hammerhead+ '\n')

F.close()
