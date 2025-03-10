import sys, re, math, subprocess, random
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

def make_deletion(seq):
    randindex = random.randrange(len(seq))
    modified_seq = seq[:randindex] + seq[randindex+1:]
    return modified_seq

def onlyACGU(RNA):
    if re.search('[^ACGU]',str(RNA)):
        return False
    else:
        return True

def make_hairpin_deletions(fasta_file):
    valid_modified_hairpins = []
    invalid_modified_hairpins = []
    for record in SeqIO.parse(fasta_file,'fasta'):
        rna_id = record.id
        seq = str(record.seq)
        dot_bracket = run_RNAfold(seq)
        if is_valid(seq, dot_bracket):
            modified_seq = make_deletion(seq)
            modified_dot_bracket = run_RNAfold(modified_seq)
            if is_valid(modified_seq, modified_dot_bracket):
                valid_modified_hairpins.append((rna_id,modified_seq,modified_dot_bracket,dot_bracket))
            else:
                invalid_modified_hairpins.append((rna_id,modified_seq,modified_dot_bracket,dot_bracket))
    return valid_modified_hairpins,invalid_modified_hairpins


##########
# MAIN   #
##########

usage = "Usage: " + sys.argv[0] + " <test fasta> <output file>"
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

fasta_file = sys.argv[1]
output_name = sys.argv[2]

valid_modified_hairpins,invalid_modified_hairpins = make_hairpin_deletions(fasta_file)

valid_count = len(valid_modified_hairpins)
invalid_count = len(invalid_modified_hairpins)

# create summary file:
output_file = output_name
F = open(output_file, 'w')
F.write("#Valid\n")
for hairpin_tuple in valid_modified_hairpins:
    F.write("%s\t%s\t%s\t%s\n" % (hairpin_tuple[0], hairpin_tuple[1], hairpin_tuple[2], hairpin_tuple[3]))
F.write("#Invalid\n")
for hairpin_tuple in invalid_modified_hairpins:
    F.write("%s\t%s\t%s\t%s\n" % (hairpin_tuple[0], hairpin_tuple[1], hairpin_tuple[2], hairpin_tuple[3]))
F.close()

print("Valid hairpins: %d/%d" % (valid_count, invalid_count+valid_count))
