import sys, re, math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy import ndimage
from Bio import SeqIO
import subprocess

#############
#Subroutines#
#############

def mfe_from_RNAfold(RNA):
    output = subprocess.check_output("echo " + str(RNA) + " | /nfs0/Hendrix_Lab/bin/RNAfold --noPS", shell=True, universal_newlines=True)
    lines = output.split('\n')
    mfe_str = lines[1].split('-')
    if(len(mfe_str) > 1):
        mfe_str = mfe_str[1].strip(')')
    else:
        return None
    mfe = float(mfe_str)
    mfe *= -1
    return mfe

def get_mfes(tabbed_file):
    hairpin_lengths = []
    hairpin_mfes = []
    with open(tabbed_file, 'r') as F:
        for line in F:
            rna_id,shuffles,seq = line.strip().split('\t')
            mfe = mfe_from_RNAfold(seq);
            length = len(seq)
            if (mfe is not None) and (length is not None):
                hairpin_lengths.append(length)
                hairpin_mfes.append(mfe)
        F.close()
    return hairpin_lengths,hairpin_mfes

def get_fasta_mfes(fasta_file):
    hairpin_lengths = []
    hairpin_mfes = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
            rna_id = record.id
            seq = str(record.seq)
            mfe = mfe_from_RNAfold(seq);
            length = len(seq)
            if (mfe is not None) and (length is not None):
                hairpin_lengths.append(length)
                hairpin_mfes.append(mfe)
    return hairpin_lengths,hairpin_mfes

######
#Main#
######

usage = "usage: " + sys.argv[0] + " <tabbed file> <figure filename>"
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

input_file = sys.argv[1]
output_file = sys.argv[2]

if input_file[-6:] == ".fasta":
    hairpin_lengths,hairpin_mfes = get_fasta_mfes(input_file)
else:
    hairpin_lengths,hairpin_mfes = get_mfes(input_file)


plt.figure()
plt.scatter(hairpin_lengths,hairpin_mfes, color = "#006A8E")
plt.xlabel('Hairpin Length')
plt.ylabel('Hairpin Minimum Free Energy calculation')
plt.savefig(output_file)
